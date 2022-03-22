import mappy as mp
import click
import sqlite3
import numpy as np
import pandas as pd
import sys


def alignment_hit_generator(args):
    aligner_preset = "map-ont"
    ref_idx = mp.Aligner(args.reference, preset=aligner_preset)
    if not ref_idx:
        click.ClickException("[ERROR] Unable load / build index of reference provided")

    test_sequences = mp.fastx_read(args.query)
    for header, seq, read_qual in test_sequences:
        for hit in ref_idx.map(seq):
            taxon = hit.ctg.split("|")[1]
            yield taxon, hit, read_qual


def sql_dict_factory(cursor, row):
    out_dict = {}
    for idx, col in enumerate(cursor.description):
        out_dict[col[0]] = row[idx]
    return out_dict


def init_db():
    db = sqlite3.connect(":memory:")
    db.row_factory = sql_dict_factory
    cursor = db.cursor()
    cursor.execute(
        "CREATE TABLE taxons(taxon_id INTEGER PRIMARY KEY AUTOINCREMENT, taxon TEXT, ref_len INTEGER, CONSTRAINT taxon_unique UNIQUE(taxon));"
    )
    cursor.execute(
        "CREATE TABLE alignments(alignment_id INTEGER PRIMARY KEY AUTOINCREMENT, taxon_id INTEGER, ref_start INTEGER, ref_end INTEGER, \
        query_start INTEGER, query_end INTEGER, strand INTEGER, read_qual TEXT, cigar TEXT, mapping_qual INTEGER, \
            FOREIGN KEY (taxon_id) REFERENCES taxons(taxon_id) ON DELETE CASCADE);"
    )
    return db, cursor


def populate_db(args):
    hits = alignment_hit_generator(args)
    db, cursor = init_db()
    for taxon, hit, read_qual in hits:
        ref_len = hit.ctg_len
        ref_start = hit.r_st
        ref_end = hit.r_en
        query_start = hit.q_st
        query_end = hit.q_en
        strand = hit.strand
        mapping_quality = hit.mapq
        cigar = hit.cigar_str
        cursor.execute(
            f"INSERT OR IGNORE INTO taxons(taxon, ref_len) VALUES('{taxon}', '{ref_len}');"
        )
        cursor.execute(f"SELECT * FROM taxons WHERE taxon='{taxon}';")
        taxon_id = cursor.fetchone()["taxon_id"]
        cursor.execute(
            f"INSERT INTO alignments(taxon_id, ref_start, ref_end, query_start, query_end, strand, cigar, mapping_qual) \
            VALUES('{taxon_id}', '{ref_start}', '{ref_end}', '{query_start}', '{query_end}', '{strand}', '{cigar}', '{mapping_quality}');"
        )
    return db, cursor


def get_ref_coverage(ref_len, alignments):
    overlaps = np.zeros(ref_len, dtype=int)
    for alignment in alignments:
        for idx in range(alignment["ref_start"] - 1, alignment["ref_end"] - 1):
            overlaps[idx] += 1
    print(overlaps)
    return np.var(overlaps)


def get_mapping_qualities(alignments):
    mapping_quals = [alignment["mapping_qual"] for alignment in alignments]
    return sum(mapping_quals) / len(mapping_quals)


def get_strand_bias_metric(alignments):
    strands = [alignment["strand"] for alignment in alignments]
    return 1 - abs(
        sum(strands) / len(strands)
    )  # Get the distance of the bias from 0 -> ideally you have an even number of forward and reverse alignments


def get_high_quality_hits_proportion(alignments, quality_threshold):
    high_quality_hits = [
        alignment
        for alignment in alignments
        if alignment["mapping_qual"] >= quality_threshold
    ]
    if len(high_quality_hits) == 0:
        return 0
    return len(high_quality_hits) / len(alignments)


def metrics_to_probability(taxon_metrics, metric_weights):
    considered_metrics = [taxon_metrics[metric] for metric in metric_weights.keys()]
    weight_array = [weight for weight in metric_weights.values()]
    return np.average(considered_metrics, weights=weight_array)


def generate_taxon_metrics(cursor, taxons_with_alignments):
    taxon_metrics = []
    weights = {
        "normalised_mean_map_qual": 0.5,
        "strand_bias": 0.25,
        "high_quality_proportion": 1,
    }
    for taxon in taxons_with_alignments:
        cursor.execute(
            f"SELECT * FROM taxons INNER JOIN alignments ON taxons.taxon_id = alignments.taxon_id WHERE alignments.taxon_id='{taxon['taxon_id']}';"
        )
        alignment_matches = cursor.fetchall()
        n_hits = len(alignment_matches)
        mean_map_qual = get_mapping_qualities(alignment_matches)
        normalised_mean_map_qual = mean_map_qual / 60
        strand_bias = get_strand_bias_metric(alignment_matches)
        high_quality_proportion = get_high_quality_hits_proportion(
            alignment_matches, 50
        )
        taxon_metric_dict = {
            "taxon": taxon["taxon"],
            "n_hits": n_hits,
            "mean_map_qual": mean_map_qual,
            "normalised_mean_map_qual": normalised_mean_map_qual,
            "strand_bias": strand_bias,
            "high_quality_proportion": high_quality_proportion,
        }
        taxon_metric_dict["probability"] = metrics_to_probability(
            taxon_metric_dict, weights
        )
        taxon_metrics.append(taxon_metric_dict)
    return taxon_metrics


def generate_report(taxon_metrics):
    metric_df = pd.DataFrame(data=taxon_metrics)
    metric_df = metric_df.set_index("taxon")

    metric_df.sort_values(by="probability", ascending=False, inplace=True)
    metric_df.to_csv(sep="\t", path_or_buf=sys.stdout, float_format="%.3f")

