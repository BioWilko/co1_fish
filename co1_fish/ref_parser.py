import pysam
import sys
import os
import click
from types import SimpleNamespace
import subprocess


def reference_filter(args):
    with open(args.big_fasta, 'r') as fh:
        filedata = fh.read()

        if " " in filedata:
            # Replace the target string
            filedata = filedata.replace(' ', '_')
            with open(args.big_fasta, 'w') as fh:
                fh.write(filedata)

    if not os.path.exists(str(args.big_fasta) + ".fai"):
        subprocess.run(["samtools", "faidx", args.big_fasta])

    in_fasta_fh = pysam.FastaFile(args.big_fasta)

    taxon_dict = {}

    taxons = (header.split("|")[1] for header in in_fasta_fh.references)

    out_entries = []

    for taxon in taxons:
        headers = (
            reference
            for reference in in_fasta_fh.references
            if reference.split("|")[1] == taxon
        )
        taxon_dict[taxon] = headers

    def header_stat(taxon, query_header, in_fasta=in_fasta_fh):
        entry = in_fasta.fetch(query_header)
        good_bases = sum(entry.count(x) for x in ["A", "C", "T", "G"])
        return (taxon, entry, good_bases)

    for taxon, headers in taxon_dict.items():
        best = ("", "", 0)
        for header in headers:
            test_stat = header_stat(taxon, header)
            if test_stat[2] > best[2]:
                best = test_stat
        out_entries.append((header, best[1]))

    for entry in out_entries:
        header = entry[0]
        sys.stdout.write(
            f">{header}\n{entry[1]}\n"
        )  # NTS -> Use ncbi taxon codes and make fasta writing less brittle

@click.command()
@click.argument("big_fasta", type=click.Path())
def main(*_, **kwargs):
    args = SimpleNamespace(**kwargs)
    reference_filter(args)


if __name__ == "__main__":
    main()
