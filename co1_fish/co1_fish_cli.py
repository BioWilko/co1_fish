import click
from types import SimpleNamespace

from . import co1_fish_funcs


@click.command()
@click.option("--hits-to-return", type=click.INT, default=10)
@click.option("--paired-data", default=False)
@click.argument("reference", type=click.Path())
@click.argument("query", type=click.Path())
def main(*_, **kwargs):
    args = SimpleNamespace(**kwargs)
    db, cursor = co1_fish_funcs.populate_db(args)
    cursor.execute("SELECT * FROM taxons;")
    aligned_taxons = cursor.fetchall()
    taxon_metric_list = co1_fish_funcs.generate_taxon_metrics(cursor, aligned_taxons)
    report_df = co1_fish_funcs.generate_report(taxon_metric_list)


if __name__ == "__main__":
    main()
