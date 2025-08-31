import typer
from ncbi_tool import gff_parser, ncbi_fetch, variant_lookup, pubmed_fetch, visualizer

app = typer.Typer(
    help="Bioinformatics Toolkit: Parse GFF files, fetch NCBI data, analyze variants, and more!"
)

app.add_typer(gff_parser.app, name="gff", help="Parse and summarize GFF/GTF files")
app.add_typer(ncbi_fetch.app, name="ncbi", help="Fetch and summarize NCBI data")
app.add_typer(variant_lookup.app, name="variant", help="Look up SNPs and ClinVar variants")
app.add_typer(pubmed_fetch.app, name="pubmed", help="Fetch and summarize PubMed papers")
app.add_typer(visualizer.app, name="viz", help="Visualize genomic data interactively")

if __name__ == "__main__":
    app()
