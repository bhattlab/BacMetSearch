import click
from bacmetsearch import *
from bacmetsearch.isolate import _isolate
from bacmetsearch.meta import _meta

@click.group()
def cli():
    """A command line tool to quickly search for closely related microbial genomes using a marker-gene based approach."""
    pass


@cli.command(short_help='Run bacmetsearch on a bacterial isolate genome.')
@click.argument('fasta', type=click.Path(exists=True))
@click.option('--outdir', '-o', default='genomesearch_isolate_uhgg_output', help='The name of the output directory.')
@click.option('--prefix', '-prefix', default='genomesearch', help='The prefix of all files in the output directory.')
@click.option('--force/--no-force', default=False, help="Force overwriting of output directory.")
@click.option('--threads', '-t', default=16, help="Number of threads to use for diamond searches.")
@click.option('--max-target-seqs', '-k', default=200, help="The maximum number of target seqs returned by the diamond search.")
@click.option('--min-percent-identity', '-mpi', default=50, help="The minimum percent identity to keep a marker match from the diamond search.")
@click.option('--keep-intermediate/--no-keep-intermediate', default=False, help="Keep intermediate files.")
@click.option('--fasta-type', '-ft', type=click.Choice(['genome', 'proteome', 'markers']), default='genome', help="Select the type of fasta input.")
def isolate(fasta, outdir, prefix, force, threads, max_target_seqs, min_percent_identity, keep_intermediate, fasta_type):
    """A click access point for the run module. This is used for creating the command line interface."""
    log_params(fasta=fasta, outdir=outdir, prefix=prefix, force=force, threads=threads,
               max_target_seqs=max_target_seqs, min_percent_identity=min_percent_identity,
               keep_intermediate=keep_intermediate, fasta_type=fasta_type)

    _isolate(fasta, outdir, prefix, force, threads, max_target_seqs, min_percent_identity, keep_intermediate, fasta_type)


@cli.command(short_help='Run bacmetsearch on a metagenomic assembly.')
@click.argument('fasta', type=click.Path(exists=True))
@click.option('--outdir', '-o', default='genomesearch_meta_refbank_output', help='The name of the output directory.')
@click.option('--prefix', '-prefix', default='genomesearch', help='The prefix of all files in the output directory.')
@click.option('--force/--no-force', default=False, help="Force overwriting of output directory.")
@click.option('--threads', '-t', default=16, help="Number of threads to use for diamond searches.")
@click.option('--max-target-seqs', '-k', default=20, help="The maximum number of target seqs returned by the diamond search.")
@click.option('--min-percent-identity', '-mpi', default=80, help="The minimum percent identity to keep a marker match from the diamond search.")
@click.option('--keep-intermediate/--no-keep-intermediate', default=False, help="Keep intermediate files.")
@click.option('--fasta-type', '-ft', type=click.Choice(['genome', 'proteome']), default='genome', help="Select the type of fasta input.")
def meta(fasta, outdir, prefix, force, threads, max_target_seqs, min_percent_identity, keep_intermediate, fasta_type):
    """A click access point for the run module. This is used for creating the command line interface."""
    log_params(fasta=fasta, outdir=outdir, prefix=prefix, force=force, threads=threads,
               max_target_seqs=max_target_seqs, min_percent_identity=min_percent_identity,
               keep_intermediate=keep_intermediate, fasta_type=fasta_type)

    _meta(fasta, outdir, prefix, force, threads, max_target_seqs, min_percent_identity, keep_intermediate, fasta_type)


def log_params(**kwargs):
    click.echo("#### PARAMETERS ####")
    click.echo('\n'.join(list(map(lambda x: ': '.join(list(map(str, x))), kwargs.items()))))
    click.echo("####################")


if __name__ == '__main__':

    cli()