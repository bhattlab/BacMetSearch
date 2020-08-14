from os import makedirs
from os.path import isdir, join
import shutil
import click
import sys
import os
from os import remove
from bacmetsearch.prodigal import run_prodigal_multithread
from bacmetsearch.diamond import run_diamond
from bacmetsearch import *
from Bio import SeqIO
from subprocess import run, DEVNULL
from collections import defaultdict
import sqlite3
from glob import glob
import numpy as np
import time
from multiprocessing import Pool
import pickle


def _meta(fasta, outdir, prefix, force, threads, max_target_seqs, min_percent_identity,
             keep_intermediate, fasta_type):

    tmpdir = join(outdir, 'tmp')

    if force and isdir(outdir):
        shutil.rmtree(outdir)
    try:
        makedirs(outdir)
        makedirs(tmpdir)
    except FileExistsError:
        click.echo("Output directory exists, please delete or overwrite with --force")
        sys.exit(1)

    prodigal_start = time.time()
    if fasta_type == 'genome':
        click.echo("Running prodigal...")
        run_prodigal_multithread(PRODIGAL_PATH, fasta, tmpdir, threads)
        proteome_path = join(tmpdir, 'prodigal.faa')
    elif fasta_type == 'proteome':
        proteome_path = fasta
    prodigal_end = time.time()

    diamond_path = join(tmpdir, 'diamond.tsv')
    run_diamond(proteome_path, diamond_path, threads, max_target_seqs, min_percent_identity)

    if not keep_intermediate:
        shutil.rmtree(tmpdir)

    print()
    print("COMPLETE.")
    print("Prodigal runtime: %f" % (prodigal_end-prodigal_start))
    print("Marker gene runtime: %f" % (marker_gene_end - marker_gene_start))
    print("Closest genome runtime: %f" % closest_genomes_time)
    print("Gene count runtime: %f" % gene_count_time)