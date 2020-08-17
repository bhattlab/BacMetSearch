from os import makedirs
from os.path import isdir, join
import shutil
import click
import sys
import os
from os import remove
from bacmetsearch.prodigal import run_prodigal_multithread
from bacmetsearch.diamond import run_diamond, parse_diamond_search
from bacmetsearch.misc import *
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


    diamond_start = time.time()
    diamond_exp_path = join(tmpdir, 'diamond.exp.tsv')
    run_diamond(proteome_path, diamond_exp_path, threads, max_target_seqs, min_percent_identity, BACMET2_EXPERIMENTAL_DMND)
    diamond_exp_results = parse_diamond_search(diamond_exp_path)

    diamond_pred_path = join(tmpdir, 'diamond.pred.tsv')
    run_diamond(proteome_path, diamond_pred_path, threads, max_target_seqs, min_percent_identity, BACMET2_PREDICTED_DMND)
    diamond_pred_results = parse_diamond_search(diamond_pred_path)
    diamond_end = time.time()


    bacmet_exp_meta = parse_bacmet_exp_metadata()
    outfile = open(join(outdir, prefix+'.results.tsv'), 'w')
    print('protein_id', 'BacMet_ID', 'Gene_name', 'Compound', sep="\t", file=outfile)
    for res in diamond_exp_results:
        bacmet_id = diamond_exp_results[res]['sseqid'].split('|')[0]
        print(res, bacmet_id, bacmet_exp_meta[bacmet_id]['Gene_name'], bacmet_exp_meta[bacmet_id]['Compound'], sep="\t", file=outfile)
        print(res, diamond_exp_results[res])
    outfile.close()


    if not keep_intermediate:
        shutil.rmtree(tmpdir)

    print()
    print("COMPLETE.")
    print("Prodigal runtime: %f" % (prodigal_end-prodigal_start))
    print("Diamond runtime: %f" % (diamond_end - diamond_start))