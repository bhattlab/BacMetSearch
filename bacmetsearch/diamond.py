from bacmetsearch import *

def run_diamond(protein_fasta_path, outfile, threads, max_target_seqs, min_percent_identity, database):
    command = '{0} blastp --id {5} -k {6} --query {1} --out {2}.dmd.tsv --outfmt 6 qseqid sseqid qlen slen pident ' \
              'length mismatch gapopen qstart qend sstart send evalue bitscore --db {3} --threads {4}'.format(
        DIAMOND_PATH, protein_fasta_path, outfile, database, threads, min_percent_identity, max_target_seqs)
    print('diamond command:', command)

    run(command.split(), stdout=DEVNULL, stderr=DEVNULL)


def parse_diamond_search(diamond_search):

    top_markers = dict()
    with open(diamond_search) as infile:
        for line in infile:
            qseqid, sseqid, qlen, slen, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore =  line.strip().split('\t')
            qlen, slen, qstart, qend, sstart, send = int(qlen), int(slen), int(qstart), int(qend), int(sstart), int(send)
            length, evalue, bitscore = int(length), float(evalue), float(bitscore)

    marker2gene = dict()
    gene2marker = dict()
    for rec in top_markers:
        marker2gene[rec] = top_markers[rec][0]
        gene2marker[top_markers[rec][0]] = rec

    records = []
    for rec in SeqIO.parse(protein_fasta_path, 'fasta'):
        if rec.id in gene2marker:
            rec.id = gene2marker[rec.id] + '__' + rec.id + '__' + prefix
            rec.description = rec.id
            records.append(rec)

    SeqIO.write(records, outfile, 'fasta')