from bacmetsearch import *
from subprocess import run, DEVNULL


def run_diamond(protein_fasta_path, outfile, threads, max_target_seqs, min_percent_identity, database):
    command = '{0} blastp --no-auto-append --id {5} -k {6} --query {1} --out {2} --outfmt 6 ' \
              'qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore --db {3} ' \
              '--threads {4}'.format(
        DIAMOND_PATH, protein_fasta_path, outfile, database, threads, min_percent_identity, max_target_seqs)
    print('diamond command:', command)

    run(command.split(), stdout=DEVNULL, stderr=DEVNULL)


def parse_diamond_search(diamond_search):

    gene2bacmet = dict()
    with open(diamond_search) as infile:
        for line in infile:
            qseqid, sseqid, qlen, slen, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore =  line.strip().split('\t')
            qlen, slen, qstart, qend, sstart, send = int(qlen), int(slen), int(qstart), int(qend), int(sstart), int(send)
            length, evalue, bitscore = int(length), float(evalue), float(bitscore)


            qaln = qend - qstart
            saln = send - sstart

            out = {
                'qseqid': qseqid, 'sseqid': sseqid, 'qlen': qlen, 'slen': slen, 'pident': pident, 'length': length,
                'mismatch': mismatch, 'gapopen': gapopen, 'qstart': qstart, 'qend': qend,
                'sstart': sstart, 'send': send, 'evalue': evalue, 'bitscore': bitscore, 'qaln':qaln, 'saln':saln,
                'size_diff_prop': min([qlen, slen]) / max([qlen, slen]), 'alnlen_prop': alnlength / min([qlen, slen])
            }

            query_smaller = True
            if slen > qlen:
                query_smaller = False

            if evalue >= 1e-6:
                continue

            if out['size_diff_prop'] < 0.85:
                continue

            if query_smaller:
                alnlength = qaln
            else:
                alnlength = saln

            if out['alnlen_prop'] < 0.85:
                continue

            if qseqid not in gene2bacmet:
                gene2bacmet[qseqid] = out
            else:
                if out['bitscore'] > gene2bacmet[qseqid]['bitscore']:
                    gene2bacmet[qseqid] = out

    print(gene2bacmet)
    #SeqIO.write(records, outfile, 'fasta')