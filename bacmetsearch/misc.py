from bacmetsearch import *


def parse_bacmet_exp_database():

    outdict = dict()
    with open(BACMET2_EXPERIMENTAL_META) as infile:
        header = infile.readline().strip().split('\t')
        for line in infile:
            line = line.strip().split('\t')
            out = {header[i]:line[i] for i in range(len(line))}
            outdict[line['BacMet_ID']] = out
    return outdict


def parse_bacmet_pred_database():

    outdict = dict()
    with open(BACMET2_PREDICTED_META) as infile:
        header = infile.readline().strip().split('\t')
        for line in infile:
            line = line.strip().split('\t')
            out = {header[i]:line[i] for i in range(len(line))}
            outdict[line['GI_number']] = out

    return outdict