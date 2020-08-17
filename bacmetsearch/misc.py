from bacmetsearch import *


def parse_bacmet_exp_metadata():

    outdict = dict()
    with open(BACMET2_EXPERIMENTAL_META) as infile:
        header = infile.readline().strip().split('\t')
        for line in infile:
            line = line.strip().split('\t')
            out = {header[i]:line[i] for i in range(len(line))}
            outdict[out['BacMet_ID']] = out
    return outdict


def parse_bacmet_pred_metadata():

    outdict = dict()
    with open(BACMET2_PREDICTED_META) as infile:
        header = infile.readline().strip().split('\t')
        for line in infile:
            line = line.strip().split('\t')
            out = {header[i]:line[i] for i in range(len(line))}
            outdict[out['GI_number']] = out

    return outdict