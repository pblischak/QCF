#!/usr/bin/env python

"""
<< qcf_boot.py >>

Bootstrap resample genes

Arguments
---------

    - infile <string>  : Name of input file with raw qcf data.
    - bootReps <int>   : Number of bootstrap replicates [default=100].
    - prefix <string>  : Prefix added to output file [default=boot].
    - quiet <bool>     : Flag to turn of printing to stdout.

Output
------

    Writes a CSV file with quartet concordance factors for each quartet,
    as well as the confidence intervals (2.5%-97.5%) estimated by bootstrap
    resampling the genes for a given quartet.
"""

from __future__ import print_function, division, unicode_literals
from sys import argv, exit
import argparse
import numpy as np

def parseRawQcfLine(line):
    """
    Take a line of raw QCF input (comma separated) and parse it into a
    key-value pair to be entered into the main dictionary storing all
    results. The format of a line is as follows:

    After splitting on commas, the first four entries are the names of the
    species in the quartet. After that, each entry corresponds to a gene,
    and the three numbers, separated by colons, represent the scores for each
    of the possible topologies.

    Return: The dictionary key returned is a four-tuple with the species' names.
    The value for the dictionary is a nested numpy arrays with dimensions
    (# genes)-by-3.
    """
    tokens = line.strip().split(",")
    quartet = (tokens[0], tokens[1], tokens[2], tokens[3])
    values = np.zeros((len(tokens[4:]),3))
    for t,tk in enumerate(tokens[4:]):
        scores = tk.split(":")
        values[t,0] = float(scores[0])
        values[t,1] = float(scores[1])
        values[t,2] = float(scores[2])
    return (quartet,values)

def combineFilesAsDict(files):
    """
    If an analysis with QCF was run for different groups of genes and has
    the output in multiple files, this will combine them all into a single
    dictionary for bootstrapping. It is basically a wrapper for calling the
    file2dict function, which does the actual conversion.
    """
    qcfDict = {}
    for file in files:
        with open(file) as f:
            for line in f.readlines():
                qcfKey, qcfValue = parseRawQcfLine(line)
                try:
                    qcfDict[qcfKey] = np.concatenate((qcfDict[qcfKey],qcfValue))
                except:
                    qcfDict[qcfKey] = qcfValue
    return qcfDict

def resample(combinedDict, reps):
    """
    This is the workhorse function that takes the combined dictionary,
    calculates the QCF values for the gene-tree quartet scores, then
    does bootstrap resampling to get the confidence intervals for each
    quartet.

    The output is a dictionary with the species quartet as the key, and
    the QCF values and their CIs as the values.

    (sp1,sp2,sp3,sp4): CF12_34, CF12_24lo, CF12_34hi, CF13_24, CF13_24lo, ...
    """
    resultsDict = {}
    for k,v in combinedDict.items():
        results = np.zeros((10))
        ngenes = v.shape[0]
        results[9] = ngenes
        sum = np.sum(v)
        cf12_34 = np.sum(v[:,0]) / sum
        cf13_24 = np.sum(v[:,1]) / sum
        cf14_23 = np.sum(v[:,2]) / sum
        assert (cf12_34 + cf13_24 + cf14_23 - 1.0) < 1e-12
        results[0] = cf12_34
        results[3] = cf13_24
        results[6] = cf14_23
        cf12_34vec = np.zeros((reps))
        cf13_24vec = np.zeros((reps))
        cf14_23vec = np.zeros((reps))
        for r in range(reps):
            bootIndices = np.random.randint(ngenes, size=ngenes)
            sum = np.sum(v[bootIndices,:])
            cf12_34vec[r] = np.sum(v[bootIndices,0]) / sum
            cf13_24vec[r] = np.sum(v[bootIndices,1]) / sum
            cf14_23vec[r] = np.sum(v[bootIndices,2]) / sum
        results[1] = np.percentile(cf12_34vec, 2.5)
        results[2] = np.percentile(cf12_34vec, 97.5)
        results[4] = np.percentile(cf13_24vec, 2.5)
        results[5] = np.percentile(cf13_24vec, 97.5)
        results[7] = np.percentile(cf14_23vec, 2.5)
        results[8] = np.percentile(cf14_23vec, 97.5)
        resultsDict[k] = results
    return resultsDict

def writeOutfile(res, pfx):
    """
    Take a results dictionary after bootstrapping and write to file.
    """
    with open(pfx+"-boot.CFs.csv", 'a') as f:
        print("taxon1", "taxon2", "taxon3", "taxon4",
              "CF12.34", "CF12.34_lo", "CF12.34_hi",
              "CF13.24", "CF13.24_lo", "CF13.24_hi",
              "CF14.23", "CF14.23_lo", "CF14.23_hi",
              "nqrts",sep=",", file=f)
        for k,v in res.items():
            print(k[0], k[1], k[2], k[3],
                  v[0], v[1], v[2], v[3], v[4],
                  v[5], v[6], v[7], v[8], int(v[9]),
                  sep=",", file=f)

if __name__ == "__main__":
    """
    Run the script.
    """
    if len(argv) < 2:
        print(__doc__)
        exit(0)

    parser = argparse.ArgumentParser(description="Options for qcf_boot.py",
                                     add_help=True)
    required = parser.add_argument_group("required arguments")
    required.add_argument('-i', '--infile', action="store", type=str, metavar='file',
                          nargs='+', help="Raw QCF files (1 or more)")
    additional = parser.add_argument_group("additional arguments")
    additional.add_argument('-b', '--bootReps', action="store", type=int, default=100,
                            metavar='\b', help="number of bootstrap reps [100]")
    additional.add_argument('--prefix', action="store", type=str, default="out",
                            metavar="\b", help="prefix for output file [out]")
    additional.add_argument('-q', '--quiet', action="store_true",
                            help="supress printing to stdout")

    args     = parser.parse_args()
    infile   = args.infile
    bootReps = args.bootReps
    prefix   = args.prefix
    quiet    = args.quiet

    print("[qcf_boot.py] Combining files across runs (if they exist)...")
    rawQcfDict = combineFilesAsDict(infile)
    print("[qcf_boot.py] Conducting bootstrap resampling...")
    qcfResultsDict = resample(rawQcfDict, bootReps)
    print("[qcf_boot.py] Writing results to file (" + prefix + "-boot.CF.csv)...")
    writeOutfile(qcfResultsDict, prefix)
    print("[qcf_boot.py] Done.")
