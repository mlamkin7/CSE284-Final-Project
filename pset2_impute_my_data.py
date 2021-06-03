#!/usr/bin/env python

"""
CSE284 - Personal Genomics for Bioinformaticians
Problem Set 2 - Ancestry: Imputation

Example usage:
./pset2_impute.py \
  ps2_impute.combined \
  ps2_impute.heldout.gen.gz \
  1000GP_Phase3_chr16.legend.gz
  
This script outputs Pearson r2 comparing true vs. imputed genotypes
"""

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import sys

try:
    imputefile = sys.argv[1]
    truthfile = sys.argv[2]
    legendfile = sys.argv[3]
except:
    sys.stderr.write(__doc__+"\n")
    sys.exit(1)

# Load imputation results
impres = pd.read_csv(imputefile, sep="\t")
impres.columns = ["chr", "position", "ref1", "ref2", "gt"]
impres['position'] = impres['position'].astype(np.int64)

# Load truth
truth = pd.read_csv(truthfile, sep=" ", names=["chrom","rsid","position","ref","alt","truth_00", "truth_01", "truth_11"])

# Merge
impres = pd.merge(impres, truth, on=["position"])

# Load legend
legend = pd.read_csv(legendfile, sep=" ")

# merge
data = pd.merge(impres, legend, on=["position"])

# Annotate best gt and score for each panel
def GetGenotype(x00, x01, x11):
    genotypes = [0, 1, 2]
    scores = [x00, x01, x11]
    ind = scores.index(max(scores))
    return genotypes[ind]

# Get r2 for each reference panel
data["truth_gt"] = data.apply(lambda x: GetGenotype(x["truth_00"], x["truth_01"], x["truth_11"]), 1)
print(" ".join(["YRI", str(pearsonr(data["truth_gt"], data["gt"])[0]**2), str(data.shape[0])])+"\n")
    
# Get r2 for different MAF thresholds
thresh = [0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
for t in thresh:
    x = data[data["ALL"] <= t]
    print(" ".join([str(t), str(pearsonr(x["truth_gt"], x["gt"])[0]**2)])+"\n")