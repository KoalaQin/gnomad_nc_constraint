#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import os

import argparse

parser = argparse.ArgumentParser()
# parser.add_argument("-pos_file_pfx", help="pos_file_pfx")
parser.add_argument("-pos_file", help="pos_file")
parser.add_argument("-obs_file", help="obs_file")
parser.add_argument("-mr_file", help="mr_file")
# parser.add_argument("-chrom", help="chrom")
args = parser.parse_args()

# -pos_file_pfx $1 -obs_file $2 -mr_file $3 -chr ${SGE_TASK_ID}

pos_file = args.pos_file
obs_file = args.obs_file
mr_file = args.mr_file
# chrom = args.chrom



po = dict([tuple(line.strip().split("\t")[:4]),float(line.strip().split("\t")[-1])] for line in open(mr_file).readlines()[1:])



d_observed = dict([line.strip().split("\t")[0],
                   int(line.strip().split("\t")[-1])] for line in 
                  open(obs_file).readlines()[1:]
                  )



df_possible = pd.read_table(pos_file,
                            header=0,
                           names = ["context","ref","alt","methylation_level","element_id","variant_count"])

df_possible = df_possible.astype({"methylation_level": str})



df_possible["po"] = df_possible.set_index(
    ["context","ref","alt","methylation_level"]).index.map(po)


df_possible["expected"] = df_possible["variant_count"]*df_possible["po"]


df_oe = df_possible.groupby(["element_id"])["variant_count","expected"].apply(sum)


df_oe["observed"] = df_oe.index.map(d_observed)
df_oe["observed"] = df_oe["observed"].fillna(value=0)

df_oe["oe"] = df_oe["observed"]/df_oe["expected"]




df_oe["chisq"] = (df_oe["observed"]-df_oe["expected"])**2 / df_oe["expected"]
df_oe["z"] = np.where(df_oe['oe'] >= 1., (-1)*np.sqrt(df_oe['chisq']), np.sqrt(df_oe['chisq']))

df_oe = df_oe.drop(columns=["chisq"]).dropna()



import csv
df_oe.reset_index().to_csv(
	obs_file.replace("observed_counts","oe_z_by_context"),
    sep="\t", quoting=csv.QUOTE_NONE,
    header=True, index=False)


# In[ ]:




# In[ ]:




