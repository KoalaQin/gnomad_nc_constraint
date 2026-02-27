#!/usr/bin/env python
# coding: utf-8

import hail as hl
hl.init()

import numpy as np
import sys
sys.path.append('/home/siwei')

from generic import *
from constraint_basics import *



met_ht = hl.read_table('gs://gnomad-nc-constraint-v31/methylation/context_prepared_methyl_CpG.ht')
met_ht = met_ht.filter(met_ht.locus.in_x_nonpar()) ### DO NOT do any other filtering


# coeff = {
#     "Sperm_raw" : 2.0018620057268968,
#     "Oocyte_raw" : -0.06265560324853964,
#     "PN_raw" : 0.2946732573212641,
#     "C2_raw" : 0.08911842730050999,
#     "C4_raw" : 0.03306105280285928,
#     "C8_raw" : 0.11128503652384455,
#     "Morula_raw" : 0.1465552252702036,
#     "ICM_raw" : 0.10121808059932674,
#     "PGC_7W_raw" : 0.9694456602679675,
#     "PGC_10W_raw" : 1.0093996893783321,
#     "PGC_11W_raw" : 0.7625750629543298,
#     "PGC_13W_raw" : 0.1681605762315797,
#     "PGC_17W_raw" : 0.4782357478205378,
#     "PGC_19W_raw" : 0.2603238455433081,
# }


# # chrX
# coeff = {
#     "Sperm_raw" : 1.2142327678822575,
#     "Oocyte_raw" : -0.10153058557280095,
#     "PN_raw" : 0.2096996493626876,
#     "C2_raw" : 0.08246182249615806,
#     "C4_raw" : 0.0803666551372927,
#     "C8_raw" : 0.040510820716452144,
#     "Morula_raw" : 0.2502153273086105,
#     "ICM_raw" : 0.0682508020601786,
#     "PGC_7W_raw" : 0.624121048631128,
#     "PGC_10W_raw" : 1.015123980247576,
#     "PGC_11W_raw" : 1.1535398138908222,
#     "PGC_13W_raw" : 0.24685249760735126,
#     "PGC_17W_raw" : 0.5309602092634527,
#     "PGC_19W_raw" : 0.37232732130631524,
# }


# # chrX (cov10x95)
# coeff = {
#     "Sperm_raw" : 1.2558393468100875,
#     "Oocyte_raw" : -0.1183727045682403,
#     "PN_raw" : 0.2175826334879446,
#     "C2_raw" : 0.10011801498417904,
#     "C4_raw" : 0.07770738171008251,
#     "C8_raw" : 0.03646455823049704,
#     "Morula_raw" : 0.26522585144187527,
#     "ICM_raw" : 0.06873108307313872,
#     "PGC_7W_raw" : 0.6383456245236468,
#     "PGC_10W_raw" : 1.0151506169174966,
#     "PGC_11W_raw" : 1.1538034133888913,
#     "PGC_13W_raw" : 0.23404050649596275,
#     "PGC_17W_raw" : 0.47873293462189914,
#     "PGC_19W_raw" : 0.4358560695003278,
# }
# chrX (exl_black_gap2_lcr_segdup_cov10x95)
# coeff = {
#     "Sperm_raw" : 1.2475789073882082,
#     "Oocyte_raw" : -0.11077940167324811,
#     "PN_raw" : 0.22166438940629254,
#     "C2_raw" : 0.10485784097878809,
#     "C4_raw" : 0.07579743825164287,
#     "C8_raw" : 0.027713251136395608,
#     "Morula_raw" : 0.2702991034740825,
#     "ICM_raw" : 0.06608753993498105,
#     "PGC_7W_raw" : 0.6410182161093382,
#     "PGC_10W_raw" : 0.9980035973409043,
#     "PGC_11W_raw" : 1.145649918124833,
#     "PGC_13W_raw" : 0.23010568906808912,
#     "PGC_17W_raw" : 0.4529948573008039,
#     "PGC_19W_raw" : 0.45913046496633564,
# }

# reqc-ed w/ exl_black_gap2_lcr_segdup_cov_10x95_med22-24
coeff = {
    "Sperm_raw" : 1.2126259621572881,
    "Oocyte_raw" : -0.11521881548437757,
    "PN_raw" : 0.21507611633339735,
    "C2_raw" : 0.10661895301648663,
    "C4_raw" : 0.06905918532650057,
    "C8_raw" : 0.03501094997618626,
    "Morula_raw" : 0.2691360920319743,
    "ICM_raw" : 0.06272120169257137,
    "PGC_7W_raw" : 0.6005496371789111,
    "PGC_10W_raw" : 0.9724133051883074,
    "PGC_11W_raw" : 1.1380968532147429,
    "PGC_13W_raw" : 0.23452365915287712,
    "PGC_17W_raw" : 0.4490618080766945,
    "PGC_19W_raw" : 0.4693400253214728,
}

weights = dict([k,np.exp(coeff[k])] for k in coeff)

met_ht = met_ht.annotate(methyl14_weighted_logit = 
                         sum([
                             met_ht.Sperm*weights["Sperm_raw"],
                             met_ht.Oocyte*weights["Oocyte_raw"],
                             met_ht.PN*weights["PN_raw"],
                             met_ht.C2*weights["C2_raw"],
                             met_ht.C4*weights["C4_raw"],
                             met_ht.C8*weights["C8_raw"],
                             met_ht.Morula*weights["Morula_raw"],
                             met_ht.ICM*weights["ICM_raw"],
                             met_ht.PGC_7W*weights["PGC_7W_raw"],
                             met_ht.PGC_10W*weights["PGC_10W_raw"],
                             met_ht.PGC_11W*weights["PGC_11W_raw"],
                             met_ht.PGC_13W*weights["PGC_13W_raw"],
                             met_ht.PGC_17W*weights["PGC_17W_raw"],
                             met_ht.PGC_19W*weights["PGC_19W_raw"]]
                         ) / sum(weights.values()) )



# met_ht.select('context', 'ref', 'alt', 'methyl14_weighted_logit').write('gs://gnomad-nc-constraint-v31/methylation/context_prepared_methyl14_weighted_logit_x.ht', overwrite=True)
# met_ht.select('context', 'ref', 'alt', 'methyl14_weighted_logit').write('gs://gnomad-nc-constraint-v31/methylation/context_prepared_methyl14_weighted_logit_exl_black_gap2_lcr_segdup_cov10x95_x.ht', overwrite=True)
met_ht.select('context', 'ref', 'alt', 'methyl14_weighted_logit').write('gs://gnomad-nc-constraint-v31/methylation/context_prepared_methyl14_weighted_logit_exl_black_gap2_lcr_segdup_cov10x95_med22-24_x.ht', overwrite=True)








