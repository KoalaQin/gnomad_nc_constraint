#!/usr/bin/env python
# coding: utf-8

import hail as hl
hl.init()


import sys
sys.path.append('/home/siwei')

from generic import *
from constraint_basics import *



# def remove_coverage_outliers(t: Union[hl.MatrixTable, hl.Table]) -> Union[hl.MatrixTable, hl.Table]:
#     """
#     Keep only loci where genome coverage was between 15 and 60
#     """
#     criteria = (t.coverage_mean >= 15) & (t.coverage_mean <= 60)
#     return t.filter_rows(criteria) if isinstance(t, hl.MatrixTable) else t.filter(criteria)
    
    
# context_ht_path = "gs://gnomad-nc-constraint-v31/context_prepared.ht" 
# genome_ht_path = "gs://gnomad-nc-constraint-v31/genome_prepared.ht" # already filtered on pass (prepare_genome_ht.py)

# full_context_ht = hl.read_table(context_ht_path)
# full_genome_ht = hl.read_table(genome_ht_path)


# context_x_ht = hl.filter_intervals(remove_coverage_outliers(full_context_ht), [hl.parse_locus_interval('chrX',reference_genome='GRCh38')])
# context_x_ht = context_x_ht.filter(context_x_ht.locus.in_x_nonpar())


# genome_x_ht = hl.filter_intervals(remove_coverage_outliers(full_genome_ht), [hl.parse_locus_interval('chrX',reference_genome='GRCh38')])
# genome_x_ht = genome_x_ht.filter(genome_x_ht.locus.in_x_nonpar())


# context_x_ht.write("gs://gnomad-nc-constraint-v31/context_prefiltered_chrX.ht")
# genome_x_ht.write("gs://gnomad-nc-constraint-v31/genome_prefiltered_chrX.ht")




context_ht_path = "gs://gnomad-nc-constraint-v31/context_prepared.ht" 
genome_ht_path = "gs://gnomad-nc-constraint-v31/genome_prepared.ht" # already filtered on pass (prepare_genome_ht.py)

full_context_ht = hl.read_table(context_ht_path)
full_genome_ht = hl.read_table(genome_ht_path)

# QIN: This is too slow, change to filter by chrX nonPAR intervals, It will get the 
# same results. 
context_x_ht = full_context_ht.filter(full_context_ht.locus.in_x_nonpar())
genome_x_ht = full_genome_ht.filter(full_genome_ht.locus.in_x_nonpar())


context_x_ht.write("gs://gnomad-nc-constraint-v31/context_prepared_chrX.ht")
genome_x_ht.write("gs://gnomad-nc-constraint-v31/genome_prepared_chrX.ht")

