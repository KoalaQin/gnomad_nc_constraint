#!/usr/bin/env python
# coding: utf-8

import hail as hl
hl.init()


# import sys
# sys.path.append('/home/siwei')

# from generic import *
# from constraint_basics import *


context_ht = hl.read_table('gs://gnomad-nc-constraint-v31/context_prefiltered_chrX.ht')
genome_ht = hl.read_table('gs://gnomad-nc-constraint-v31/genome_prefiltered_chrX.ht')



grouping_variables = (("cpg","variant_type","variant_type_model",
                       "methylation_ENCFF918PML","methylation_GSM429321",
                       "methylation_level_ENCFF918PML","methylation_level_GSM429321"))
context_ht = context_ht.select('context', 'ref', 'alt', *grouping_variables)
genome_ht = genome_ht.select('context', 'ref', 'alt','freq', 'pass_filters', *grouping_variables)

genome_join = genome_ht[context_ht.key]




ac_cutoff = 5
downsampling_level = 1000
freq_index = hl.eval(
    hl.zip_with_index(genome_ht.freq_meta).find(lambda f:
                                                (f[1].get('downsampling') == str(downsampling_level)) &
                                                (f[1].get('pop') == 'global') &
                                                (f[1].get('group') == 'adj') &
                                                (f[1].size() == 3)))[0]
# In the denominator, only keep variants not in the genome dataset, or with AC <= ac_cutoff and passing filters
context_ht = context_ht.filter(
    hl.is_missing(genome_join) | ((genome_join.freq[freq_index].AC <= ac_cutoff) & genome_join.pass_filters)
)
# Keep only AC <= ac_cutoff in numerator
genome_ht = genome_ht.filter((genome_ht.freq[freq_index].AC <= ac_cutoff) & genome_ht.pass_filters)


genome_ht.write("gs://gnomad-nc-constraint-v31/genome_downsampled_1000_chrX.ht", overwrite = True)
context_ht.write("gs://gnomad-nc-constraint-v31/context_downsampled_1000_chrX.ht", overwrite = True)



