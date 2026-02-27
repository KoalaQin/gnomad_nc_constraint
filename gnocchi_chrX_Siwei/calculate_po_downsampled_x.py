#!/usr/bin/env python
# coding: utf-8

import hail as hl
hl.init()


import sys
sys.path.append('/home/siwei')

from generic import *
from constraint_basics import *




context_ht = hl.read_table('gs://gnomad-nc-constraint-v31/context_downsampled_1000_chrX.ht')
genome_ht = hl.read_table('gs://gnomad-nc-constraint-v31/genome_downsampled_1000_chrX.ht')

# # ##
# coverage_ht = hl.read_table("gs://gcp-public-data--gnomad/release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.ht")
# context_ht = context_ht.filter(
#     # (coverage_ht[context_ht.locus].over_10>=0.95) &
#     (coverage_ht[context_ht.locus].median_approx>=22) 
#     & (coverage_ht[context_ht.locus].median_approx<=24) 
#     )
# # # context_ht = context_ht.filter(
# # #     (coverage_ht[context_ht.locus].median_approx>=15) & 
# # #     (coverage_ht[context_ht.locus].median_approx<=25) )
# # ##

# def filter_bad_regions(ht: Union[hl.MatrixTable, hl.Table], 
#                      bed_path: str, bed_reference_genome: str) -> Union[hl.MatrixTable, hl.Table]:
#     bed = hl.import_bed(bed_path, bed_reference_genome, skip_invalid_intervals=True)
#     ht = ht.filter(hl.is_defined(bed[ht.locus]), keep=False)
#     return ht

# # path = "gs://gnomad-nc-constraint-v31/mutation_rate/probalematic_regions_ENCFF356LFX_repeatMasker.bed"
# # path = "gs://gnomad-nc-constraint-v31/mutation_rate/blacklist_gap2.bed"
# path = "gs://gnomad-nc-constraint-v31/mutation_rate/blacklist_gap2_lcr_segdup.merged.bed"
# # path = "gs://gnomad-nc-constraint-v31/mutation_rate/cpgIslandExt.bed"
# reference_genome = 'GRCh38'
# context_ht = filter_bad_regions(context_ht,path,reference_genome)  
# ###

# # # filter_for_mu
# # ht = hl.read_table('gs://gnomad-nc-constraint-v31/context_prefiltered_chrX.ht')
# # context_ht =context_ht.filter( (ht[context_ht.key].vep.most_severe_consequence == 'intron_variant') | (ht[context_ht.key].vep.most_severe_consequence == 'intergenic_variant') )

# ##########
# genome_ht = genome_ht.semi_join(context_ht)
# ##########




# met_ht = hl.read_table('gs://gnomad-nc-constraint-v31/methylation/context_prepared_methyl14_weighted_logit.ht')
# annotation = {
#     'methyl14_weighted_logit_level': hl.case().when(
#         met_ht.methyl14_weighted_logit > 0.9, 15
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.85, 14
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.8, 13
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.75, 12
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.7, 11
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.65, 10
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.6, 9
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.55, 8
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.5, 7
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.3, 6
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.25, 5
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.2, 4
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.15, 3
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.1, 2
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.05, 1
#     ).default(0)
#     }

# met_ht = hl.read_table('gs://gnomad-nc-constraint-v31/methylation/context_prepared_methyl14_weighted_logit_x.ht')
# annotation = {
#     'methyl14_weighted_logit_level': hl.case().when(
#         met_ht.methyl14_weighted_logit > 0.85, 14
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.8, 13
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.75, 12
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.7, 11
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.65, 10
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.6, 9
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.55, 8
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.5, 7
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.45, 6
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.4, 5
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.2, 4
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.15, 3
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.1, 2
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.05, 1
#     ).default(0)
#     }

# met_ht = hl.read_table('gs://gnomad-nc-constraint-v31/methylation/context_prepared_methyl14_weighted_logit_exl_black_gap2_lcr_segdup_cov10x95_x.ht')
# # logit20-x_exl_black_gap2_lcr_segdup_cov10x95_x
# annotation = {
#     'methyl14_weighted_logit_level': hl.case().when(
#         met_ht.methyl14_weighted_logit > 0.85, 13
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.8, 12
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.75, 11
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.7, 10
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.65, 9
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.6, 8
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.55, 7
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.5, 6
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.45, 5
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.2, 4
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.15, 3
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.1, 2
#     ).when(
#         met_ht.methyl14_weighted_logit > 0.05, 1
#     ).default(0)
#     }



# logit20_exl_black_gap2_lcr_segdup_cov10x95_med22-24_x
# met_ht = hl.read_table('gs://gnomad-nc-constraint-v31/methylation/context_prepared_methyl14_weighted_logit_exl_black_gap2_lcr_segdup_cov10x95_med22-24_x.ht')
met_ht = hl.read_table('gs://gnomad-nc-constraint-v31/methylation/context_prepared_methyl14_weighted_logit_x.ht')

annotation = {
    'methyl14_weighted_logit_level': hl.case().when(
        met_ht.methyl14_weighted_logit > 0.8, 12
    ).when(
        met_ht.methyl14_weighted_logit > 0.75, 11
    ).when(
        met_ht.methyl14_weighted_logit > 0.7, 10
    ).when(
        met_ht.methyl14_weighted_logit > 0.65, 9
    ).when(
        met_ht.methyl14_weighted_logit > 0.6, 8
    ).when(
        met_ht.methyl14_weighted_logit > 0.55, 7
    ).when(
        met_ht.methyl14_weighted_logit > 0.5, 6
    ).when(
        met_ht.methyl14_weighted_logit > 0.45, 5
    ).when(
        met_ht.methyl14_weighted_logit > 0.2, 4
    ).when(
        met_ht.methyl14_weighted_logit > 0.15, 3
    ).when(
        met_ht.methyl14_weighted_logit > 0.1, 2
    ).when(
        met_ht.methyl14_weighted_logit > 0.05, 1
    ).default(0)
    }

met_ht = met_ht.annotate(**annotation)

context_ht = context_ht.annotate(methyl_level=hl.if_else(hl.is_missing(met_ht[context_ht.key]), 0, met_ht[context_ht.key].methyl14_weighted_logit_level))
genome_ht = genome_ht.annotate(methyl_level=hl.if_else(hl.is_missing(met_ht[genome_ht.key]), 0, met_ht[genome_ht.key].methyl14_weighted_logit_level))




observed_ht = count_variants(genome_ht, count_downsamplings=['global', 'nfe', 'afr'],
                             additional_grouping=('methyl_level',), omit_methylation=True)


# observed_ht.write("gs://gnomad-nc-constraint-v31/mutation_rate/observed_counts_by_context_methyl_logit15_exl_black_gap2_downsampled_pop_x.ht",True)
# observed_ht.write("gs://gnomad-nc-constraint-v31/mutation_rate/observed_counts_by_context_methyl_logit14-x_exl_black_gap2_downsampled_pop_x.ht",True)
# observed_ht.write("gs://gnomad-nc-constraint-v31/mutation_rate/observed_counts_by_context_methyl_logit13-x_exl_black_gap2_lcr_segdup_cov10x95_downsampled_pop_x.ht",True)
# observed_ht.write("gs://gnomad-nc-constraint-v31/mutation_rate/observed_counts_by_context_methyl_logit12_exl_black_gap2_lcr_segdup_cov10x95_med22-24_downsampled_pop_x.ht",True)
# observed_ht.write("gs://gnomad-nc-constraint-v31/mutation_rate/observed_counts_by_context_methyl_logit12_exl_black_gap2_lcr_segdup_downsampled_pop_x.ht",True)
# observed_ht.write("gs://gnomad-nc-constraint-v31/mutation_rate/observed_counts_by_context_methyl_logit12_csq_exl_black_gap2_lcr_segdup_downsampled_pop_x.ht",True)
# observed_ht.write("gs://gnomad-nc-constraint-v31/mutation_rate/observed_counts_by_context_methyl_logit12_exl_black_gap2_lcr_segdup_cov22-24_downsampled_pop_x.ht",True)
observed_ht.write("gs://gnomad-nc-constraint-v31/mutation_rate/observed_counts_by_context_methyl_logit12_downsampled_pop_x.ht",True)


grouping = hl.struct(context=context_ht.context, ref=context_ht.ref, alt=context_ht.alt, 
                     methylation_level=context_ht.methyl_level)

output = {'variant_count': hl.agg.count()}
possible_ht = context_ht.group_by(**grouping).aggregate(**output)

possible_ht.export('gs://gnomad-nc-constraint-v31/mutation_rate/possible_counts_by_context_methyl_logit12_downsampled_x.txt')











