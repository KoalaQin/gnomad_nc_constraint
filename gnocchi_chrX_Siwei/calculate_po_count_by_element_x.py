#!/usr/bin/env python
# coding: utf-8

import hail as hl
hl.init()


import sys
sys.path.append('/home/siwei')

from generic import *
from constraint_basics import *



def filter_for_element(ht: Union[hl.MatrixTable, hl.Table], 
                     bed_path: str, bed_reference_genome: str) -> Union[hl.MatrixTable, hl.Table]:

    # bed = hl.import_bed(bed_path, bed_reference_genome)
    bed = hl.import_bed(bed_path, bed_reference_genome, skip_invalid_intervals=True)
    
    ht = ht.filter(hl.is_defined(bed[ht.locus]))
    # return ht.annotate(element_id = bed[ht.locus].target)
    return ht.annotate(element_id = bed.index(ht.locus,all_matches=True).target).explode('element_id')


 


context_ht = hl.read_table("gs://gnomad-nc-constraint-v31/context_prefiltered_chrX.ht")
genome_ht = hl.read_table("gs://gnomad-nc-constraint-v31/genome_prefiltered_chrX.ht")

###
# coverage_ht = hl.read_table("gs://gcp-public-data--gnomad/release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.ht")
# context_ht = context_ht.filter(
#     (coverage_ht[context_ht.locus].over_10>=0.95) & 
#     (coverage_ht[context_ht.locus].median_approx>=22) & 
#     (coverage_ht[context_ht.locus].median_approx<=24) )


def filter_bad_regions(ht: Union[hl.MatrixTable, hl.Table], 
                     bed_path: str, bed_reference_genome: str) -> Union[hl.MatrixTable, hl.Table]:
    bed = hl.import_bed(bed_path, bed_reference_genome, skip_invalid_intervals=True)
    ht = ht.filter(hl.is_defined(bed[ht.locus]), keep=False)
    return ht

# path = "gs://gnomad-nc-constraint-v31/mutation_rate/probalematic_regions_ENCFF356LFX_repeatMasker.bed"
path = "gs://gnomad-nc-constraint-v31/mutation_rate/blacklist_gap2.bed"
# path = "gs://gnomad-nc-constraint-v31/mutation_rate/blacklist_gap2_lcr_segdup.merged.bed"
# path = "gs://gnomad-nc-constraint-v31/mutation_rate/cpgIslandExt.bed"
reference_genome = 'GRCh38'
context_ht = filter_bad_regions(context_ht,path,reference_genome)  

###
# experimental - restrict coverage
context_ht = context_ht.filter((context_ht.coverage_mean >= 23) & (context_ht.coverage_mean <= 24))

##########
genome_ht = genome_ht.semi_join(context_ht)
##########



# filter for element

path = 'gs://gnomad-nc-constraint-v31/genome_element/hg38.chromX.1kb.bed'
reference_genome = 'GRCh38'

# filter for element
context_ht = filter_for_element(context_ht,path,reference_genome)
genome_ht = filter_for_element(genome_ht,path,reference_genome)    


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





# select & join
context_ht = context_ht.select('context', 'ref', 'alt', 'methyl_level','element_id','variant_type')

# genome_ht = genome_ht.select('context', 'ref', 'alt', 'methyl_level', 'element_id','variant_type',
#                              'freq', 'pass_filters')
genome_ht = genome_ht.select('element_id','variant_type',
                             'freq', 'pass_filters')


# filter pass_filters

genome_join = genome_ht[context_ht.key]

# context_ht = context_ht.filter(
#     hl.is_missing(genome_join) | (genome_join.pass_filters)
# )

# genome_ht = genome_ht.filter(genome_ht.pass_filters)


# filter no_common (remove_common_ordinary)
af_cutoff = 0.001
# ac_cutoff = 1
context_ht = context_ht.filter(
    hl.is_missing(genome_join) | (
        (genome_join.freq[0].AF <= af_cutoff) & genome_join.pass_filters
        # (genome_join.freq[0].AC <= ac_cutoff) & genome_join.pass_filters
    )
)

genome_ht = genome_ht.filter(
    (genome_ht.freq[0].AF <= af_cutoff) & genome_ht.pass_filters
    # (genome_ht.freq[0].AC <= ac_cutoff) & genome_ht.pass_filters
)




# grouping criteria
grouping = hl.struct(context=context_ht.context, ref=context_ht.ref, alt=context_ht.alt, 
                     methylation_level=context_ht.methyl_level, element_id = context_ht.element_id,
                     # variant_type = context_ht.variant_type, ### better to include if for random forest
                     )


# count possible (context_ht)
output = {'variant_count': hl.agg.count()}
possible_ht = context_ht.group_by(**grouping).aggregate(**output)

# possible_ht.export('gs://gnomad-nc-constraint-v31/genome_element/possible_counts_by_context_logit13-x.genome_1kb_exl_black_gap2_lcr_segdup_cov10x95_chrX.txt') ### current model
# possible_ht.export('gs://gnomad-nc-constraint-v31/genome_element/possible_counts_by_context_methyl_logit12.genome_1kb_exl_black_gap2_lcr_segdup_cov10x95_med22-24_chrX.txt') # reqc
possible_ht.export('gs://gnomad-nc-constraint-v31/genome_element/possible_counts_by_context_methyl_logit12.genome_1kb_exl_black_gap2_cov23-24_chrX.txt') # reqc





# compute observed (genome_ht)
# grouping = hl.struct(context=genome_ht.context, ref=genome_ht.ref, alt=genome_ht.alt, 
#                      methylation_level=genome_ht.methyl_level, element_id = genome_ht.element_id,
#                      variant_type = genome_ht.variant_type,
#                      ) 
grouping = hl.struct(element_id = genome_ht.element_id,
                     # variant_type = genome_ht.variant_type, ### needed if for random forest
                     )
output = {'variant_count': hl.agg.count()}
observed_ht = genome_ht.group_by(**grouping).aggregate(**output)

observed_ht.export('gs://gnomad-nc-constraint-v31/genome_element/observed_counts.genome_1kb_exl_black_gap2_cov23-24_chrX.txt')





