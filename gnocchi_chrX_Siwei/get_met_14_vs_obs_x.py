import hail as hl
hl.init()


import sys
sys.path.append('/home/siwei')

from generic import *
from constraint_basics import *


context_ht = hl.read_table('gs://gnomad-nc-constraint-v31/methylation/context_prepared_methyl_CpG.ht')

grouping_variables = ('Sperm_raw','Oocyte_raw','PN_raw','C2_raw','C4_raw','C8_raw','Morula_raw','ICM_raw',
                      'PGC_7W_raw','PGC_10W_raw','PGC_11W_raw','PGC_13W_raw','PGC_17W_raw','PGC_19W_raw')
grouping_variables2 = list(grouping_variables) + ['Sperm_methyl_level','Oocyte_methyl_level','PN_methyl_level',
                                                 'C2_methyl_level','C4_methyl_level','C8_methyl_level',
                                                 'Morula_methyl_level','ICM_methyl_level','PGC_7W_methyl_level',
                                                 'PGC_10W_methyl_level','PGC_11W_methyl_level','PGC_13W_methyl_level',
                                                 'PGC_17W_methyl_level','PGC_19W_methyl_level','context']
grouping_variables2 = tuple(grouping_variables2)
# context_ht = context_ht.filter( (hl.all(lambda x: x, [hl.is_defined(context_ht[x]) for x in grouping_variables])) &
#                                 (context_ht.alleles[0] == "C")).select(*grouping_variables2)
context_ht = context_ht.filter( hl.all(lambda x: x, [hl.is_defined(context_ht[x]) for x in grouping_variables])).select(*grouping_variables2)


# chrX
context_ht_x = hl.read_table("gs://gnomad-nc-constraint-v31/context_prefiltered_chrX.ht")
###
coverage_ht = hl.read_table("gs://gcp-public-data--gnomad/release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.ht")
context_ht_x = context_ht_x.filter(
    (coverage_ht[context_ht_x.locus].over_10>=0.95) & 
    (coverage_ht[context_ht_x.locus].median_approx>=22) & 
    (coverage_ht[context_ht_x.locus].median_approx<=24) )


def filter_bad_regions(ht: Union[hl.MatrixTable, hl.Table], 
                     bed_path: str, bed_reference_genome: str) -> Union[hl.MatrixTable, hl.Table]:
    bed = hl.import_bed(bed_path, bed_reference_genome, skip_invalid_intervals=True)
    ht = ht.filter(hl.is_defined(bed[ht.locus]), keep=False)
    return ht

# path = "gs://gnomad-nc-constraint-v31/mutation_rate/probalematic_regions_ENCFF356LFX_repeatMasker.bed"
# path = "gs://gnomad-nc-constraint-v31/mutation_rate/blacklist_gap2.bed"
path = "gs://gnomad-nc-constraint-v31/mutation_rate/blacklist_gap2_lcr_segdup.merged.bed"
# path = "gs://gnomad-nc-constraint-v31/mutation_rate/cpgIslandExt.bed"
reference_genome = 'GRCh38'
context_ht_x = filter_bad_regions(context_ht_x,path,reference_genome)  
###



context_ht = context_ht.semi_join(context_ht_x)


genome_ht = hl.read_table("gs://gnomad-nc-constraint-v31/genome_prefiltered_chrX.ht")
genome_join = genome_ht[context_ht.key]

af_cutoff = 0.001
# af_cutoff = 0.00001
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


context_ht = context_ht.annotate(obs = hl.is_defined(genome_ht[context_ht.key]))


# context_ht.export('gs://gnomad-nc-constraint-v31/methylation/context_methyl_CpG_14_obs01_passaf01_x.txt')
# context_ht.export('gs://gnomad-nc-constraint-v31/methylation/context_methyl_CpG_14_obs01_passaf01_cov10x95_x.txt')
# context_ht.export('gs://gnomad-nc-constraint-v31/methylation/context_methyl_CpG_14_obs01_passaf01_exl_black_gap2_lcr_segdup_cov10x95_x.txt')
context_ht.export('gs://gnomad-nc-constraint-v31/methylation/context_methyl_CpG_14_obs01_af01_exl_black_gap2_lcr_segdup_cov10x95_med22-24_x.txt')







