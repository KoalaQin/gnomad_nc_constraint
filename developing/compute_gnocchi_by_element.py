#!/usr/bin/env python
# coding: utf-8

import hail as hl
hl.init()


import sys
sys.path.append('/home/qinhe')

from constraint_utils.generic import *
from constraint_utils.constraint_basics import *

# gcloud compute scp /Users/siweichen/Desktop/Ben/gnomAD/lof/constraint_utils/*.py siwei@nc-l-m:/home/siwei/

def filter_black_regions(ht: Union[hl.MatrixTable, hl.Table]) -> Union[hl.MatrixTable, hl.Table]:
    bed = hl.import_bed('gs://gnomad-nc-constraint-v31-paper/misc/blacklist_gap2.bed', reference_genome='GRCh38', skip_invalid_intervals=True)
    ht = ht.filter(hl.is_defined(bed[ht.locus]), keep=False)
    return ht

def annotate_genome_element(ht: Union[hl.MatrixTable, hl.Table], 
                            bed_path: str) -> Union[hl.MatrixTable, hl.Table]:
    bed = hl.import_bed(bed_path, reference_genome='GRCh38', skip_invalid_intervals=True)
    ht = ht.filter(hl.is_defined(bed[ht.locus]))
    return ht.annotate(element_id = bed.index(ht.locus,all_matches=True).target).explode('element_id')  

def calculate_z(input_ht: hl.Table, obs: hl.expr.NumericExpression, exp: hl.expr.NumericExpression, output: str = 'z_raw') -> hl.Table:
    """
    Compute the signed raw z score using observed and expected variant counts.
    
    The raw z scores are positive when the transcript had fewer variants than expected, and are negative when transcripts had more variants than expected.
    
    The following annotation is included in the output Table in addition to the `input_ht` keys:
        - `output` - the raw z score

    :param input_ht: Input Table.
    :param obs: Observed variant count expression.
    :param exp: Expected variant count expression.
    :param output: The annotation label to use for the raw z score output, defaults to 'z_raw'.
    :return: Table with raw z scores.
    """
    ht = input_ht.select(_obs=obs, _exp=exp)
    ht = ht.annotate(_chisq=(ht._obs - ht._exp) ** 2 / ht._exp)
    return ht.select(**{output: hl.sqrt(ht._chisq) * hl.cond(ht._obs > ht._exp, -1, 1)})

public_bucket = 'gs://gnomad-nc-constraint-v31-paper'
output_bucket = 'gs://qin-gnocchi/tmp-30day'

pfx="hg38.chrom.1kb.chr22"
sfx="genome_1kb_chr22"


context_ht = hl.read_table('{0}/context_prefiltered.ht'.format(public_bucket))
genome_ht = hl.read_table('{0}/genome_prefiltered.ht'.format(public_bucket))
context_rr_ht = hl.read_table('{0}/context_prepared_mutation_rate_po_rr_1kb.tmp.ht'.format(public_bucket))
context_ht = filter_black_regions(context_ht.filter((context_ht.coverage_mean >= 30) & (context_ht.coverage_mean <= 32)))  
# genome_ht = genome_ht.semi_join(context_ht)
context_rr_ht = context_rr_ht.semi_join(context_ht)

# annotate by 1kb element
context_rr_ht = annotate_genome_element(context_rr_ht,'{0}/misc/{1}.bed'.format(public_bucket, pfx))
genome_ht = annotate_genome_element(genome_ht,'{0}/misc/{1}.bed'.format(public_bucket, pfx))    

# select & join
# context_rr_ht = context_rr_ht.select('context', 'ref', 'alt', 'methyl_level', 'element_id')
genome_ht = genome_ht.select('context', 'ref', 'alt', 'methyl_level', 'element_id' ,'freq', 'pass_filters')
genome_join = genome_ht[context_rr_ht.key]

# filter for rare & pass qc variants
af_cutoff = 0.001
context_rr_ht = context_rr_ht.filter(hl.is_missing(genome_join) | ((genome_join.freq[0].AF <= af_cutoff) & genome_join.pass_filters))
genome_ht = genome_ht.filter((genome_ht.freq[0].AF <= af_cutoff) & genome_ht.pass_filters)


# annotate mutation rates (for computing expected) observed (1 or 0)
context_rr_ht = context_rr_ht.annotate(
    fitted_po_adj = context_rr_ht.fitted_po * context_rr_ht.rr,
    is_observed = hl.is_defined(genome_ht[context_rr_ht.key]) 
)

# compute observed, expected, and possible variants per element
by_element_ht = (
    context_rr_ht
    .group_by(
        element_id = context_rr_ht.element_id
    )
    .aggregate(
        observed=hl.agg.sum(hl.int(context_rr_ht.is_observed)),
        expected = hl.agg.sum(context_rr_ht.fitted_po),
        expected_adj = hl.agg.sum(context_rr_ht.fitted_po_adj),
        possible = hl.agg.count(),
    )
)

by_element_ht.write('{0}/oe_by_{1}.ht'.format(output_bucket,sfx), overwrite=True)

# compute gnocchi (context-only and adj)
oe_ht = hl.read_table('{0}/oe_by_{1}.ht'.format(output_bucket,sfx))

z1 = calculate_z(oe_ht, obs=oe_ht.observed, exp=oe_ht.expected, output="gnocchi")
z2 = calculate_z(oe_ht, obs=oe_ht.observed, exp=oe_ht.expected_adj, output="gnocchi_adj")

oe_ht = oe_ht.annotate(
    gnocchi = z1[oe_ht.key].gnocchi,
    gnocchi_adj = z2[oe_ht.key].gnocchi_adj,
)

oe_ht.export('{0}/gnocchi_by_{1}.txt'.format(output_bucket,sfx))



