# This code is adapted from Siwei's code `compute_gnocchi_by_element.py`
from typing import Optional, Sequence, Union, Dict, List

import hail as hl
import argparse
from gnomad.utils.filtering import filter_by_intervals, filter_low_conf_regions
from constraint_utils.nc_constraint_utils import remove_coverage_outliers

# ----------------------------
# Config
# ----------------------------
PUBLIC_BUCKET = "gs://gnomad-nc-constraint-v31-paper"

REFGENOME_NAME = "GRCh38"
REFGENOME = hl.get_reference(REFGENOME_NAME)

BED_1KB_PATH = f"{PUBLIC_BUCKET}/misc/hg38.chrom.1kb.bed"
BLACKLIST_BED_PATH = f"{PUBLIC_BUCKET}/misc/blacklist_gap2.bed"
# This blacklist covers 20kb of chrX PAR and 8.6Mb of the chrX non-PAR region.

REGION_INTERVALS = {
    "chrX_PAR": [
        "[chrX:10000-2781479]",
        "[chrX:155701383-156030895]",
    ],
    "chrX_nonPAR": [
        "[chrX:1-9999]",
        "[chrX:2781480-155701382]",
        "[chrX:156030896-156040895]",
    ],
}


def parse_intervals(interval_strs: Sequence[str]) -> Sequence[hl.Interval]:
    return [
        hl.parse_locus_interval(x, reference_genome=REFGENOME) for x in interval_strs
    ]


def filter_black_regions(ht: hl.Table, blacklist_bed: hl.Table) -> hl.Table:
    return ht.filter(hl.is_defined(blacklist_bed[ht.locus]), keep=False)


def annotate_genome_element(
    ht: hl.Table,
    element_bed: hl.Table,
    element_field: str = "element_id",
) -> hl.Table:
    """
    Filter to loci present in element_bed and annotate element IDs.
    If multiple overlapping intervals exist, explode to one row per match.
    """
    ht = ht.filter(hl.is_defined(element_bed[ht.locus]))
    # all_matches=True returns an array of structs (interval records)
    # "target" is BED's 4th column (name) if present.
    ht = ht.annotate(
        **{element_field: element_bed.index(ht.locus, all_matches=True).target}
    )
    return ht.explode(element_field)


def calculate_z(
    input_ht: hl.Table,
    obs: hl.expr.NumericExpression,
    exp: hl.expr.NumericExpression,
    output: str = "z_raw",
) -> hl.Table:
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
    return ht.select(
        **{output: hl.sqrt(ht._chisq) * hl.if_else(ht._obs > ht._exp, -1, 1)}
    )


def prefilter_input_ht(
    ht: hl.Table,
    blacklist_ht: Optional[hl.Table] = None,
    intervals: Optional[Sequence[hl.Interval]] = None,
) -> hl.Table:
    """
    1) interval filter (optional)
    2) blacklist filter (optional)
    3) remove coverage outliers if coverage_mean exists
    4) lcr and segdup regions
    """
    if intervals is not None:
        ht = filter_by_intervals(ht, intervals, reference_genome=REFGENOME_NAME)

    if blacklist_ht is not None:
        ht = filter_black_regions(ht, blacklist_ht)

    if "coverage_mean" in ht.row.dtype.fields:
        ht = remove_coverage_outliers(ht)

    ht = filter_low_conf_regions(
        ht, filter_lcr=True, filter_decoy=False, filter_segdup=True
    )

    return ht


def main(args):
    intervals = resolve_intervals(args)

    if args.dry_run:
        print("=== DRY RUN ===")
        print(f"Region:        {args.region}")
        print(f"Intervals:     {[str(i) for i in hl.eval(intervals)]}")
        print(f"AF cutoff:     {args.af_cutoff}")
        print(f"Output bucket: {args.output_bucket}")
        print(f"Output suffix: {args.output_suffix}")
        print("==============")
        return

    hl.default_reference(REFGENOME_NAME)

    AF_CUTOFF = args.af_cutoff
    OUTPUT_BUCKET = args.output_bucket
    SFX = args.output_suffix

    # ----------------------------
    # Load reference BEDs
    # ----------------------------
    blacklist_bed = hl.import_bed(
        BLACKLIST_BED_PATH,
        reference_genome=REFGENOME_NAME,
        skip_invalid_intervals=True,
    )

    element_bed = hl.import_bed(
        BED_1KB_PATH,
        reference_genome=REFGENOME_NAME,
        skip_invalid_intervals=True,
    )

    # ----------------------------
    # Read base tables
    # ----------------------------
    if args.region == "chrX_nonPAR":
        context_rr_ht = hl.read_table(
            f"{OUTPUT_BUCKET}/context_prepared_chrX_nonPAR_methyl_po_1kb.ht"
        )
        context_ht = hl.read_table(
            f"{OUTPUT_BUCKET}/context_prepared_chrX_nonPAR.ht"
        )
        genome_ht = hl.read_table(f"{OUTPUT_BUCKET}/genome_prepared_chrX_nonPAR.ht")
    else:
        context_rr_ht = hl.read_table(
            f"{PUBLIC_BUCKET}/context_prepared_mutation_rate_po_rr_1kb.tmp.ht"
        )
        context_ht = hl.read_table(f"{PUBLIC_BUCKET}/context_prepared.ht")
        genome_ht = hl.read_table(f"{PUBLIC_BUCKET}/genome_prepared.ht")
    
    # ----------------------------
    # Prefilter
    # ----------------------------

    context_rr_ht = prefilter_input_ht(
        context_rr_ht, intervals=intervals, blacklist_ht=blacklist_bed
    )
    genome_ht = prefilter_input_ht(
        genome_ht, intervals=intervals, blacklist_ht=blacklist_bed
    )
    context_ht = prefilter_input_ht(
        context_ht, intervals=intervals, blacklist_ht=blacklist_bed
    )

    if args.region == "chrX_nonPAR":
        # Based on Siwei's notes
        coverage_mean_lower = 23
        coverage_mean_upper = 24
    else:
        coverage_mean_lower = 30
        coverage_mean_upper = 32

    context_ht = context_ht.filter(
        (context_ht.coverage_mean >= coverage_mean_lower)
        & (context_ht.coverage_mean <= coverage_mean_upper)
    )
    context_rr_ht = context_rr_ht.semi_join(context_ht)

    # ----------------------------
    # Checkpoint
    # ----------------------------
    context_rr_ht = context_rr_ht.checkpoint(
        f"{OUTPUT_BUCKET}/context_rr_{SFX}.ht", overwrite=True
    )
    genome_ht = genome_ht.checkpoint(f"{OUTPUT_BUCKET}/genome_{SFX}.ht", overwrite=True)

    # ----------------------------
    # Annotate by 1kb elements
    # ----------------------------
    context_rr_ht = annotate_genome_element(context_rr_ht, element_bed)
    genome_ht = annotate_genome_element(genome_ht, element_bed)

    # ----------------------------
    # Select & join
    # ----------------------------
    genome_ht = genome_ht.select(
        "context", "ref", "alt", "methyl_level", "element_id", "freq", "pass_filters"
    )
    genome_join = genome_ht[context_rr_ht.key]

    # ----------------------------
    # filter for rare & pass qc variants
    # ----------------------------
    context_rr_ht = context_rr_ht.filter(
        hl.is_missing(genome_join)
        | ((genome_join.freq[0].AF <= AF_CUTOFF) & genome_join.pass_filters)
    )
    genome_ht = genome_ht.filter(
        (genome_ht.freq[0].AF <= AF_CUTOFF) & genome_ht.pass_filters
    )

    # ----------------------------
    # annotate mutation rates (for computing expected) observed (1 or 0)
    # ----------------------------
    context_rr_ht = context_rr_ht.annotate(
        fitted_po_adj=context_rr_ht.fitted_po * context_rr_ht.rr,
        is_observed=hl.is_defined(genome_ht[context_rr_ht.key]),
    )

    # ----------------------------
    # Aggregate
    # ----------------------------
    by_element_ht = context_rr_ht.group_by(
        element_id=context_rr_ht.element_id
    ).aggregate(
        observed=hl.agg.sum(hl.int32(context_rr_ht.is_observed)),
        expected=hl.agg.sum(context_rr_ht.fitted_po),
        expected_adj=hl.agg.sum(context_rr_ht.fitted_po_adj),
        possible=hl.agg.count(),
    )

    out_ht = f"{OUTPUT_BUCKET}/oe_by_{SFX}.ht"
    by_element_ht.write(out_ht, overwrite=True)

    # ----------------------------
    # compute Z scores
    # ----------------------------
    oe_ht = hl.read_table(out_ht)

    z = calculate_z(oe_ht, oe_ht.observed, oe_ht.expected, "gnocchi")
    z_adj = calculate_z(oe_ht, oe_ht.observed, oe_ht.expected_adj, "gnocchi_adj")

    oe_ht = oe_ht.annotate(
        gnocchi=z[oe_ht.key].gnocchi,
        gnocchi_adj=z_adj[oe_ht.key].gnocchi_adj,
    )

    # ----------------------------
    # QC elements by pass_filters and coverage
    # ----------------------------
    pass_ht = (
        hl.import_table(
            f"{PUBLIC_BUCKET}/misc/genome_1kb_gnomad_v31_pass.txt",
            no_header=True,
            delimiter="\t",
            impute=True,
        )
        .rename({"f0": "element_id", "f1": "pct_PASS"})
        .key_by("element_id")
    )

    cov_ht = (
        hl.import_table(
            f"{PUBLIC_BUCKET}/misc/genome_1kb_gnomad_v31_coverage.txt",
            no_header=True,
            delimiter="\t",
            impute=True,
        )
        .rename({"f0": "element_id", "f1": "mean_coverage"})
        .key_by("element_id")
    )

    oe_ht = oe_ht.annotate(
        pct_pass=pass_ht[oe_ht.key].pct_PASS,
        mean_coverage=cov_ht[oe_ht.key].mean_coverage,
    )

    if args.region == "chrX_nonPAR":
        cov_lower = 20
        cov_upper = 25
    else:
        cov_lower = 25
        cov_upper = 35

    oe_ht = oe_ht.annotate(
        pass_qc=hl.if_else(
            (oe_ht.pct_pass >= 0.8)
            & (oe_ht.mean_coverage >= cov_lower)
            & (oe_ht.mean_coverage <= cov_upper)
            & (oe_ht.possible >= 1000),
            True,
            False,
        )
    )

    oe_ht.export(f"{OUTPUT_BUCKET}/gnocchi_by_{SFX}.txt")


def resolve_intervals(args) -> Sequence[hl.Interval]:
    if args.interval:
        interval_strs = args.interval
    elif args.region:
        interval_strs = REGION_INTERVALS[args.region]
    else:
        raise ValueError("Must provide either --region or --interval")

    return parse_intervals(interval_strs)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run gnocchi for chrX PAR / nonPAR or custom intervals"
    )

    parser.add_argument(
        "--region",
        choices=REGION_INTERVALS.keys(),
        help="Predefined region to run (chrX_PAR or chrX_nonPAR)",
    )

    parser.add_argument(
        "--interval",
        action="append",
        help="Custom interval(s), e.g. '[chrX:10000-2781479]'. Can be passed multiple times.",
    )

    parser.add_argument(
        "--af-cutoff",
        type=float,
        default=0.001,
        help="AF cutoff for observed variants (default: 0.001)",
    )

    parser.add_argument(
        "--output-bucket",
        required=True,
        help="GCS bucket/path for outputs (e.g. gs://qin-gnocchi/tmp-30day)",
    )

    parser.add_argument(
        "--output-suffix",
        required=True,
        help="Suffix for output files (e.g. genome_1kb_chrX_PAR)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Resolve intervals and print config, but do not run Hail",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(args)
