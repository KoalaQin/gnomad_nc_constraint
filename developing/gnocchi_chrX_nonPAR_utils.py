
"""
Shared utilities for the chrX nonPAR gnocchi pipeline.

Pipeline steps (run via run_pipeline):
  1  --prefilter        : filter context/genome HTs to chrX nonPAR intervals
  2  --downsampling     : filter to 1000-genome AC cutoff for mu estimation
  3a --annotate-methyl  : annotate CpG loci with raw per-stage methylation values
  3b --annotate-obs     : apply config filters; annotate obs (0/1) for logit training
  3c --compute-coeff    : fit logistic regression → weighted logit level HT
  4  --compute-mu       : count obs/possible, estimate mu, fit po~mu → fitted_po table
  5 --compute-element-z : annotate fitted_po/element_id per locus; compute per-element oe and Z scores (Hail)
"""

from typing import List, Optional, Sequence, Tuple
from types import SimpleNamespace
import logging
import time

import hail as hl
import argparse
from gnomad.utils.filtering import filter_by_intervals
from constraint_utils.nc_constraint_utils import remove_coverage_outliers, get_downsamplings
from constraint_utils.generic import count_variants, annotate_variant_types

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────
# Bucket / path constants
# ──────────────────────────────────────────────────────────────────
PUBLIC_BUCKET = "gs://gnomad-nc-constraint-v31-paper"
OUTPUT_BUCKET = "gs://qin-gnocchi/gnocchi_files"

# External gnomAD genome coverage HT — provides the over_10 field used by Config B
COVERAGE_HT_PATH = (
    "gs://gcp-public-data--gnomad/release/3.0.1/coverage/genomes/"
    "gnomad.genomes.r3.0.1.coverage.ht"
)

REFGENOME_NAME = "GRCh38"

BED_1KB_PATH = f"{PUBLIC_BUCKET}/misc/hg38.chrom.1kb.bed"

# Region-exclusion BEDs
BLACKLIST_BED_PATH        = f"{PUBLIC_BUCKET}/misc/blacklist_gap2.bed"
BLACKLIST_LCR_SEGDUP_PATH = f"{OUTPUT_BUCKET}/blacklist_gap2_lcr_segdup.merged.bed"

# Per-element QC tables (used for pass_qc annotation in the element Z step)
ELEMENT_PASS_TXT  = f"{PUBLIC_BUCKET}/misc/genome_1kb_gnomad_v31_pass.txt"
ELEMENT_COV_TXT   = f"{PUBLIC_BUCKET}/misc/genome_1kb_gnomad_v31_coverage.txt"

# Coverage range for pass_qc on chrX nonPAR (matches gnocchi_autosome_par.py)
CHRX_NONPAR_COV_QC_RANGE = (20.0, 25.0)

# chrX nonPAR intervals (excludes PARs on both ends)
CHRX_NONPAR_INTERVALS = [
    "[chrX:1-9999]",
    "[chrX:2781480-155701382]",
    "[chrX:156030896-156040895]",
]

# 14 germline methylation stages: display name → field in meth_ht
STAGES = {
    "Sperm":   "methyl_level",
    "Oocyte":  "methyl_level_1",
    "PN":      "methyl_level_2",
    "C2":      "methyl_level_3",
    "C4":      "methyl_level_4",
    "C8":      "methyl_level_5",
    "Morula":  "methyl_level_6",
    "ICM":     "methyl_level_7",
    "PGC_7W":  "methyl_level_8",
    "PGC_10W": "methyl_level_9",
    "PGC_11W": "methyl_level_10",
    "PGC_13W": "methyl_level_11",
    "PGC_17W": "methyl_level_12",
    "PGC_19W": "methyl_level_13",
}

# ──────────────────────────────────────────────────────────────────
# Methylation discretization thresholds  (descending order)
# Each entry: (score_threshold, level_assigned_if_score_exceeds_threshold)
# ──────────────────────────────────────────────────────────────────

# Config A — logit12: top bin >0.80, levels 0–12
METHYL_THRESHOLDS_LOGIT12: List[Tuple[float, int]] = [
    (0.80, 12), (0.75, 11), (0.70, 10), (0.65, 9), (0.60, 8),
    (0.55, 7),  (0.50, 6),  (0.45, 5),  (0.20, 4), (0.15, 3),
    (0.10, 2),  (0.05, 1),
]

# Config B — logit13: top bin >0.85, levels 0–13
METHYL_THRESHOLDS_LOGIT13: List[Tuple[float, int]] = [
    (0.85, 13),
    (0.80, 12), (0.75, 11), (0.70, 10), (0.65, 9), (0.60, 8),
    (0.55, 7),  (0.50, 6),  (0.45, 5),  (0.20, 4), (0.15, 3),
    (0.10, 2),  (0.05, 1),
]

# logit14: top bin >0.85, levels 0–14 (adds >0.40→5 between logit13's 0.45 and 0.20)
METHYL_THRESHOLDS_LOGIT14: List[Tuple[float, int]] = [
    (0.85, 14),
    (0.80, 13), (0.75, 12), (0.70, 11), (0.65, 10), (0.60, 9),
    (0.55, 8),  (0.50, 7),  (0.45, 6),  (0.40, 5),  (0.20, 4),
    (0.15, 3),  (0.10, 2),  (0.05, 1),
]


# ──────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────

def parse_intervals(interval_strs: Sequence[str]) -> Sequence[hl.Interval]:
    return [hl.parse_locus_interval(x, reference_genome=REFGENOME_NAME) for x in interval_strs]


def filter_by_bed(ht: hl.Table, bed_path: str) -> hl.Table:
    """Remove loci overlapping a BED file."""
    bed = hl.import_bed(bed_path, reference_genome=REFGENOME_NAME, skip_invalid_intervals=True)
    return ht.filter(hl.is_defined(bed[ht.locus]), keep=False)


def annotate_genome_element(
        ht: hl.Table,
        element_bed: hl.Table,
        element_field: str = "element_id",
) -> hl.Table:
    """Filter to loci covered by element_bed and annotate element_id.
    Explodes rows when a locus overlaps multiple intervals."""
    ht = ht.filter(hl.is_defined(element_bed[ht.locus]))
    ht = ht.annotate(**{element_field: element_bed.index(ht.locus, all_matches=True).target})
    return ht.explode(element_field)


def signed_chisq_z(obs, exp):
    chisq = hl.if_else(exp > 0, (obs - exp) ** 2 / exp, hl.missing(hl.tfloat64))
    z = hl.if_else(chisq == 0, 0.0, hl.sqrt(chisq))
    return hl.if_else(obs > exp, -z, z)


def _effective_cfg(cfg: SimpleNamespace) -> SimpleNamespace:
    """Return cfg.mu_cfg if defined (used by mu/element steps), otherwise cfg itself."""
    return getattr(cfg, "mu_cfg", cfg)


def apply_coverage_filter(
        ht: hl.Table,
        cfg: SimpleNamespace,
        coverage_ht: Optional[hl.Table] = None,
) -> hl.Table:
    """Apply the config-appropriate coverage filter.

    Config A: filter on coverage_mean embedded in the prefiltered context HT.
    Config B: filter on over_10 from the external gnomAD coverage HT.
    """
    if cfg.use_external_coverage:
        if coverage_ht is None:
            raise ValueError("coverage_ht must be provided when use_external_coverage=True (Config B)")
        return ht.filter(coverage_ht[ht.locus].over_10 >= cfg.min_over_10)
    elif cfg.coverage_mean_range is not None:
        lo, hi = cfg.coverage_mean_range
        return ht.filter((ht.coverage_mean >= lo) & (ht.coverage_mean <= hi))
    return ht


def apply_region_filters(ht: hl.Table, cfg: SimpleNamespace) -> hl.Table:
    """Apply region exclusion filters.

    When both lcr and segdup are excluded (Config B), uses the pre-merged
    blacklist_gap2_lcr_segdup.merged.bed (single BED import).
    Otherwise applies blacklist_gap2.bed only (Config A / Basic).
    """
    if cfg.exclude_blacklist_lcr_segdup:
        return filter_by_bed(ht, BLACKLIST_LCR_SEGDUP_PATH)
    return filter_by_bed(ht, BLACKLIST_BED_PATH)


def discretize_methyl(
        ht: hl.Table,
        logit_field: str,
        thresholds: List[Tuple[float, int]],
) -> hl.expr.Expression:
    """Map a continuous logit score to integer methyl levels via a case expression."""
    expr = ht[logit_field]
    case = hl.case()
    for thr, lvl in thresholds:
        case = case.when(expr > thr, lvl)
    return case.default(0)


# ──────────────────────────────────────────────────────────────────
# Step 1: Prefilter to chrX nonPAR
# ──────────────────────────────────────────────────────────────────

def prefilter_input_ht(
        ht: hl.Table,
        intervals: Optional[Sequence[hl.Interval]] = None,
        coverage_outliers: bool = False,
) -> hl.Table:
    """Restrict to chrX nonPAR intervals and optionally remove coverage outliers.

    Config-specific coverage and region filters are applied separately via
    apply_coverage_filter / apply_region_filters at each step that needs them.
    """
    if intervals is not None:
        ht = filter_by_intervals(ht, intervals, reference_genome=REFGENOME_NAME)
    if coverage_outliers and "coverage_mean" in ht.row.dtype.fields:
        ht = remove_coverage_outliers(ht)
    return ht

# ──────────────────────────────────────────────────────────────────
# Step 2: Downsampling
# ──────────────────────────────────────────────────────────────────
def downsample_prefilter_input_ht(
        context_ht: hl.Table,
        genome_ht: hl.Table,
):
    """Filter context and genome HTs to the 1000-sample downsampling AC cutoff.

    Locates the global adj downsampling frequency index for 1000 samples in freq_meta,
    then applies AC <= 5 as a proxy for rare variants within that downsampling.

    - context_ht (denominator): keeps loci absent from genome_ht OR with AC <= 5 and
      passing QC filters at the 1000-sample downsampling level.
    - genome_ht (numerator): keeps only AC <= 5, QC-passing variants at the same level.

    :param context_ht: Prefiltered context HT (locus, alleles keyed).
    :param genome_ht: Prefiltered genome HT (must have freq and freq_meta fields).
    :return: Tuple of (filtered context_ht, filtered genome_ht).
    """
    grouping_variables = (("cpg", "variant_type", "variant_type_model",
                           )
    )
    context_ht = context_ht.select('context', 'ref', 'alt', *grouping_variables)
    genome_ht = genome_ht.select('context', 'ref', 'alt', 'freq', 'pass_filters', *grouping_variables)

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

    return context_ht, genome_ht


# ──────────────────────────────────────────────────────────────────
# Step 3a: Annotate raw per-stage methylation values on CpG loci
# ──────────────────────────────────────────────────────────────────

def annotate_methyl_step(
        context_ht: hl.Table,
        meth_ht: hl.Table,
        output_path: str,
        overwrite: bool = False,
) -> hl.Table:
    """Filter context_ht to CpG loci and annotate raw per-stage methylation values.
    Missing values are filled with per-stage means. Config-independent."""
    context_ht = context_ht.filter(context_ht.cpg).key_by("locus")
    context_ht = context_ht.select("alleles", "context", "ref", "alt")

    m = meth_ht[context_ht.locus]
    context_ht = context_ht.annotate(
        **{f"{stage}_raw": m[field] for stage, field in STAGES.items()}
    )

    means = meth_ht.aggregate(
        hl.struct(**{stage: hl.agg.mean(meth_ht[field]) for stage, field in STAGES.items()})
    )
    context_ht = context_ht.annotate(
        **{stage: hl.or_else(context_ht[f"{stage}_raw"], means[stage]) for stage in STAGES}
    )

    # Discretize into 0/1/2 levels
    def to_level(x):
        return (
            hl.case()
            .when(x > 0.6, 2)
            .when(x > 0.2, 1)
            .default(0)
        )

    level_exprs = {
        f"{stage}_methyl_level": to_level(context_ht[stage])
        for stage in STAGES.keys()
    }
    context_ht = context_ht.annotate(**level_exprs)

    # Summarize max + sum
    levels = [context_ht[f"{stage}_methyl_level"] for stage in STAGES.keys()]

    context_ht = context_ht.annotate(
        methyl_level_max=hl.max(levels),
        methyl_level_sum=hl.sum(levels),
    )

    return context_ht.key_by("locus", "alleles").checkpoint(output_path, overwrite=overwrite)


# ──────────────────────────────────────────────────────────────────
# Step 3b: Build logistic regression training dataset
# ──────────────────────────────────────────────────────────────────

def annotate_obs_step(
        cpg_ht: hl.Table,
        context_ht: hl.Table,
        genome_ht: hl.Table,
        output_path: str,
        cfg: SimpleNamespace,
        coverage_ht: Optional[hl.Table] = None,
        overwrite: bool = False,
        af_cutoff: float = 0.001,
) -> hl.Table:
    """Build training data: CpG loci passing config-specific coverage/region filters,
    annotated with obs=True if a rare QC-passing variant exists in genome_ht.

    For Config B, applies an additional median_approx filter on the external coverage HT
    (matching get_met_14_vs_obs_x.py active state: over_10>=0.95 AND median_approx 22-24).
    This extra filter is training-only — the element step uses over_10 alone.
    """
    RAW_VARS = tuple(f"{stage}_raw" for stage in STAGES)
    LEVEL_VARS = tuple(f"{stage}_methyl_level" for stage in STAGES)

    grouping_variables = RAW_VARS
    grouping_variables2 = RAW_VARS + LEVEL_VARS + ("context",)

    # Drops any row where any of the *_raw fields is missing
    # Keeps only the columns listed in grouping_variables2
    cpg_ht = cpg_ht.filter(
        hl.all(lambda x: x,
               [hl.is_defined(cpg_ht[x]) for x in grouping_variables])
    ).select(*grouping_variables2)

    context_ht = apply_coverage_filter(context_ht, cfg, coverage_ht)
    # Training-only: also filter by median_approx when using external coverage (Config B)
    if (
        getattr(cfg, "training_median_approx_range", None) is not None
        and coverage_ht is not None
    ):
        lo, hi = cfg.training_median_approx_range
        context_ht = context_ht.filter(
            (coverage_ht[context_ht.locus].median_approx >= lo)
            & (coverage_ht[context_ht.locus].median_approx <= hi)
        )
    context_ht = apply_region_filters(context_ht, cfg)
    cpg_ht = cpg_ht.semi_join(context_ht)

    genome_join = genome_ht[cpg_ht.key]

    # Keep CpG loci that either:
    #  - have no variant observed, OR
    #  - have only rare (AF ≤ 0.001), QC-passing variants
    cpg_ht = cpg_ht.filter(
        hl.is_missing(genome_join) | (
                (genome_join.freq[0].AF <= af_cutoff) &
                genome_join.pass_filters
        )
    )

    # Filter genome variants to rare, QC-passing variants only
    genome_ht = genome_ht.filter(
        (genome_ht.freq[0].AF <= af_cutoff) &
        genome_ht.pass_filters
    )

    # Annotate each CpG locus with whether a qualifying variant is observed
    cpg_ht = cpg_ht.annotate(
        obs=hl.is_defined(genome_ht[cpg_ht.key])
    )

    return cpg_ht.checkpoint(output_path, overwrite=overwrite)

# ──────────────────────────────────────────────────────────────────
# Step 3c: Fit logistic regression → weighted logit level HT
# ──────────────────────────────────────────────────────────────────

def compute_coeff_step(
        obs_ht: hl.Table,
        cpg_ht: hl.Table,
        output_path: str,
        cfg: SimpleNamespace,
        overwrite: bool = False,
) -> hl.Table:
    """Compute per-locus weighted logit score and discretize using cfg.methyl_thresholds.

    If cfg.hardcoded_coeff is set, uses those log-odds coefficients directly instead of
    fitting a logistic regression. Otherwise fits obs ~ stage_raw fields via statsmodels.
    """
    import numpy as np

    x_vars = [f"{stage}_raw" for stage in STAGES]

    hardcoded = getattr(cfg, "hardcoded_coeff", None)
    if hardcoded is not None:
        print("Using hardcoded coefficients (skipping logistic regression):")
        for k, v in hardcoded.items():
            print(f"  {k}: {v}")
        weights = {k: float(np.exp(hardcoded[k])) for k in x_vars}
    else:
        import statsmodels.api as sm

        df = obs_ht.select("obs", *x_vars).to_pandas()
        df["obs"] = df["obs"].astype(int)

        X = sm.add_constant(df[x_vars].astype(float), has_constant="add")
        fit = sm.Logit(df["obs"], X).fit(disp=False)
        coeff = fit.params

        print("Logistic regression coefficients:")
        print(fit.summary2().tables[1][["Coef.", "Std.Err.", "P>|z|"]])

        weights = {k: float(np.exp(coeff[k])) for k in x_vars}

    w_arr = hl.array([hl.float64(weights[f"{s}_raw"]) for s in STAGES])
    x_arr = hl.array([hl.float64(cpg_ht[s]) for s in STAGES])
    cpg_ht = cpg_ht.annotate(
        methyl14_weighted_logit=hl.sum(x_arr * w_arr) / hl.sum(w_arr)
    )
    cpg_ht = cpg_ht.annotate(
        methyl14_weighted_logit_level=discretize_methyl(
            cpg_ht, "methyl14_weighted_logit", cfg.methyl_thresholds
        )
    )
    cpg_ht = cpg_ht.select("context", "ref", "alt","methyl14_weighted_logit", "methyl14_weighted_logit_level")

    return cpg_ht.checkpoint(output_path, overwrite=overwrite)


# ──────────────────────────────────────────────────────────────────
# Step 4a: Count observed/possible for mu estimation
# ──────────────────────────────────────────────────────────────────

def calculate_chrx_nonpar_po(
        context_ht: hl.Table,
        genome_ht: hl.Table,
        met_ht: hl.Table,
        cfg: SimpleNamespace,
        coverage_ht: Optional[hl.Table] = None,
        af_cutoff: float = 0.001,
):
    """Count observed and possible variants grouped by (context, ref, alt, methyl_level)
    after applying config-specific coverage and region filters."""
    mcfg = _effective_cfg(cfg)
    context_ht = apply_coverage_filter(context_ht, mcfg, coverage_ht)
    context_ht = apply_region_filters(context_ht, mcfg)
    genome_ht = genome_ht.semi_join(context_ht)

    context_ht = context_ht.annotate(
        methyl_level=hl.or_else(met_ht[context_ht.key].methyl14_weighted_logit_level, 0)
    )
    genome_ht = genome_ht.annotate(
        methyl_level=hl.or_else(met_ht[genome_ht.key].methyl14_weighted_logit_level, 0)
    )

    context_ht = context_ht.select("context", "ref", "alt", "methyl_level")
    genome_ht = genome_ht.select("context", "ref", "alt", "methyl_level", "freq", "pass_filters")

    genome_join = genome_ht[context_ht.key]
    context_ht = context_ht.filter(
        hl.is_missing(genome_join)
        | ((genome_join.freq[0].AF <= af_cutoff) & genome_join.pass_filters)
    )
    genome_ht = genome_ht.filter((genome_ht.freq[0].AF <= af_cutoff) & genome_ht.pass_filters)

    observed_ht = (
        genome_ht
        .group_by(context=genome_ht.context, ref=genome_ht.ref,
                  alt=genome_ht.alt, methyl_level=genome_ht.methyl_level)
        .aggregate(variant_count=hl.agg.count())
        .naive_coalesce(10)
    )
    possible_ht = (
        context_ht
        .group_by(context=context_ht.context, ref=context_ht.ref,
                  alt=context_ht.alt, methyl_level=context_ht.methyl_level)
        .aggregate(variant_count=hl.agg.count())
        .naive_coalesce(10)
    )

    return observed_ht, possible_ht


# ──────────────────────────────────────────────────────────────────
# Step 4b: Estimate mu from downsampled data
# ──────────────────────────────────────────────────────────────────

def compute_mu_step(
        context_ds_ht: hl.Table,
        genome_ds_ht: hl.Table,
        met_ht: hl.Table,
        output_bucket: str,
        cfg: SimpleNamespace,
        coverage_ht: Optional[hl.Table] = None,
) -> str:
    """Estimate per-(context, ref, alt, methyl_level) mutation rates from downsampled data.
    Config B applies the same spatial filters as the element step to the downsampled HTs."""
    if cfg.filter_downsampled:
        context_ds_ht = apply_coverage_filter(context_ds_ht, cfg, coverage_ht)
        context_ds_ht = apply_region_filters(context_ds_ht, cfg)
        genome_ds_ht = genome_ds_ht.semi_join(context_ds_ht)

    context_ds_ht = context_ds_ht.annotate(
        methyl_level=hl.if_else(
            hl.is_missing(met_ht[context_ds_ht.key]), 0,
            met_ht[context_ds_ht.key].methyl14_weighted_logit_level,
        )
    )
    genome_ds_ht = genome_ds_ht.annotate(
        methyl_level=hl.if_else(
            hl.is_missing(met_ht[genome_ds_ht.key]), 0,
            met_ht[genome_ds_ht.key].methyl14_weighted_logit_level,
        )
    )

    observed_ht = count_variants(
        genome_ds_ht,
        count_downsamplings=["global", "nfe", "afr"],
        additional_grouping=("methyl_level",),
        omit_methylation=True,
    )

    # key field name must match count_variants output (methyl_level, not methylation_level)
    possible_ht = context_ds_ht.group_by(
        context=context_ds_ht.context,
        ref=context_ds_ht.ref,
        alt=context_ds_ht.alt,
        methyl_level=context_ds_ht.methyl_level,
    ).aggregate(variant_count=hl.agg.count())

    observed_ht = observed_ht.annotate(
        possible_variants=possible_ht[observed_ht.key].variant_count
    )

    total_bases = observed_ht.aggregate(hl.agg.sum(observed_ht.possible_variants)) // 3
    total_mu = 1.2e-08
    correction_factors = observed_ht.aggregate(
        total_mu / (hl.agg.array_sum(observed_ht.downsampling_counts_global) / total_bases)
    )

    observed_ht = annotate_variant_types(observed_ht.annotate(
        downsamplings_frac_observed=(
            observed_ht.downsampling_counts_global / observed_ht.possible_variants
        ),
        downsamplings_mu_snp=(
            hl.literal(correction_factors)
            * observed_ht.downsampling_counts_global
            / observed_ht.possible_variants
        ),
    ))

    downsamplings = list(map(lambda x: x[1], get_downsamplings(observed_ht)))
    index_1kg = downsamplings.index(1000)

    observed_ht = observed_ht.annotate(
        observed_1kg=observed_ht.downsampling_counts_global[index_1kg],
        proportion_observed_1kg=(
            observed_ht.downsampling_counts_global[index_1kg] / observed_ht.possible_variants
        ),
        mu=observed_ht.downsamplings_mu_snp[index_1kg],
    )

    mu_path = f"{output_bucket}/mu_by_context_methyl_downsampled_1000_chrX_nonpar_config{cfg.name}.txt"
    observed_ht.select(
        "transition", "cpg", "variant_type", "variant_type_model",
        "possible_variants", "observed_1kg", "proportion_observed_1kg", "mu",
    ).export(mu_path)

    return mu_path


# ──────────────────────────────────────────────────────────────────
# Step 4c: Fit po ~ mu → fitted_po per context/methyl group
# ──────────────────────────────────────────────────────────────────

def fit_mutation_rate(
        pos_txt: str,
        obs_txt: str,
        mu_txt: str,
        output_path: str,
) -> str:
    """Fit log(1 - proportion_observed) ~ mu with inverse-variance weights.
    Produces fitted_po per (context, ref, alt, methyl_level) group."""
    import csv
    import numpy as np
    import pandas as pd
    from sklearn.metrics import r2_score

    def _sem(x, n):
        p = float(x) / float(n)
        return np.sqrt(p * (1 - p) / float(n))

    index_cols = ["context", "ref", "alt", "methyl_level"]

    df_possible = pd.read_csv(pos_txt, sep="\t", index_col=index_cols).rename(
        columns={"variant_count": "possible"}
    )
    df_observed = pd.read_csv(obs_txt, sep="\t", index_col=index_cols).rename(
        columns={"variant_count": "observed"}
    )

    df_po = df_possible.join(df_observed)
    df_po["proportion_observed"] = df_po["observed"] / df_po["possible"]

    df_mu = pd.read_csv(mu_txt, sep="\t").set_index(index_cols)
    df_mu["sem"] = df_mu.apply(
        lambda row: _sem(row["observed_1kg"], row["possible_variants"]), axis=1
    )

    df_po = df_po.join(df_mu[["mu", "sem"]])
    df_po = df_po.replace([float("inf"), float("-inf")], float("nan")).dropna(
        subset=["mu", "sem", "proportion_observed"]
    )
    df_po["sem"] = df_po["sem"].clip(lower=1e-12)

    eps = 1e-12
    p = df_po["proportion_observed"].clip(eps, 1 - eps)
    y = np.log(1 - p)
    w = 1.0 / df_po["sem"]  # matches Siwei: w=1/sem (not 1/sem^2 as recommended by chatGPT for inverse-variance weighting )

    A, B = np.polyfit(df_po["mu"], y, 1, w=w)
    df_po["fitted_po"] = (1 - np.exp(B + A * df_po["mu"])).clip(0, 1)

    r2 = r2_score(df_po["proportion_observed"], df_po["fitted_po"], sample_weight=w)
    print(f"R² = {r2:.4f},  A = {A:.6e},  B = {B:.6e}")

    df_out = df_po.drop(columns=["sem"]).reset_index()
    df_out.to_csv(output_path, sep="\t", quoting=csv.QUOTE_NONE, header=True, index=False)

    # Also write a Hail HT keyed by (context, ref, alt, methyl_level) for fast lookup
    # in compute_element_z_hail_step via hl.read_table instead of hl.import_table.
    ht = hl.Table.from_pandas(df_out).key_by("context", "ref", "alt", "methyl_level")
    ht.write(output_path.replace(".txt", ".ht"), overwrite=True)

    return output_path


# ──────────────────────────────────────────────────────────────────
# Steps 5: per-element obs/expected/Z in a single Hail pass
# ──────────────────────────────────────────────────────────────────

def compute_element_z_hail_step(
        context_ht: hl.Table,
        genome_ht: hl.Table,
        met_ht: hl.Table,
        element_bed: hl.Table,
        mr_txt: str,
        output_txt: str,
        cfg: SimpleNamespace,
        coverage_ht: Optional[hl.Table] = None,
        af_cutoff: float = 0.001,
        overwrite: bool = False,
) -> None:
    """Compute per-element observed, expected, oe and Z score in a single Hail pass.

    Annotates fitted_po per locus from the mutation-rate table, explodes by element_id,
    then aggregates:
      observed = sum(is_observed)  — number of loci with a rare QC-passing variant
      expected = sum(fitted_po)    — sum of per-locus mutation probabilities
      possible = count()           — total loci considered

    rr = 1 for chrX nonPAR, so fitted_po_adj = fitted_po (no relative-risk correction).

    Z formula (signed chi-square, matches Siwei's calculate_oe_z_from_po_x.submit.py):
        chisq = (obs - exp)^2 / exp
        z = -sqrt(chisq) if obs > exp  else  sqrt(chisq)
    """
    mcfg = _effective_cfg(cfg)

    context_ht = apply_coverage_filter(context_ht, mcfg, coverage_ht)
    context_ht = apply_region_filters(context_ht, mcfg)

    # Keep genome minimal; just for observed lookups
    genome_ht = genome_ht.semi_join(context_ht)
    genome_ht = genome_ht.select("freq", "pass_filters")  # drop element_id here

    # methyl
    context_ht = context_ht.annotate(
        methyl_level=hl.or_else(met_ht[context_ht.key].methyl14_weighted_logit_level, 0)
    )

    # mutation rate (ideally read_table from .ht)
    mr_ht = hl.read_table(mr_txt.replace(".txt", ".ht"))  # best
    # or fallback: import_table + checkpoint once
    context_ht = context_ht.annotate(
        fitted_po=mr_ht[context_ht.context, context_ht.ref, context_ht.alt, context_ht.methyl_level].fitted_po,
        rr=1.0,
    )

    # Filter context_ht: remove loci where a common or QC-failing variant exists
    genome_join = genome_ht[context_ht.key]
    context_ht = context_ht.filter(
        hl.is_missing(genome_join)
        | ((genome_join.freq[0].AF <= af_cutoff) & genome_join.pass_filters)
    )

    # Build observed ht once (genome already semi_joined to original context above)
    obs_ht = genome_ht.filter((genome_ht.freq[0].AF <= af_cutoff) & genome_ht.pass_filters)

    # annotate observed BEFORE element annotation
    context_ht = context_ht.annotate(
        is_observed=hl.is_defined(obs_ht[context_ht.key]),
        fitted_po_adj=context_ht.fitted_po * context_ht.rr,
    )

    # checkpoint here (optional but often good)
    context_ht = context_ht.select("fitted_po", "is_observed", "fitted_po_adj").checkpoint(
        output_txt.replace(".txt", "_tmp_locus.ht"), overwrite=True
    )

    # NOW annotate elements (only once, only on context)
    context_ht = annotate_genome_element(context_ht, element_bed)

    # Aggregate by element
    by_element_ht = context_ht.group_by("element_id").aggregate(
        observed=hl.agg.sum(hl.int(context_ht.is_observed)),
        expected=hl.agg.sum(context_ht.fitted_po),
        expected_adj=hl.agg.sum(context_ht.fitted_po_adj),
        possible=hl.agg.count(),
    )

    # compute z inline
    by_element_ht = by_element_ht.annotate(
        gnocchi=signed_chisq_z(by_element_ht.observed, by_element_ht.expected),
        gnocchi_adj=signed_chisq_z(by_element_ht.observed, by_element_ht.expected_adj),
    )

    # QC annotation: pct_pass, mean_coverage, pass_qc (matches gnocchi_autosome_par.py)
    pass_ht = (
        hl.import_table(ELEMENT_PASS_TXT, no_header=True, delimiter="\t", impute=True)
        .rename({"f0": "element_id", "f1": "pct_pass"})
        .key_by("element_id")
    )
    cov_ht = (
        hl.import_table(ELEMENT_COV_TXT, no_header=True, delimiter="\t", impute=True)
        .rename({"f0": "element_id", "f1": "mean_coverage"})
        .key_by("element_id")
    )
    cov_lo, cov_hi = CHRX_NONPAR_COV_QC_RANGE
    pos = by_element_ht.element_id.split("-")
    by_element_ht = by_element_ht.annotate(
        pct_pass=pass_ht[by_element_ht.key].pct_pass,
        mean_coverage=cov_ht[by_element_ht.key].mean_coverage,
        chrom=pos[0],
        start=hl.int(pos[1]),
        end=hl.int(pos[2]),
    )
    by_element_ht = by_element_ht.annotate(
        pass_qc=hl.if_else(
            (by_element_ht.pct_pass >= 0.8)
            & (by_element_ht.mean_coverage >= cov_lo)
            & (by_element_ht.mean_coverage <= cov_hi)
            & (by_element_ht.possible >= 1000),
            True,
            False,
        )
    )
    by_element_ht = by_element_ht.select(
        "chrom", "start", "end",
        "observed", "possible", "expected", "expected_adj",
        "gnocchi", "gnocchi_adj", "pct_pass", "mean_coverage", "pass_qc",
    )

    by_element_ht.checkpoint(output_txt.replace(".txt", ".ht"), overwrite=overwrite).export(output_txt)

# ──────────────────────────────────────────────────────────────────
# Resource paths
# ──────────────────────────────────────────────────────────────────

def get_pipeline_resources(output_bucket: str, cfg: SimpleNamespace) -> SimpleNamespace:
    """Return all input/output GCS paths for the chrX nonPAR gnocchi pipeline."""
    sfx = f"config{cfg.name}"
    return SimpleNamespace(
        # --- Inputs ---
        context_ht=f"{PUBLIC_BUCKET}/context_prepared.ht",
        genome_ht=f"{PUBLIC_BUCKET}/genome_prepared.ht",
        # Raw methylation HT (14 germline stages per CpG locus).
        # Original: gs://gnomad-nc-constraint-v31/methylation/methylation_dev_stages_liftover.ht
        # Update this path to match your copy:
        meth_ht=f"{PUBLIC_BUCKET}/methylation",
        coverage_ht=COVERAGE_HT_PATH,

        # --- Step 1: Prefilter (config-independent, shared between both configs) ---
        prefiltered_context_ht=f"{output_bucket}/context_prefiltered_chrx_nonpar.ht",
        prefiltered_genome_ht=f"{output_bucket}/genome_prefiltered_chrx_nonpar.ht",

        # --- Step 2: Downsampling (config-independent) ---
        downsampled_context_ht=f"{output_bucket}/context_prefiltered_chrx_nonpar_downsampled_1000.ht",
        downsampled_genome_ht=f"{output_bucket}/genome_prefiltered_chrx_nonpar_downsampled_1000.ht",

        # --- Step 3a: Raw methylation annotation (config-independent) ---
        methyl_cpg_ht=f"{output_bucket}/context_prepared_chrX_nonPAR_methyl_CpG.ht",

        # --- Steps 3b–6: config-specific (suffixed with configA or configB) ---
        methyl_cpg_obs_ht=f"{output_bucket}/context_prepared_chrX_nonPAR_methyl_CpG_obs_{sfx}.ht",
        methyl_level_ht=f"{output_bucket}/context_prepared_chrX_nonPAR_methyl_logit_level_{sfx}.ht",
        obs_txt=f"{output_bucket}/observed_counts_by_context_methyl_chrX_nonpar_{sfx}.txt",
        pos_txt=f"{output_bucket}/possible_counts_by_context_methyl_chrX_nonpar_{sfx}.txt",
        mu_txt=f"{output_bucket}/mu_by_context_methyl_downsampled_1000_chrX_nonpar_{sfx}.txt",
        mutation_rate_txt=f"{output_bucket}/mutation_rate_by_context_methyl_chrX_nonpar_{sfx}.txt",
        element_z_txt=f"{output_bucket}/gnocchi_1kb_chrX_nonpar_{sfx}.txt",
    )


# ──────────────────────────────────────────────────────────────────
# Shared pipeline runner and CLI
# ──────────────────────────────────────────────────────────────────

def run_pipeline(args, cfg: SimpleNamespace) -> None:
    """Execute the requested pipeline steps using the provided config."""
    intervals = parse_intervals(CHRX_NONPAR_INTERVALS)

    if args.dry_run:
        cov_desc = (
            f"over_10 >= {cfg.min_over_10} (external coverage HT)"
            if cfg.use_external_coverage
            else f"coverage_mean {cfg.coverage_mean_range} (prefiltered HT field)"
        )
        bl_desc = "blacklist_gap2 + lcr + segdup (merged)" if cfg.exclude_blacklist_lcr_segdup else "blacklist_gap2"
        lvl_desc = f"logit{'13' if len(cfg.methyl_thresholds) == 13 else '12'}"
        print("=== DRY RUN ===")
        print(f"Config:        {cfg.name}")
        print(f"Coverage:      {cov_desc}")
        print(f"Blacklist:     {bl_desc}")
        print(f"Methyl levels: {lvl_desc}")
        print(f"AF cutoff:     {args.af_cutoff}")
        print(f"Output bucket: {args.output_bucket}")
        print("==============")
        return

    hl.default_reference(REFGENOME_NAME)

    af_cutoff = args.af_cutoff
    output_bucket = args.output_bucket
    overwrite = args.overwrite

    res = get_pipeline_resources(output_bucket, cfg)

    # Load external coverage HT once (only needed for Config B)
    coverage_ht = hl.read_table(res.coverage_ht) if cfg.use_external_coverage else None

    # Step 1: Prefilter
    if args.prefilter:
        logger.info("Step 1 [prefilter]: starting")
        t0 = time.time()
        context_ht = hl.read_table(res.context_ht)
        genome_ht = hl.read_table(res.genome_ht)

        context_ht = prefilter_input_ht(context_ht, intervals=intervals, coverage_outliers=True)
        genome_ht = prefilter_input_ht(genome_ht, intervals=intervals, coverage_outliers=True)

        ctx_fields = ("context", "ref", "alt", "coverage_mean", "cpg", "variant_type", "variant_type_model")
        gnm_fields = ctx_fields + ("freq", "pass_filters")

        # Steps 1-3a outputs are config-independent; never overwrite once written
        context_ht.select(*ctx_fields).naive_coalesce(1000).write(
            res.prefiltered_context_ht, overwrite=False
        )
        genome_ht.select(*gnm_fields).naive_coalesce(1000).write(
            res.prefiltered_genome_ht, overwrite=False
        )
        logger.info("Step 1 [prefilter]: done in %.1fs", time.time() - t0)

    # Step 2: Downsampling
    if args.downsampling:
        logger.info("Step 2 [downsampling]: starting")
        t0 = time.time()
        context_ht = hl.read_table(res.prefiltered_context_ht)
        genome_ht = hl.read_table(res.prefiltered_genome_ht)

        context_ht, genome_ht = downsample_prefilter_input_ht(context_ht, genome_ht)

        context_ht.write(res.downsampled_context_ht, overwrite=False)
        genome_ht.write(res.downsampled_genome_ht, overwrite=False)
        logger.info("Step 2 [downsampling]: done in %.1fs", time.time() - t0)

    # Step 3a: Annotate raw methylation
    if args.annotate_methyl:
        logger.info("Step 3a [annotate-methyl]: starting")
        t0 = time.time()
        context_ht = hl.read_table(res.context_ht)
        context_ht = prefilter_input_ht(context_ht, intervals=intervals)
        meth_ht = hl.read_table(res.meth_ht)
        annotate_methyl_step(context_ht, meth_ht, res.methyl_cpg_ht, overwrite=False)
        logger.info("Step 3a [annotate-methyl]: done in %.1fs", time.time() - t0)

    # Step 3b: Build logit training dataset
    if args.annotate_obs:
        logger.info("Step 3b [annotate-obs] config=%s: starting", cfg.name)
        t0 = time.time()
        cpg_ht = hl.read_table(res.methyl_cpg_ht)
        context_ht = hl.read_table(res.prefiltered_context_ht)
        genome_ht = hl.read_table(res.prefiltered_genome_ht)
        annotate_obs_step(
            cpg_ht, context_ht, genome_ht,
            output_path=res.methyl_cpg_obs_ht,
            cfg=cfg, coverage_ht=coverage_ht,
            overwrite=overwrite, af_cutoff=af_cutoff,
        )
        logger.info("Step 3b [annotate-obs]: done in %.1fs", time.time() - t0)

    # Step 3c: Compute weighted logit level HT
    # If cfg.hardcoded_coeff is set, obs_ht is not needed (coefficients are fixed).
    if args.compute_coeff:
        logger.info("Step 3c [compute-coeff] config=%s: starting", cfg.name)
        t0 = time.time()
        obs_ht = (
            None if getattr(cfg, "hardcoded_coeff", None) is not None
            else hl.read_table(res.methyl_cpg_obs_ht)
        )
        cpg_ht = hl.read_table(res.methyl_cpg_ht)
        compute_coeff_step(obs_ht, cpg_ht, res.methyl_level_ht, cfg=cfg, overwrite=overwrite)
        logger.info("Step 3c [compute-coeff]: done in %.1fs", time.time() - t0)

    # Step 4: Count obs/possible, estimate mu, fit po~mu
    if args.compute_mu:
        logger.info("Step 4 [compute-mu] config=%s: starting", cfg.name)
        t0 = time.time()
        context_ht = hl.read_table(res.prefiltered_context_ht)
        genome_ht = hl.read_table(res.prefiltered_genome_ht)
        context_ds_ht = hl.read_table(res.downsampled_context_ht)
        genome_ds_ht = hl.read_table(res.downsampled_genome_ht)
        met_ht = hl.read_table(res.methyl_level_ht)

        obs_ht, pos_ht = calculate_chrx_nonpar_po(
            context_ht, genome_ht, met_ht,
            cfg=cfg, coverage_ht=coverage_ht, af_cutoff=af_cutoff,
        )
        obs_ht.checkpoint(res.obs_txt.replace(".txt", ".ht"), overwrite=overwrite).export(res.obs_txt)
        pos_ht.checkpoint(res.pos_txt.replace(".txt", ".ht"), overwrite=overwrite).export(res.pos_txt)
        logger.info("Step 4 [compute-mu]: obs/possible exported (%.1fs so far)", time.time() - t0)

        mu_txt = compute_mu_step(
            context_ds_ht, genome_ds_ht, met_ht,
            output_bucket=output_bucket, cfg=cfg, coverage_ht=coverage_ht,
        )
        fit_mutation_rate(
            pos_txt=res.pos_txt, obs_txt=res.obs_txt,
            mu_txt=mu_txt, output_path=res.mutation_rate_txt,
        )
        logger.info("Step 4 [compute-mu]: done in %.1fs", time.time() - t0)

    # Step 5: per-element obs/expected/Z in a single Hail pass
    if args.compute_element_z:
        logger.info("Step 5 [compute-element-z] config=%s: starting", cfg.name)
        t0 = time.time()
        context_ht = hl.read_table(res.prefiltered_context_ht)
        genome_ht = hl.read_table(res.prefiltered_genome_ht)
        met_ht = hl.read_table(res.methyl_level_ht)
        element_bed = hl.import_bed(
            BED_1KB_PATH, reference_genome=REFGENOME_NAME, skip_invalid_intervals=True
        )
        compute_element_z_hail_step(
            context_ht, genome_ht, met_ht, element_bed,
            mr_txt=res.mutation_rate_txt,
            output_txt=res.element_z_txt,
            cfg=cfg, coverage_ht=coverage_ht,
            af_cutoff=af_cutoff, overwrite=overwrite,
        )
        logger.info("Step 5 [compute-element-z]: done in %.1fs → %s", time.time() - t0, res.element_z_txt)


def parse_args(description: str) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--af-cutoff", type=float, default=0.001,
        help="AF cutoff for observed variants (default: 0.001)",
    )
    parser.add_argument(
        "--output-bucket", required=True,
        help="GCS bucket/path for outputs (e.g. gs://qin-gnocchi/gnocchi_files)",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print config and exit without running Hail")

    parser.add_argument("--prefilter", action="store_true",
                        help="Step 1: filter context/genome to chrX nonPAR intervals")
    parser.add_argument("--downsampling", action="store_true",
                        help="Step 2: filter to downsampled 1000 AC cutoff")
    parser.add_argument("--annotate-methyl", action="store_true",
                        help="Step 3a: annotate CpG loci with raw per-stage methyl values")
    parser.add_argument("--annotate-obs", action="store_true",
                        help="Step 3b: apply config filters and annotate obs (0/1)")
    parser.add_argument("--compute-coeff", action="store_true",
                        help="Step 3c: fit logistic regression, compute weighted logit level HT")
    parser.add_argument("--compute-mu", action="store_true",
                        help="Step 4: count obs/possible, estimate mu, fit po~mu")
    parser.add_argument("--compute-element-z", action="store_true",
                        help="Steps 5+6: annotate fitted_po/element_id and compute per-element oe and Z scores (Hail)")

    return parser.parse_args()