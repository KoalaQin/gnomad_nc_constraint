"""
Config A — chrX nonPAR gnocchi pipeline
Matches the most recently active code state in Siwei's scripts
(calculate_po_x.py + calculate_po_count_by_element_x.py active state).

Parameters:
  Coverage:      coverage_mean 23–24 (field in prefiltered context HT)
  Blacklist:     blacklist_gap2.bed only
  Methyl levels: logit12 (top bin >0.80, levels 0–12)
  Downsampled:   no spatial filters on mu estimation

Usage:
  python gnocchi_chrX_nonPAR_configA.py --output-bucket gs://... [step flags]

Step flags (can combine multiple):
  --prefilter          Step 1:  filter context/genome to chrX nonPAR intervals
  --downsampling       Step 2:  filter to downsampled 1000 AC cutoff
  --annotate-methyl    Step 3a: annotate CpG loci with raw per-stage methyl values
  --annotate-obs       Step 3b: apply Config A filters; annotate obs (0/1)
  --compute-coeff      Step 3c: fit logistic regression → logit12 level HT
  --compute-mu         Step 4:  count obs/possible, estimate mu, fit po~mu
  --compute-element-z  Steps 5+6: annotate fitted_po/element_id per locus; compute per-element oe and Z scores
"""

from types import SimpleNamespace
from gnocchi_chrX_nonPAR_utils import (
    METHYL_THRESHOLDS_LOGIT12,
    run_pipeline,
    parse_args,
)

CONFIG_A = SimpleNamespace(
    name="A",
    # Coverage: use coverageean field embedded in prefiltered context HT
    use_external_coverage=False,
    coverage_mean_range=(23.0, 24.0),
    min_over_10=None,
    # Region exclusion: blacklist_gap2 only
    exclude_blacklist_lcr_segdup=False,
    # Methylation discretization: logit12 (top threshold >0.80, 12 levels)
    methyl_thresholds=METHYL_THRESHOLDS_LOGIT12,
    # Downsampled mu: no spatial filters (consistent with Siwei's active calculate_po_downsampled_x.py)
    filter_downsampled=False,
    # No median_approx filter in training (Config A uses coverage_mean from prefiltered HT, not external HT)
    training_median_approx_range=None,
)

if __name__ == "__main__":
    args = parse_args(description="Gnocchi chrX nonPAR pipeline — Config A (logit12, cov_mean 23-24, gap2 only)")
    run_pipeline(args, CONFIG_A)