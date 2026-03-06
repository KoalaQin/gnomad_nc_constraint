"""
Config Basic — chrX nonPAR gnocchi pipeline.
Training (steps 3b/3c): no coverage or blacklist filters.
Mu/element steps (4–6): coverage_mean 23–24, blacklist_gap2 only (same spatial filters as Config A).

Parameters (training):
  Coverage:      none
  Blacklist:     none
  Methyl levels: logit12

Parameters (mu/element):
  Coverage:      coverage_mean 23–24 (prefiltered HT field)
  Blacklist:     blacklist_gap2 only

Usage:
  python gnocchi_chrX_nonPAR_configBasic.py --output-bucket gs://... [step flags]
"""

from types import SimpleNamespace
from gnocchi_chrX_nonPAR_utils import (
    METHYL_THRESHOLDS_LOGIT12,
    run_pipeline,
    parse_args,
)

CONFIG_BASIC = SimpleNamespace(
    name="Basic",
    # Training filters (steps 3b/3c): no coverage or region filters
    use_external_coverage=False,
    coverage_mean_range=None,
    min_over_10=None,
    exclude_blacklist_lcr_segdup=False,
    methyl_thresholds=METHYL_THRESHOLDS_LOGIT12,
    filter_downsampled=False,
    training_median_approx_range=None,
    # Mu/element filters (steps 4–6): coverage_mean 23–24, blacklist_gap2 only
    mu_cfg=SimpleNamespace(
        use_external_coverage=False,
        coverage_mean_range=(23.0, 24.0),
        min_over_10=None,
        exclude_blacklist_lcr_segdup=False,
        filter_downsampled=False,
    )
)

if __name__ == "__main__":
    args = parse_args(description="Gnocchi chrX nonPAR — Config Basic (no coverage/blacklist filters, logit12)")
    run_pipeline(args, CONFIG_BASIC)