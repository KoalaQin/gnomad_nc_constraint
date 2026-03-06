"""
Config Siwei14 — chrX nonPAR gnocchi pipeline.
Uses Siwei's exact hardcoded "chrX basic" logit coefficients (same as configSiwei),
but with logit14 methylation discretization (calculate_po_x.py lines 121–151):
  top bin >0.85 → 14 levels (adds >0.40→5 split vs logit13).

Spatial filters match configSiwei/configBasic (= Candidate A from Siwei's active code):
  Mu/element: coverage_mean 23–24, blacklist_gap2 only

Usage:
  python gnocchi_chrX_nonPAR_configSiwei14.py --output-bucket gs://... [step flags]

Note: --annotate-obs and --compute-coeff can be skipped (hardcoded_coeff bypasses obs step).
Run directly: --prefilter --downsampling --annotate-methyl --compute-coeff --compute-mu --compute-element-z
"""

from types import SimpleNamespace
from gnocchi_chrX_nonPAR_utils import (
    METHYL_THRESHOLDS_LOGIT14,
    run_pipeline,
    parse_args,
)

# Siwei's "chrX basic" coefficients (log-odds scale) from calculate_methyl_score_x.py lines 38-54
SIWEI_CHRX_BASIC_COEFF = {
    "Sperm_raw":   1.2142327678822575,
    "Oocyte_raw":  -0.10153058557280095,
    "PN_raw":      0.2096996493626876,
    "C2_raw":      0.08246182249615806,
    "C4_raw":      0.0803666551372927,
    "C8_raw":      0.040510820716452144,
    "Morula_raw":  0.2502153273086105,
    "ICM_raw":     0.0682508020601786,
    "PGC_7W_raw":  0.624121048631128,
    "PGC_10W_raw": 1.015123980247576,
    "PGC_11W_raw": 1.1535398138908222,
    "PGC_13W_raw": 0.24685249760735126,
    "PGC_17W_raw": 0.5309602092634527,
    "PGC_19W_raw": 0.37232732130631524,
}

CONFIG_SIWEI14 = SimpleNamespace(
    name="Siwei14",
    # Training: hardcoded coefficients — annotate-obs and logistic regression are skipped
    hardcoded_coeff=SIWEI_CHRX_BASIC_COEFF,
    use_external_coverage=False,
    coverage_mean_range=None,
    min_over_10=None,
    exclude_blacklist_lcr_segdup=False,
    methyl_thresholds=METHYL_THRESHOLDS_LOGIT14,
    filter_downsampled=False,
    training_median_approx_range=None,
    # Mu/element filters: coverage_mean 23–24, blacklist_gap2 only
    mu_cfg=SimpleNamespace(
        use_external_coverage=False,
        coverage_mean_range=(23.0, 24.0),
        min_over_10=None,
        exclude_blacklist_lcr_segdup=False,
        filter_downsampled=False,
    ),
)

if __name__ == "__main__":
    args = parse_args(description="Gnocchi chrX nonPAR — Config Siwei14 (hardcoded chrX basic coeff + logit14)")
    run_pipeline(args, CONFIG_SIWEI14)