"""
Config Final — chrX nonPAR gnocchi pipeline.
Mirrors the autosome/PAR pipeline in run_nc_constraint_gnomad_v31_main.py:
  - coverage_mean filter (23–25 for chrX nonPAR, vs 30–32 for autosomes/PAR)
  - blacklist_gap2 only (no lcr/segdup exclusion)
  - same spatial filters applied uniformly to training, mu, and element steps
  - logit12 methylation discretization
  - fresh logistic regression fit

Usage:
  python gnocchi_chrX_nonPAR_configFinal.py --output-bucket gs://... [step flags]
"""

from types import SimpleNamespace
from gnocchi_chrX_nonPAR_utils import (
    METHYL_THRESHOLDS_LOGIT12,
    run_pipeline,
    parse_args,
)

CONFIG_FINAL = SimpleNamespace(
    name="Final",
    # Training filters: coverage_mean 23–25, blacklist_gap2 only
    use_external_coverage=False,
    coverage_mean_range=(23.0, 25.0),
    min_over_10=None,
    exclude_blacklist_lcr_segdup=False,
    methyl_thresholds=METHYL_THRESHOLDS_LOGIT12,
    filter_downsampled=False,
    training_median_approx_range=None,
    # No mu_cfg: mu and element steps use the same filters as training
)

if __name__ == "__main__":
    args = parse_args(description="Gnocchi chrX nonPAR — Config Final (coverage_mean 23–25, blacklist_gap2, logit12)")
    run_pipeline(args, CONFIG_FINAL)

#                 Coef.  Std.Err.          P>|z|
# const       -0.548092  0.008387   0.000000e+00
# Sperm_raw    1.206945  0.010400   0.000000e+00
# Oocyte_raw  -0.117388  0.017722   3.495737e-11
# PN_raw       0.214849  0.017027   1.673425e-36
# C2_raw       0.090357  0.015553   6.256902e-09
# C4_raw       0.090604  0.018276   7.134252e-07
# C8_raw       0.040094  0.013477   2.930844e-03
# Morula_raw   0.246798  0.018460   9.153623e-41
# ICM_raw      0.062240  0.010038   5.625071e-10
# PGC_7W_raw   0.573498  0.023537  3.962800e-131
# PGC_10W_raw  0.950499  0.049618   8.593124e-82
# PGC_11W_raw  1.122677  0.047043  7.132102e-126
# PGC_13W_raw  0.236625  0.027243   3.758510e-18
# PGC_17W_raw  0.500412  0.035258   1.014944e-45
# PGC_19W_raw  0.417720  0.030436   7.242764e-43