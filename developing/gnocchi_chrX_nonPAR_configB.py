"""
Config B — chrX nonPAR gnocchi pipeline
Matches Siwei's final/updated result, explicitly cited in her step-7 command
(prefilter_1kb_ht_x.py → calculate_genome_context_vt_count_by_element_x.py path).

Parameters:
  Coverage:      over_10 >= 0.95 from external gnomAD coverage HT
  Blacklist:     blacklist_gap2.bed + lcr.bed + segdup.bed
  Methyl levels: logit13 (top bin >0.85, levels 0–13)
  Downsampled:   same spatial filters applied to mu estimation

Usage:
  python gnocchi_chrX_nonPAR_configB.py --output-bucket gs://... [step flags]

Step flags (can combine multiple):
  --prefilter          Step 1:  filter context/genome to chrX nonPAR intervals
  --downsampling       Step 2:  filter to downsampled 1000 AC cutoff
  --annotate-methyl    Step 3a: annotate CpG loci with raw per-stage methyl values
  --annotate-obs       Step 3b: apply Config B filters; annotate obs (0/1)
  --compute-coeff      Step 3c: fit logistic regression → logit13 level HT
  --compute-mu         Step 4:  count obs/possible, estimate mu, fit po~mu
  --compute-element-z  Steps 5+6: annotate fitted_po/element_id per locus; compute per-element oe and Z scores

Note: Steps 1 and 2 (--prefilter, --downsampling) write to the same output paths
as Config A. If you already ran those steps with Config A, you can skip them here.
"""

from types import SimpleNamespace
from gnocchi_chrX_nonPAR_utils import (
    METHYL_THRESHOLDS_LOGIT13,
    run_pipeline,
    parse_args,
)

CONFIG_B = SimpleNamespace(
    name="B",
    # Coverage: use over_10 field from external gnomAD coverage HT
    use_external_coverage=True,
    coverage_mean_range=None,
    min_over_10=0.95,
    # Region exclusion: blacklist_gap2 + LCR + segdup (merged BED)
    exclude_blacklist_lcr_segdup=True,
    # Methylation discretization: logit13 (top threshold >0.85, 13 levels)
    methyl_thresholds=METHYL_THRESHOLDS_LOGIT13,
    # Downsampled mu: apply same spatial filters as the element step
    filter_downsampled=True,
    # Training-only extra filter: median_approx 22-24 from external coverage HT
    # (matches get_met_14_vs_obs_x.py active state; NOT applied to element or mu steps)
    training_median_approx_range=(22, 24),
)

if __name__ == "__main__":
    args = parse_args(description="Gnocchi chrX nonPAR pipeline — Config B (logit13, over_10>=0.95, gap2+lcr+segdup)")
    run_pipeline(args, CONFIG_B)