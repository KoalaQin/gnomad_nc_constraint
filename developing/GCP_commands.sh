hailctl dataproc start qh1 \
  --master-machine-type n1-highmem-16 \
  --zone us-central1-c \
  --requester-pays-allow-all \
  --pkgs="git+https://github.com/koalaqin/gnomad_lof.git@main,\
          git+https://github.com/broadinstitute/gnomad_methods.git@main,\
          git+https://github.com/broadinstitute/gnomad_qc.git@main" \
  --autoscaling-policy=max-10 \
  --labels owner=qin,workload=gnocchi \
  --max-idle 60m

gcloud compute scp /Users/qinhe/PycharmProjects/gnomad_lof/constraint_utils/*.py heqin27@qh1-m:/home/heqin27/

hailctl dataproc submit qh1 developing/gnocchi_autosome_par.py --\
    --region chrX_PAR \
    --output-bucket gs://qin-gnocchi/gnocchi_files \
    --output-suffix chrX_par \
    --af-cutoff 0.001

hailctl dataproc submit qh3 developing/gnocchi_autosome_par.py \
    --output-bucket gs://qin-gnocchi/gnocchi_files \
    --output-suffix chr1-22 \
    --af-cutoff 0.001

hailctl dataproc submit qh1 gnocchi_chrX_nonPAR_configFinal.py \
    --output-bucket gs://qin-gnocchi/gnocchi_files \
    --annotate-obs --compute-coeff --compute-mu --compute-element-z \
    --pyfiles /Users/qinhe/PyCharmMiscProject/gnomad_nc_constraint/developing/gnocchi_chrX_nonPAR_utils.py \
    --overwrite

hailctl dataproc submit qh2 gnocchi_chrX_nonPAR_configFinal.py \
    --output-bucket gs://qin-gnocchi/gnocchi_files \
    --compute-element-z \
    --pyfiles /Users/qinhe/PyCharmMiscProject/gnomad_nc_constraint/developing/gnocchi_chrX_nonPAR_utils.py \
    --overwrite

hailctl dataproc submit qh1 developing/postprocess_gnocchi_outputs.py