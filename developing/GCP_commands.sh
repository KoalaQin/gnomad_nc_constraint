hailctl dataproc start qh1 \
  --zone us-central1-c \
  --requester-pays-allow-all \
  --pkgs="git+https://github.com/koalaqin/gnomad_lof.git@main,\
          git+https://github.com/broadinstitute/gnomad_methods.git@main,\
          git+https://github.com/broadinstitute/gnomad_qc.git" \
  --autoscaling-policy=max-10 \
  --labels owner=qin,workload=gnocchi \
  --max-idle 60m

gcloud compute scp /Users/qinhe/PycharmProjects/gnomad_lof/constraint_utils/*.py heqin27@qh1-m:/home/heqin27/


hailctl dataproc submit qh1 gnocchi_1kb_by_interval.py -- \
  --region chrX_nonPAR \
  --output-bucket gs://qin-gnocchi/gnocchi_files \
  --output-suffix genome_1kb_chrX_nonPAR \
  --dry-run

hailctl dataproc submit qh1 gnocchi_1kb_by_interval.py -- \
  --region chrX_PAR \
  --output-bucket gs://qin-gnocchi/gnocchi_files \
  --output-suffix genome_1kb_chrX_PAR

hailctl dataproc submit qh1 gnocchi_1kb_by_interval.py -- \
  --interval '[chr22:start-end]' \
  --output-bucket gs://qin-gnocchi/gnocchi_files \
  --output-suffix genome_1kb_chr22


oe_ht = hl.read_table(out_ht)

    z = calculate_z(oe_ht, oe_ht.observed, oe_ht.expected, "gnocchi")
    z_adj = calculate_z(oe_ht, oe_ht.observed, oe_ht.expected_adj, "gnocchi_adj")

    oe_ht = oe_ht.annotate(
        gnocchi=z[oe_ht.key].gnocchi,
        gnocchi_adj=z_adj[oe_ht.key].gnocchi_adj,
    )

