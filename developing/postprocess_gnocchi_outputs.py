"""
Post-process gnocchi output TXTs:
  - In-house copy : rename expected->expected_raw, expected_adj->expected,
                    gnocchi->gnocchi_raw, gnocchi_adj->gnocchi
  - Sharing copy  : same, but drop expected_raw and gnocchi_raw
"""
import pandas as pd

OUTPUT_BUCKET = "gs://qin-gnocchi/gnocchi_files"

FILES = {
    "chr1-22":     f"{OUTPUT_BUCKET}/gnocchi_1kb_chr1-22.txt",
    "chrX_par":    f"{OUTPUT_BUCKET}/gnocchi_1kb_chrX_par.txt",
    "chrX_nonpar": f"{OUTPUT_BUCKET}/gnocchi_1kb_chrX_nonpar_configFinal.txt",
}

RENAME = {
    "expected":     "expected_raw",
    "expected_adj": "expected",
    "gnocchi":      "gnocchi_raw",
    "gnocchi_adj":  "gnocchi",
}

PUBLIC_DROP = ["expected_raw", "gnocchi_raw"]

for label, src in FILES.items():
    print(f"Processing {label} ...")
    df = pd.read_csv(src, sep="\t")
    df = df.rename(columns=RENAME)

    # In-house: all columns with renamed fields
    inhouse_path = f"{OUTPUT_BUCKET}/gnocchi_1kb_{label}_inhouse.txt"
    df.to_csv(inhouse_path, sep="\t", index=False)
    print(f"  in-house -> {inhouse_path}")

    # public: drop raw columns
    df_share = df.drop(columns=PUBLIC_DROP)
    share_path = f"{OUTPUT_BUCKET}/gnocchi_1kb_{label}_public.txt"
    df_share.to_csv(share_path, sep="\t", index=False)
    print(f"  sharing  -> {share_path}")

print("Done.")
