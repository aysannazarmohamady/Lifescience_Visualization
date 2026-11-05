"""Build manifest.csv (image, plot, truth) for the seven new plot types.
Run from the repo root after Phase1_Benchmark/new7/run_all.sh.
DRAFT NOTE: drafted with LLM assistance."""
import csv, os
PREFIX = {"Enrollment over time": ("EN", "enrollment_bench"), "Spider plot": ("SP", "spider_bench"),
          "Volcano": ("VP", "volcano_bench"), "AE bar chart": ("TX", "ae_bench"),
          "OncoPrint": ("OP", "onco_bench"), "PK concentration-time": ("PK", "pk_bench"),
          "Cumulative incidence": ("CI", "cif_bench")}
tax = list(csv.DictReader(open("Phase1_Benchmark/new7/taxonomy_new7.csv")))
rows = []
for plot, (pre, d) in PREFIX.items():
    rows.append({"image": f"Outputs/{d}/{pre}_clean.png", "plot": plot, "truth": "NONE"})
    for t in tax:
        if t["plot_type"] == plot:
            rows.append({"image": f"Outputs/{d}/{pre}_{t['error_id']}.png", "plot": plot, "truth": t["error_id"]})
missing = [r["image"] for r in rows if not os.path.exists(r["image"])]
assert not missing, missing
with open("Phase1_Benchmark/eval/manifest.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["image", "plot", "truth"]); w.writeheader(); w.writerows(rows)
print(len(rows), "images")
