"""Build the self-contained evaluation package (60 images, 12 plot types) from rendered outputs.
Run from the repo root after:  bash Phase1_Benchmark/new7/run_all.sh  and the five
Phase1_Benchmark/variant2/run_*_v2.R runners.  Writes Outputs/oncviz_eval/ and checks it.
Images are flattened onto a white background (some renderers write transparent PNGs).
DRAFT NOTE: drafted with LLM assistance."""
import csv, hashlib, os, shutil
from PIL import Image

OUT = "Outputs/oncviz_eval"
shutil.rmtree(OUT, ignore_errors=True); os.makedirs(f"{OUT}/images")
NEW7 = {"Enrollment over time": ("EN", "enrollment_bench"), "Spider plot": ("SP", "spider_bench"),
        "Volcano": ("VP", "volcano_bench"), "AE bar chart": ("TX", "ae_bench"),
        "OncoPrint": ("OP", "onco_bench"), "PK concentration-time": ("PK", "pk_bench"),
        "Cumulative incidence": ("CI", "cif_bench")}
OLD5 = {"kaplan_meier": ("Kaplan-Meier curve", "KM", "km_v2"), "forest": ("Forest plot", "FP", "forest_v2"),
        "waterfall": ("Waterfall plot", "WF", "waterfall_v2"), "swimmer": ("Swimmer plot", "SW", "swimmer_v2"),
        "consort": ("CONSORT diagram", "CN", "consort_v2")}
tax, man = [], []

def add(src, plot, truth):
    m = Image.open(src)
    if "A" in m.getbands():
        bg = Image.new("RGB", m.size, "white"); bg.paste(m, mask=m.getchannel("A")); m = bg
    m = m.convert("RGB"); name = os.path.basename(src); m.save(f"{OUT}/images/{name}")
    man.append({"image": f"images/{name}", "plot": plot, "truth": truth})

for r in csv.DictReader(open("Phase1_Benchmark/new7/taxonomy_new7.csv")):
    tax.append({"plot_type": r["plot_type"], "error_id": r["error_id"], "description": r["description"]})
for plot, (pre, d) in NEW7.items():
    add(f"Outputs/{d}/{pre}_clean.png", plot, "NONE")
    for t in tax:
        if t["plot_type"] == plot: add(f"Outputs/{d}/{pre}_{t['error_id']}.png", plot, t["error_id"])
for t in csv.DictReader(open("Phase1_Benchmark/taxonomy/taxonomy.csv")):
    plot, pre, d = OLD5[t["plot_type"]]
    tax.append({"plot_type": plot, "error_id": t["error_id"], "description": t["description"]})
for key, (plot, pre, d) in OLD5.items():
    add(f"Outputs/{d}/{pre}_v2_clean.png", plot, "NONE")
    for t in tax:
        if t["plot_type"] == plot: add(f"Outputs/{d}/{pre}_v2_{t['error_id']}.png", plot, t["error_id"])

for fn, rows, cols in [("manifest.csv", man, ["image", "plot", "truth"]),
                       ("taxonomy.csv", tax, ["plot_type", "error_id", "description"])]:
    with open(f"{OUT}/{fn}", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
shutil.copy("Phase1_Benchmark/eval/run_eval.py", f"{OUT}/run_eval.py")
shutil.copy("Phase1_Benchmark/eval/PACKAGE_README.txt", f"{OUT}/README.txt")

# ---- consistency checks ----
plots = sorted({r["plot"] for r in man}); assert len(plots) == 12, plots
for p in plots:
    rs = [r for r in man if r["plot"] == p]
    ids = sorted(t["error_id"] for t in tax if t["plot_type"] == p)
    assert len(rs) == 5 and sum(r["truth"] == "NONE" for r in rs) == 1, p
    assert sorted(r["truth"] for r in rs if r["truth"] != "NONE") == ids and len(ids) == 4, p
md5 = {r["image"]: hashlib.md5(open(f"{OUT}/{r['image']}", "rb").read()).hexdigest() for r in man}
assert len(set(md5.values())) == len(man), "duplicate images"
for r in man:
    im = Image.open(f"{OUT}/{r['image']}"); assert im.mode == "RGB" and max(im.size) <= 8000
    assert os.path.getsize(f"{OUT}/{r['image']}") < 5_000_000
print(f"OK: {len(man)} images, {len(tax)} errors, {len(plots)} plot types, all RGB, unique, <5 MB")
