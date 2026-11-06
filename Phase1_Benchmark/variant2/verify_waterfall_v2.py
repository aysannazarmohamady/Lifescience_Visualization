"""Verify the variant-2 waterfall images: the clean image must pass every
checker and each flawed image must fail exactly the checker of its own error.
Run from the repository root after run_waterfall_v2.R.
DRAFT NOTE: drafted with LLM assistance."""
import hashlib, sys
import pandas as pd
sys.path.insert(0, "Phase1_Benchmark/checkers")
import checkers as c

OUT = "Outputs/waterfall_v2/"

def category_ok(order, pr_th=-30, pd_th=20):
    exp = order.pct.apply(lambda p: "CR" if p <= -100 else "PR" if p <= pr_th
                          else "PD" if p >= pd_th else "SD")
    return bool((order.resp == exp).all())

def run_checks(tag):
    order = pd.read_csv(f"{OUT}{tag}_order.csv")
    hl = [float(x) for x in open(f"{OUT}{tag}_hlines.txt").read().split()]
    return {
        "WF01/WF04 colour matches % change": category_ok(order),
        "WF02 PD line present": c.check_waterfall_pd_line_present(hl),
        "WF03 bars sorted": c.check_waterfall_bars_sorted(order.pct.tolist()),
    }

def md5(f):
    return hashlib.md5(open(f, "rb").read()).hexdigest()

clean = run_checks("WF_v2_clean")
print("CLEAN:", clean)
assert all(clean.values()), "clean image fails a checker"
expected_fail = {"WF01": "WF01/WF04 colour matches % change",
                 "WF04": "WF01/WF04 colour matches % change",
                 "WF02": "WF02 PD line present",
                 "WF03": "WF03 bars sorted"}
ok = True
for err, name in expected_fail.items():
    res = run_checks(f"WF_v2_{err}")
    failed = [k for k, v in res.items() if not v]
    differs = md5(f"{OUT}WF_v2_{err}.png") != md5(f"{OUT}WF_v2_clean.png")
    good = failed == [name] and differs
    ok &= good
    print(f"{err}: fails={failed} image_differs_from_clean={differs} -> {'OK' if good else 'PROBLEM'}")
print("ALL OK" if ok else "SOME PROBLEM")
