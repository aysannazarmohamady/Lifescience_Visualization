"""Verify the enrollment benchmark images from the data written next to each
figure. Ground truth comes from the CLEAN ADRAND/ADSL. Clean must pass every
check; each flawed image must fail exactly its own check.
DRAFT NOTE: drafted with LLM assistance."""
import hashlib, re
import pandas as pd

OUT = "Outputs/enrollment_bench/"
adrand = pd.read_csv("Data/V1/ADRAND.csv"); adsl = pd.read_csv("Data/V1/ADSL.csv")
true_arm = adsl[adsl.USUBJID.isin(adrand[adrand.RANDFL == "Y"].USUBJID)].ARM.value_counts().to_dict()
true_n = int(sum(true_arm.values()))

def run(tag):
    curve = pd.read_csv(f"{OUT}{tag}_curve.csv", parse_dates=["date"])
    meta = {l.split(" ", 1)[0]: l.split(" ", 1)[1].strip() for l in open(f"{OUT}{tag}_meta.txt")}
    closure = pd.Timestamp(meta["closure"])
    m = re.search(r"overall n=(\d+) \(Treatment (\d+) \+ Control (\d+)\)", meta["title"])
    overall, t, c = (int(x) for x in m.groups())
    return {"EN01 cumulative curve non-decreasing": bool((curve.cum.diff().dropna() >= 0).all()),
            "EN02 arm totals sum to overall": (t + c == overall) and {"TREATMENT": t, "CONTROL": c} == true_arm and overall == true_n,
            "EN03 no enrollment after closure": bool((curve.date <= closure).all()),
            "EN04 target line present": meta["target_line"] == "TRUE"}

md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest()
clean = run("EN_clean"); print("CLEAN:", clean)
assert all(clean.values()), "clean image fails a checker"
ok = True
for e, name in {"EN01": "EN01 cumulative curve non-decreasing", "EN02": "EN02 arm totals sum to overall",
                "EN03": "EN03 no enrollment after closure", "EN04": "EN04 target line present"}.items():
    res = run(f"EN_{e}"); failed = [k for k, v in res.items() if not v]
    diff = md5(f"{OUT}EN_{e}.png") != md5(f"{OUT}EN_clean.png")
    good = failed == [name] and diff; ok &= good
    print(f"{e}: fails={failed} image_differs_from_clean={diff} -> {'OK' if good else 'PROBLEM'}")
print("ALL OK" if ok else "SOME PROBLEM")
