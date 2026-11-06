"""Verify the variant-2 swimmer images against ground truth computed from the
CLEAN ADSL/ADRS: the clean image must pass every check and each flawed image
must fail exactly the check of its own error. Run from the repository root
after run_swimmer_v2.R. DRAFT NOTE: drafted with LLM assistance."""
import hashlib
import pandas as pd

OUT = "Outputs/swimmer_v2/"
adsl = pd.read_csv("Data/V1/ADSL.csv"); adrs = pd.read_csv("Data/V1/ADRS.csv")
trt = adsl[adsl.ARM == "TREATMENT"].copy()
trt["bar"] = (pd.to_datetime(trt.TRTEDT) - pd.to_datetime(trt.TRTSDT)).dt.days
true_bar = trt.set_index("USUBJID").bar.to_dict()
true_ong = set(trt.USUBJID[trt.EOSSTT == "ONGOING"])
start = pd.to_datetime(trt.set_index("USUBJID").TRTSDT).to_dict()
ov = adrs[(adrs.PARAMCD == "OVRLRESP") & adrs.AVALC.isin(["CR", "PR", "SD", "PD"]) & adrs.ADT.notna()]
true_day = {}
for r in ov.itertuples():          # earliest assessment per (patient, response)
    if r.USUBJID in start:
        d = (pd.to_datetime(r.ADT) - start[r.USUBJID]).days
        k = (r.USUBJID, r.AVALC); true_day[k] = min(true_day.get(k, d), d)
TUMORS = ["NSCLC", "BRCA", "HCC", "CRC", "PDAC"]

def checks(tag):
    rows = pd.read_csv(f"{OUT}{tag}_rows.csv"); mk = pd.read_csv(f"{OUT}{tag}_markers.csv")
    ong = set(l.strip() for l in open(f"{OUT}{tag}_ongoing.txt") if l.strip())
    marker_ok = all(abs(true_day[(r.USUBJID, r.resp)] - r.day) < 1e-6 for r in mk.itertuples())
    bar_ok = all(abs(true_bar[r.USUBJID] - r.bar_days) < 1e-6 for r in rows.itertuples())
    arrow_ok = true_ong.issubset(ong)
    ordered = True
    for t in TUMORS:
        s = rows[rows.tumor == t].bar_days.tolist()
        ordered &= s == sorted(s, reverse=True)
    return {"SW01 marker at assessment date": marker_ok, "SW02 bar length = dates": bar_ok,
            "SW03 ongoing arrows present": arrow_ok, "SW04 rows sorted": ordered}

md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest()
clean = checks("SW_v2_clean"); print("CLEAN:", clean)
assert all(clean.values()), "clean image fails a checker"
exp = {"SW01": "SW01 marker at assessment date", "SW02": "SW02 bar length = dates",
       "SW03": "SW03 ongoing arrows present", "SW04": "SW04 rows sorted"}
ok = True
for e, name in exp.items():
    res = checks(f"SW_v2_{e}"); failed = [k for k, v in res.items() if not v]
    diff = md5(f"{OUT}SW_v2_{e}.png") != md5(f"{OUT}SW_v2_clean.png")
    good = failed == [name] and diff; ok &= good
    print(f"{e}: fails={failed} image_differs_from_clean={diff} -> {'OK' if good else 'PROBLEM'}")
print("ALL OK" if ok else "SOME PROBLEM")
