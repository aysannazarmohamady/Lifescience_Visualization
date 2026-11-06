"""Score results_*.jsonl: detection (flawed vs clean) and, in closed mode, identification.
  python3 Phase1_Benchmark/eval/score_eval.py Phase1_Benchmark/eval/results_*.jsonl
DRAFT NOTE: drafted with LLM assistance."""
import json, sys, collections
for fn in sys.argv[1:]:
    rs = [json.loads(l) for l in open(fn)]
    flawed = [r for r in rs if r["truth"] != "NONE"]; clean = [r for r in rs if r["truth"] == "NONE"]
    rec = sum(r["has_error"] is True for r in flawed) / max(1, len(flawed))
    fpr = sum(r["has_error"] is True for r in clean) / max(1, len(clean))
    ident = sum(r["error_id"] == r["truth"] for r in flawed) / max(1, len(flawed))
    pf = sum(r["error_id"] == "PARSE_FAIL" for r in rs)
    print(f"{fn}\n  n={len(rs)} (flawed {len(flawed)}, clean {len(clean)})  detection recall={rec:.2f}  "
          f"false-alarm on clean={fpr:.2f}  exact error identified={ident:.2f}  parse failures={pf}")
    per = collections.defaultdict(list)
    for r in flawed: per[r["plot"]].append(r["error_id"] == r["truth"])
    for p, v in per.items(): print(f"    {p}: {sum(v)}/{len(v)} identified")
