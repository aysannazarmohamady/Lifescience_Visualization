#!/usr/bin/env bash
# run all 12 benchmarks (render + verify) from repo root; prints one line per plot
cd "$(dirname "$0")"
for k in enrollment spider volcano ae onco pk cif; do
  Rscript Phase1_Benchmark/new7/run_${k}_bench.R >/tmp/run_$k.log 2>&1 || echo "RENDER FAIL $k"
  echo "$k: $(python3 -I Phase1_Benchmark/new7/verify_${k}_bench.py 2>&1 | tail -1)"
done
for k in waterfall swimmer forest consort km; do
  Rscript Phase1_Benchmark/variant2/run_${k}_v2.R >/tmp/run_$k.log 2>&1 || echo "RENDER FAIL $k"
  echo "$k: $(python3 -I Phase1_Benchmark/variant2/verify_${k}_v2.py 2>&1 | tail -1)"
done
