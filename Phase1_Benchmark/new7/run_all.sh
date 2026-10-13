#!/usr/bin/env bash
# Run all seven new benchmarks (render + verify). From the repository root:
#   bash Phase1_Benchmark/new7/run_all.sh
set -e
for k in enrollment spider volcano ae onco pk cif; do
  echo "== $k"; Rscript Phase1_Benchmark/new7/run_${k}_bench.R >/dev/null 2>&1
  python3 Phase1_Benchmark/new7/verify_${k}_bench.py | tail -1
done
