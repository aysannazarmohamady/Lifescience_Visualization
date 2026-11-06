ONCVIZ VLM evaluation package: 60 images = 12 plot types x [1 clean + 4 flawed] (48 injected errors).
7 newer plot types + 5 original plot types (second implementation). Python 3 standard library only.
Keys are read from the environment (ANTHROPIC_API_KEY, OPENAI_API_KEY); optional
ANTHROPIC_BASE_URL / OPENAI_BASE_URL for a gateway.
  python3 run_eval.py --mock      offline self-test
  python3 run_eval.py --limit 3   cheap real trial
  python3 run_eval.py             full run (writes report.md)
  options: --claude-model <id>  --openai-model <id>
Ground truth: manifest.csv (truth column). Error definitions: taxonomy.csv.
