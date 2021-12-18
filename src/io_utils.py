import os
import re
import pandas as pd

RUNS_DIR = os.path.join("..", "runs")
EVAL_DIR = os.path.join("..", "eval")


def create_run_dir(meth, n_genes):
    runs = next(os.walk(RUNS_DIR))[1]
    prev_run_id = int(max(map(lambda x: int(re.search(r"\d+", x).group()), runs)))
    run_id = prev_run_id + 1 if runs else 0
    run_dir = os.path.join(RUNS_DIR, f"{run_id}_{meth}_{n_genes}")
    os.mkdir(run_dir)
    return run_dir


def sort_genes(sols):
    vals = sols.values
    vals.sort(axis=1)
    sols = pd.DataFrame(vals, sols.index, sols.columns).drop_duplicates()
    return sols
