import pandas as pd
import os

from io_utils import RUNS_DIR, EVAL_DIR


def load_df(path):
    df = pd.read_csv(path, dtype={"disease": str})
    df = df.fillna(value={"disease": "control"})
    df = df.iloc[:, 2:]
    df = df[df.source != "hOB"]
    df = df[
        (df.tissue != "Cells - EBV-transformed lymphocytes")
        & (df.tissue != "Cells - Cultured fibroblasts")
        ]

    return df


def load_solutions():
    dfs = []
    for n in [5, 10, 20]:
        temp1 = pd.read_csv(f"{EVAL_DIR}/{n}_rf.csv")
        temp1["approach"] = "rf"
        temp2 = pd.read_csv(
            f"{EVAL_DIR}/{n}_single_median_average_quartile_quantile.csv "
        )
        temp2["approach"] = "ga"
        df = pd.concat([temp1, temp2]).iloc[:, 1:].sort_values("single", ascending=False).reset_index(drop=True)
        dfs.append(df)
    return dfs


class Results:
    def __init__(self, scoring):
        self.scoring = scoring

    def _clean_and_eval(self, dfs, num_genes):
        sols = pd.concat(dfs)
        sols = sols[[str(i) for i in range(num_genes)]]
        sols = sols.drop_duplicates()
        metrics = sols.apply(self.scoring.score, axis=1)
        sols = pd.concat([sols, metrics], axis=1)
        sols = sols.sort_values(by="single", ascending=False)
        return sols

    def collect_and_save(self, num_genes, metrics):
        runs = next(os.walk(RUNS_DIR))[1]
        runs = list(
            filter(
                lambda x: f"_{num_genes}" in x and any(metric in x for metric in metrics),
                runs,
            )
        )
        print(runs)

        runs = [os.path.join(RUNS_DIR, run_dir, "best_sols.csv") for run_dir in runs]
        dfs = []
        for run in runs:
            dfs.append(pd.read_csv(run))

        sols = self._clean_and_eval(dfs, num_genes)
        sols.to_csv(os.path.join(EVAL_DIR, f'{num_genes}_{"_".join(metrics)}.csv'))
        return sols
