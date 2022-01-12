import pandas as pd
import os

from io_utils import RUNS_DIR, EVAL_DIR


def load_solutions(ns):
    dfs = []
    for n in ns:
        # try:
        #     temp1 = pd.read_csv(f"{EVAL_DIR}/{n}_rf.csv")
        #     temp1["approach"] = "rf"
        # except FileNotFoundError:
        #     temp1 = pd.DataFrame()
        #     print("Didn't find rf solutions")
        #
        # try:
        #     temp2 = pd.read_csv(
        #         f"{EVAL_DIR}/{n}_single_median_average_quartile_quantile.csv"
        #     )
        #     temp2["approach"] = "ga"
        # except FileNotFoundError:
        #     temp2 = pd.DataFrame()
        #     print("Didn't find ga solutions")

        try:
            temp3 = pd.read_csv(
                f"{EVAL_DIR}/{n}_single_log_median_log_average_log_quartile_log_quantile_log.csv"  # TODO
            )
            temp3["approach"] = "new"
        except FileNotFoundError:
            temp3 = pd.DataFrame()
            print("Didn't find new solutions")
        # df = pd.concat([temp1, temp2, temp3]).sort_values("single", ascending=False).reset_index(drop=True)
        df = temp3.sort_values("single_log", ascending=False).iloc[:, 1:].reset_index(drop=True)
        dfs.append(df)
    return dfs


class Results:
    def __init__(self, scoring):
        self.scoring = scoring

    def _clean_and_eval(self, dfs, num_genes, runs):
        for df, run in zip(dfs, runs):
            df['run_id'] = run.split('/')[-1]
        conc = pd.concat(dfs)
        # print(conc)
        sols = conc[[str(i) for i in range(num_genes)]]
        sols = sols.drop_duplicates()
        metrics = sols.apply(lambda x: self.scoring.score(x, log=True), axis=1)
        # metrics_log = sols.apply(lambda x: self.scoring.score(x, log=True), axis=1)
        sols = pd.concat([sols, metrics], axis=1)
        sols = sols.join(conc['run_id'])
        print(list(sols.columns))
        sols = sols.sort_values(by="single_log", ascending=False)
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

        sols = self._clean_and_eval(dfs, num_genes, runs)
        sols.to_csv(os.path.join(EVAL_DIR, f'{num_genes}_{"_".join(metrics)}.csv'))
        return sols
