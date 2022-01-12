import pandas as pd

import numpy as np

_linkage_map_abs = {
    "single": lambda h, s: s.min() - h.max(),
    # "complete": lambda h, s: s.max() - h.min(),
    "median": lambda h, s: s.median() - h.median(),
    "average": lambda h, s: s.mean() - h.mean(),
    "quartile": lambda h, s: s.quantile(0.25) - h.quantile(0.75),
    "quantile": lambda h, s: s.quantile(0.05) - h.quantile(0.95),

}

_linkage_map_log = {
    "single_log": lambda h, s: np.log2(s.min()) - np.log2(h.max()),  # TODO +1?
    # "complete_log": lambda h, s: s.max() - h.min(),
    "median_log": lambda h, s: np.log2(s.median()) - np.log2(h.median()),
    "average_log": lambda h, s: np.log2(s.mean()) - np.log2(h.mean()),
    "quartile_log": lambda h, s: np.log2(s.quantile(0.25)) - np.log2(h.quantile(0.75)),
    "quantile_log": lambda h, s: np.log2(s.quantile(0.05)) - np.log2(h.quantile(0.95)),
}

log_linkages = list(_linkage_map_log.keys())
abs_linkages = list(_linkage_map_abs.keys())

_linkage_map = dict(_linkage_map_abs)
_linkage_map.update(_linkage_map_log)


def _split_control(data: pd.DataFrame):
    df_healthy = data.loc[~data["disease"]]
    df_sick = data.loc[data["disease"]]
    return df_healthy, df_sick


class Scoring:
    def __init__(self, df: pd.DataFrame):
        self.df = df
        self.df_control, self.df_sick = _split_control(df)
        self.h_tissues = self.df_control.tissue.unique()

    def summary(self, log: bool):
        df_summary = self.df.groupby(["disease", "tissue"]).agg(
            {"tpm_sum": ["min", "max", "median", "mean", ("3rd_quartile", lambda x: x.quantile(0.75)),
                         ("1st_quartile", lambda x: x.quantile(0.25))]}) \
            .tpm_sum
        return df_summary.apply(np.log2) if log else df_summary

    def score(self, selection: [str], log: bool, key: str = None):
        if key in ['single', 'single_log', 'complete', 'complete_log']:
            return self.score_without_group(selection, log, key)
        scores = self.score_by_group(selection, log, key)
        return scores.min()

    def score_without_group(self, selection: [str], log: bool, key: str):
        hs = self.df_control[selection].sum(axis=1)
        ss = self.df_sick[selection].sum(axis=1)

        if log:
            key = key + '_log'
        fun = _linkage_map[key]
        return fun(hs, ss)

    def score_by_group(self, selection: [str], log: bool, key: str = None):
        hs = pd.DataFrame(index=self.df_control.index)
        hs['tissue'] = self.df_control.tissue
        hs['tpm_sum'] = self.df_control[selection].sum(axis=1, numeric_only=True)
        hs = hs.groupby('tissue')
        ss = self.df_sick[selection].sum(axis=1)

        if key and log:
            key = key + '_log'

        if key:
            fun = _linkage_map[key]
            return fun(hs, ss)

        df_scores = pd.DataFrame(data=None, index=self.h_tissues)
        for method, fun in _linkage_map.items():
            df_scores[method] = fun(hs, ss)

        return df_scores
