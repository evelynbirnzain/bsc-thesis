import pandas as pd

import numpy as np

linkage_map = {
    "single": lambda h, s: s.min() - h.max(),
    # "complete": lambda h, s: s.max() - h.min(),
    "median": lambda h, s: s.median() - h.median(),
    "average": lambda h, s: s.mean() - h.mean(),
    "quartile": lambda h, s: s.quantile(0.25) - h.quantile(0.75),
    "quantile": lambda h, s: s.quantile(0.05) - h.quantile(0.95),
}

linkage_map_log = {
    "single_log": lambda h, s: np.log2(s.min()) - np.log2(h.max()),
    # "complete_log": lambda h, s: s.max() - h.min(),
    "median_log": lambda h, s: np.log2(s.median()) - np.log2(h.median()),
    "average_log": lambda h, s: np.log2(s.mean()) - np.log2(h.mean()),
    "quartile_log": lambda h, s: np.log2(s.quantile(0.25)) - np.log2(h.quantile(0.75)),
    "quantileLog": lambda h, s: np.log2(s.quantile(0.05)) - np.log2(h.quantile(0.95)),
}


def _split_control(data):
    is_healthy = data["disease"] == "control"
    df_healthy = data.loc[is_healthy]
    df_sick = data.loc[~is_healthy]
    return df_healthy, df_sick


class Scoring:
    def __init__(self, df):
        self.df_control, self.df_sick = _split_control(df)
        self.h_tissues = self.df_control.tissue.unique()

    def score_abs(self, selection, key, log=False):
        hs = self.df_control[selection].sum(axis=1)
        ss = self.df_sick[selection].sum(axis=1)
        fun = linkage_map[key] if not log else linkage_map_log[key]
        return fun(hs, ss)

    def score_by_group(self, selection, key=None, log=False):
        hs = pd.DataFrame(index=self.df_control.index)
        hs['tissue'] = self.df_control.tissue
        hs["tpm_sum"] = self.df_control[selection].sum(axis=1, numeric_only=True)
        hs = hs.groupby('tissue')
        ss = self.df_sick[selection].sum(axis=1)

        if key:
            fun = linkage_map[key] if not log else linkage_map_log[key]
            return fun(hs, ss)

        df_scores = pd.DataFrame(data=None, index=self.h_tissues)
        for method, fun in linkage_map.items() if not log else linkage_map_log.items():
            df_scores[method] = fun(hs, ss)

        return df_scores

    def score(self, selection, key=None, log=False):
        if key in ['single', 'single_log', 'complete']:
            return self.score_abs(selection, key, log)
        return self.score_by_group(selection, key, log).min()
