import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

from scoring import Scoring


def _y_minor_ticks(plot):
    plot.tick_params(axis="y", which="minor", bottom=True)
    plot.minorticks_on()
    plot.xaxis.set_tick_params(which="minor", bottom=False)


def _prettify_axes(plot):
    sns.set_style("ticks")
    plot.set_xticklabels(plot.get_xticklabels(), rotation=40, ha="right")
    plot.get_figure().set_size_inches(14, 8)
    _y_minor_ticks(plot)
    sns.despine()
    plt.tight_layout()


class Plotting:
    def __init__(self, scoring: Scoring):
        self.df_unscaled = scoring.df
        self.tissues = list(self.df_unscaled.tissue.unique())
        self.scoring = scoring

    def boxplot(self, selection: [str], log: bool, path: str = None):
        plt.clf()
        self.df_unscaled['tpm_sum'] = self.df_unscaled[selection].sum(axis=1)
        self.df_unscaled['tpm_sum_log'] = np.log2(self.df_unscaled['tpm_sum'])
        pal = sns.color_palette()
        pal = [pal[0], pal[2], pal[1]]
        key = "tpm_sum_log" if log else "tpm_sum"
        plot = sns.boxplot(
            data=self.df_unscaled,
            x="tissue",
            y=key,
            hue="group",
            linewidth=1,
            flierprops=dict(markersize=2),
            dodge=False,
            palette=pal
        )

        plt.title(f"{len(selection)} selected antigens: {', '.join(selection)}"
                  if len(selection) <= 20 else f"{len(selection)} selected antigens")
        plt.ylabel("log2(TPM sum)" if log else "TPM sum")

        _prettify_axes(plot)
        self._plot_summary_stats(selection, log)

        if path:
            plot.get_figure().savefig(path, format="svg", bbox_inches="tight")

        return plot

    def _plot_summary_stats(self, selection: [str], log: bool):
        scores = self.scoring.score(selection, log=log)
        summary = self.scoring.summary(log)
        key1, key2, key3 = "single", "iqr", "median"
        if log:
            key1, key2, key3 = key1 + "_log", key2 + "_log", key3 + "_log"

        x = self._plot_lg(summary, scores, key1, "max", "min", 1)
        x = self._plot_lg(summary, scores, key2, "outlier_upper", "outlier_lower", 0.50) # TODO
        self._plot_lg(summary, scores, key3, "median", "median", 0.3, x)
        # self._plot_lg(summary, scores, "quartile_log", "3rd_quartile", "1st_quartile", 0.5, x)

    lw = 0.7


    def _l_hline(self, y, c, a):
        plt.axhline(y, color=c, alpha=a, linewidth=self.lw)

    def _plot_lg(self, sm, scores, ln, hk, sk, alpha, prev_x=None):
        sm = sm.reset_index()
        sm_h = sm.loc[~sm.disease][hk].max()
        sm_s = sm.loc[sm.disease][sk].min()

        col = "red" if scores[ln] < 0 else "green"
        self._l_hline(sm_h, col, alpha)
        self._l_hline(sm_s, col, alpha)
        grp = sm[sm[hk] == sm_h].tissue.iloc[0]

        x = self.tissues.index(grp)
        # print(x)

        ymin = min(sm_h, sm_s)
        ymax = max(sm_h, sm_s)
        plt.vlines(
            x=x,
            ymin=ymin,
            ymax=ymax,
            alpha=alpha,
            color=col,
            linewidth=self.lw,
        )
        if prev_x and prev_x == x:
            rotation = -90
            offset = 0.1
        else:
            rotation = 90
            offset = -0.6
        plt.text(
            x + offset,
            (ymin + ymax) / 2,
            "{:.2f}".format(scores[ln]),
            rotation=rotation,
        )
        return x
