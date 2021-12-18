import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


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


def _summary(data):
    return data.groupby(["disease", "tissue"]).agg(
        {"tpm_sum": ["min", "max", "median", "mean", ("3rd_quartile", lambda x: x.quantile(0.75))]}) \
        .tpm_sum.apply(np.log2)


class Plotting:
    def __init__(self, df_unscaled, scoring):
        self.df_unscaled = df_unscaled
        self.tissues = list(df_unscaled.tissue.unique())
        self.scoring = scoring

    def boxplot(self, selection, path=None):
        plt.clf()
        self.df_unscaled['tpm_sum'] = self.df_unscaled[selection].sum(axis=1)
        self.df_unscaled['tpm_sum_log'] = np.log2(self.df_unscaled['tpm_sum'])

        plot = sns.boxplot(
            data=self.df_unscaled,
            x="tissue",
            y="tpm_sum_log",
            hue="disease",
            linewidth=1,
            flierprops=dict(markersize=2),
            dodge=False,
        )

        plt.title(f"{len(selection)} selected antigens: {', '.join(selection)}")
        plt.ylabel("log2(TPM sum)")

        _prettify_axes(plot)
        self._plot_summary_stats(selection)
        if path:
            plot.get_figure().savefig(path, format="svg", bbox_inches="tight")

        return plot

    lw = 0.7

    def _l_hline(self, y, c, a):
        plt.axhline(y, color=c, alpha=a, linewidth=self.lw)

    def _plot_lg(self, sm, scores, ln, hk, sk, alpha, prev_x=None):
        sm_h = sm.loc['control'][hk].max()
        sm_s = sm.loc['neuroblastoma'][sk].min()

        col = "red" if scores[ln] < 0 else "green"
        self._l_hline(sm_h, col, alpha)
        self._l_hline(sm_s, col, alpha)
        grp = sm[sm[hk] == sm_h].index[0][1]
        x = self.tissues.index(grp)
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

    def _plot_summary_stats(self, selection):
        scores = self.scoring.score(selection, log=True)
        sm = _summary(self.df_unscaled)

        x = self._plot_lg(sm, scores, "single_log", "max", "min", 1)
        self._plot_lg(sm, scores, "median_log", "median", "median", 0.5, x)
