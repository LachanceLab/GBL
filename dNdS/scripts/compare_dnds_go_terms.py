#!/usr/bin/env python
import pandas as pd
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from matplotlib.lines import Line2D


def split_go_terms(df):
    """
    Some entries may contain multiple GO terms
    @param df: pd.DataFrame,
    return: pd.DataFrame, GO term, dNdS
    """
    new_df = pd.DataFrame()
    n = 0
    for i in range(df.shape[0]):
        go_terms = [go.split('~')[0] for go in df.iloc[i, 1].split(',')]
        for go in go_terms:
            if not go.startswith('GO'):
                continue
            new_df.loc[n, 'GO'] = go
            new_df.loc[n, 'dnds'] = df.iloc[i, 0]
            n += 1
    return new_df


def find_overlap(gbl, gl):
    """
    Find GO terms that are present on branch to GBL and Glass lizard
    @param gbl: pd.DataFrame, GO terms and dNdS values for GBL
    @param gl: pd.DataFrame, GO terms and dNdS values for Glass lizard
    return: dict, mapping GO terms to dNdS values in background set and set that is possible under selection
    """
    map_dict = {}
    for i in range(gl.shape[0]):
        try:
            gbl.loc[gl.index[i]]
        except KeyError:
            continue
        map_dict[gl.index[i]] = {}
        map_dict[gl.index[i]]['gl'] = gl.iloc[i]
        map_dict[gl.index[i]]['gbl'] = gbl.loc[gl.index[i]]
    return map_dict


def plot_dnds(gbl, gl, go_dnds_mapping, alpha=0.05, min_genes=3, go_column='GO_BP_DIRECT'):
    """
    Create Box plot of dNdS for GO terms that have more than min_genes genes and perform Mann-Whitney U test to test
    if dNdS of genes that evolved at different rates on the branch leading to the GBL and Asian glass lizard, respectively.
    Bonferroni correction is performed.
    @param go_dnds_mapping: dict
    @param alpha: float, significance level
    @param min_genes: int, minimum number of genes to test significance
    @param go_column: str, which GO term is represented
    """
    fig, ax = plt.subplots(figsize=(9, 4.8))#1, 2, sharey=True, gridspec_kw={"width_ratios": [2, 18]})
    color_gbl = 'cornflowerblue'
    color_gl ="yellowgreen"
    pval = mannwhitneyu(gbl, gl)[1]
    print(f"Median dN/dS GBL: {np.median(gbl)}\nMedian dN/dS GL: {np.median(gl)}")
    print('Pval general dN/dS GBL vs. GL: {:.5f}'.format(pval))
    medianprops = dict(color='black')
    deviations_gbl = np.random.uniform(-0.05, .05, size=gbl.shape[0])
    deviations_gl = np.random.uniform(-0.05, .05, size=gl.shape[0])

    ax.scatter(-1 + deviations_gbl - 0.25, gbl, color=color_gbl, s=8, alpha=0.8)
    ax.scatter(-1 + deviations_gl + 0.25, gl, color=color_gl, s=8, alpha=0.8)
    ax.boxplot([gbl, gl], positions=[-1 - 0.25, -1 + 0.25], showfliers=False, medianprops=medianprops, widths=0.2)

    # ax.set_xticks([1, 2])
    # ax.set_xticklabels(['GBL', 'Glass lizard'], rotation=90)
    ax.set_ylabel('dN/dS')
    ax.set_ylim([0, 1.25])
    ax.axvline(0, 0, 1.25, color='black')
    ticklabels = ['Genome wide']
    ticks = [-1]
    pvals = [pval]
    i = 1
    for n, (go, dnds) in enumerate(go_dnds_mapping.items()):
        gbl = dnds['gbl']
        gl = dnds['gl']
        if len(gbl) >= min_genes and len(gl) >= min_genes:
            pvals.append(mannwhitneyu(gl, gbl)[1])
            ticklabels.append(go)
            ticks.append(i)
            ax.boxplot([gbl, gl], positions=[i - 0.25, i + 0.25], widths=0.2, showfliers=False,
                       medianprops=medianprops)

            for dnds in gbl:
                deviation = np.random.uniform(-0.05, 0.05)
                if dnds > 1.25:
                    dnds = 1.25
                ax.scatter(i - 0.25 + deviation, dnds, color=color_gbl, alpha=0.8, s=8)
            for dnds in gl:
                deviation = np.random.uniform(-0.05, 0.05)
                if dnds > 1.25:
                    dnds = 1.25
                ax.scatter(i + 0.25 + deviation, dnds, color=color_gl, alpha=0.8, s=8)
            i += 1
    ax.set_ylim([0, 1.25])
    ax.set_xticks(ticks)
    # n is the number of go terms
    bonferroni = alpha / n
    signficant_labels = ['*' + l if p < bonferroni else l for l, p in zip(ticklabels, pvals)]
    ax.set_xticklabels(signficant_labels, rotation=90)
    [l.set_fontweight('bold') for l, p in zip(ax.get_xticklabels(), pvals) if p < bonferroni]
    [print('{}: {:.4f}'.format(l, p)) for l, p in zip(ticklabels, pvals) if p < bonferroni]
    legend_elements = [Line2D([0], [0], marker='o', color=color_gbl, label='GBL', ls=''),
                       Line2D([0], [0], marker='o', color=color_gl, label='Glass lizard', ls='')]
    ax.legend(handles=legend_elements, ncol=2, loc='upper center', bbox_to_anchor=(0.5, -.32))
    fig.savefig(f'go_dnds_mapping_{go_column.lower()}.png', bbox_inches='tight')


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--summary_selection', help='Summary selection file [summary_dnds_diff_rate.tab]',
                        default='summary_dnds_diff_rate.tab')
    parser.add_argument('--annotations', help='GO term annotations obtained from DAVID. '
                                              '[functional_annotation_background_david.tab]',
                        default='functional_annotation_background_david.tab')
    parser.add_argument('--alpha', help='Signficance level when testing if GO terms of genes that have evolved at '
                                        'different rate have de-/increased dNdS. [0.05]',
                        type=float, default=0.05)
    parser.add_argument('--min_genes', help='Minimum number of genes in each group to perform Mann-Whitney U test. [30]',
                        default=30, type=int)
    parser.add_argument('--go_column', help='Which go column to use from DAVID [GOTERM_BP_5]',
                        default='GOTERM_BP_5')
    args = parser.parse_args()
    dnds_selection = pd.read_csv(args.summary_selection, sep='\t', index_col=0)
    annotations = pd.read_csv(args.annotations, sep='\t')
    dnds_selection = dnds_selection.set_index('gene').join(annotations.set_index('ID'))
    # glass lizard
    dnds_gl = dnds_selection.loc[:, ['dnds_gl_selection', args.go_column]]
    # guatemalan beaded lizard
    dnds_gbl = dnds_selection.loc[:, ['dnds_gbl_selection', args.go_column]]
    dnds_gl.dropna(inplace=True)
    dnds_gbl.dropna(inplace=True)
    go_dnds_mapping_gl = split_go_terms(dnds_gl)
    go_dnds_mapping_gbl = split_go_terms(dnds_gbl)
    go_terms_gl = go_dnds_mapping_gl.groupby('GO')['dnds'].apply(list)
    go_terms_gbl = go_dnds_mapping_gbl.groupby('GO')['dnds'].apply(list)
    go_dnds_mapping = find_overlap(go_terms_gbl, go_terms_gl)
    plot_dnds(dnds_gbl.dnds_gbl_selection.values, dnds_gl.dnds_gl_selection.values, go_dnds_mapping, args.alpha, args.min_genes, args.go_column)


if __name__ == '__main__':
    main(sys.argv[1:])