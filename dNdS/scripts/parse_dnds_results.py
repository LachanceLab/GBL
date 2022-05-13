#!/usr/bin/env python
import sys
import argparse
import pandas as pd
from scipy.stats import chi2
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import os
from collections import defaultdict


def parse_codeml_file(input_file):
    """
    @param input_file: str, codeml output file
    @return: float, float, dN/dS and likelihood of model
    """
    # open file
    with open(input_file, "r") as infile:
        for line in infile:
            # line with dN/dS ratio if branch model (one ratio for each branch) was used
            if line.startswith("w ratios as labels for TreeView:"):
                results_line = infile.readline().strip()
                # GBL
                dnds_gbl = float(results_line.strip().split(',')[0].split('#')[1])
                #glass lizard
                dnds_gl = float(results_line.strip().split(',')[1].split('#')[1].split(' )')[0])
                # komodo dragon
                dnds_komodo = float(results_line.strip().split(',')[2].split('#')[1].split(' )')[0])
                # (GBL, GL) - Komodo
                dnds_gbl_and_gl_komodo = float(results_line.strip().split(',')[1].split('#')[-1])
                # (GBL, GL, Komodo) - anole)
                dnds_gbl_and_gl_and_komodo_anole = float(results_line.strip().split(',')[2].split('#')[-1])
                # Anole
                dnds_anole = float(results_line.strip().split(',')[3].split('#')[1].split(')')[0])
                break
            # line with dN/dS ratio if simple model (one rate) was used
            elif line.startswith("omega (dN/dS) =  "):
                dnds_gbl = float(line.strip().split('omega (dN/dS) =  ')[-1])
                dnds_gl = float(line.strip().split('omega (dN/dS) =  ')[-1])
                dnds_komodo = float(line.strip().split('omega (dN/dS) =  ')[-1])
                dnds_anole = float(line.strip().split('omega (dN/dS) =  ')[-1])
                dnds_gbl_and_gl_komodo = float(line.strip().split('omega (dN/dS) =  ')[-1])
                dnds_gbl_and_gl_and_komodo_anole = float(line.strip().split('omega (dN/dS) =  ')[-1])
                break
            # line with likelihood
            elif line.startswith('lnL'):
                likelihood = float(line.split('):')[1].split("+")[0])
    infile.close()
    return dnds_gbl, dnds_gl, dnds_komodo, dnds_anole, dnds_gbl_and_gl_komodo, dnds_gbl_and_gl_and_komodo_anole,\
           likelihood


def get_single_copy_orthologs(input_file):
    """
    Extract single copy orthologs from OrthoFinder results (N0.tsv file)
    @params input_file: str, path to N0.tsv file generated by OrthoFinder
    @return: pd.DataFrame, with gene ids of single copy orthologs for each species
    """
    orthologs = pd.read_csv(input_file, sep='\t', header=0)
    orthologs.dropna(inplace=True)
    # get species names
    species = [c.split('_cds_aa')[0] for c in orthologs.columns if "cds_aa" in c]
    # identify single copy orthologs and extract
    single_copy_orthologs = orthologs[orthologs.applymap(lambda x: len(x.split(', '))).iloc[:, -len(species):].sum(axis=1) == len(species)]
    # set index
    single_copy_orthologs.set_index('HOG', drop=True, inplace=True)
    # only keep columns with gene ids for each species
    single_copy_orthologs = single_copy_orthologs.iloc[:, -len(species):]
    single_copy_orthologs.columns = species
    return single_copy_orthologs


def parse_results(input_files):
    """
    Helper function to parse multiple codeml output files
    @param input_files: list, codeml result files for different HOGs
    @return: pd.DataFrame, dNdS ratio and likelihood for each HOG
    """
    hog_dnds = defaultdict(list)
    for infile in input_files:
        hog = infile.split('/')[-1].split('_')[0]
        # parse file
        dnds_gbl, dnds_gl, dnds_komodo, dnds_anole, dnds_gbl_and_gl_komodo, \
        dnds_gbl_and_gl_and_komodo_anole, likelihood = parse_codeml_file(infile)
        # save data
        hog_dnds['hog'].append(hog)
        hog_dnds['dnds_gbl'].append(dnds_gbl)
        hog_dnds['dnds_gl'].append(dnds_gl)
        hog_dnds['dnds_komodo'].append(dnds_komodo)
        hog_dnds['dnds_anole'].append(dnds_anole)
        hog_dnds['dnds_(gbl_gl)_komodo'].append(dnds_gbl_and_gl_komodo)
        hog_dnds['dnds_(gbl_gl_komodo)_anole'].append(dnds_gbl_and_gl_and_komodo_anole)
        hog_dnds['likelihood'].append(likelihood)
    # create df
    dnds_df = pd.DataFrame.from_dict(hog_dnds).set_index('hog')
    return dnds_df


def plot_dnds(dnds_bg, dnds_selection, labels=['Background', 'Evolving at different rate'],
              filename='histogram_dnds_bg_and_selection.png'):
    fig, ax = plt.subplots()
    bins = np.arange(0, 1.35, 0.05).tolist()
    hist_bg, bin_edges_bg = np.histogram(dnds_bg, bins=bins, weights=np.repeat(1 / dnds_bg.shape[0], dnds_bg.shape[0]))
    hist_sel, bin_edges_sel = np.histogram(dnds_selection, bins=bins,
                                           weights=np.repeat(1 / dnds_selection.shape[0], dnds_selection.shape[0]))
    ax.bar(bin_edges_bg[:-1] + 0.02, hist_bg, label=labels[0], width=0.02)
    ax.bar(bin_edges_sel[:-1] + 0.04, hist_sel, label=labels[1], width=0.02)
    ax.set_xticks(np.arange(0.03, 1.33, 0.05).tolist())
    tick_labels = [f"{round(l, 2)} - {round(u, 2)}" for l, u in zip(np.arange(0, 1.25, 0.05), np.arange(0.05, 1.3, 0.05))]
    tick_labels.append('>1.25')
    ax.set_xticklabels(tick_labels,
                       rotation=90)
    ax.set_xlabel('dN/dS')
    ax.set_ylabel('Density')
    ax.legend(bbox_to_anchor=(0.5, -.26), loc='upper center', ncol=2)
    fig.savefig(filename, bbox_inches='tight')


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--branch_model', nargs='+', help='Result files generated by codeml using branch model',
                        required=True)
    parser.add_argument('-s', '--simple_model', nargs='+', help='Result files generated by codeml using simple model',
                        required=True)
    parser.add_argument('--orthologs', help='N0.tsv file from OrthoFinder', required=True)
    parser.add_argument('-a', '--alpha', help='Significance level [0.05]', default=0.05, type=float)
    parser.add_argument('-k', '--degrees_of_freedom', help='Degrees of freedom of chi2 function [2]', default=2,
                        type=int)
    parser.add_argument('--dnds_thresh', help='Extract genes whose dN/dS is in bottom and top dnds_thresh percentiles. [10]',
                        default=10, type=float)
    parser.add_argument('--anole_id_gene_mapping', help='ID to gene symbol mapping anole. [anole_id_gene_mapping.tab]',
                        default='data/anole_id_gene_mapping.tab')
    # parser.add_argument('--output_selection', help='Summary output file (Only genes that pass likelihood ratio test and '
    #                                                'dNdS > dNdS thresh', required=True)
    parser.add_argument('--output_diff_rate', help='Summary output file (Only genes that pass likelihood ratio test',
                        required=True)
    parser.add_argument('--output_complete', help='Output file (All single-copy orthologs', required=True)
    parser.add_argument('--gene_list_background', help='Output file containing anole gene names of background',
                        required=True)
    parser.add_argument('--gene_list_diff_rate',
                        help='Output file containing anole gene names of genes evolving at different rate',
                        required=True)
    parser.add_argument('--gene_list_top',
                        help='Output file containing anole gene names of genes possibly whose dN/dS ratio on '
                             'the branch leading to the GBL is among the top dnds_thresh% ', required=True)

    parser.add_argument('--gene_list_bottom',
                        help='Output file containing anole gene names of genes possibly whose dN/dS ratio on '
                             'the branch leading to the GBL is among the bottom dnds_thresh% ', required=True)

    args = parser.parse_args()
    id_gene_mapping = pd.read_csv(args.anole_id_gene_mapping, header=None, sep='\t', names=['transcript_id',
                                                                                            'official_gene_id'])
    branch_model = args.branch_model
    simple_model = args.simple_model
    # parse results
    dnds_selection = parse_results(branch_model)
    dnds_neutral = parse_results(simple_model)
    # join dataframes
    dnds_df = dnds_selection.join(dnds_neutral, lsuffix='_selection', rsuffix='_neutral')
    # compute likelihood ratio
    dnds_df['lrt'] = 2 * (dnds_df.loc[:, 'likelihood_selection'] - dnds_df.loc[:, 'likelihood_neutral'])
    single_copy_orthologs = get_single_copy_orthologs(args.orthologs)
    # rename columns
    single_copy_orthologs.columns = ['transcript_id_' + col for col in single_copy_orthologs.columns.values]
    # join df with single copy orthologs to get original gene ids for each species
    df = dnds_df.join(single_copy_orthologs)
    df = df.set_index('transcript_id_anole').join(id_gene_mapping.set_index('transcript_id'))
    df.loc[:, 'official_gene_id'].to_csv(args.gene_list_background, header=False, index=False)
    df.to_csv(args.output_complete, sep='\t', header=True, index=True)

    # bonferroni correct significance level
    alpha = args.alpha / dnds_df.shape[0]
    # test if likelihood ratio is greater than critical chi2 statistic
    critical_chi2 = chi2.ppf(1 - alpha, df=args.degrees_of_freedom)
    # if yes --> gene evolving at different rate
    df = df[df.lrt > critical_chi2]
    df.to_csv(args.output_diff_rate, sep='\t', header=True, index=True)

    df.loc[:, 'official_gene_id'].to_csv(args.gene_list_diff_rate, header=False, index=False)
    # get top percentiles: if dnds > thresh possibly less conserved
    print('dN/dS cutoff for top {:.0f}%: {:.3f}'.format(args.dnds_thresh,
                                                        np.percentile(df.dnds_gbl_selection.values, 100 - args.dnds_thresh)))
    df_top = df[df.dnds_gbl_selection > np.percentile(df.dnds_gbl_selection.values, 100 - args.dnds_thresh)].copy()

    # sort by dNdS ratio
    df_top.sort_values("dnds_gbl_selection", inplace=True)
    # save
    # df_top.to_csv(args.output_selection, sep='\t', header=True, index=True)
    df_top.loc[:, 'official_gene_id'].to_csv(args.gene_list_top, header=False, index=False)

    # get bottom percentiles: if dnds < thresh possibly conserved
    print('dN/dS cutoff for bottom {:.0f}%: {:.3f}'.format(args.dnds_thresh,
                                                        np.percentile(df.dnds_gbl_selection.values, args.dnds_thresh)))
    df_bottom = df[df.dnds_gbl_selection < np.percentile(df.dnds_gbl_selection.values, args.dnds_thresh)].copy()

    # sort by dNdS ratio
    df_bottom.sort_values("dnds_gbl_selection", inplace=True)
    # save
    # df_bottom.to_csv(args.output_selection, sep='\t', header=True, index=True)
    df_bottom.loc[:, 'official_gene_id'].to_csv(args.gene_list_bottom, header=False, index=False)


if __name__ == '__main__':
    main(sys.argv[1:])