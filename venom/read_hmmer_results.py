#!/usr/bin/env python
import sys
import argparse
import pandas as pd
from collections import defaultdict
import re
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Diamond output file (fmt 6)')
    parser.add_argument('-a', '--annotations', help='Venom annotations')
    parser.add_argument('-n', '--nr_queries', help='Number of queries. [7370]', default=7370, type=int)
    parser.add_argument('--dnds', help='Path to summary file with dNdS ratios. [../dNdS/summary_dnds.tab]',
                        default='../dNdS/summary_dnds.tab')
    parser.add_argument('--output', help='Path to write output file to.', required=False)
    args = parser.parse_args()
    if args.output is None:
        output = args.input.replace('.tab', '_dnds.tab')
    else:
        output = args.output
    columns = ['target', 't_acc', 'query','q_acc', 'evalue_full', 'score_full', 'bias_full',
               'evalue_domain', 'score_domain', 'bias_domain', 'exp', 'reg', 'clu', 'ov',
               'env', 'dom', 'rep', 'inc', 'description']

    hits = defaultdict(list)

    with open(args.input) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            hit = re.split("\s+", line.strip())
            for i, col in enumerate(columns):
                if i < len(columns) - 1:
                   hits[col].append(hit[i])
                elif i == len(columns) - 1:
                    hits[col].append(' '.join(hit[i:]))

    results = pd.DataFrame.from_dict(hits)

    results.evalue_full = results.evalue_full.astype(float)
    results.evalue_domain = results.evalue_domain.astype(float)
    results.evalue_full *= args.nr_queries
    results = results[results.evalue_full < 0.0001]
    results = results[results.evalue_domain < 0.01]

    results.sort_values(['evalue_full', 'evalue_domain'], inplace=True)
    results.target = [t.replace('tr|', 'sp|') for t in results.target.values]
    annotations = pd.read_csv(args.annotations, sep='\t')
    annotations.loc[:, 'seqid'] = 'sp|' + annotations.Entry + '|' + annotations['Entry name']
    annotations.loc[:, "Protein names"] = [pn.split(' (')[0].split(' [')[0] if not isinstance(pn, float) else pn for pn in annotations["Protein names"].values]
    merged = results.set_index("target").join(annotations.set_index('seqid'))
    # merged['Protein families'].fillna(merged['Protein names'], inplace=True)
    protein_families_grouped = merged.dropna(subset=['Protein names'], axis=0)[['query', 'Protein names']].groupby('Protein names')
    protein_family_gene_mapping = []
    for cdd, group in protein_families_grouped:
        print(f"{cdd}; CN={group['query'].unique().shape[0]}")
        protein_family_gene_mapping.append(group.set_index("query"))
    protein_family_gene_mapping = pd.concat(protein_family_gene_mapping)
    dnds = pd.read_csv(args.dnds, sep='\t')
    dnds_venomous_genes = dnds.loc[np.isin(dnds.transcript_id_gbl.values, protein_family_gene_mapping.index.values), 'dnds_gbl_selection'].values
    dnds_background = dnds.loc[~np.isin(dnds.transcript_id_gbl.values, protein_family_gene_mapping.index.values), 'dnds_gbl_selection'].values
    pval = mannwhitneyu(dnds_venomous_genes, dnds_background)[1]
    print('Pval dN/dS venomous genes vs background: {:.6f}'.format(pval))
    fig, ax = plt.subplots()
    ax.boxplot([dnds_venomous_genes, dnds_background], positions=[0, 1], showfliers=False)
    ax.set_xticklabels(['Potentially venomous genes', 'GBL background'], rotation=60)
    ax.set_ylabel('dN/dS')
    ax.set_ylim([0, 1.2])
    fig.savefig('dnds_venom_vs_non_venom.png', bbox_inches='tight')

    df = merged.set_index('query').join(dnds.loc[:, ["dnds_gbl_selection", "transcript_id_gbl"]].set_index("transcript_id_gbl"))
    df.fillna('-', inplace=True)
    df.to_csv(output, sep='\t', index=True, header=True)
    breakpoint()


if __name__ == '__main__':
    main(sys.argv[1:])