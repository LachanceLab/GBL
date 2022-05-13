#!/usr/bin/bash

#conda activate snakemake
cd data

# find orthologs
snakemake -c16 --use-conda --conda-frontend conda
# orthofinder doesn't like to be executed from within snakemake
if [[ ! -f proteome_dir/orthologs/Results_1/Phylogenetic_Hierarchical_Orthogroups/N0.tsv ]]
then
    orthofinder.py -t 16 -f proteome_dir -o proteome_dir/orthologs -n 1
fi
cd ..
# test for selection
snakemake -c16 --use-conda --conda-frontend conda
#awk -F '\t' '{print $7}' summary_dnds_selection.tab | xargs -I {} grep {} data/anole_annotation.gff | awk -F 'gene=' '{print $2}' | cut -d ';' -f1 | sort | uniq > genes_under_selection.txt;
#awk -F '\t' '{print $7}' summary_dnds.tab | xargs -I {} grep {} data/anole_annotation.gff | awk -F 'gene=' '{print $2}' | cut -d ';' -f1 | sort | uniq > genes_background.txt
echo "To perform GSEA, upload 'genes_background.txt to DAVID, select 'OFFICIAL_GENE_SYMBOL' as identifier, convert the gene IDs to Entrez IDs, and use this list as background. Then upload 'genes_dnds_[top, bottom]_[1-9]+.txt' to DAVID as a gene list, and start the analysis."