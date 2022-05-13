#!/usr/bin/bash
conda activate hmmer
jackhmmer -o venom_alignments_hmmer.out --cpu 16 --tblout venom_alignments_hmmer.tab ../genome_annotation/gene_annotation/augustus.hints.aa \
databases/heloderma_and_viper_venom_proteins.fasta
conda deactivate
./read_hmmer_results.py -i venom_alignments_hmmer.tab -a databases/heloderma_and_viper_venom_proteins.tab
