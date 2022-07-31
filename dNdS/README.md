## Inferring dN/dS ratios for predicted genes in a draft genome of the Guatemalan Beaded lizard
For this analysis, gene annotations for the Asian Glass Lizard<sup>1</sup>, the Green Anole<sup>2</sup>, and the Komodo Dragon<sup>3</sup> are required, 
which will be automatically retrieved by the script.<br>

### Input files
1. Genomes and annotations in gff format of Green Anole, Asian Glass Lizard, and Komodo Dragon will be retrieved automatically into `data` subdirectory
2. Softmasked genome of Guatamalan Beaded Lizard (../genome_annotation/GBL_genome_assembly.fasta.combined.masked)
3. Annotation of protein coding regions in the Guatemalan Beaded Lizard (../genome_annotation/gene_annotation/augustus.hints.gff3)

### Output files
1. genes_background.txt: Green Anole gene IDs of identified single-copy orthologs.
2. genes_diff_rate.txt: Green Anole gene IDs of genes for which the more complex branch model was accepted, i.e., the gene evolved at a different rate on different branches in the phylogenetic tree
3. gene_dnds_bottom/top_10.txt: Green Anole gene IDs of genes whose dN/dS is among the bottom/top 10% on the branch leading to Guatemalan Beaded Lizard
4. summary_dnds.txt: Summary of dN/dS ratios of all single-copy orthologs
5. summary_dnds_diff_rate.txt: Summary of dN/dS ratios of genes that passed likelihood ratio test, i.e., they evolved at a different rate on different branches in the phylogenetic tree
6. tree.nwk: phylogenetic tree
7. data/proteome_dir/orthologs/Results_1/: OrthoFinder output

### Do the analysis
The analysis is then carried out in three steps:
1. First, the data is retrieved, and annotated coding regions are extract using [gffread](https://github.com/gpertea/gffread)<sup>4</sup>.<br>
This part is implemented in a separated Snakemake workflow<sup>5</sup> in the `data` directory (This was necessary because OrthoFinder can't be easily called from within Snakemake). To execture this part do:<br>
```
cd data
snakemake -c16 --use-conda --conda-frontend conda
```
2. Then, orthologous genes are identified using [OrthoFinder](https://github.com/davidemms/OrthoFinder)<sup>6</sup>. From within the `data` directory call:<br>
`orthofinder.py -t 16 -f proteome_dir -o proteome_dir/orthologs -n 1`

3. Finally, the dN/dS ratios are determined for single-copy orthologs using the [codeml](http://abacus.gene.ucl.ac.uk/software/paml.html) program included in the PAML package<sup>7</sup>. First, a multiple sequence alignment of protein sequences is performed, which is then reverse translated into nucleotide alignments. The alignments are concatenated to build a phylogenetic tree to support the analysis. dN/dS for each gene are computed based on nucleotide MSA using two codeml models: First, the null model that computes one ratio for the entire tree. Second, a branch specific model that calculates a ratios for each branch in the phylogenetic tree. If the branch-specific model passes the likelihood-ratio test, a gene is assumed to have evolved at a differential rate in the GBL. To perform this analysis, run:<br>
```
cd ..
# test for selection
snakemake -c16 --use-conda --conda-frontend conda
```
4. Upload gene list "genes_background.txt" to [DAVID](https://david.ncifcrf.gov/)<sup>8</sup>, select "OFFICIAL_GENE_SYMBOL" as identifier, Green Anole as species. Create and download a functional annotation table with the desired specificity (recommended 5) to "functional_annotation_background_david.tab". Then, run `scripts/compare_dnds_go_terms.py`. To perform GSEA, upload "genes_background.txt" to DAVID, select 'OFFICIAL_GENE_SYMBOL' as identifier, convert the gene IDs to Entrez IDs, and use this list as background. Then, upload "genes_under_selection" to DAVID as a gene list, and start the analysis.

Steps 1-3 can also be executed by running `./execute_workflow.sh'.

### Refrences
1. Bo Song, Shifeng Cheng, Yanbo Sun, Xiao Zhong, Jieqiong Jin, Rui Guan, Robert W Murphy, Jing Che, Yaping Zhang, Xin Liu, A genome draft of the legless anguid lizard, Ophisaurus gracilis, GigaScience, Volume 4, Issue 1, December 2015, s13742–015–0056–7, https://doi.org/10.1186/s13742-015-0056-7
2. Alföldi J, Di Palma F, Grabherr M, Williams C, Kong L, Mauceli E, Russell P, Lowe CB, Glor RE, Jaffe JD, Ray DA, Boissinot S, Shedlock AM, Botka C, Castoe TA, Colbourne JK, Fujita MK, Moreno RG, ten Hallers BF, Haussler D, Heger A, Heiman D, Janes DE, Johnson J, de Jong PJ, Koriabine MY, Lara M, Novick PA, Organ CL, Peach SE, Poe S, Pollock DD, de Queiroz K, Sanger T, Searle S, Smith JD, Smith Z, Swofford R, Turner-Maier J, Wade J, Young S, Zadissa A, Edwards SV, Glenn TC, Schneider CJ, Losos JB, Lander ES, Breen M, Ponting CP, Lindblad-Toh K. The genome of the green anole lizard and a comparative analysis with birds and mammals. Nature. 2011 Aug 31;477(7366):587-91. doi: 10.1038/nature10390. PMID: 21881562; PMCID: PMC3184186.
3. Lind, A.L., Lai, Y.Y.Y., Mostovoy, Y. et al. Genome of the Komodo dragon reveals adaptations in the cardiovascular and chemosensory systems of monitor lizards. Nat Ecol Evol 3, 1241–1252 (2019). https://doi.org/10.1038/s41559-019-0945-8
4. Pertea G and Pertea M. GFF Utilities: GffRead and GffCompare [version 1; peer review: 3 approved]. F1000Research 2020, 9:304 (https://doi.org/10.12688/f1000research.23297.1)
5. Mölder F, Jablonski KP, Letcher B et al. Sustainable data analysis with Snakemake [version 1; peer review: 1 approved, 1 approved with reservations]. F1000Research 2021, 10:33 (https://doi.org/10.12688/f1000research.29032.1)
6. Emms, D.M., Kelly, S. OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biol 20, 238 (2019). https://doi.org/10.1186/s13059-019-1832-y
7. Ziheng Yang, PAML 4: Phylogenetic Analysis by Maximum Likelihood, Molecular Biology and Evolution, Volume 24, Issue 8, August 2007, Pages 1586–1591, https://doi.org/10.1093/molbev/msm088
8. Brad T Sherman, Ming Hao, Ju Qiu, Xiaoli Jiao, Michael W Baseler, H Clifford Lane, Tomozumi Imamichi, Weizhong Chang, DAVID: a web server for functional enrichment analysis and functional annotation of gene lists (2021 update), Nucleic Acids Research, 2022;, gkac194, https://doi.org/10.1093/nar/gkac194
