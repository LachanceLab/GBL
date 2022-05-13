## Genome Annotation of draft assembly of GBL
Here, we provide a Snakemake<sup>1</sup> workflow for the annotation of the GBL draft genome. We use RepeatModeler<sup>2</sup>, RepeatMasker<sup>3</sup>, and TRF<sup>4</sup> to softmask repeats in the genome, and subsequentyl annotate genes with BRAKER2<sup>5</sup>. tRNA genes are predicted with tRNAscan-SE 2.0<sup>6</sup>.

Because other squamate genomes (e.g., X. tropicalis) have long tandem repeats (i.e., tandem repeats with period higher > 10) which are not masked by RepeatMasker we run TRF separately and combined the two maskings. For this, we follow the workflow by Brůna et al. (2021), which is described in the Supplemental Information (1.4.2).<sup>4</sup>
The script [parseTrfOutput.py](https://github.com/gatech-genemark/BRAKER2-exp/blob/master/bin/trf-scripts/parseTrfOutput.py) is a custom script by Tomáš Brůna, and needs be retrieved from his GitHub and added to the `scripts` directory.<br>

We recommend to run TRF before executing the Snakemake workflow since for some unknown reason the rule calling TRF fails without error message, although correctly generating the desired output files. TRF is called with the following command:<br>
```trf <genome.fasta> 2 7 7 80 10 50 500 -d -m -h```<br>
where `<genome.fasta>` is the path to the assembly to be analyzed.

To assess the sensitivity of the gene predictions, we assess the overlap of predicted protein-coding regions with single-copy orthologs that are expected to be in every vertebrate genome using BUSCO.<sup>7</sup> Furthermore, the predicted proteins were functionally annotated using InterProScan version 5.55-88.0.<sup>8,9</sup>

Because tRNAscan-SE predicts more than 3000 tRNA genes, we filter tRNA gene predictions based on overlaps with annotated repeat regions. Only tRNA genes that overlap at most to 50% with a repeat are retained.

### Steps:
1. Build custom repeat library using RepeatModeler<sup>5</sup>
2. Run RepeatMasker and TRF, and combined the softmasking
3. Run BRAKER2
4. Run BUSCO in protein-mode with annotated proteins (uses a modified version of the `generate_plot.py` script to plot the results)
5. Functioannly annotated predicted proteins with InterProScan
6. Run tRNAscan-SE
7. Filter predicted tRNAs based on overlaps with annotated repeats (>50%)

### Input files
1. A copy of the draft genome assembly must be stored in this directory. The BioProject accession is PRJNA834834. Note that you may have to update the `assembly` path in `config/config.yaml`<br>
2. OrthoDB protein database. Prepare as described in the [ProtHint](https://github.com/gatech-genemark/ProtHint)<sup>10</sup> GitHub repository and update `odb` path in `config/config.yaml`

### Output files
1. RepeatModeler<br>
  a) gbl_repeats-families.fa<br>
  b) gbl_repeats-families.stk<br>
  c) gbL_repeats.*<br>
2. RepeatMasker<br>
  a) GBL_genome_assembly.fasta.masked<br>
  b) GBL_genome_assembly.fasta.out<br>
  c) GBL_genome_assembly.fasta.tbl<br>
  d) GBL_genome_assembly.fasta.cat.gz<br>   
3. TRF<br>
  a) GBL_genome_assembly.fasta..2.7.7.80.10.50.500.dat<br>
  b) GBL_genome_assembly.fasta.2.7.7.80.10.50.500.mask<br>
  c) GBL_genome_assembly.fasta.2.7.7.80.10.50.500.raw.gff<br>
  d) GBL_genome_assembly.fasta.2.7.7.80.10.50.500.sorted.gff<br>
  e) GBL_genome_assembly.fasta.2.7.7.80.10.50.500.merged.gff<br>
4. Merged repeat annotation from RepeatMasker and TRF<br>
  a) GBL_genome_assembly.fasta.combined.masked<br>
  b) GBL_genome_assembly.fasta.combined.masked.gff<br>
5. BRAKER2<br>
  a) gene_annotation/augustus.hints.gff3 (GFF file of annotated coding regions)<br>
  b) gene_annotation/augustus.hints.aa (Protein sequences of annotated coding regions)<br>
  c) See [BRAKER2 repository](https://github.com/Gaius-Augustus/BRAKER) for the description of the other output files generated in `gene_annotation`<br>
6. InterProScan<br>
  a) functional_annotation/augustus.hints.aa.cleaned.tsv<br>
  b) See [InterProScan docs](https://interproscan-docs.readthedocs.io/en/latest/) for a description of the other output files<br>
 7. BUSCO (proteome)<br>
  a) plots/busco_figure.png<br>
  b) short_summary_specific.vertebrate_odb10.<species>.txt<br>
  c) BUSCO will also create subdirectories for each analyzed species<br>
8. tRNAscan-SE<br>
  a) filtered_trnas.gff (removed tRNAs that are overlapping with annotated repeat region)<br>
  b) raw tRNAscan-SE output files (trnas.*)<br>
9. Annotated genes and (filtered) tRNAs merged into one gff file (sorted according to chrom sweep algortihm)
  a) predicted_genes_and_trnas_merged.gff
  
### Do the analysis
The Snakemake workflow can be executed by the following command using 16 cores:<br>
`snakemake -c16 --use-conda --conda-frontend conda`


### References
1. Mölder F, Jablonski KP, Letcher B et al. Sustainable data analysis with Snakemake [version 1; peer review: 1 approved, 1 approved with reservations]. F1000Research 2021, 10:33 (https://doi.org/10.12688/f1000research.29032.1)
2. Smit, AFA, Hubley, R. RepeatModeler Open-1.0. 2008-2015 <http://www.repeatmasker.org>
3. Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0. 2013-2015 <http://www.repeatmasker.org>.
4. G. Benson, "Tandem repeats finder: a program to analyze DNA sequences" Nucleic Acids Research (1999) Vol. 27, No. 2, pp. 573-580.
5. Tomáš Brůna, Katharina J Hoff, Alexandre Lomsadze, Mario Stanke, Mark Borodovsky, BRAKER2: automatic eukaryotic genome annotation with GeneMark-EP+ and AUGUSTUS supported by a protein database, NAR Genomics and Bioinformatics, Volume 3, Issue 1, March 2021, lqaa108, https://doi.org/10.1093/nargab/lqaa108
6. Chan, Patricia P., Brian Y. Lin, Allysia J. Mak, and Todd M. Lowe. 2021. “TRNAscan-SE 2.0: Improved Detection and Functional Classification of Transfer RNA Genes.” Nucleic Acids Research 49 (16): 9077–96. https://doi.org/10.1093/NAR/GKAB688.
7. Manni, Mosè, Matthew R. Berkeley, Mathieu Seppey, Felipe A. Simão, and Evgeny M. Zdobnov. 2021. “BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes.” Molecular Biology and Evolution 38 (10): 4647–54. https://doi.org/10.1093/MOLBEV/MSAB199.
8. The InterPro protein families and domains database: 20 years on Matthias Blum, Hsin-Yu Chang, Sara Chuguransky, Tiago Grego, Swaathi Kandasaamy, Alex Mitchell, Gift Nuka, Typhaine Paysan-Lafosse, Matloob Qureshi, Shriya Raj, Lorna Richardson, Gustavo A Salazar, Lowri Williams, Peer Bork, Alan Bridge, Julian Gough, Daniel H Haft, Ivica Letunic, Aron Marchler-Bauer, Huaiyu Mi, Darren A Natale, Marco Necci, Christine A Orengo, Arun P Pandurangan, Catherine Rivoire, Christian J A Sigrist, Ian Sillitoe, Narmada Thanki, Paul D Thomas, Silvio C E Tosatto, Cathy H Wu, Alex Bateman, Robert D Finn Nucleic Acids Research (2020), gkaa977, PMID: 33156333
9. InterProScan 5: genome-scale protein function classification Philip Jones, David Binns, Hsin-Yu Chang, Matthew Fraser, Weizhong Li, Craig McAnulla, Hamish McWilliam, John Maslen, Alex Mitchell, Gift Nuka, Sebastien Pesseat, Antony F. Quinn, Amaia Sangrador-Vegas, Maxim Scheremetjew, Siew-Yit Yong, Rodrigo Lopez, Sarah Hunter Bioinformatics (2014), PMID: 24451626
10. Tomáš Brůna, Alexandre Lomsadze, Mark Borodovsky, GeneMark-EP+: eukaryotic gene prediction with self-training in the space of genes and proteins, NAR Genomics and Bioinformatics, Volume 2, Issue 2, June 2020, lqaa026, https://doi.org/10.1093/nargab/lqaa026
