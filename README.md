## Compelementary analysis of a draft genome of the Guatemalan Beaded Lizard (Heloderma charlesbogerti) 

The first step will be to annotate the draft genome, as some of the subsequent analysis requires information about protein coding regions. Then, the effective population size over time is estimated using a pairwise-sequential markov chain as well as the ratio of non-synonymous to synonymous substitution (dN/dS) for each gene using. Additionally, we perform a crude search for potentially venomous genes. Detailed description of each step are provided in the corresponding subdirectories.


### Software requirements

For annotating the genome, estimating the effective population size over time, and the dN/dS ratios you will need to have snakemake<sup>1</sup> installed. The easiest way to do that and to ensure that all the other required package will be installed along the way, is by creating a conda environment from the provide .yaml file This can be done by running `conda env create -f snakemake.yaml`. Otherwise snakemake can also be installed from scratch or already installed version can be used. Additional package that will be required are: `biopython v1.79`, `matplotlibv3.5.0`, `numpy v1.21.4`, `pandas v1.3.4`, and `scipy v1.7.3`<br>

Additonal software requirements :
- [gffread v0.12.7](https://github.com/gpertea/gffread)<sup>2</sup> (needs to be added to the PATH variable)
- [OrthoFinder v2.5.4](https://github.com/davidemms/OrthoFinder)<sup>3</sup> (needs to be added to the PATH variable)
- [codeml v4.9](http://abacus.gene.ucl.ac.uk/software/paml.html)<sup>4</sup> (assumed to be in dN/dS subdirectory)
- [HMMER3 v3.1b2](http://hmmer.org/)<sup>5</sup> (can be installed from provided .yaml file in `venom/envs` subdirectory)
- [RepeatModeler v2.0.3 & RepeatMasker v4.1.2](https://www.repeatmasker.org/)<sup>6,7</sup> (can be installed from provided .yaml in `genome_annotation/envs` subdirectory)
- [TRF v4.09.1](https://tandem.bu.edu/trf/trf.html)<sup>8</sup> (can be installed from provided .yaml in `genome_annotation` subdirectory)
- [BUSCO v5.3.2](https://busco.ezlab.org/)<sup>9</sup> (can be installed from provided .yaml in `genome_annotation/envs` subdirectory)
- [BRAKER2 v2.1.6](https://github.com/Gaius-Augustus/BRAKER)<sup>10</sup> (needs to be manually installed by following the instruction provided in the linked GitHub. We provide a .yaml to install some of the dependencies in the `genome_annotation/envs` subdirectory)
- [tRNAscan-SE v2.0.9](https://github.com/UCSC-LoweLab/tRNAscan-SE)<sup>11</sup> (can be install from provided .yaml in `genome_annotation/envs` subdirectory)
- [InterProScan v55.5-88.0](https://interproscan-docs.readthedocs.io/en/latest/Introduction.html)<sup>12,13</sup> (needs to be added to the PATH variable)
- [Exonerate v2.4.0](https://github.com/nathanweeks/exonerate)<sup>14</sup> (can be installed from provided .yaml in `venom` subdirectory)

### Data requirements
The specific input and output files are described in detail in the respective subdirectories. However, you will need to retrieve the draft genome assembly of the Guatemalan Beaded Lizard from and the bam files with the raw long-reads from NCBI. The BioProject accession is PRJNA834834. A copy of the draft genome assembly should be short in the subdirectory `genome_annotation`, and a copy of the draft genome assembly and the bam files should be added to the subdirectory `effective_population_size`.

To perform the individual analyses, `cd` into the corresponding directories and follow the instructions in the README for execution of the commands.
For example, to annotate genes, do:
```
# change directory
cd genome_annotation
conda activate snakemake
# Do analysis
snakemake --use-conda --conda-frontend conda --latency-wait 20 --cores 36
```
## References

1. Mölder F, Jablonski KP, Letcher B et al. Sustainable data analysis with Snakemake [version 1; peer review: 1 approved, 1 approved with reservations]. F1000Research 2021, 10:33 (https://doi.org/10.12688/f1000research.29032.1)
2. Pertea G and Pertea M. GFF Utilities: GffRead and GffCompare [version 1; peer review: 3 approved]. F1000Research 2020, 9:304 (https://doi.org/10.12688/f1000research.23297.1)
3. Emms, D.M., Kelly, S. OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biol 20, 238 (2019). https://doi.org/10.1186/s13059-019-1832-y
4. Ziheng Yang, PAML 4: Phylogenetic Analysis by Maximum Likelihood, Molecular Biology and Evolution, Volume 24, Issue 8, August 2007, Pages 1586–1591, https://doi.org/10.1093/molbev/msm088
5. Jaina Mistry, Robert D. Finn, Sean R. Eddy, Alex Bateman, Marco Punta, Challenges in homology search: HMMER3 and convergent evolution of coiled-coil regions, Nucleic Acids Research, Volume 41, Issue 12, 1 July 2013, Page e121, https://doi.org/10.1093/nar/gkt263
6. Smit, AFA, Hubley, R. RepeatModeler Open-1.0. 2008-2015 <http://www.repeatmasker.org>
7. Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0. 2013-2015 <http://www.repeatmasker.org>.
8. G. Benson, "Tandem repeats finder: a program to analyze DNA sequences" Nucleic Acids Research (1999) Vol. 27, No. 2, pp. 573-580.
9. Manni, Mosè, Matthew R. Berkeley, Mathieu Seppey, Felipe A. Simão, and Evgeny M. Zdobnov. 2021. “BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes.” Molecular Biology and Evolution 38 (10): 4647–54. https://doi.org/10.1093/MOLBEV/MSAB199.
10. Tomáš Brůna, Katharina J Hoff, Alexandre Lomsadze, Mario Stanke, Mark Borodovsky, BRAKER2: automatic eukaryotic genome annotation with GeneMark-EP+ and AUGUSTUS supported by a protein database, NAR Genomics and Bioinformatics, Volume 3, Issue 1, March 2021, lqaa108, https://doi.org/10.1093/nargab/lqaa108
11. Chan, Patricia P., Brian Y. Lin, Allysia J. Mak, and Todd M. Lowe. 2021. “TRNAscan-SE 2.0: Improved Detection and Functional Classification of Transfer RNA Genes.” Nucleic Acids Research 49 (16): 9077–96. https://doi.org/10.1093/NAR/GKAB688.
12. The InterPro protein families and domains database: 20 years on Matthias Blum, Hsin-Yu Chang, Sara Chuguransky, Tiago Grego, Swaathi Kandasaamy, Alex Mitchell, Gift Nuka, Typhaine Paysan-Lafosse, Matloob Qureshi, Shriya Raj, Lorna Richardson, Gustavo A Salazar, Lowri Williams, Peer Bork, Alan Bridge, Julian Gough, Daniel H Haft, Ivica Letunic, Aron Marchler-Bauer, Huaiyu Mi, Darren A Natale, Marco Necci, Christine A Orengo, Arun P Pandurangan, Catherine Rivoire, Christian J A Sigrist, Ian Sillitoe, Narmada Thanki, Paul D Thomas, Silvio C E Tosatto, Cathy H Wu, Alex Bateman, Robert D Finn Nucleic Acids Research (2020), gkaa977, PMID: 33156333
13. InterProScan 5: genome-scale protein function classification Philip Jones, David Binns, Hsin-Yu Chang, Matthew Fraser, Weizhong Li, Craig McAnulla, Hamish McWilliam, John Maslen, Alex Mitchell, Gift Nuka, Sebastien Pesseat, Antony F. Quinn, Amaia Sangrador-Vegas, Maxim Scheremetjew, Siew-Yit Yong, Rodrigo Lopez, Sarah Hunter Bioinformatics (2014), PMID: 24451626
14. Slater, G.S.C., Birney, E. Automated generation of heuristics for biological sequence comparison. BMC Bioinformatics 6, 31 (2005). https://doi.org/10.1186/1471-2105-6-31

