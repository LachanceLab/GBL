## Identifying potentially venomous genes in the draft genome assembly of the Guatemalan Beaded Lizard

Predicted protein-coding regions are searched against known venom proteins of heldodermatid lizards and vipers, which can be convineintly accessed through [VenomZone](https://venomzone.expasy.org/)<sup>1-3</sup>. The protein search is done using `jackhmmer` from the [HMMER3](http://hmmer.org/) package<sup>4</sup>. We use a Bonferroni corrected e-value cut-off of 0.0001 / n - where n is the number of queries - for the full sequence and 0.01 / n for the single-best scoring domain. Furthermore, only the best for each query is retained.

### Install HMMER3 and required Python packages
`conda env create -f hmmer3.yaml`

### Input files
1. Protein sequences of annotated protein-coding regions<br>
  a) ../genome_annotation/gene_annotations/augustus.hints.aa<br>
2. Database of venome proteins from helodermatid lizards and vipers (included in repository)
3. Summary file from dN/dS analysis (../dNdS/summary_dnds.tab)
4. Softmasked draft genome assembly<br>
  a) We use the draft genome assembly in the `genome_annotation` subdirectory (`../genome_annotation/GBL_genome_assembly.fasta.combined.masked`) 

### Output files
1. `jackhmmer` output:<br>
  a) venom_alignments_hmmer.out<br>
  b) venom_alignments_hmmer.tab<br>
2. `jackhmmer` table output with add dN/dS ratios of hits putative venom proteins<br>
  a) venom_alignments_hmmer_dnds.tab<br>
3. Boxplot comparing the dN/dS ratios of putative venome proteins and non-venom proteins<br>
  a) dnds_venom_vs_non_venom.png<br>
4. Exonerate results from Helofensin search<br>
  a) exonerate_results.gff<br>

### Do the analysis
```
conda activate hmmer
jackhmmer -o venom_alignments_hmmer.out --cpu 16 --tblout venom_alignments_hmmer.tab ../genome_annotation/gene_annotation/augustus.hints.aa \
databases/heloderma_and_viper_venom_proteins.fasta
conda deactivate
./read_hmmer_results.py -i venom_alignments_hmmer.tab -a databases/heloderma_and_viper_venom_proteins.tab
```

Alternatively, you can also execute `./search_venomous_proteins_hmmer.sh`, which executes the above two commands.

## Searching for instances of proteins belonging to Helofensin type class
After the first iteration of searching for putatively venomous proteins that have similarity to known helodermatid venoms, representative genes of the Helofensin class were missing. Therefore, we aligned these protein sequences to genome using Exonerate<sup>5</sup> with the following command:

```
exonerate --model protein2genome --percent 50 --softmasktarget TRUE --showvulgar FALSE --showalignment FALSE --showtargetgff TRUE databases/heloderma_helofensin_class.fasta ../genome_annotation/GBL_genome_assembly.fasta.combined.masked | grep '^contig' > exonerate_results_helofensin.gff
```

Exonerate can be installed using conda from the provided .yaml file:<br>
`conda env create -f exonerate.yaml`

### References

1. “VenomZone.” n.d. Accessed April 12, 2022. https://venomzone.expasy.org/.
2. Mackessy, Stephen P. 2022. “Venom Production and Secretion in Reptiles.” Journal of Experimental Biology 225 (7). https://doi.org/10.1242/jeb.227348.
3. Koludarov, Ivan, Timothy N.W. Jackson, Kartik Sunagar, Amanda Nouwens, Iwan Hendrikx, and Bryan G. Fry. 2014. “Fossilized Venom: The Unusually Conserved Venom Profiles of Heloderma Species (Beaded Lizards and Gila Monsters).” Toxins 2014, Vol. 6, Pages 3582-3595 6 (12): 3582–95. https://doi.org/10.3390/TOXINS6123582.
4. Jaina Mistry, Robert D. Finn, Sean R. Eddy, Alex Bateman, Marco Punta, Challenges in homology search: HMMER3 and convergent evolution of coiled-coil regions, Nucleic Acids Research, Volume 41, Issue 12, 1 July 2013, Page e121, https://doi.org/10.1093/nar/gkt263
5. Slater, G.S.C., Birney, E. Automated generation of heuristics for biological sequence comparison. BMC Bioinformatics 6, 31 (2005). https://doi.org/10.1186/1471-2105-6-31
