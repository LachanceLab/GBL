## Inferring the effective population size of the Guatemalan Beaded Lizard (GBL)

### Input files
1. Draft genome of Guatemalan Beaded Lizard
2. Raw long-reads

The BioProject accession is PRJNA834834. The files need to be added to this directory. Note that you may need to update the `subreads_A`, `subreads_B`, and `reference` paths in `config/config.yaml`.

### Output files
1. psmc_results/:PSMC results including bootstraps
2. gbl.png: Visualization of PSMC results
3. Other files generated along the way include:<br>
  a) gbl.subreads.bam: Merged raw reads <br>
  b) GBL_genome_assembly.[fasta.fai, mmi]: Indexed draft genome<br>
  c) gbl_ref.alignmentset.bam[.bai]: Alignment of raw reads to draft genome<br>
  d) gbl_ref.alignmentset.ref.collapsed.fasta[.fai]: Homopolymer collapsed draft genome<br>
  e) gbl_ref.alignmentset.ref.collapsed_filtered_roh.fasta[.fai]: Homopolymer and ROH masked draft genome<br>
  f) vcfs_per_scaff/: Contains variant calls per contig<br>
  g) gbl.variants.vcf.gz[.tbi]: merged called variants<br>
  h) gbl.variants_filtered.vcf.gz[.tbi]: merged & filtered variants<br>
  i) contig_stats.tab: Contig statistics (size and heterozygosity)<br>
  j) selected_contigs.txt: Contigs that fulfill size and heterozygosity requirements<br>
  k) gbl.consensus.[fasta, fastq, pmscfa]: IUPAC consensus sequence prepared for PSMC input<br>
  l) gbl.consensus_split.psmcfa: PSMC input prepared for bootstrapping<br>


### Do the analysis

The analysis includes the following steps:<br>

1. Align raw reads to assembled genome using [PacBio's minimap2](https://github.com/PacificBiosciences/pbmm2) for long-reads<br>

2. Call variants with [Longshot](https://github.com/pjedge/longshot)<sup>1</sup> for long-reads<br>

3. Build consensus sequence using IUPAC encoding (applied filter: DP >=40; AC >= 10; GQ >=20) using bcftools<sup>2</sup><br>

4. Run [PSMC](https://github.com/lh3/psmc)<sup>3</sup> with bootstrapping to get confident intervals<br>

Parameters for the different analysis steps can be modified in the `config/config.yaml` file.

The long-reads and draft-assembly need to be downloaded separately. The entire workflow is implemented using [Snakemake](https://snakemake.readthedocs.io/en/stable/)<sup>4</sup>, and can be executed using 16 cores by running `snakemake -c16 --use-conda --conda-frontend conda`. The script `execute_snakefile.pbs` is a PBS script to run the workflow on a cluster.<br>

The results can visualized using the provided `./plot_psmc.py` script. For example, run:
```
./plot_psmc.py -i psmc_results/gbl_diploid.psmc -d psmc_results/ -o ./ -p gbl_psmc -m 6.69e-9 -g 8 -l GBL -x 1e4 -X 1e6
``` 
To get more options, do `./plot_psmc.py --help`.

### References
1. Edge, P., Bansal, V. Longshot enables accurate variant calling in diploid genomes from single-molecule long read sequencing. Nat Commun 10, 4660 (2019). https://doi.org/10.1038/s41467-019-12493-y
2. Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li, Twelve years of SAMtools and BCFtools, GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008
3. Li, Heng, and Richard Durbin. 2011. “Inference of Human Population History from Individual Whole-Genome Sequences.” https://doi.org/10.1038/nature10231.
4. Mölder F, Jablonski KP, Letcher B et al. Sustainable data analysis with Snakemake [version 1; peer review: 1 approved, 1 approved with reservations]. F1000Research 2021, 10:33 (https://doi.org/10.12688/f1000research.29032.1)


