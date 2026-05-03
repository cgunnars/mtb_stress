### Motivation
This repo reflects data cleaning, collation, and exploration used to prioritize transcriptional reporter genes in Mycobacterium tuberculosis. The goal was to identify genes that were acid-associated, rather than generic stress markers, across multiple Mtb isolates (CDC1551, Erdman, H37Rv).

Prior literature on acid stress used the strain CDC1551 and *aprA* as a reporter gene, but *aprA*-based reporters are not acid-induced in H37Rv. This is associated with, but not entirely explained by, [a SNP in the two-component system *phoR*](https://doi.org/10.1099/mic.0.000036).

At the time, the [iModulonDB](https://doi.org/10.1128/msphere.00033-22) project was not published for Mtb, so I needed to manually curate gene expression datasets over Mtb strains and stress-relevant conditions. I included microarray and RNAseq data from axenic stress conditions and intracellular infections, but excluded datasets that were predominantly focused on antibiotic stress.

### Overview
1. Align raw reads. I use a Makefile target `raw_reads`, which requires a text file containing SRA accessions to download and align `data/SRR_Acc_List.txt`. Each line specifies a valid SRA accession. It additionally requires:
- fastq-dump from sra-tools, which handles downloading from SRA.
- repair.sh from bbmap, which handles fixing paired-end fastqs
- bowtie2 for alignment and a valid pre-built index fastq/H37Rv. `bowtie2-build reference_sequence.fasta index_name` 
- samtools
- featureCounts from subRead and a valid genome annotation file. I used the function `src/gffToSaf.py` to parse my GFF into a simple annotation format. This is because the gff for Mtb from Mycobrowser was not in valid GFF3 format.

2. Download external data that isn't in count tables. I use the Makefile targets `RNAseq-CDC1551_Abramovitch(2)` and `microarray-CDC1551_Russell`, or the target `externalDE` to get everything.

3. Download raw counts using the makefile target `loadgse`. 

4. Standardize gene nomenclature in counts tables, using the target `data/GSEXXXX_raw_count_std.txt`. Individual counts files need to be treated differently because gene names/locus names and standards for capitalization differ -- this would be resolved by re-processing everything with the same GFF to generate uniformly annotated counts tables. Here, I use awk to avoid computational cost associated with re-processing everything.

5. Convert genes to filtered, log tpm data, using the target `data/GSEXXXX_tpm.txt`, `data/GSEXXXX_filt-tpm.txt`, or `data/GSEXXXX_ln-filt-tpm.txt`.

6. Join standardized counts or filtered, log-normalized data. We do this for everything (`data/all_raw_count_std.txt`, `data/all_tpm.txt`) and for pH datasets specifically (`data/ph_tpm.txt`).

7. Batch-correct data using combat applied to log TPM data `data/all-combat_ln-filt-tpm.txt`. This also requires sample and batch information `data/all_pheno.txt`. We treat each accession as a batch. Running `var-genes` will generate a list of variable genes with or without batch correction.

8. Perform differential expression analysis on raw counts using the rule `DE`.

9. Get consensus pH, iron, and SDS-regulated genes using the rule `data/consensus_pH-DE.csv`

10. Plot consensus, most variable, and canonical stress genes on a heatmap using `consensus_fig`, `fig/%_stress-tx.pdf`,`fig/%_var-genes.pdf`. You can also use `fig` to generate all figures, which include visualizations of PCA (`fig_PCA`), correlation (`fig_corr`), and DE analyses (`fig_DE`).

### Installation
This repo was made in 2019 on an x86 Mac, so some aspects of this repo may not reproduce well. Instead, this repo is intended to be a resource for general data cleaning operations.

That said, see 1. for general commandline requirements for fastq alignment. The `packrat.lock` file specifies R packages used and `requirements.txt` specifies Python packages used.
