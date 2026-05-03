SHELL := /bin/bash
# get table of gene lengths for tpm
gene_info: ./data/H37Rv_gene-lengths.csv
./data/H37Rv_gene-lengths.csv: src/makeGeneLengths.py	
	python $< -i data/Mycobacterium_tuberculosis_H37Rv_txt_v2.txt -o $@

#ALIGN RAW READS (PAIRED END) 
raw_srr := $(foreach n, $(cat data/SRR_Acc_List.txt), $(addprefix $(addprefix fastq/, n), -counts.txt)); 
raw_reads: $(raw_srr)
fastq/%-counts.txt: 
	parallel -j 1 fastq-dump --skip-technical -F --split-files -O fastq {} ::: $*; repair.sh in1=fastq/$*_1.fastq in2=fastq/$*_2.fastq out1=fastq/$*_fixed_1.fastq out2=fastq/$*_fixed_2.fastq; rm fastq/$*_1.fastq fastq/$*_2.fastq; bowtie2 -x fastq/H37Rv -1 fastq/$*_fixed_1.fastq -2 fastq/$*_fixed_2.fastq | samtools view -Su | samtools sort | featureCounts -p -F SAF -O -a fastq/H37Rv_saf.gff -o fastq/$*-counts.txt; rm fastq/$*_fixed_1.fastq fastq/$*_fixed_2.fastq


#DOWNLOAD EXTERNAL DATA (i.e. not in count tables)
cdc_gse := 52020 63917 89106 
RNAseq-CDC1551_Abramovitch_raw := data/GSE89106_CDC1551_filt-tpm.txt 
RNAseq-CDC1551_Abramovitch_gene-list := data/GSE89106_CDC1551_pH-DE_DMSO5-vs-DMSO7.txt
RNAseq-CDC1551_Abramovitch: $(RNAseq-CDC1551_Abramovitch_raw) $(RNAseq-CDC1551_Abramovitch_gene-list)
data/GSE89106_CDC1551_filt-tpm.txt: src/proc/GSE89106.R data/GSE89106/GSE89106_Table_S3.xlsx
	Rscript $< -i $(word 2, $^) -sheet 3 -o data/GSE89106_CDC1551
data/GSE89106_CDC1551_pH-DE_DMSO5-vs-DMSO7.txt: data/GSE89106_CDC1551_filt-tpm.txt
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi

microarray-CDC1551_Russell_full-DE := data/GSE8867_CDC1551_pH-DE_6.5-vs-7.csv data/GSE8867_CDC1551_pH-DE_5.5-vs-7.csv 
microarray-CDC1551_Russell_gene-list := data/GSE8867_CDC1551_pH-DE_6.5-vs-7.txt data/GSE8867_CDC1551_pH-DE_5.5-vs-7.txt
microarray-CDC1551_Russell: $(microarray-CDC1551_Russell_full-DE) $(microarray-CDC1551_Russell_full-DE)
data/GSE8867_CDC1551_pH-DE_6.5-vs-7.csv: src/proc/GSE8867.R data/GSE8867-external-DE.xls
	Rscript $< -i $(word 2, $^) -sheet 19 20 21 22 -o data/GSE8867_CDC1551_pH-DE 
data/GSE8867_CDC1551_pH-DE_5.5-vs-7.csv: data/GSE8867_CDC1551_pH-DE_6.5-vs-7.csv
	@if test -f $@; then :; else \
		rm -f $<; \
		make $<; \
	fi
data/GSE8867_CDC1551_pH-DE_6.5-vs-7.txt: data/GSE8867_CDC1551_pH-DE_6.5-vs-7.csv
	@if test -f $@; then :; else \
		rm -f $<; \
		make $<; \
	fi
data/GSE8867_CDC1551_pH-DE_6.5-vs-7.txt: data/GSE8867_CDC1551_pH-DE_6.5-vs-7.csv
	@if test -f $@; then :; else \
		rm -f $<; \
		make $<; \
	fi

RNAseq-CDC1551_Abramovitch2: data/GSE63917_CDC1551_pH-DE_5.7-vs-7
data/GSE63917_CDC1551_pH-DE_5.7-vs-7.txt: src/proc/GSE63917.R
	Rscript $< -i data/GSE63917/GSE63917_SupplementalTable1A_DMSO_pH5.7_v_pH7_induction.txt \
	data/GSE63917/GSE63917_SupplementalTable1B_DMSO_pH5.7_v_pH7_repression.txt -o data/GSE63917_CDC1551

data/GSE52020_CDC1551-glycerol_pH-DE_5.7-vs-7.txt: src/proc/GSE52020.R
	Rscript $< -i data/GSE52020/GSE52020_SupplementalTable2.xls -sheet 1 -o data/GSE52020_CDC1551-glycerol
data/GSE52020_CDC1551-pyruvate_pH-DE_5.7-vs-7.txt: src/proc/GSE52020.R
	Rscript $< -i data/GSE52020/GSE52020_SupplementalTable3.xls -sheet 1 -o data/GSE52020_CDC1551-pyruvate

externalDE: $(microarray-CDC1551_Russell) $(RNAseq-CDC1551_Abramovitch) $(RNAseq-CDC1551_Abramovitch2)

#DOWNLOAD RAW COUNT DATA
raw_data: ./data/GSE132354/GSE132354_Mtb_raw_counts_matrix.txt \
	  ./data/GSE67035/GSE67035_U19_expr_counts.txt \
	  ./data/GSE123267/GSE123267_raw_data_counts_pH.txt \
	  ./data/GSE123267/GSE123267_raw_data_counts.txt
gse := 67035 123267 132354 
loadgse: 
	for n in $(gse); do\
		Rscript src/downloadGSE.R GSE$$n; \
		yes n | gunzip data/GSE$$n/*.gz; \
		rm -rf data/GSE$$n/*.gz; \
	done

data_prefix := $(foreach n, $(gse), $(addprefix data/GSE, $n))
#GET INTO FORMAT WHERE GENE/RV FIRST COL, REST IS SAMPLE COUNTS
std_data := $(foreach n, $(data_prefix), $(addprefix $n, _raw_count_std.txt))
std-srp_data := data/SRP142541_H37Rv_raw_count_std.txt $(std_data)

#fix capitalization "DnaA" > "dnaA" except for Rvs, PPE, PE
data/GSE132354_Erdman_raw_count_std.txt: data/GSE132354/GSE132354_Mtb_raw_counts_matrix.txt
	awk 'BEGIN{FS=OFS="\t"}{if (substr($$1,1,3)!="PPE" && substr($$1,1,2)!="PE" && substr($$1,1,2)!="Rv"){$$1=(tolower(substr($$1,1,1))substr($$1,2))}; print}' $< > $@
# concatenating data, remove unnecessary quotation marks
data/GSE123267_H37Rv_raw_count_std.txt: data/GSE123267/GSE123267_raw_data_counts.txt data/GSE123267/GSE123267_raw_data_counts_pH.txt
	cut -f 2- $< | paste -d '	' $(word 2, $^) - | sed 's/\"//g' > $@ 
# process by locus, remove extra metadata about gene
# fix manB, lysS, sdhA duplicate issues manually
# Rv3257c = pmmA
# Rv1640c = lysX
# Rv0248c = sdhA_1
data/GSE67035_H37Rv_raw_count_std.txt: data/GSE67035/GSE67035_U19_expr_counts.txt
	awk 'BEGIN{OFS="\t"}{if ($$2=="-"){$$2=$$1} else if ($$1=="Rv3257c"){$$2="pmmA"} else if ($$1=="Rv1640c"){$$2="lysX"} else if ($$1=="Rv0248c"){$$2="sdhA_1"}; print}' $< | cut -d $$'\t' -f 2,5-26 > $@;

#CONCATENATE FOR BATCH CORR
data/all_raw_count_std.txt: $(std_data)
	./src/xjoin.sh $^ > $@
data/all-srp_raw_count_std.txt: 
	./src/xjoin.sh data/GSE*_raw_count_std.txt data/SRP*_raw_count_std.txt > $@

#CONVERT TO TPM
tpm_data := $(foreach n, $(data_prefix), $(addprefix $n, _tpm.txt))
tpm-srp_data := $(tpm_data) data/SRP142541_H37Rv_tpm.txt
data/%_tpm.txt: src/countToTpm.py data/%_raw_count_std.txt data/H37Rv_gene-lengths.csv
	python $< -i $(word 2, $^) -g $(word 3, $^) -o $@
data/GSE67035_H37Rv_tpm.txt: src/countToTpm.py data/GSE67035_H37Rv_raw_count_std.txt data/H37Rv_gene-lengths.csv
	python $< -i $(word 2, $^) -g $(word 3, $^) -o $@ # --locus

#MAKE FIGURES
fig_prefix := $(foreach n, $(gse), $(addprefix fig/GSE, $n))

#plot readct distribution
fig_read-dist := $(foreach n, $(fig_prefix), $(addprefix $n, _read-dist.pdf))
fig/%_read-dist.pdf: src/plotReadDist.py data/%_tpm.txt src/util.py
	python $< -i $(word 2, $^) -o $@

#compile gse files
data/ph_tpm.txt: data/GSE67035_H37Rv_tpm.txt data/GSE123267_H37Rv_tpm.txt
	./src/xjoin.sh $^ > $@
data/all_tpm.txt: $(tpm_data)
	./src/xjoin.sh $^ > $@
data/all-srp_tpm.txt: $(tpm-srp_data)
	./src/xjoin.sh data/GSE*_tpm.txt data/SRP*_tpm.txt > $@

filt-tpm := $(for each n, $(data_prefix), $(addprefix $n, _filt-tpm.txt)) data/SRP142541_H37Rv_filt-tpm.txt data/all_filt-tpm.txt data/ph_filt-tpm.txt data/all-srp_filt-tpm.txt
ln-filt-tpm := $(for each n, $(data_prefix), $(addprefix $n, _ln-filt-tpm.txt)) data/SRP142541_H37Rv_ln-filt-tpm.txt data/all_ln-filt-tpm.txt data/ph_ln-filt-tpm.txt \
	       data/all-srp_ln-filt-tpm.txt
data/%_filt-tpm.txt: data/%_tpm.txt src/util.py
	python -c "from src.util import filter_genes; import pandas as pd; df = filter_genes(pd.read_csv('$<', sep='\t', index_col=0)); df.to_csv('$@', sep='\t')"
data/%_ln-filt-tpm.txt: data/%_filt-tpm.txt src/util.py
	python -c "import numpy as np; import pandas as pd; df = np.log(1 + pd.read_csv('$<', sep='\t', index_col=0).astype('float'));print(df.head()); df.to_csv('$@', sep='\t')"
#compile phenotypic data
pheno: data/SRP142541_H37Rv_pheno.txt
data/all_pheno.txt: $(tpm_data)
	echo "name	condition	batch" > $@; i=0; for f in $^; do \
		((i=i+1)); \
		head -1 $$f | cut -d '	' -f2- | tr '\t' '\n' | awk -v awkvar="$$i" 'BEGIN{FS=OFS="\t"}{print $$1, substr($$1, 1, length($$1)-1), awkvar}' >> $@; \
	done
data/SRP142541_H37Rv_pheno.txt: data/SRP142541_H37Rv_raw_count_std.txt
	echo "name	condition	batch" > $@; head -1 $^ | cut -d '	' -f2- | tr '\t' '\n' | awk 'BEGIN{FS=OFS="\t"}{print $$1, substr($$1, 1, length($$1)-5), 9}' >> $@
data/all-srp_pheno.txt: data/all_pheno.txt data/SRP142541_H37Rv_pheno.txt 
	(head -1 data/all_pheno.txt && tail -n +2 -q data/*_pheno.txt) > $@

data/all-combat_ln-filt-tpm.txt: src/runCombat.py data/all_ln-filt-tpm.txt data/all_pheno.txt
	python $< -i $(word 2, $^) -f $(word 3, $^) -o $@
data/all-DESeq2_filt-tpm.txt: src/makeDE.R data/all_raw_count_std.txt data/all_pheno.txt
	Rscript $< -i $(word 2, $^) -f $(word 3, $^) -o $@


#(3) get most variable genes
fig/all-combat_cv.pdf: src/getVarGenes.py data/all-combat_ln-filt-tpm.txt
	python $< -i $(word 2, $^) -o $@

var-genes: $(foreach n, $(data_prefix), $(addprefix $n, _var-genes.txt)) data/all_var-genes.txt data/ph_var-genes.txt data/all-combat_var-genes.txt data/all-srp-combat_var-genes.txt
data/%_var-genes.txt: data/%_filt-tpm.txt src/util.py
	python -c "from src.util import get_var_genes; import pandas as pd; genes = get_var_genes(pd.read_csv('$<', sep='\t', index_col=0), thresh=200); f = open('$@', 'w'); f.writelines(gene + '\n' for gene in genes); f.close()"

data/all-genes.txt: data/all_tpm.txt
	tail -n +2 $^ | cut -f 1  > $@

DE: data/GSE67035_H37Rv_pH-DE.csv data/GSE67035_H37Rv_Fe-DE.csv data/GSE123267_H37Rv_pH-DE.csv data/SRP142541_H37Rv_Fe-DE.csv
data/GSE67035_H37Rv_pH-DE.csv: src/makeDE.R data/GSE67035_H37Rv_raw_count_std.txt 
	Rscript $< -i data/all_raw_count_std.txt -f data/all_pheno.txt -o $(basename $@) -conds 37Rv_TOL_ 37Rv_TOL_LpH_ 
data/GSE67035_H37Rv_Fe-DE.csv: src/makeDE.R data/GSE67035_H37Rv_raw_count_std.txt
	Rscript $< -i data/all_raw_count_std.txt -f data/all_pheno.txt -o $(basename $@) -conds 37Rv_HiFe_ 37Rv_LoFe_1D_ 37Rv_LoFe_1W_
data/GSE123267_H37Rv_pH-DE.csv: src/makeDE.R data/GSE123267_H37Rv_raw_count_std.txt
	Rscript $< -i data/all_raw_count_std.txt -f data/all_pheno.txt -o $(basename $@) -conds wt_6.6_ wt_4.5_
data/SRP142541_H37Rv_Fe-DE.csv: src/makeDE.R data/SRP142541_H37Rv_raw_count_std.txt
	Rscript $< -i $(word 2, $^) -f data/SRP142541_H37Rv_pheno.txt -o $(basename $@) -conds nt_lofe_minAtc nt_ox_minAtc nt_sds_minAtc
data/SRP142541_H37Rv_sds-DE.csv: src/makeDE.R data/SRP142541_H37Rv_raw_count_std.txt
	Rscript $< -i $(word 2, $^) -f data/SRP142541_H37Rv_pheno.txt -o $(basename $@) -conds nt_sds_minAtc nt_lofe_minAtc nt_ox_minAtc

pH_sig := data/GSE67035_H37Rv_pH-DE_37Rv_TOL_LpH_-vs-37Rv_TOL_.txt data/GSE123267_H37Rv_pH-DE_wt_4.5_-vs-wt_6.6_.txt \
	data/GSE8867_CDC1551_pH-DE_5.5-vs-7.txt data/GSE89106_CDC1551_pH-DE_DMSO5-vs-DMSO7.txt data/GSE63917_CDC1551_pH-DE_DMSO5.7-vs-DMSO7.txt
Fe_sig := data/SRP142541_H37Rv_Fe-DE_nt_ox_minAtc-vs-nt_lofe_minAtc.txt data/SRP142541_H37Rv_Fe-DE_nt_sds_minAtc-vs-nt_lofe_minAtc.txt \
	data/GSE67035_H37Rv_Fe-DE_37Rv_LoFe_1D_-vs-37Rv_HiFe_.txt data/GSE67035_H37Rv_Fe-DE_37Rv_LoFe_1W_-vs-37Rv_HiFe_.txt
sds_sig := data/SRP142541_H37Rv_sds-DE_nt_lofe_minAtc-vs-nt_sds_minAtc.txt data/SRP142541_H37Rv_sds-DE_nt_ox_minAtc-vs-nt_sds_minAtc.txt

data/consensus_pH-DE.csv: src/consensusGenes.py $(pH_sig)
	python $< -i $(pH_sig) -o $(basename $@)
data/consensus_Fe-DE.csv: src/consensusGenes.py $(Fe_sig)
	python $< -i $(Fe_sig) -o $(basename $@)
data/consensus_sds-DE.csv: src/consensusGenes.py #$(sds_sig)
	python $< -i $(sds_sig) -o $(basename $@)	

fig_var-genes := $(foreach n, $(fig_prefix), $(addprefix $n, _var-genes.pdf))
fig_stress-tx := $(foreach n, $(fig_prefix), $(addprefix $n, _stress-tx.pdf))
fig_PCA       := $(foreach n, $(fig_prefix), $(addprefix $n, _PCA.pdf))
fig_corr      := $(foreach n, $(fig_prefix), $(addprefix $n, _corr.pdf))
fig_DE        := fig/GSE67035_H37Rv_pH-DE.pdf

#plot txs of interest
fig: $(fig_var-genes) $(fig_stress-tx) $(fig_read-dist) $(fig_PCA) $(fig_corr) $(fig_DE) fig/all_var-genes.pdf fig/all_stress-tx.pdf fig/all_PCA.pdf fig/all-combat_var-genes.pdf fig/all-srp-combat_var-genes.pdf fig/all-combat_stress-tx.pdf fig/all-combat_PCA.pdf fig/all-combat_corr.pdf fig/all_corr.pdf \
     fig/ph_var-genes.pdf fig/ph_stress-tx.pdf fig/ph_PCA.pdf fig/ph_corr.pdf
fig/%_stress-tx.pdf: src/plotTx.py data/%_ln-filt-tpm.txt data/stress-tx.txt src/util.py 	
	python $< -i $(word 2, $^) -f $(word 3, $^) -o $@
fig/%_var-genes.pdf: src/plotTx.py data/%_ln-filt-tpm.txt data/%_var-genes.txt src/util.py
	python $< -i $(word 2, $^) -f $(word 3, $^) -o $@

consensus_fig: fig/all-combat_consensus_pH-DE.pdf fig/all-combat_consensus_pH-DE.pdf
fig/%_consensus_pH-DE.pdf: src/plotTx.py data/%_ln-filt-tpm.txt data/consensus_pH-DE.txt src/util.py
	python $< -i $(word 2, $^) -f $(word 3, $^) -o $@
#perform PCA
fig/%_PCA.pdf: src/performPCA.py data/%_ln-filt-tpm.txt data/%_var-genes.txt src/util.py
	python $< -i $(word 2, $^) -f $(word 3, $^) -o fig/$*_PCA
fig/all_PCA-all-genes.pdf: src/performPCA.py data/all_ln-filt-tpm.txt data/all-genes.txt src/util.py
	python $< -i $(word 2, $^) -f $(word 3, $^) -o fig/all_PCA-all-genes
fig/all_PCA-var-genes.pdf: src/performPCA.py data/all_ln-filt-tpm.txt data/all_var-genes.txt src/util.py
	python $< -i $(word 2, $^) -f $(word 3, $^) -o fig/all_PCA-var-genes
fig/all_PCA-GSE67035-var-genes.pdf: src/performPCA.py data/all_ln-filt-tpm.txt data/GSE67035_H37Rv_var-genes.txt src/util.py
	python $< -i $(word 2, $^) -f $(word 3, $^) -o fig/all_PCA-var-genes
fig/all_GSE67035-var-genes.pdf: src/plotTx.py data/all_ln-filt-tpm.txt data/GSE67035_H37Rv_var-genes.txt src/util.py
	python $< -i $(word 2, $^) -f $(word 3, $^) -o fig/all_GSE67035-var-genes
fig/ph_PCA-all-genes.pdf: src/performPCA.py data/ph_ln-filt-tpm.txt data/all-genes.txt src/util.py
	python $< -i $(word 2, $^) -f $(word 3, $^) -o $@
fig/ph_PCA-var-genes.pdf: src/performPCA.py data/ph_ln-filt-tpm.txt data/ph_var-genes.txt src/util.py
	python $< -i $(word 2, $^) -f $(word 3, $^) -o $@
fig/ph_PCA-GSE67035-var-genes.pdf: src/performPCA.py data/ph_ln-filt-tpm.txt data/GSE67035_H37Rv_var-genes.txt src/util.py
	python $< -i $(word 2, $^) -f $(word 3, $^) -o $@
fig/%_corr.pdf: src/corr.py data/%_ln-filt-tpm.txt data/%_var-genes.txt src/util.py
	python $< -i $(word 2, $^) -f $(word 3, $^) -o $@

fig/%-DE.pdf: src/plotVolcano.R data/%-DE-37Rv_TOL_LpH_vs_37Rv_TOL_.csv
	Rscript $< -i $(word 2, $^) -o $@
