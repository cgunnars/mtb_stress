library('DESeq2')
library(argparse)
library(tidyverse)

padj_cutoff = 0.05
log2FC_cutoff = 1.5

parser <- ArgumentParser()
parser$add_argument('-i', type="character", nargs=1)
parser$add_argument('-f', type="character", nargs=1)
parser$add_argument('-o', type="character", nargs=1)
parser$add_argument('-conds', type="character", nargs='+')
parser$add_argument('--down', action='store_true')
args <- parser$parse_args()


conds = args$conds

pheno = read.delim(args$f, stringsAsFactors=TRUE)
pheno$batch <- factor(pheno$batch) #convert to factor for compatibility
pheno_cond <- subset(pheno, condition %in% conds)
dat = read.delim(args$i, check.names=FALSE)
rownames(dat) <- dat[,1]
dat_cond <- select(dat, pheno_cond$name)

ddsFullCountTable <- DESeqDataSetFromMatrix(countData = dat_cond, colData = pheno_cond, design= ~ condition)
keep <- rowSums(counts(ddsFullCountTable)) >= 10 #magic number rn
dds <- ddsFullCountTable[keep,]
dds$condition <- relevel(dds$condition, ref=conds[1])

dds <- DESeq(dds)
for (cond in conds[-1]) {
	res <- results(dds, alpha = padj_cutoff, contrast=c('condition', cond, conds[1]))
	res_sig <- res %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble() %>%	
		            filter(padj < padj_cutoff)
	if (args$down) {
	       res_sig <- res_sig %>% filter(log2FoldChange < -log2FC_cutoff)
	} else {
		res_sig <- res_sig %>% filter(log2FoldChange > log2FC_cutoff)
	}
	res_sig_file <- file(paste0(args$o, '_', cond, '-vs-', conds[1], '.txt'))
	writeLines(res_sig$gene, res_sig_file)
	close(res_sig_file)	
	write.csv(as.data.frame(res), file=paste0(args$o, '-', cond, 'vs_', conds[1], '.csv'))
	
}

