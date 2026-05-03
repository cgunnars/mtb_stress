library(ggplot2)
library(readxl)
library(argparse)
library(ggrepel)
source('../mcsf-gmcsf/src/plotutils.R')
source('src/util.R')

parser <-ArgumentParser()
parser$add_argument('-i', type='character', nargs=1)
parser$add_argument('-o', type='character', nargs=1)
args <- parser$parse_args()

res <- read.csv(args$i, stringsAsFactors=FALSE)
cdc_sig_file <- file('data/GSE89106_CDC1551_pH-DE_DMSO5-vs-DMSO7.txt', open='r')
cdc_sig <- readLines(cdc_sig_file)

res_sig <- subset(res, abs(log2FoldChange) > 1.5 & padj < 0.05)

head(res)
volcano <- volcano_helper(ggplot(res, aes(log2FoldChange, -log10(padj)))) + 
	   geom_point(data = res_sig, color = 'deepskyblue1', size = 0.5, stroke = 0, shape = 16) +
   	   geom_point(data = subset(res_sig, X %in% cdc_sig), color = 'deeppink1', size = 0.5, stroke = 0, shape = 16) +
	   geom_text_repel(size = (6 * 3 / 14), min.segment.length = 0.1, segment.size = 0.1 * 5 /14,
		           box.padding = unit(3, 'pt'), force = 100, 
			   aes(x = log2FoldChange, y = -log10(padj), label = ifelse(X %in% c('Rv2390c', 'Rv0575c', 'Rv3054c'), res$X, "")))

ggsave(args$o, plot=volcano, height=1.7, width=1.7, units='in')
