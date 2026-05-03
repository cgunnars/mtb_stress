library(GEOquery)
args = commandArgs(trailingOnly = TRUE)
gse_num = args[1]
print(gse_num)
getGEOSuppFiles(gse_num, baseDir='./data/')
gse <- getGEO(gse_num, destdir=paste0('./data/', gse_num))
