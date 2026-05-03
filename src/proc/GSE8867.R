source('src/util.R')
library(tidyverse)
args <- io_sheet()
df <- lapply(args$sheet, function(x) read_sheet(args$i, x, header_row = 4))

#df has duplicate columns, first represent pH_5.5, second are pH_6.5
df <- lapply(df, function (x) { colnames(x) <- c('Gene_name', 'Log2_Fold_Change_pH_5.5', 'Adjusted_p_value_pH_5.5', 'Log2_Fold_Change_pH_6.5', 'Adjusted_p_value_pH_6.5', 'Function'); x})

#sheets are Up_6.5, Down_6.5, Up_5.5, Down_5.5
df_6.5 <- rbind(df[[1]], df[[2]])
df_5.5 <- rbind(df[[3]], df[[4]])

write.csv(df_6.5, file=paste0(args$o, '_6.5-vs-7.csv'))
write.csv(df_5.5, file=paste0(args$o, '_5.5-vs-7.csv'))

sig_file_6.5 <- file(paste0(args$o, '_6.5-vs-7.txt'))
sig_file_5.5 <- file(paste0(args$o, '_5.5-vs-7.txt'))

sig_6.5 <- get_sig(df_6.5, 'Adjusted_p_value_pH_6.5', 'Log2_Fold_Change_pH_6.5')
sig_5.5 <- get_sig(df_5.5, 'Adjusted_p_value_pH_5.5', 'Log2_Fold_Change_pH_5.5')

writeLines(sig_6.5$Gene_name, sig_file_6.5)
writeLines(sig_5.5$Gene_name, sig_file_5.5)

close(sig_file_6.5)
close(sig_file_5.5)

