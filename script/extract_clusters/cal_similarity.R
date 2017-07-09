#! /usr/bin/Rscript
library(seqinr)
library(gplots)

args <- commandArgs(trailingOnly = TRUE)

seqs <- read.alignment(args[1], format = "fasta")
amino_seq <- dist.alignment(seqs, matrix = "similarity")
similarity_matrix <- as.matrix(amino_seq)
correct_similarity_matrix <- 1-similarity_matrix
colnames(correct_similarity_matrix)
data_matrix <- correct_similarity_matrix[1:as.numeric(args[2]), (as.numeric(args[2])+1):as.numeric(args[3])]
data_matrix
colnames(data_matrix)
row.names(data_matrix)
png(args[4],  width = 8, height = 8, units = 'in', res = 600)
heatmap.2(data_matrix , Rowv = FALSE, Colv = FALSE, dendrogram = "none", cexRow  = 0.8, cexCol = 0.75,margins =c(17,15), keysize = 1, key.xlab="sequence similarity",density = "none", trace="none",col = colorpanel(256,"yellow","mediumseagreen","midnightblue"))
dev.off()




