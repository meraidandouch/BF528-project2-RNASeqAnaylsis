# Load libraries
library(tidyverse)

# Load in data
filepath <- "/projectnb/bf528/users/hedgehog_2022/project2/Programmer-Project2/cuffdiff_out/gene_exp.diff"
diff_data <- read_tsv(filepath)

# Order data by q_value ascending
diff_data <- diff_data[order(diff_data$q_value),]

# Take top ten differentially expressed genes by q value
top10_diff <- diff_data[1:10, c(3, 8, 9, 10, 12, 13)]

# Calculate the number of genes with the minimum q-value
min_q_genes <- nrow(subset(diff_data, q_value == 0.00106929))
print(paste("the number of genes with the minimum q-value is:", min_q_genes))

# Parameter that allows two plots to be displayed on the same image
par(mfrow=c(1,2))

# Create histogram of log2.foldchange for all genes
hist(diff_data$`log2(fold_change)`, breaks = 40, xlim = c(-6, 6),
     col = 'skyblue', xlab = 'log2(fold change)',
     main = "log2(fold change)\n for All Genes")
# Adds figure label
mtext("(A)", font = 2, side = 3, adj = -0.35)

# Subset data that is significant
sig_check <- which(diff_data$significant == "yes")
sig_diff_data <- diff_data[sig_check,]

# Create histogram of log2.foldchange for all genes
hist(sig_diff_data$`log2(fold_change)`, breaks = 40, xlim = c(-6, 6),
     col = 'skyblue', xlab = 'log2(fold change)',
     main = "log2(fold change)\n for Significant Genes")
# Adds figure label
mtext("(B)", font = 2, side = 3, adj = -0.35)

# Separate down-regulated and up-regulated genes
negative_sig_diff <- subset(sig_diff_data, `log2(fold_change)` < 0)
positive_sig_diff <- subset(sig_diff_data, `log2(fold_change)` > 0)

# Output number of up and down regulated genes
downregulated <- nrow(negative_sig_diff)
print(paste("Number of down-regulated genes: ", downregulated, sep = ""))
upregulated <- nrow(positive_sig_diff)
print(paste("Number of up-regulated genes: ", upregulated, sep = ""))

# Extract gene names
downregulated_genes <- negative_sig_diff$gene
upregulated_genes <- positive_sig_diff$gene

# Looking for presence of specific genes
# Not required but used in the report, some manual operation of the
# following lines is required
"Stat3" %in% sig_diff_data$gene
"Stat3" %in% upregulated_genes # Stat3 significant + upregulated

"Postn" %in% sig_diff_data$gene 
"Postn" %in% downregulated_genes # Postn significant + downregulated

"Il13" %in% diff_data$gene # present in original set
"Il13" %in% sig_diff_data$gene # not significant

# Check if user wants to write tables to file
write_check <- readline(prompt="Write tables to file in hedgehog project 2 folder Y/N?")
if (write_check == "Y" | write_check == "y") {
  write.table(downregulated_genes, "/projectnb/bf528/users/hedgehog_2022/project2/Analyst-Project2/downregulated_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(upregulated_genes, "/projectnb/bf528/users/hedgehog_2022/project2/Analyst-Project2/upregulated_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  readr::write_csv(top10_diff, "/projectnb/bf528/users/hedgehog_2022/project2/Analyst-Project2/top10_genes.csv")
}
