# Load required libraries
library(DESeq2)
library(dplyr)

# Step 1: Read the expression matrix
exprs <- read.csv("reordered.csv", header = TRUE, check.names = FALSE)

# Step 2: Group by gene_name and average duplicate genes
deduplicated_exprs <- exprs %>%
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

# Step 3: Convert to data frame and set rownames
counts <- as.data.frame(deduplicated_exprs)
rownames(counts) <- counts$gene_name
counts$gene_name <- NULL

# Step 4: Read metadata
coldata <- read.table("metadata.txt", header = TRUE, row.names = 1)

# Step 5: Ensure counts columns match metadata rows
counts <- counts[, rownames(coldata)]

# Step 6: Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = coldata,
                              design = ~ Condition)

# Step 7: Prefilter
dds <- dds[rowSums(counts(dds)) > 1, ]

# Step 8: Run DESeq
dds <- DESeq(dds)

# Step 9: Extract results
res <- results(dds)
res <- lfcShrink(dds, coef = "Condition_vo_vs_sham", type = "apeglm")

# Step 10: Save results
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Add threshold label
res_df$threshold <- as.factor(
  ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
         ifelse(res_df$log2FoldChange > 1, "Up", "Down"),
         "NotSig")
)

write.csv(res_df, file = "sig_DESeq2_results.csv")

# Volcano Plot
library(ggplot2)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Up" = "blue", "NotSig" = "grey", "Down" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 adjusted p-value") +
  theme(legend.title = element_blank())
