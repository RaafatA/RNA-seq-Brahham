# Load packages -----
library(tidyverse) # already know about this from Step 1 script
library(matrixStats) # let's us easily calculate stats on rows or columns of a data matrix
library(cowplot) # allows you to combine multiple plots in one figure



# Examine your data up to this point ----
myTPM <- txi_gene$abundance
myCounts <- txi_gene$counts
colSums(myTPM)
colSums(myCounts)


# Generate summary stats for your data ----
# 1st, calculate summary stats for each transcript or gene, and add these to your data matrix
# then use the base R function 'transform' to modify the data matrix (equivalent of Excel's '=')
# then we use the 'rowSds', 'rowMeans' and 'rowMedians' functions from the matrixStats package
myTPM.stats <- transform(myTPM, 
                         SD=rowSds(myTPM), 
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM))

# look at what you created
head(myTPM.stats)

# Create your first plot using ggplot2 ----
# produce a scatter plot of the transformed data
ggplot(myTPM.stats) + 
  aes(x = SD, y = MED) +
  geom_point( size=3)
# Experiment with point shape and size in the plot above
# Experiment with other plot types (e.g. 'geom_hex' instead of 'geom_point')
# Add a theme to your ggplot code above.  Try 'theme_bw()'
# How would these graphs change if you log2 converted the data?

# Let's expand on the plot above a bit more and take a look at each 'layer' of the ggplot code
ggplot(myTPM.stats) + 
  aes(x = SD, y = MED) +
  geom_point(shape=16, size=2) +
  geom_smooth(method=lm) +
  geom_hex(show.legend = FALSE) +
  labs(y="Median", x = "Standard deviation",
       title="Transcripts per million (TPM)",
       subtitle="unfiltered, non-normalized data",) +
  theme_classic() +
  theme_dark() + 
  theme_bw()

sampleLabels <- sample$sample


colnames(myCounts) <- sample$sample






library(DESeq2)
# Assuming you have count data and sample metadata
countData <- round(myCounts)

sample$group <- as.factor(sample$group)

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = sample,
                              design = ~ group) # Replace with 



# Initial transformation for visualization (equivalent to log2 CPM)
rld <- rlog(dds, blind=TRUE)
rlog_counts <- assay(rld)

# Convert to tibble for plotting
rlog_counts.df <- as_tibble(rlog_counts, rownames = "geneID")
# colnames(rlog_counts.df) <- c("geneID", sampleLabels)
rlog_counts.df.pivot <- pivot_longer(rlog_counts.df,
                                     cols = -1,
                                     names_to = "samples",
                                     values_to = "expression")

p1 <- ggplot(rlog_counts.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="rlog expression", x = "sample",
       title="Regularized Log Transformation",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# Pre-filtering using DESeq2 approach
dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds, normalized=TRUE) >= 5) >= 3
dds_filtered <- dds[keep,]

# Alternative: more stringent filtering equivalent to original
# keep <- rowSums(counts(dds, normalized=TRUE) >= 5) >= 3
# dds_filtered <- dds[keep,]

# Apply rlog transformation to filtered data
rld_filtered <- rlog(dds_filtered, blind=TRUE)
rlog_filtered_counts <- assay(rld_filtered)

rlog_filtered_counts.df <- as_tibble(rlog_filtered_counts, rownames = "geneID")
colnames(rlog_filtered_counts.df) <- c("geneID", sampleLabels)
rlog_filtered_counts.df.pivot <- pivot_longer(rlog_filtered_counts.df,
                                              cols = -1,
                                              names_to = "samples",
                                              values_to = "expression")

p2 <- ggplot(rlog_filtered_counts.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="rlog expression", x = "sample",
       title="Regularized Log Transformation",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()



# Run full DESeq2 pipeline (includes normalization via size factors)
dds_filtered <- DESeq(dds_filtered)

# Apply rlog transformation after full pipeline
rld_norm <- rlog(dds_filtered, blind=FALSE)
rlog_norm_counts <- assay(rld_norm)

rlog_norm_counts.df <- as_tibble(rlog_norm_counts, rownames = "geneID")
colnames(rlog_norm_counts.df) <- c("geneID", sampleLabels)
rlog_norm_counts.df.pivot <- pivot_longer(rlog_norm_counts.df,
                                          cols = -1,
                                          names_to = "samples",
                                          values_to = "expression")

p3 <- ggplot(rlog_norm_counts.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="rlog expression", x = "sample",
       title="Regularized Log Transformation",
       subtitle="filtered, DESeq2 normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# Create combined plot
plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)


#-------------------------------------------------------------------------------------------------------

# dIFFERENTIAL GENE EXPRESSION


dds_filtered <- DESeq(dds_filtered)


res <- results(dds_filtered, contrast = c("group", "Tumor", "Normal"))
res

DESeq2::plotMA(res,main = "DESeq2 MA-Plot: Tumor vs Normal")


# Apply the variance-stabilizing transform instead of rlog:
vsd <- varianceStabilizingTransformation(dds_filtered, blind = FALSE)

# Extract the transformed matrix for downstream analyses (PCA, clustering, etc.):
plotPCA(vsd, intgroup = c("group"))


# Select top 30 TRUE         AND     TRUE 
resSig <- res[res$padj < 0.05 & res$log2FoldChange < -2, ]






mat <- assay(vsd)
top30 <- rownames(head(resSig[order(-abs(resSig$log2FoldChange)), ], 30))
mat_top <- mat[top30, ]



# 3c. Prepare annotation for columns
annotation_col <- as.data.frame(colData(vsd)[, "group", drop=FALSE])
library(pheatmap)
# 3d. Draw heatmap
p_heat <- pheatmap(mat_top,
                   scale          = "row",
                   annotation_col = annotation_col,
                   show_rownames  = TRUE,
                   show_colnames      = FALSE,    # <� turn off sample names
                   fontsize_row   = 8,
                   clustering_method = "complete",
                   main           = "Top 30 DE Genes")




# Remove NA values
res <- na.omit(res)

# Create a data frame for plotting
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)  # Ensure gene names are stored correctly
res_df$significant <- "Not Significant"
res_df$significant[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1] <- "Significant"

# Load ggplot2 for plotting
library(ggplot2)


# Create volcano plot with gene labels for significant genes
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text(data = subset(res_df, significant == "Significant"), 
            aes(label = gene), size = 2, vjust = 1.5, hjust = 0.5, check_overlap = TRUE, inherit.aes = FALSE) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Tumor vs Normal",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")



# Create volcano plot with gene labels
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant, label = gene)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text(data = subset(res_df, significant == "Significant"), 
            aes(label = gene), size = 2, vjust = 1.5, hjust = 0.5, check_overlap = TRUE) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Tumor vs Normal",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")
