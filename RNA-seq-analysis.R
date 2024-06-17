# Load libraries
library(DESeq2)
library(readr)

# Define file paths
normal_dir <- "C:/Users/user/Desktop/sCienTificaRticlE/normal"
tumor_dir <- "C:/Users/user/Desktop/sCienTificaRticlE/tumor"

# List all TSV files in the directories
normal_files <- list.files(normal_dir, full.names = TRUE, recursive = TRUE, pattern = "\\.tsv$")
tumor_files <- list.files(tumor_dir, full.names = TRUE, recursive = TRUE, pattern = "\\.tsv$")

# Combine into a single list
files <- c(normal_files, tumor_files)

# Create sample names and conditions
sample_names <- c(paste0("normal_", 1:length(normal_files)), paste0("tumor_", 1:length(tumor_files)))
conditions <- c(rep("normal", length(normal_files)), rep("tumor", length(tumor_files)))

# Create metadata dataframe
metadata <- data.frame(
  sampleName = sample_names,
  fileName = files,
  condition = conditions
)

# Function to read a file and skip initial metadata lines with correct column names
read_tsv_correctly <- function(file) {
  df <- readr::read_tsv(file, skip = 1, col_names = c("gene_id", "gene_name", "gene_type", "unstranded", "stranded_first", "stranded_second", "tpm_unstranded", "fpkm_unstranded", "fpkm_uq_unstranded"))
  return(df)
}

# Read all files
data_list <- lapply(metadata$fileName, read_tsv_correctly)

# Define a safe conversion function
safe_as_numeric <- function(x) {
  as.numeric(as.character(x))
}

# Extract counts and convert to numeric safely, then combine
counts_matrix <- do.call(cbind, lapply(data_list, function(df) safe_as_numeric(df$unstranded)))
colnames(counts_matrix) <- metadata$sampleName
rownames(counts_matrix) <- data_list[[1]]$gene_id

# Check for NA values
anyNA(counts_matrix)

# Replace NA values with zero
counts_matrix[is.na(counts_matrix)] <- 0

# Convert conditions to factors
metadata$condition <- as.factor(metadata$condition)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = metadata, design = ~condition)
dds <- DESeq(dds)

res <- results(dds)
significant_genes <- subset(res, padj < 0.05)
summary(res)



library(ggplot2)

# Convert results to a data frame
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Create a volcano plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.4) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Gene Expression", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  scale_color_manual(values = c("grey", "red")) +
  theme(legend.position = "none")

print(volcano_plot)




annotation_col <- as.data.frame(colData(dds)$condition)
colnames(annotation_col) <- "Condition"
print(annotation_col)




library(pheatmap)

# Ensure the conditions are correctly labeled in the metadata
metadata$condition <- factor(metadata$condition, levels = c("normal", "tumor"))

# Check the annotation data
annotation_col <- data.frame(Condition = metadata$condition)
rownames(annotation_col) <- metadata$sampleName

# Select top 20 significant genes and remove non-gene rows
top_genes <- head(rownames(significant_genes), 20)
top_genes <- top_genes[!grepl("N_", top_genes)]
mat <- assay(dds)[top_genes, ]
mat <- mat - rowMeans(mat)

# Ensure all values are non-negative
mat[mat < 0] <- 0

# Log-transform the data to reduce the range of values
mat_log <- log2(mat + 1)  # Adding 1 to avoid log(0)

# Replace any NaN values with zeros (if any)
mat_log[is.na(mat_log)] <- 0

# Create a heatmap with log-transformed data
pheatmap(mat_log, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE,
         annotation_col = annotation_col,
         annotation_colors = list(Condition = c(normal = "blue", tumor = "red")),
         fontsize_row = 5,
         fontsize_col = 6,
         height = 12,  # Adjust plot height
         width = 15,  # Adjust plot width
         cellheight = 6,  # Adjust cell height
         cellwidth = 4)   # Adjust cell width






library(pheatmap)

# Select top 20 significant genes and remove non-gene rows
top_genes <- head(rownames(significant_genes), 20)
top_genes <- top_genes[!grepl("N_", top_genes)]
mat <- assay(dds)[top_genes, ]
mat <- mat - rowMeans(mat)

# Ensure all values are non-negative
mat[mat < 0] <- 0

# Log-transform the data to reduce the range of values
mat_log <- log2(mat + 1)  # Adding 1 to avoid log(0)

# Replace any NaN values with zeros (if any)
mat_log[is.na(mat_log)] <- 0

# Create a heatmap with log-transformed data
pheatmap(mat_log, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE,
         annotation_col = annotation_col,
         annotation_colors = list(Condition = c(normal = "blue", tumor = "red")),
         fontsize_row = 4,  # Increase font size for better readability
         fontsize_col = 3,  # Increase font size for better readability
         height = 8,  # Increase plot height
         width = 40) 





library(DESeq2)
library(ggplot2)

vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot of RNA-seq Data") +
  theme_minimal()



summary_stats <- summary(res)
print(summary_stats)



# Convert results to data frame
results_table <- as.data.frame(res)

# Display the first few rows of the table
head(results_table)

# Save the results to a CSV file (optional)
write.csv(results_table, "deseq2_results.csv", row.names = TRUE)

getwd()











