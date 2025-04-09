# -------------------------------------------------------------------
# Part 1: Install Required Libraries
# -------------------------------------------------------------------
# Check if BiocManager is installed, and install it if necessary
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install necessary libraries from Bioconductor
BiocManager::install("MuSiC")
BiocManager::install("SingleCellExperiment")
BiocManager::install("Biobase")
BiocManager::install("SummarizedExperiment")

# Load MuSiC library after installation
library(MuSiC)


# -------------------------------------------------------------------
# Part 2: Load and Prepare Data
# -------------------------------------------------------------------

## A) Load Bulk RNA-seq Matrix (GSE125554 dataset)
# Read the CSV file containing the bulk RNA-seq expression data
expr <- read.csv("GSE125554_zikv_counts.csv", row.names = 1)

# Check the first few rows to ensure correct loading
head(expr[, 1:4])
dim(expr)  # Check dimensions

# -------------------------------------------------------------------
## B) Load Brain Single-Cell RNA-seq Reference (GSE104276 dataset)
# -------------------------------------------------------------------
# Load necessary libraries
library(SingleCellExperiment)
library(readxl)

# Step 1: Load TPM Matrix from .xls File
tpm_file <- "GSE104276_all_pfc_2394_UMI_TPM_NOERCC.xls"

# Read the .xls file (cell x gene expression matrix)
expr_matrix <- read.delim(tpm_file, row.names = 1, check.names = FALSE)

# Convert to matrix format
expr_matrix <- as.matrix(expr_matrix)

# Check dimensions and headers
dim(expr_matrix)
head(colnames(expr_matrix))
head(rownames(expr_matrix))

# Step 2: Load Metadata for Single-Cell Data
meta_data <- read_xlsx("GSE104276_readme_sample_barcode.xlsx", sheet = "CellNameBarcodeID")

# View the metadata columns
head(meta_data)
colnames(meta_data)

# Reshape metadata from wide to long format (cell names only)
library(tidyverse)

cell_name_cols <- meta_data %>% select(-starts_with("Barcode_ID"))

meta_long <- cell_name_cols %>%
  pivot_longer(cols = everything(),
               names_to = "donor_id",
               values_to = "cell_name") %>%
  filter(!is.na(cell_name))

# Step 3: Add Cell Type Labels
meta_types <- read_xlsx("GSE104276_readme_sample_barcode.xlsx", sheet = "PFC_NPCneuglia_monocle_ident")
head(meta_types)
colnames(meta_types)

# Rename ...1 column to cell_name for merging
meta_types <- meta_types %>%
  rename(cell_name = `...1`)

# Merge metadata (donor and cell type info)
metadata_final <- meta_long %>%
  left_join(meta_types, by = "cell_name")

# Check the cell types in the merged data
table(metadata_final$cell_type)

# Clean metadata: Filter and remove NAs
metadata_final <- metadata_final %>%
  filter(cell_name %in% colnames(expr_matrix)) %>%
  distinct(cell_name, .keep_all = TRUE) %>%
  filter(!is.na(cell_name)) %>%
  arrange(match(cell_name, colnames(expr_matrix)))

# Final check to ensure no mismatch between cell names and expression matrix columns
stopifnot(all(metadata_final$cell_name == colnames(expr_matrix)))

# -------------------------------------------------------------------
# Part 3: Build a SingleCellExperiment Object
# -------------------------------------------------------------------
# Create SingleCellExperiment (SCE) object
sce <- SingleCellExperiment(
  assays = list(counts = expr_matrix),
  colData = data.frame(
    cell_type = metadata_final$cell_type,
    donor_id = metadata_final$donor_id,
    row.names = metadata_final$cell_name
  )
)

# Save the SingleCellExperiment object to disk
saveRDS(sce, "GSE104276_sce.rds")


# -------------------------------------------------------------------
# Part 4: Deconvolution using MuSiC
# -------------------------------------------------------------------
# Load the SCE object for MuSiC
sce <- readRDS("GSE104276_sce.rds")

# Load Bulk RNA-seq Expression Data
expr <- read.csv("GSE125554_zikv_counts.csv", row.names = 1)

# Intersect genes between single-cell and bulk datasets
common_genes <- intersect(rownames(sce), rownames(expr))

# Filter both matrices to include only common genes
expr <- expr[common_genes, ]
sce <- sce[common_genes, ]

# -------------------------------------------------------------------
# Step 1: Convert Ensembl IDs to Gene Symbols (if needed)
# -------------------------------------------------------------------
# If your bulk dataset uses Ensembl IDs while the reference uses gene symbols,
# we need to map the Ensembl IDs to gene symbols using the biomaRt package.

library(biomaRt)

# Connect to Ensembl database
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve mapping between Ensembl IDs and gene symbols
id_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(expr),
  mart = mart
)

# Clean up the mapping: Remove duplicates and empty gene symbols
id_mapping <- id_mapping[!duplicated(id_mapping$ensembl_gene_id), ]
id_mapping <- id_mapping[id_mapping$hgnc_symbol != "", ]

# Add gene symbols to the expr matrix
expr$gene_symbol <- id_mapping$hgnc_symbol[match(rownames(expr), id_mapping$ensembl_gene_id)]

# Filter out rows with missing or invalid gene symbols
expr <- expr[!is.na(expr$gene_symbol) & expr$gene_symbol != "", ]

# Collapse duplicates by summing counts for the same gene symbol
expr <- expr %>%
  group_by(gene_symbol) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  as.data.frame()

# Set the gene symbols as rownames and remove the gene_symbol column
rownames(expr) <- expr$gene_symbol
expr$gene_symbol <- NULL


# -------------------------------------------------------------------
# Part 5: Run MuSiC for Deconvolution
# -------------------------------------------------------------------
# Load MuSiC after installation
library(MuSiC)

# Rename the assay to "counts" for consistency
assayNames(sce) <- "counts"

# Run MuSiC Deconvolution
music_results <- music_prop(
  bulk.mtx = as.matrix(expr),
  sc.sce = sce,
  clusters = 'cell_type',  # Column containing cell types
  samples = 'donor_id',    # Column containing sample/donor IDs
  verbose = TRUE
)

# Inspect the deconvolution results (cell type proportions)
head(music_results$Est.prop.weighted)
str(music_results$Est.prop.weighted)


# -------------------------------------------------------------------
# Part 6: Data Visualization
# -------------------------------------------------------------------
# Replace NAs in cell type proportions with "Unknown"
colnames(music_results$Est.prop.weighted)[is.na(colnames(music_results$Est.prop.weighted))] <- "Unknown"

# Convert results into a long format for plotting
library(ggplot2)
df <- as.data.frame(music_results$Est.prop.weighted)
df$Sample <- rownames(df)

df_long <- pivot_longer(df, -Sample, names_to = "CellType", values_to = "Proportion")

# Create bar plot of estimated cell type proportions
ggplot(df_long, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Estimated Cell Type Proportions from MuSiC")
