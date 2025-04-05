# Part 1: Installation 

#Download devtools

install.packages("devtools")

install.packages("tibble")  

library(tibble)

# Download CIBERSORT.R and LM22.txt

dir.create("C:/Users/WA/Downloads/cibersort/CIBESORT.R", showWarnings = FALSE)

# Download R script from a public GitHub repo

download.file(
  
  url = "https://raw.githubusercontent.com/JingYangSciBio/Immu-Mela/master/CIBERSORT.R",
  
  destfile = "C:/Users/WA/Downloads/cibersort/CIBERSORT.R",
  
  mode = "wb"
  
)

# Download LM22.txt

download.file(
  
  url = "https://raw.githubusercontent.com/mdozmorov/Immuno_notes/master/data/Cibersoft/LM22.txt",
  
  destfile = "C:/Users/WA/Downloads/cibersort/LM22.txt",
  
  mode = "wb"
  
)

# set the paths in immunedeconv

library(immunedeconv)

set_cibersort_binary("C:/Users/WA/Downloads/cibersort/CIBERSORT.R")
set_cibersort_mat("C:/Users/WA/Downloads/cibersort/LM22.txt")


# Inside R

Sys.setenv(GITHUB_PAT = "ghp_personal token")

# Install Immunedeconv

devtools::install_local(".")

# Load library(immunedeconv)

library(immunedeconv)

# List the methods

names(deconvolution_methods)

# Part 2: Download dataset from GSE93385 and GSE87750


# Load GSE87750.raw.tar

setwd("C:/Users/WA/Downloads/GSE87750_RAW/")

# Part 3: Build a gene expression matrix from multiple sample-level FPKM files

library(dplyr)

# List all .fpkm_tracking files
files <- list.files("C:/Users/WA/Downloads/GSE87750_RAW/")
files
############################Define a function to extract FPKM per sample#####################

read_fpkm <- function(file) {
  
  dat <- read.delim(file)
  sample_name <- tools::file_path_sans_ext(basename(file))
  
  sample_name <- sub(".fpkm_tracking", "", sample_name)
  
  dat <- dat[, c("gene_id", "FPKM")]
  
  dat <- dat[!is.na(dat$gene_id) & dat$gene_id != "", ]
  
  dat <- dat[!duplicated(dat$gene_id), ]
  
  tibble(gene = dat$gene_id, !!sample_name := dat$FPKM)
  
}


##################################Read and merge all samples into one matrix#######################

expr_list <- lapply(files, read_fpkm)

expr_matrix <- Reduce(function(x, y) full_join(x, y, by = "gene"), expr_list) %>%
  
  filter(!is.na(gene) & gene != "") %>%
  
  distinct(gene, .keep_all = TRUE) %>%
  
  as.data.frame()



rownames(expr_matrix) <- expr_matrix$gene


expr_matrix$gene <- NULL

head(rownames(expr_matrix), 10)


##############################################Normalization########################################################

install.packages("preprocessCore")

library(preprocessCore)

# Convert to matrix
mat <- as.matrix(expr_matrix)

# Quantile normalize
mat_norm <- normalize.quantiles(mat)

# Preserve row and column names
rownames(mat_norm) <- rownames(expr_matrix)
colnames(mat_norm) <- colnames(expr_matrix)

# Convert back to data.frame if needed
expr_matrix_normalized <- as.data.frame(mat_norm)

# View a snippet
head(expr_matrix_normalized[, 1:5])

##################### comparison between the expression distribution before and after normalization
boxplot(expr_matrix, 
        las = 2, 
        main = "Expression Distributions per Sample", 
        outline = FALSE)


boxplot(expr_matrix_normalized, 
        las = 2, 
        main = "Expression Distributions per Normalized Sample", 
        outline = FALSE)
############################################################ mapping #######################################3
# Install biomaRt if you don't have it
install.packages("biomaRt")

# Load biomaRt library
library(biomaRt)

# Connect to Ensembl database (use human genes)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Map Ensembl IDs to HUGO Gene Symbols (if your gene IDs are Ensembl IDs)
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id", 
                  values = rownames(expr_matrix_normalized),
                  mart = ensembl)

# Merge the gene map with expr_matrix
expr_matrix$gene_id <- rownames(expr_matrix)

expr_matrix_with_symbols <- merge(expr_matrix, gene_map, by.x = "gene_id", by.y = "ensembl_gene_id", all.x = TRUE)



rownames(expr_matrix_with_symbols) <- expr_matrix_with_symbols$gene_id
expr_matrix_with_symbols$gene_id <- NULL  # Remove the temporary gene_id column

expr_matrix <- expr_matrix_with_symbols

########################################################## lm22 quality check ####################


# Reload original LM22 to get back the gene names
lm22 <- read.delim("C:/Users/WA/Downloads/cibersort/LM22_fixed.txt", check.names = FALSE)

# Optional: Rename first column to "GeneSymbol"
colnames(lm22)[1] <- "GeneSymbol" 

# Convert all other columns to numeric safely
lm22[, -1] <- sapply(lm22[, -1], as.numeric)

# Replace NA values with 0
lm22[is.na(lm22)] <- 0

# Sanity check – make sure GeneSymbol column still has text
head(lm22$GeneSymbol)


# Save it
write.table(lm22,
            file = "C:/Users/WA/Downloads/cibersort/LM22_fixed.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

########################################### Filter to LM22 genes ####################################
lm22_genes <- lm22[1]

sum(rownames(expr_matrix) %in% lm22_genes)


expr_matrix_filtered <- expr_matrix[rownames(expr_matrix) %in% lm22_genes, ]


# Ensure numeric matrix and clean formatting
expr_matrix_filtered <- expr_matrix_filtered[complete.cases(expr_matrix_filtered), ]
expr_matrix_filtered <- expr_matrix_filtered[!duplicated(rownames(expr_matrix_filtered)), ]
expr_matrix_filtered <- as.matrix(expr_matrix_filtered)
storage.mode(expr_matrix_filtered) <- "numeric"

# Write filtered expression matrix to a text file for CIBERSORT
write.table(expr_matrix_filtered,
            file = "C:/Users/WA/Downloads/cibersort/expr_for_cibersort.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

######################################### Run CIBERSORT ####################################

# Load the CIBERSORT R script
source("C:/Users/WA/Downloads/cibersort/CIBERSORT.R")

# Run CIBERSORT
results <- CIBERSORT(
  sig_matrix = "C:/Users/WA/Downloads/cibersort/LM22_fixed.txt",
  mixture_file = "C:/Users/WA/Downloads/cibersort/expr_for_cibersort.txt",
  perm = 100,
  QN = FALSE  # for RNA-seq
)

# View results
head(results)





