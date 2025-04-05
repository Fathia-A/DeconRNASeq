# Part 1: Installation 

#Download devtools

install.packages("devtools")

install.packages("tibble")  

install.

library(tibble)

# Download CIBERSORT.R and LM22.txt

dir.create("~/Documents/cibersort", showWarnings = FALSE)

# Download R script from a public GitHub repo

download.file(
  
  url = "https://raw.githubusercontent.com/JingYangSciBio/Immu-Mela/master/CIBERSORT.R",
  
  destfile = "~/Documents/cibersort/CIBERSORT.R",
  
  mode = "wb"
  
)

# Download LM22.txt

download.file(
  
  url = "https://raw.githubusercontent.com/mdozmorov/Immuno_notes/master/data/Cibersoft/LM22.txt",
  
  destfile = "~/Documents/cibersort/LM22.txt",
  
  mode = "wb"
  
)

# set the paths in immunedeconv

library(immunedeconv)

set_cibersort_binary("~/Documents/cibersort/CIBERSORT.R")

set_cibersort_mat("~/Documents/cibersort/LM22.txt")



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
# Define a function to extract FPKM per sample

read_fpkm <- function(file) {
  
  dat <- read.delim(file)
  
  sample_name <- tools::file_path_sans_ext(basename(file))
  
  sample_name <- sub(".fpkm_tracking", "", sample_name)
  
  dat <- dat[, c("gene_id", "FPKM")]
  
  dat <- dat[!is.na(dat$gene_id) & dat$gene_id != "", ]
  
  dat <- dat[!duplicated(dat$gene_id), ]
  
  tibble(gene = dat$gene_id, !!sample_name := dat$FPKM)
  
}


# Read and merge all samples into one matrix

expr_list <- lapply(files, read_fpkm)

expr_matrix <- Reduce(function(x, y) full_join(x, y, by = "gene"), expr_list) %>%
  
  filter(!is.na(gene) & gene != "") %>%
  
  distinct(gene, .keep_all = TRUE) %>%
  
  as.data.frame()

rownames(expr_matrix) <- expr_matrix$gene

expr_matrix$gene <- NULL

head(rownames(expr_matrix), 10)

# Map Ensembl IDs to HUGO symbols

if (!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt")

library(biomaRt)

mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

id_map <- getBM(
  
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  
  filters = "ensembl_gene_id",
  
  values = rownames(expr_matrix),
  
  mart = mart
  
)

# Clean mapping

id_map <- id_map[id_map$hgnc_symbol != "", ]

id_map <- id_map[!duplicated(id_map$ensembl_gene_id), ]

# Apply mapped symbols to matrix

expr_matrix <- expr_matrix[rownames(expr_matrix) %in% id_map$ensembl_gene_id, ]

rownames(expr_matrix) <- id_map$hgnc_symbol[match(rownames(expr_matrix), id_map$ensembl_gene_id)]

# Filter to LM22 genes

lm22 <- read.delim("~/Documents/cibersort/LM22.txt", check.names = FALSE)

lm22_genes <- rownames(lm22)

expr_matrix_filtered <- expr_matrix[rownames(expr_matrix) %in% lm22_genes, ]

# Run Cibersort

res <- deconvolute(as.matrix(expr_matrix_filtered), method = "cibersort")

head(res)

This is what I have but I couldn't run the res <- deconvolute(as.matrix(expr_matrix_filtered), method = "cibersort") because it keeps giving me Error in if (dups > 0) { : argument is of length zero
 