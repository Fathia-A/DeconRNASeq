---
title: "README.md"
output: html_document
date: "2025-04-05"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## README

# DeconRNASeq: Zika Virus Transcriptomic Deconvolution

This repository contains R scripts and resources for performing transcriptomic deconvolution analysis on RNA-seq data from Zika virus (ZIKV)-infected cells. The aim is to disentangle mixed cellular signals and better understand the impact of ZIKV infection on gene expression at the cell-type level.

## 🧪 Project Overview

**Main Script:** `zikv.R`  
This script performs the following:
- Preprocessing and normalization of RNA-seq count data
- Deconvolution of bulk RNA-seq using reference cell type signatures
- Visualization of cell-type proportions and transcriptional changes

## 📁 Repository Structure

. ├── zikv.R # Main R script for deconvolution ├── data/ # Folder for input RNA-seq and reference data ├── output/ # Folder for saving plots and results └── README.md # Project description and usage guide


## 📊 Requirements

Make sure you have the following R packages installed:
```r
install.packages(c("DeconRNASeq", "ggplot2", "dplyr", "readr"))

## 🚀 How to Run

#1-Clone the repository:

git clone https://github.com/Fathia-A/DeconRNASeq.git
cd DeconRNASeq

#2-Open zikv.R in RStudio or run from the terminal:

source("zikv.R")




















