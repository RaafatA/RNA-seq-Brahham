# Introduction to the Step 1 script----
# Step 1 Learning Objectives:
# 1 - Step 1 serves as your gateway to R scripts and, as such, you will learn the proper 'anatomy' for any R script.
# 2 - Learn how to install packages and load libraries into your R environment
# 3 - Understand the various file types that describe RNAseq data and how to import these files (e.g. kallisto read mapping data) into R
# 4 - Learn basic tools for annotation

# Notes:
# This script is organized into 'chunks' of code, and the final chunk (called 'the essentials') is a minimal representation of this script.

# load packages----
library(rhdf5) #provides functions for handling hdf5 file formats (kallisto outputs bootstraps in this format)
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl
library(EnsDb.Hsapiens.v86) #replace with your organism-specific database package
library(beepr) #just for fun

install.packages(c("rhdf5", "tidyverse", "tximport"))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ensembldb")

install.packages("ensembldb")



# read in your study design ----
#there are LOTS of ways to read data into R, but the readr package (from tidyverse) is one of the simplest
sample <- read_tsv("study_design.txt")


path <- file.path(sample$sample, "abundance.tsv") # set file paths to your mapped data
# now check to make sure this path is correct by seeing if the files exist
all(file.exists(path))
########################################################################################
# Annotation 

Tx <- transcripts(EnsDb.Hsapiens.v86, columns = c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")


# import Kallisto transcript counts into R using Tximport ----
# copy the abundance files to the working directory and rename so that each sample has a unique name


txi_gene <- tximport(path,
                     type = "kallisto",
                     tx2gene = Tx,
                     txOut = FALSE,
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE
                     )

beep(sound =2 )

txi_gene$counts







