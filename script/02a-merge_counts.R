setwd("/Volumes/TucosHDD/Bioinformatics/workspace/TCGA-HNSCC_retrieve_mRNA_and_miRNA/")
library('data.table')
library('tidyverse')
library('miRBaseConverter')


# RNAseq unificataion -----------------------------------------------------------------------
# Samplesheet with filename to patient/sample-id matching
Samplesheet <- read_delim("data/metadata/gdc_sample_sheet.2023-01-18.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
# Generate a named vector for easily converting filename to sample-id 
filename_to_sampleID <- set_names(Samplesheet$`Sample ID`, Samplesheet$`File Name`)

# Parse RNAseq ------------------------------------------------
# Function for parsing read counts
parse_RNAseq_reads <- function(path) {
  # extract file name
  file_name <- str_extract(path, "(?<=/)[^/]*$")
  # replace file name with sample id
  sample_name <- filename_to_sampleID[file_name]
  
  counts <- fread(path) %>%
    # remove stats rows
    filter(!gene_id %in% c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")) %>%
    # remove ENSG version but preserve PAR_Y annotation (y chromosome allele)
    mutate(gene_id = str_remove(gene_id, "\\.[0-9_]+$")) %>%
    select(gene_id, gene_name, gene_type, unstranded) %>%
    # rename the counts column with the sample name
    rename(!!as.symbol(sample_name) := unstranded)
}

# List files to parse
RNAseq_files <- list.files("data/RNAseq_counts/", pattern = "star_gene_counts.tsv$", recursive = T, full.names = T)

# parse all read counts and join them in a single table
RNAseq <- map(RNAseq_files, parse_RNAseq_reads) %>%
  reduce(full_join, by=c('gene_id', 'gene_name', 'gene_type'))

# export
fwrite(RNAseq, file = "result/TCGA-HNSC.RNAseq.counts.tsv", quote = F, sep = "\t")


# Parse miRNAseq -----------------------------------------------------------------------

# Function for parsing read counts
parse_miRNA_reads <- function(path) {
  # extract file name
  file_name <- str_extract(path, "(?<=/)[^/]*$")
  # replace file name with sample id
  sample_name <- filename_to_sampleID[file_name]
  counts <- fread(path) %>%
    select(miRNA = miRNA_ID, read_count) %>%
    # rename the counts column with the sample name
    rename(!!as.symbol(sample_name) := read_count)
}

# List files to parse
miRNAseq_files <- list.files("data/mirnaRNAseq_quant", pattern = "mirnas.quantification.txt$", recursive = T, full.names = T)

# parse all read counts and join them in a single table
miRNAseq <- map(miRNAseq_files, parse_miRNA_reads) %>%
  reduce(full_join, by='miRNA')

# export
fwrite(miRNAseq, file = "result/TCGA-HNSC.miRNAseq.counts.tsv", quote = F, sep = "\t")


# Parse miRNAseq isoform -----------------------------------------------------------------------

mirbase_version <- "v21"

# Function for parsing read counts
parse_miRNA_isoform_reads <- function(path) {
  # extract file name
  file_name <- str_extract(path, "(?<=/)[^/]*$")
  # replace file name with sample id
  sample_name <- filename_to_sampleID[file_name]
  
  counts <- fread(path, data.table = F) %>%
    filter(str_detect(miRNA_region, "mature")) %>%
    # extract MIMAT ID
    mutate(MIMAT_ID = str_extract_all(miRNA_region, "MIMAT.*")) %>%
    group_by(MIMAT_ID) %>% 
    # sum counts for different isoform mapping to the same MIMAT ID
    summarize(read_count = sum(read_count)) %>%
    # rename the counts column with the sample name
    rename(!!as.symbol(sample_name) := read_count)
}

# List files to parse
miRNAseq_isoform_files <- list.files("data/mirnaRNAseq_isoform", pattern = "isoforms.quantification.txt$", recursive = T, full.names = T)

# parse all read counts and join them in a single table
miRNAseq_isoform <- map(miRNAseq_isoform_files, parse_miRNA_isoform_reads) %>%
  reduce(full_join, by='MIMAT_ID')

# convert MIMAT ID to stranded miRNA name
miRNAseq_isoform["miRNA"] <- miRNA_AccessionToName(miRNAseq_isoform$MIMAT_ID,
                                                   targetVersion = mirbase_version)$TargetName

# replace NAs with 0 counts
miRNAseq_isoform[is.na(miRNAseq_isoform)] <- 0

# reorder columns
miRNAseq_isoform <- miRNAseq_isoform %>% select(MIMAT_ID, miRNA, everything())

# export
fwrite(miRNAseq_isoform, file = "result/TCGA-HNSC.miRNAseq_isoform.counts.tsv", quote = F, sep = "\t")
