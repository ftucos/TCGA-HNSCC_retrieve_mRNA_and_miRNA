#!/bin/bash
set ue -o pipefail

WD="/Volumes/TucosHDD/Bioinformatics/workspace/TCGA-HNSCC_retrieve_mRNA_and_miRNA"
cd $WD

mkdir -p $WD/data/RNAseq_counts
cd $WD/data/RNAseq_counts
gdc-client download -m  $WD/data/manifest/gdc_manifest.2023-01-18.RNAseq_counts.txt

mkdir -p $WD/data/mirnaRNAseq_quant
cd $WD/data/mirnaRNAseq_quant
gdc-client download -m  $WD/data/manifest/gdc_manifest.2023-01-18.miRNAseq_quantification.txt

mkdir -p $WD/data/mirnaRNAseq_isoform
cd $WD/data/mirnaRNAseq_isoform
gdc-client download -m  $WD/data/manifest/gdc_manifest.2023-01-18.miRNAseq_isoform.txt