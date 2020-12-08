# Mise en place des librairies

# Supervised

# RNAseq

## DeconvSeq
library(deconvSeq)

## DeconRNASeq
library(DeconRNASeq)

## CellMix
# ?

## MuSic
library(MuSiC)

# methylation arrays

## EpiDISH
library(EpiDISH)

# BOTH

## DIABLO
library('mixOmics')

# Unsupervised

# RNA-seq

## deconica
library(deconica)

## fastICA
library(fastICA)

# Methylation arrays

## RefFreeEWAS
library(RefFreeEWAS)

## medepir
library(medepir)

## DecompPipeline
library(DecompPipeline)

## MeDeCom
#???

##  FactorViz
library(FactorViz)

# ANY

## Decoder
library(decoder)

# BOTH

## intNMF
library(IntNMF)

# Packages annexes nécessaires
library("DESeq2")

# Mise en place des donnees

# Localisation Data
setwd(dir = "/home/manie/Documents/INSA_5BIM/GenomMedProjet/Projet_4")
# A ajuster selon les ordis (localisation des donnees)

## RNASeq Data

### Mesothelioma
#We have the raw reads counts for the mesothelioma samples, therefore this matrix must to be normalised. To do so we will use the method developped in DeSeq2. Furthermore, we will remove the genes associated with the sexual chromosomes for this dataset. Doing those two operations the preprocessing steps applied to the mesothelioma data set will be different from those applied to the others samples (LNEN + Carcinoids), for which have the FPKM matrices including the sexual chromosomes (See: the folllowing section).

gene_counts_mesothelioma = read.csv('Data/gene_count_matrix_1passMeso.csv',row.names = 1)
print(gene_counts_mesothelioma[1:3,1:15])
print(dim(gene_counts_mesothelioma))

# Remove sex chromosomes genes before doing the normalisation by using the annotations from:
expr_annot = read.table("Data/TCGA-3H-AB3K-01A_pass1_gene_abund.tab",header = T,sep = "\t")[,1:6]

# Normalisation of raw reads counts through DeSeq2:
# creation object deseq2
deseqexpr = DESeqDataSetFromMatrix(gene_counts_mesothelioma, colData = data.frame(colnames(gene_counts_mesothelioma)), design = ~1, tidy = F)
# Accroche la matrice annot à notre matrice count
annot_ordered = expr_annot[sapply( rownames(deseqexpr), function(x) which(expr_annot$Gene.ID==x)[1] ),]
# Enleve les chromosomes sexuels et Mitochondrial (garde les differents on a une perte de 10)
deseqexpr_nosex = deseqexpr[!annot_ordered$Reference %in% c("chrM","chrX","chrY"),]

vstexpr_nosex = varianceStabilizingTransformation(deseqexpr_nosex,blind = T)
vstexpr_nosex_meso= assay(vstexpr_nosex) ## Final matrix to work with

### LNEN
#Transcriptomic data from LNEN samples are already normalized (FPKM) but the sexual chromosomes should be removed.**
gene_counts_FPKM_LNEN = read.csv('Data/gene_FPKM_matrix_LnenSamples.csv',row.names = 1)
print(dim(gene_counts_FPKM_LNEN))

### Carcinoids
#As before the normalized are already avaible.
gene_counts_FPKM_carcinoids = read.csv('Data/gene_FPKM_matrix_CarcinoidSamples.csv',row.names = 1)
print(dim(gene_counts_FPKM_carcinoids))

# Combine both datasets
gene_counts_FPKM_LNNEN_carcinoids = cbind(gene_counts_FPKM_LNEN, gene_counts_FPKM_carcinoids)

## Methylation arrays
# For methylation arrays are only available for the LNENs sample. we propose you to work Wwith the the beta values matrix, nevertheless the M values matrix is also available if needed. Both have already been normalised. Finally these matrices exclude the sexual chromosomes.
BetalValNormLNENs <- read.csv('Data/NormalisedFilteredBetaTable_LnenSamples.csv', header = T, row.names = 1)
print(head(BetalValNormLNENs))
print(dim(BetalValNormLNENs ))

# Les fichiers à utiliser sont : vstexpr_nosex_meso, gene_counts_FPKM_LNNEN_carcinoids et BetalValNormLNENs pour la méthylation
