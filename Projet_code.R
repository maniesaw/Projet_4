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
#??

# BOTH

## intNMF
library(IntNMF)

# Mise en place des donnees

# Localisation Data
setwd(dir = "/Data")
# A ajuster selon les ordis (localisation des donnees)

GCM1Car = read.csv("gene_count_matrix_1pass_CarcinoidSamples.csv", row.names = "gene_id")
GCM1Lnen = read.csv("gene_count_matrix_1pass_LnenSamples.csv", row.names = "gene_id")

GFPKMCar = read.csv("gene_FPKM_matrix_CarcinoidSamples.csv", row.names = "gene_id")
GFPKMLnen = read.csv("gene_FPKM_matrix_LnenSamples.csv", row.names = "gene_id")

NFB = read.csv("NormalisedFilteredBetaTable_CarcinoidSamples.csv", row.names = "X")
NFM = read.csv("NormalisedFilteredMTable_noInf_CarcinoidSamples.csv", row.names = "X")

# Analyse
