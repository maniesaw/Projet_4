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

# Mise en place des donnees

# Localisation Data
setwd(dir = "/home/manie/Documents/INSA_5BIM/GenomMedProjet/Projet_4/Data")
# A ajuster selon les ordis (localisation des donnees)

GCM1Car = read.csv("gene_count_matrix_1pass_CarcinoidSamples.csv", row.names = "gene_id")
GCM1Lnen = read.csv("gene_count_matrix_1pass_LnenSamples.csv", row.names = "gene_id")

GFPKMCar = read.csv("gene_FPKM_matrix_CarcinoidSamples.csv", row.names = "gene_id")
GFPKMLnen = read.csv("gene_FPKM_matrix_LnenSamples.csv", row.names = "gene_id")

NFB = read.csv("NormalisedFilteredBetaTable_CarcinoidSamples.csv", row.names = "X")
NFM = read.csv("NormalisedFilteredMTable_noInf_CarcinoidSamples.csv", row.names = "X")



print(gene_counts_mesothelioma[1:3,1:15])

# Analyse
