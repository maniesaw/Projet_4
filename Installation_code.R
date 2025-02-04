# Mise en place des packages

# Prerequis
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")

# Supervised

# RNAseq

## DeconvSeq
install.packages("devtools")
library(devtools) # Si pb prendre le lien dans drive

install_github("rosedu1/deconvSeq", dependencies=TRUE) # (all si demande: 1)
library(deconvSeq)

## DeconRNASeq


BiocManager::install("DeconRNASeq")
library(DeconRNASeq)

browseVignettes("DeconRNASeq")

## CellMix
if( !require(BiocInstaller) ){
   # enable Bioconductor repositories
   # -> add Bioc-software
   setRepositories() 
   
   install.packages('BiocInstaller')
   library(BiocInstaller)
} # mETTRE chiffre?

source('https://www.bioconductor.org/biocLite.R')

# install (NB: this might ask you to update some of your packages)
biocLite('CellMix', siteRepos = 'http://web.cbio.uct.ac.za/~renaud/CRAN', type='both')

BiocManager::install("CellMix")

## MuSic
devtools::install_github('xuranw/MuSiC')
library(MuSiC)

# methylation arrays

## EpiDISH
BiocManager::install("EpiDISH")
library(EpiDISH)

# BOTH

## DIABLO
BiocManager::install('mixOmics')
library('mixOmics')

# Unsupervised

# RNA-seq

## deconica
devtools::install_github("UrszulaCzerwinska/DeconICA", build_vignettes = TRUE, dependencies = TRUE)
library(deconica)

## fastICA
install.packages("fastICA")
library(fastICA)
   
# Methylation arrays

## RefFreeEWAS
install.packages("RefFreeEWAS")
library(RefFreeEWAS)

## medepir
remotes::install_github("bcm-uga/medepir")
library(medepir)

## DecompPipeline
devtools::install_github("CompEpigen/DecompPipeline")
library(DecompPipeline)

## MeDeCom
install.packages("remotes")
remotes::install_github("lutsik/MeDeCom")
library("MedeCom")

##  FactorViz
devtools::install_github("CompEpigen/FactorViz")
library(FactorViz)

   # ANY

## Decoder
install.packages("decoder")
library(decoder)   

# BOTH

## intNMF
install.packages("IntNMF")
library(IntNMF)


# Annexe
BiocManager::install("DESeq2")
library(DeSeq2)

BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

BiocManager::install("biomaRt") # au besoin pour la conversion
