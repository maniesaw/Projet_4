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
setwd(dir = "/home/anne/Melanie/Projet_4")
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


# Sauvegarde des fichiers
write.csv2(vstexpr_nosex_meso,"data_vstexpr_nosex_meso.csv")
write.csv2(gene_counts_FPKM_LNNEN_carcinoids, "data_gene_counts_FPKM_LNNEN_carcinoids.csv")

#####################################################################################################

# Ouvertures des donnees
setwd(dir = "/home/anne/Melanie/Projet_4")
vstexpr_nosex_meso = read.csv("data_vstexpr_nosex_meso.csv", row.names = 1, sep=';', header=TRUE, dec=",")
gene_counts_FPKM_LNNEN_carcinoids = read.csv("data_gene_counts_FPKM_LNNEN_carcinoids.csv", row.names = 1, sep=';', header=TRUE)
BetalValNormLNENs <- read.csv('Data/NormalisedFilteredBetaTable_LnenSamples.csv', header = T, row.names = 1)

# Les fichiers à utiliser sont : vstexpr_nosex_meso, gene_counts_FPKM_LNNEN_carcinoids et BetalValNormLNENs pour la méthylation

# On fait correspondre les ID (ens avec biomaRt): (SI BESOIN EST)
library(biomaRt)

Names = rownames(vstexpr_nosex_meso)

ensembl<-  useMart("ensembl",dataset="hsapiens_gene_ensembl") # Get Ensembl data for human
# From the "mart" object get the attribute of interest, here we extract the ensembl Id such as ENS000 and REf RNASeq ID (something like NM_...)
ensembl_ID_refseq_mrna <- getBM(attributes=c('ensembl_gene_id_version','refseq_mrna'), mart = ensembl)

#vstexpr_nosex_meso_merge = merge(vstexpr_nosex_meso, ensembl_ID_refseq_mrna)
 
#####################################################################################

# Tests des packages

#Epidish
# https://bioconductor.org/packages/devel/bioc/vignettes/EpiDISH/inst/doc/EpiDISH.html
library(EpiDISH)

out.l <- epidish(beta.m = BetalValNormLNENs, ref.m = centDHSbloodDMC.m, method = "RPC")
boxplot(out.l$estF)

### ICI on a des resultats

# Music

# https://cran.r-project.org/web/packages/BisqueRNA/vignettes/bisque.html
# https://xuranw.github.io/MuSiC/reference/music_prop.html
# https://www.synapse.org/#!Synapse:syn21041850/wiki/600865
# https://xuranw.github.io/MuSiC/articles/MuSiC.html

library(Biobase)
library(MuSiC)
library(BisqueRNA)

Matrix_bulk = as.matrix(vstexpr_nosex_meso)
bulk.eset <- Biobase::ExpressionSet(assayData = Matrix_bulk) # Transformation donnees bulk

hlca = read.csv("Data2/hlca_counts.csv", row.names = 1) # recuperation des donnees single cell


# recuperation des genes et des ind
Col = colnames(hlca)
Ind <- c()
Gen <- c()

for (i in 1:length(Col)){
Indi = unlist(strsplit(Col[i], split = ".", fixed=TRUE))[1]
Geni = unlist(strsplit(Col[i], split = ".", fixed=TRUE))[3]

Ind <- c(Ind,Indi)
Gen <- c(Gen,Geni)
}

individual.labels = Ind # recuperation nom pour la suite
cell.type.labels = Gen
  
sample.ids <- colnames(hlca)

# individual.ids and cell.types should be in the same order as in sample.ids
sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=sample.ids,
                       SubjectName=individual.labels,
                       cellType=cell.type.labels)

sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))

sc.pdata <- new("AnnotatedDataFrame",
                data=sc.pheno,
                varMetadata=sc.meta)

sc.eset <- Biobase::ExpressionSet(assayData=as.matrix(hlca),
                                  phenoData=sc.pdata)

# Music mise en place
library(xbioc)
Est.prop = music_prop(bulk.eset, sc.eset, clusters = 'cellType',
                      samples = 'sampleID')
names(Est.prop)


# Exemple TUTO

GSE50244.bulk.eset = readRDS("Data2/GSE50244bulkeset.rds")
EMTAB.eset = readRDS("Data2/EMTABesethealthy.rds")

Est.prop.GSE50244 = music_prop(bulk.eset = GSE50244.bulk.eset, sc.eset = EMTAB.eset, clusters = 'cellType',
                               samples = 'sampleID', select.ct = c('alpha', 'beta', 'delta', 'gamma',
                                                                   'acinar', 'ductal'), verbose = F)
names(Est.prop.GSE50244)
      