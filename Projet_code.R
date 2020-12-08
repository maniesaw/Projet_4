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

# On fait correspondre les ID (ens avec biomaRt):
library(biomaRt)

Names = rownames(vstexpr_nosex_meso)

ensembl<-  useMart("ensembl",dataset="hsapiens_gene_ensembl") # Get Ensembl data for human
# From the "mart" object get the attribute of interest, here we extract the ensembl Id such as ENS000 and REf RNASeq ID (something like NM_...)
ensembl_ID_refseq_mrna <- getBM(attributes=c('ensembl_gene_id_version','refseq_mrna'), mart = ensembl)

#vstexpr_nosex_meso_merge = merge(vstexpr_nosex_meso, ensembl_ID_refseq_mrna)
 
#####################################################################################

# Tests des packages

# Music
library(ggplot2)
library(MuSiC)
library(xbioc)

GSE50244.bulk.eset = readRDS('Data2/GSE50244bulkeset.rds')
GSE50244.bulk.eset    

EMTAB.eset = readRDS('Data2/EMTABesethealthy.rds')
EMTAB.eset


Est.prop.GSE50244 = music_prop(bulk.eset = GSE50244.bulk.eset, sc.eset = EMTAB.eset, clusters = 'cellType',
                               samples = 'sampleID', select.ct = c('alpha', 'beta', 'delta', 'gamma',
                                                                   'acinar', 'ductal'), verbose = F)
names(Est.prop.GSE50244)

jitter.fig = Jitter_Est(list(data.matrix(Est.prop.GSE50244$Est.prop.weighted),
                             data.matrix(Est.prop.GSE50244$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')

m.prop.GSE50244 = rbind(melt(GSE50244.EMTAB.prop$Est.prop.weighted), 
                        melt(GSE50244.EMTAB.prop$Est.prop.allgene))

colnames(m.prop.GSE50244) = c('Sub', 'CellType', 'Prop')
m.prop.GSE50244$CellType = factor(m.prop.GSE50244$CellType, levels = c('alpha', 'beta', 'delta', 'gamma', 'acinar', 'ductal'))
m.prop.GSE50244$Method = factor(rep(c('MuSiC', 'NNLS'), each = 89*6), levels = c('MuSiC', 'NNLS'))
m.prop.GSE50244$HbA1c = rep(GSE50244.bulk.eset$hba1c, 2*6)
m.prop.GSE50244 = m.prop.GSE50244[!is.na(m.prop.GSE50244$HbA1c), ]
m.prop.GSE50244$Disease = factor(c('Normal', 'T2D')[(m.prop.GSE50244$HbA1c > 6.5)+1], levels = c('Normal', 'T2D'))
m.prop.GSE50244$D = (m.prop.GSE50244$Disease == 'T2D')/5
m.prop.GSE50244 = rbind(subset(m.prop.GSE50244, Disease == 'Normal'), subset(m.prop.GSE50244, Disease != 'Normal'))

jitter.new = ggplot(m.prop.GSE50244, aes(Method, Prop)) + 
  geom_point(aes(fill = Method, color = Disease, stroke = D, shape = Disease), 
             size = 2, alpha = 0.7, position = position_jitter(width=0.25, height=0)) +
  facet_wrap(~ CellType, scales = 'free') + scale_colour_manual( values = c('white', "gray20")) +
  scale_shape_manual(values = c(21, 24))+ theme_minimal()

plot_grid(jitter.fig, jitter.new, labels = 'auto')


#Epidish
library(EpiDISH)
data(centEpiFibIC.m)
data(DummyBeta.m)

out.l <- epidish(beta.m = BetalValNormLNENs, ref.m = centDHSbloodDMC.m, method = "RPC") 
out.l$estF
boxplot(out.l$estF)

