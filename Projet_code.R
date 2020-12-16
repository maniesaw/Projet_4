# Data Preprocessing

## Localisation Data
setwd(dir = "/home/anne/Melanie/Projet_4")
# Adjust localisation

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

################################################################################

# Data opening
setwd(dir = "/home/anne/Melanie/Projet_4")

vstexpr_nosex_meso = read.csv("preprocessing/data_vstexpr_nosex_meso.csv", row.names = 1, sep=';', header=TRUE, dec=",")
gene_counts_FPKM_LNNEN_carcinoids = read.csv("preprocessing/data_gene_counts_FPKM_LNNEN_carcinoids.csv", row.names = 1, sep=';', header=TRUE, dec=",")
BetalValNormLNENs <- read.csv('Data/NormalisedFilteredBetaTable_LnenSamples.csv', header = T, row.names = 1)

# Files : vstexpr_nosex_meso, gene_counts_FPKM_LNNEN_carcinoids et BetalValNormLNENs for méthylation


################################################################################
################################################################################

# Tests of packages


## Epidish

# https://bioconductor.org/packages/devel/bioc/vignettes/EpiDISH/inst/doc/EpiDISH.html
library(EpiDISH)

out.l <- epidish(beta.m = BetalValNormLNENs, ref.m = centDHSbloodDMC.m, method = "RPC")
boxplot(out.l$estF)

################################################################################

## Music

# https://cran.r-project.org/web/packages/BisqueRNA/vignettes/bisque.html
# https://xuranw.github.io/MuSiC/reference/music_prop.html
# https://www.synapse.org/#!Synapse:syn21041850/wiki/600865
# https://xuranw.github.io/MuSiC/articles/MuSiC.html


### Preprocessing eset

#### bulk eset MESO
Namebulk <- row.names(vstexpr_nosex_meso)

Namebulk2 = Namebulk
for (i in 1:length(Namebulk)){
  A = Namebulk[i]
  B = unlist(strsplit(A, split = ".", fixed=TRUE))[1]
  Namebulk2[i]=B
}

vstexpr_nosex_meso2 = vstexpr_nosex_meso
vstexpr_nosex_meso2$Namebulk = Namebulk2
row.names(vstexpr_nosex_meso2) <- vstexpr_nosex_meso2$Namebulk
vstexpr_nosex_meso2 <- vstexpr_nosex_meso2[,-87]
write.csv2(vstexpr_nosex_meso2,"vstexpr_nosex_meso2.csv")

bulk.eset<- Biobase::ExpressionSet(assayData = as.matrix(vstexpr_nosex_meso2))
saveRDS(object = bulk.eset, "bulk.eset4.rds")
remove(vstexpr_nosex_meso2, bulk.eset)


#### bulk eset Carcinoid
Namebulk <- row.names(gene_counts_FPKM_LNNEN_carcinoids)

Namebulk2 = Namebulk
for (i in 1:length(Namebulk)){
  A = Namebulk[i]
  B = unlist(strsplit(A, split = ".", fixed=TRUE))[1]
  Namebulk2[i]=B
}

gene_counts_FPKM_LNNEN_carcinoids2 = gene_counts_FPKM_LNNEN_carcinoids
gene_counts_FPKM_LNNEN_carcinoids2$Namebulk = Namebulk2

Dup =  rev(which(duplicated(Namebulk2)))
for (i in Dup){
  gene_counts_FPKM_LNNEN_carcinoids2 <- gene_counts_FPKM_LNNEN_carcinoids2[-i,]
}

row.names(gene_counts_FPKM_LNNEN_carcinoids2) <- gene_counts_FPKM_LNNEN_carcinoids2$Namebulk
gene_counts_FPKM_LNNEN_carcinoids2 <- gene_counts_FPKM_LNNEN_carcinoids2[,-240]
write.csv2(gene_counts_FPKM_LNNEN_carcinoids2,"gene_counts_FPKM_LNNEN_carcinoids2.csv")

bulk.eset<- Biobase::ExpressionSet(assayData = as.matrix(gene_counts_FPKM_LNNEN_carcinoids2))
saveRDS(object = bulk.eset, "bulk.esetCarci.rds")
remove(gene_counts_FPKM_LNNEN_carcinoids2, bulk.eset)


#### single RNA esest
hlca = read.csv("Data2/hlca_counts.csv", row.names = 1) # recuperation single cell data

# Conversion gene name
library(biomaRt)
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)

annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'hgnc_symbol',
    'ensembl_gene_id'),
  uniqueRows = TRUE)

head(annotLookup)


HGCN =row.names(hlca)
a= 0
for (i in 1:length(HGCN)){
  if (!(is.element(HGCN[i], annotLookup$hgnc_symbol))){
    
  }
  else{
    Elmt = which(annotLookup$hgnc_symbol==HGCN[i])
    HGCN[i] = annotLookup$ensembl_gene_id[Elmt]
    a = a+1}
}

hlca$HGCN <- HGCN

# Supression duplicate 
Dup =  rev(which(duplicated(HGCN)))
for (i in Dup){
  hlca <- hlca[-i,]
}

row.names(hlca) <- hlca$HGCN
hlca <- hlca[,-9410]

remove(annotLookup, mart, a, Dup, Elmt,HGCN,i)

write.csv2(hlca,"hlca2.csv")
remove(hlca2,hlca)


# Remove not "ENS"
hlca2 = read.csv("hlca2.csv", row.names = 1, sep=';', header=TRUE, dec=",")

namehlca <-rownames(hlca2)
Ind <- c()
library(stringr)
for (i in 1:length(namehlca)){
  val = length(namehlca)-i+1
  A = namehlca[val][1][1]
  if (str_detect(A, "ENSG")){
    
  }
  else{
   Ind <- c(Ind,val)
  }
  print(val)
}


hlca3 <- hlca2[-Ind,]

write.csv2(hlca3,"hlca3.csv")
remove(hlca2,A,i,Ind,namehlca,val)


# Processing to have same ID
hlca3 = read.csv("hlca3.csv", row.names = 1, sep=';', header=TRUE, dec=",")
GeneSc <-row.names(hlca3)

GeneSc2 = GeneSc
for (i in 1:length(GeneSc)){
  A = GeneSc[i]
  B = unlist(strsplit(A, split = ".", fixed=TRUE))[1]
  GeneSc2[i]=B
}
hlca3$GeneSc = GeneSc2
row.names(hlca3) <- hlca3$GeneSc
hlca3 <- hlca3[,-9410]
write.csv2(hlca3,"hlca4.csv")

hlca4 = hlca3
remove(hlca3,GeneSc,GeneSc2,i,A,B)
hlca4 = read.csv("hlca/hlca4.csv", row.names = 1, sep=';', header=TRUE, dec=",")

# recuperation des genes et des ind
Col = colnames(hlca4)
Ind_id <- c()
Gen_id <- c()
Ind <- c()
Gen <- c()
Ind_name <- c()
Gen_name <- c()



for (i in 1:length(Col)){
A = unlist(strsplit(Col[i], split = ".", fixed=TRUE))[1]
Indi = unlist(strsplit(A, split = "_", fixed=TRUE))[1]
Geni = unlist(strsplit(A, split = "_", fixed=TRUE))[2]

if (is.element(Indi,Ind)){
  R <- which(Ind==Indi)
  Ind_id <- c(Ind_id,R)
}
else{
  Ind <-c(Ind,Indi)
  R <- which(Ind==Indi)
  Ind_id <- c(Ind_id,R)
}

Ind_name <- c(Ind_name,Indi)

if (is.element(Geni,Gen)){
  R <- which(Gen==Geni)
  Gen_id <- c(Gen_id,R)
}
else{
  Gen <-c(Gen,Geni)
  R <- which(Gen==Geni)
  Gen_id <- c(Gen_id,R)
}

Gen_name <- c(Gen_name,Geni)
}

sample.ids <- colnames(hlca4)

remove( Indi,Geni,i,A,R)
remove(Ind,Gen, Col)

# individual.ids and cell.types should be in the same order as in sample.ids
sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=sample.ids,
                       sampleID=as.factor(Ind_name),
                       SubjectName=Ind_id,
                       cellTypeID=Gen_id,
                       cellType=as.factor(Gen_name))

remove(Ind_id,Gen_id,Gen_name, Ind_name)

sc.meta <- data.frame(labelDescription=c("sampleID","SubjectName",
                                         "cellTypeID","cellType"),
                      row.names=c("sampleID","SubjectName",
                                  "cellTypeID","cellType"))

sc.pdata <- new("AnnotatedDataFrame",
                data=sc.pheno,
                varMetadata=sc.meta)

remove(sc.meta,sc.pheno)

sc.eset <- Biobase::ExpressionSet(assayData=as.matrix(hlca4),
                                  phenoData=sc.pdata)


# Cleaning
remove(hlca3,individual.labels, cell.type.labels, sc.pdata, sample.ids)
saveRDS(object = sc.eset, "sc.eset4.rds")
remove(hlca4)
remove( sc.eset)

######################################################################

## MuSiC
library(xbioc)
library(Biobase)
library(MuSiC)
library(BisqueRNA)
library(cowplot)
library(reshape2)

setwd(dir = "/home/anne/Melanie/Projet_4")
bulk.eset <- readRDS(file = "bulk.eset4.rds")
sc.eset <- readRDS(file = "sc.eset4.rds")

# En parametre on donne ce qui sert a faire la clusterisation et les echantillons,verbose c'est si on veut afficher les resultats
Est.prop = music_prop(bulk.eset, sc.eset, clusters = 'cellType', samples = 'sampleID')
names(Est.prop)

saveRDS(object = Est.prop, "Est.prop.rds")
#1 B001222" "B001224" "B001227" "B001235" "B001237" "B001239" "B002578" "B002579"
#2 "BP12"    "BP7"     "B001223" "B001225" "B001229" "B001234" "BP4"     "BP9"    
#3 "B002580" "B002592" "B002593" "B001226" "B001228" "B001238" "B002586" "B002591"
#4 "BP10"    "BP5"     "BP6"     "BP11"    "BP1"     "BP8"     "BP2"     "BP3"    
#5 "B002014" "B002382" "B003138" "B003140" "B003141" "B002383" "B003135" "B003136"
#6 "B003137" "B002386" "B003139" "B002385" "B003142" "B001773" "B002460" "B003269"
#7 "B002461" "B001772" "B001774" "B002464" "B001769" "B001771"

jitter.fig = Jitter_Est(list(data.matrix(Est.prop$Est.prop.weighted),
                             data.matrix(Est.prop$Est.prop.allgene)),
                        method.name = c('MuSiC','NNLS'), title = 'Jitter plot of Est Proportions')


plot_grid(jitter.fig, labels = 'auto')


# Part 1
Est.prop1 = music_prop(bulk.eset, sc.eset, clusters = 'cellType', samples = 'sampleID', 
                       select.ct  = c("B001222","B001224", "B001227" ,"B001235" ,"B001237", "B001239","B002578","B002579"))

jitter.fig = Jitter_Est(list(data.matrix(Est.prop1$Est.prop.weighted),
                             data.matrix(Est.prop1$Est.prop.allgene)),
                        method.name = c('MuSiC','NNLS'), title = 'Jitter plot of Est Proportions')


plot_grid(jitter.fig, labels = 'auto')

# Part 2
Est.prop2 = music_prop(bulk.eset, sc.eset, clusters = 'cellType', samples = 'sampleID', 
                       select.ct  = c("BP12","BP7", "B001223", "B001225", "B001229", "B001234" ,"BP4","BP9"   ))

jitter.fig = Jitter_Est(list(data.matrix(Est.prop2$Est.prop.weighted),
                             data.matrix(Est.prop2$Est.prop.allgene)),
                        method.name = c('MuSiC','NNLS'), title = 'Jitter plot of Est Proportions')


plot_grid(jitter.fig, labels = 'auto')

# Part 3
Est.prop3 = music_prop(bulk.eset, sc.eset, clusters = 'cellType', samples = 'sampleID', 
                       select.ct  = c("B002580","B002592", "B002593", "B001226" ,"B001228", "B001238" ,"B002586", "B002591"))

jitter.fig = Jitter_Est(list(data.matrix(Est.prop3$Est.prop.weighted),
                             data.matrix(Est.prop3$Est.prop.allgene)),
                        method.name = c('MuSiC','NNLS'), title = 'Jitter plot of Est Proportions')


plot_grid(jitter.fig, labels = 'auto')

# Part 4
Est.prop4 = music_prop(bulk.eset, sc.eset, clusters = 'cellType', samples = 'sampleID', 
                       select.ct  = c( "BP10",    "BP5" ,    "BP6",     "BP11" ,   "BP1" ,    "BP8"   ,  "BP2" ,    "BP3"))

jitter.fig = Jitter_Est(list(data.matrix(Est.prop4$Est.prop.weighted),
                             data.matrix(Est.prop4$Est.prop.allgene)),
                        method.name = c('MuSiC','NNLS'), title = 'Jitter plot of Est Proportions')


plot_grid(jitter.fig, labels = 'auto')

# Part 5
Est.prop5 = music_prop(bulk.eset, sc.eset, clusters = 'cellType', samples = 'sampleID', 
                       select.ct  = c("B002014","B002382","B003138","B003140", "B003141" ,"B002383", "B003135", "B003136"))

jitter.fig = Jitter_Est(list(data.matrix(Est.prop5$Est.prop.weighted),
                             data.matrix(Est.prop5$Est.prop.allgene)),
                        method.name = c('MuSiC','NNLS'), title = 'Jitter plot of Est Proportions')


plot_grid(jitter.fig, labels = 'auto')

# Part 6
Est.prop6 = music_prop(bulk.eset, sc.eset, clusters = 'cellType', samples = 'sampleID', 
                       select.ct  = c("B003137", "B002386", "B003139", "B002385", "B003142", "B001773", "B002460", "B003269"))

jitter.fig = Jitter_Est(list(data.matrix(Est.prop6$Est.prop.weighted),
                             data.matrix(Est.prop6$Est.prop.allgene)),
                        method.name = c('MuSiC','NNLS'), title = 'Jitter plot of Est Proportions')


plot_grid(jitter.fig, labels = 'auto')

# Part 7
Est.prop7 = music_prop(bulk.eset, sc.eset, clusters = 'cellType', samples = 'sampleID', 
                       select.ct  = c("B002461", "B001772", "B001774","B002464", "B001769", "B001771"))

jitter.fig = Jitter_Est(list(data.matrix(Est.prop7$Est.prop.weighted),
                             data.matrix(Est.prop7$Est.prop.allgene)),
                        method.name = c('MuSiC','NNLS'), title = 'Jitter plot of Est Proportions')


plot_grid(jitter.fig, labels = 'auto')

#############################################################################################

# test Carcinoid

bulk.eset <- readRDS(file = "bulk.esetCarci.rds")
sc.eset <- readRDS(file = "sc.eset4.rds")

Est.prop = music_prop(bulk.eset, sc.eset, clusters = 'cellType', samples = 'sampleID')
names(Est.prop)

saveRDS(object = Est.prop, "Est.propCarci.rds")

jitter.fig = Jitter_Est(list(data.matrix(Est.prop$Est.prop.weighted),
                             data.matrix(Est.prop$Est.prop.allgene)),
                        method.name = c('MuSiC','NNLS'), title = 'Jitter plot of Est Proportions')


# MuSiC only
jitter.fig = Jitter_Est(list(data.matrix(Est.prop$Est.prop.weighted)),
                        method.name = c('MuSiC'), title = 'Jitter plot of Est Proportions')

plot_grid(jitter.fig, labels = 'auto')

## Clusterisation

Mousesub.eset = readRDS(file = "sc.eset/sc.eset4.rds")
Mousesub.basis = music_basis(Mousesub.eset, clusters = 'cellType', samples = 'sampleID')

# Plot the dendrogram of design matrix and cross-subject mean of realtive abundance
par(mfrow = c(1, 2))
d <- dist(t(log(Mousesub.basis$Disgn.mtx + 1e-6)), method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )
# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1, main = 'Cluster log(Design Matrix)')

d <- dist(t(log(Mousesub.basis$M.theta + 1e-8)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
hc2 <- hclust(d, method = "complete")
# Plot the obtained dendrogram
plot(hc2, cex = 0.6, hang = -1, main = 'Cluster log(Mean of RA)')


#########################################################################################################

## Annexe de MuSiC

### Verification des concordances des genes
GeneSc2 = GeneSc
for (i in 1:length(GeneSc)){
A = GeneSc[i]
B = unlist(strsplit(A, split = ".", fixed=TRUE))[1]
GeneSc2[i]=B
}

Namebulk2 = Namebulk
for (i in 1:length(Namebulk)){
  A = Namebulk[i]
  B = unlist(strsplit(A, split = ".", fixed=TRUE))[1]
  Namebulk2[i]=B
}

### Recuperation des genes en commun 
comp = 0
for (i in Namebulk2){
  if (is.element(i,GeneSc2)){
    comp = comp+1
  }
   print(i) 
}

## Exemple TUTO MuSiC
GSE50244.bulk.eset = readRDS("Data2/GSE50244bulkeset.rds")
EMTAB.eset = readRDS("Data2/EMTABesethealthy.rds")

Est.prop.GSE50244 = music_prop(bulk.eset = GSE50244.bulk.eset, sc.eset = EMTAB.eset, clusters = 'cellType',
                               samples = 'sampleID', select.ct = c('alpha', 'beta', 'delta', 'gamma',
                                                                   'acinar', 'ductal'))
names(Est.prop.GSE50244)
      
jitter.fig = Jitter_Est(list(data.matrix(Est.prop.GSE50244$Est.prop.weighted),
                             data.matrix(Est.prop.GSE50244$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')

GSE50244.EMTAB.prop=Est.prop.GSE50244
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


## Clusterisation

Mousesub.eset = readRDS('Data2/Mousesubeset.rds')
Mousesub.basis = music_basis(Mousesub.eset, clusters = 'cellType', samples = 'sampleID')

# Plot the dendrogram of design matrix and cross-subject mean of realtive abundance
par(mfrow = c(1, 2))
d <- dist(t(log(Mousesub.basis$Disgn.mtx + 1e-6)), method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )
# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1, main = 'Cluster log(Design Matrix)')

d <- dist(t(log(Mousesub.basis$M.theta + 1e-8)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
hc2 <- hclust(d, method = "complete")
# Plot the obtained dendrogram
plot(hc2, cex = 0.6, hang = -1, main = 'Cluster log(Mean of RA)')

############################################################################################
