
library(SCENIC)
library(SeuratObject)
library(Seurat)

library(devtools)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(devtools)
library(RcisTarget)

library(dplyr)
library(Matrix)
library(GENIE3)
library(tidyr)
library(AUCell)

#library(SCopeLoomR)
#library(hdf5r)
#library(arrow)
#library(KernSmooth)
#library(RColorBrewer)
#library(GSEABase)
#library(remotes)



# Extract command-line arguments
#args <- commandArgs(trailingOnly = TRUE)


### Load data
#seurat_obj <- args[1]
#seurat_obj <- get(seurat_obj)
#assign(scenic_obj, readRDS(args[1]))
#readRDS(args[1])

seurat_obj <- readRDS('../Desktop/P1_seurat_M.RDS')

# Update cell type information
#seurat_obj <- get(scenic_obj)
Idents(seurat_obj) <- 'celltype'


#cellInfo
seurat_obj_cellInfo <- data.frame(cellType = Idents(seurat_obj))
nGene <- seurat_obj@meta.data$nGene
nGene <- as.data.frame(nGene)
nUMI <- seurat_obj@meta.data$nUMI
nUMI <- as.data.frame(nUMI)
seurat_obj_cellInfo["nGene"] <- nGene
seurat_obj_cellInfo["nUMI"] <- nUMI

#exprMatL2 <- GSE130148_Lung_small_seuratobject@assays$RNA@data


#dim(exprMatL22)

exprMat <- as.matrix(seurat_obj[["RNA"]]@counts)
dim(exprMat)

#DefaultAssay(combined) <- 'RNA'
#combined <- NormalizeData(combined)


#### change motifAnnotations with motifAnnotations_mgi to make it work ####
#### for human change it to motifAnnotations_hgnc???
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9


### Initialize settings
org <- "hgnc" # or hgnc, or dmel
dbDir <- "../Desktop" # RcisTarget databases location
myDatasetTitle <- "P1" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=5)

#scenicOptions@inputDatasetInfo$GSE130148_Lung_cellInfo <- "int/GSE130148_Lung_cellInfo.Rds"
#saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# Modify if needed
scenicOptions@inputDatasetInfo$scenic_obj_cellInfo <- "int/scenic_obj_cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/scenic_obj_colVars.Rds"
# Databases:
# scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr"="mm9-tss-centered-5kb-10species.mc8nr.feather")
# scenicOptions@settings$db_mcVersion <- "v8"

# Save to use at a later time...
####saveRDS(scenicOptions, file="int/scenicOptions.Rds")


### Co-expression network
scenicOptions <- readRDS("int/scenicOptions.Rds")
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions, minCountsPerGene=3*.01*ncol(exprMat), minSamples=ncol(exprMat)*.01)

#interestingGenes <- c("Sox9", "Sox10", "Dlx5")
# any missing?
#interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)


logMat <- log2(exprMat+1)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#exprMatL2_filtered <- log2(exprMatL2_filtered+1)
#dim(exprMatL2_filtered)
options(error = recover)
set.seed(123)


runCorrelation(exprMat_filtered, scenicOptions)

runGenie3(exprMat_filtered, scenicOptions, nParts = 5, resumePreviousRun = FALSE, allTFs = getDbTfs(scenicOptions))


### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top10perTarget"))
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

saveRDS(scenicOptions, file="int/scenicOptions.Rds")


# Optional: Binarize activity
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
savedSelections <- shiny::runApp(aucellApp)
print(tsneFileName(scenicOptions))

# Save the modified thresholds:
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# scenicOptions@settings$devType="png"
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

#`AUCell_plotTSNE()` to save static plots:
#tsneTfExpression, fig.height=6, fig.width=8}
print(tsneFileName(scenicOptions))
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

# Show TF expression:
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat_filtered, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Dlx5")],], plots="Expression")
