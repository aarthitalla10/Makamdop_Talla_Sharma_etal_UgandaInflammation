######################################
#      Load required libraries
######################################
suppressPackageStartupMessages(library(package = "gdata"))
suppressPackageStartupMessages(library(package = "EDASeq"))
suppressPackageStartupMessages(library(package = "edgeR"))
suppressPackageStartupMessages(library(package = "ggplot2"))
suppressPackageStartupMessages(library(package = "pheatmap"))
suppressPackageStartupMessages(library(package = "grid"))
suppressPackageStartupMessages(library(package = "WriteXLS"))
suppressPackageStartupMessages(library(package = "dplyr"))
suppressPackageStartupMessages(library(package = "tidyr"))
suppressPackageStartupMessages(library(package = "tibble"))
suppressPackageStartupMessages(library(package = "readxl"))
suppressPackageStartupMessages(library(package = "sva"))
suppressPackageStartupMessages(library(package = "igraph"))
suppressPackageStartupMessages(library(package = "EpiDISH"))
suppressPackageStartupMessages(library(package = "Biobase"))
suppressPackageStartupMessages(library(package = "GEOquery"))
suppressPackageStartupMessages(library(package = "biomaRt"))
suppressPackageStartupMessages(library(package = "mixOmics"))


######################################
# Define the global options of the script
######################################
options(stringsAsFactors = FALSE)
options(useFancyQuotes = FALSE)

######################################
# Initializing directories
######################################
rawDir <- "raw"
dataDir <- "data"
diagnosticDir <- "diagnostic_plot"
exploDir <- "exploratory_plot"
degDir <- "deg_plot"
gseaDir <- "gsea"
regDir <- "reg"
figDir <- "Manuscript/Figures_20190223"

########################################################################
#  Create Count Matrix / Samplennotation / Feature annotation
########################################################################
# Read counts files
fileLS <- list.files(path = rawDir,
                     pattern = "genecounts$",
                     full.names = TRUE,
                     recursive = TRUE)
countMat <- lapply(fileLS, FUN = function(file){
                return(value = read.table(file = file,
                                          sep = "\t",
                                          col.names = c("id", "value")))
                            })
# Verify that all the files have the same  ids
id <- countMat[[1]][, "id"]
flag <- sapply(countMat, FUN = function(file) {
                return(value = all(id == file[, "id"]))})
if (any(!flag)) {
    print("warning some tag id missing in some of the count files")
    }
# Merge all the count files
countMat <- sapply(countMat, FUN = function(file) {
             file[, "value"] })
rownames(countMat) <- id
colnames(countMat) <- gsub(pattern = "_genecounts",
                           replacement = "",
                           basename(fileLS))

## Read samplesheet and add in Batch information
sampleSheet <- read_excel("data/Samplesheet.xlsx", sheet = 1) %>% as.data.frame()
batch2Samples <- read_excel("data/Samplesheet.xlsx", sheet = 1) %>%
                 as.data.frame() %>%
                 filter(timePoint %in% "M12") %>%
                .$SampleID
batch3Samples <- read_excel("data/Samplesheet.xlsx", sheet = 2) %>%
                 as.data.frame() %>%
                 .$SampleID
sampleSheet <- sampleSheet %>%
               mutate(Batch = ifelse(SampleID %in% batch3Samples, "B3",
                              ifelse(SampleID %in% batch2Samples, "B2", "B1")))
rownames(sampleSheet) <- sampleSheet$SampleID

## Subset count matrix on samples present in the samplesheet
countMat <- countMat[, rownames(sampleSheet)]
table(rownames(sampleSheet) == colnames(countMat))
countDF <- as.data.frame(countMat)

## Feature Annotation
# name columns for Features Annotation file from 'GTF'
cNames <- c("seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attributes")
featuresAnnotationFile <- "data/genes.gtf"
featuresAnnotation <- read.table(file = featuresAnnotationFile,
                                 sep = "\t",
                                 col.names = cNames)
featuresAnnotation$"gene_id" <- gsub(pattern = ".*gene_id ([^;]*);.*",
                                     replacement = "\\1",
                                     featuresAnnotation$"attributes")
featuresAnnotation$"gene_name" <- gsub(pattern = ".*gene_name ([^;]*);.*",
                                       replacement = "\\1",
                                       featuresAnnotation$"attributes")
featuresAnnotation <- unique(featuresAnnotation[, c("seqname",
                                                    "strand",
                                                    "gene_id",
                                                    "gene_name")])
rownames(featuresAnnotation) <- featuresAnnotation$"gene_id"
featuresAnnotation <- featuresAnnotation[rownames(countDF), ]


######################################
#       Build Expression Set
######################################
# Raw Counts
esetRaw <- newSeqExpressionSet(counts = as.matrix(countDF),
                               featureData = featuresAnnotation,
                               phenoData = sampleSheet)

## Average technical replicates
pData(esetRaw)$repID <- interaction(esetRaw$subjectID,
                                    esetRaw$timePoint,
                                    esetRaw$cellSubset,
                                    esetRaw$Batch) %>% as.character()
replicates <- names(which(table(pData(esetRaw)[["repID"]]) > 1))
replicates <- replicates[-length(replicates)]
for(RID in replicates) {
    dupSample <- sampleNames(esetRaw)[pData(esetRaw)[["repID"]] == RID]
    counts(esetRaw)[, sampleNames(esetRaw) %in% dupSample[1]] <-
            apply(counts(esetRaw[, sampleNames(esetRaw) %in% dupSample]),
                  MARGIN = 1,
                  FUN = mean)
    esetRaw <- esetRaw[, !(sampleNames(esetRaw) %in% dupSample[2:length(dupSample)])]
 }


# Normalized counts
dge <- DGEList(counts = counts(esetRaw), remove.zeros = TRUE)
# note : Removed 842 out of 29881 rows with all zero counts
dge <- calcNormFactors(object = dge, method = "TMM")
normCounts <- cpm(dge, normalized.lib.sizes = TRUE)

# Build SeqExpressionSet with normalized counts
eset <- newSeqExpressionSet(counts = as.matrix(normCounts),
                            featureData = featuresAnnotation[rownames(normCounts), ],
                            phenoData = pData(esetRaw))


############################################################
# FIG S5A and B
############################################################
#  PCA
esetTemp <- eset
mat <- log2(counts(esetTemp) + 0.25)
matpca <- prcomp(mat, center=TRUE, scale. = TRUE)
pcaDF <- as.data.frame(matpca[[2]]) %>% dplyr::select(PC1, PC2) %>%
         rownames_to_column() %>%
         dplyr::rename(SampleID = rowname) %>%
         mutate(cellSubset = esetTemp$type[match(SampleID,
                                                 table = esetTemp$SampleID)],
                timePoint = esetTemp$timePoint[match(SampleID,
                                                     table = esetTemp$SampleID)],
                Batch = esetTemp$Batch[match(SampleID,table = esetTemp$SampleID)])
Scols <- c("Tcells" = "royalblue3",
           "mDC" = "springgreen",
           "Monocytes" = "orchid",
           "NK" = "gold",
           "WATERCONTROL" = "black")
plotScat <- ggplot(data = pcaDF,
                   mapping = aes(x = PC1,
                                 y = PC2,
                                 color = cellSubset)) +
            scale_color_manual(values=Scols) +
            geom_point(size = 2) +
            theme_bw() +
            theme(axis.text.x = element_text(size=8, color = "black"),
                  axis.text.y = element_text(size=8, color = "black"),
                  axis.line = element_line(color = "grey"),
                  panel.background = element_blank(),
                  panel.grid = element_blank())
pdf(file = file.path(figDir, "PCA_withWATERCONTROL_V2.pdf"), width = 4, height = 3)
print(plotScat)
dev.off()

# PCA without water control
sampleSheet2 <- pData(eset) %>% filter(!cellSubset %in% "WATERCONTROL")
rownames(sampleSheet2) <- sampleSheet2$SampleID
countDF2 <- countDF[, rownames(sampleSheet2)]
table(colnames(countDF2) == rownames(sampleSheet2))

# re-normalize counts (post removal of water control samples)
dge <- DGEList(counts = countDF2, remove.zeros = TRUE)
dge <- calcNormFactors(object = dge, method = "TMM")
normCounts <- cpm(dge, normalized.lib.sizes = TRUE)

# PCA
mat <- log2(normCounts + 0.25)
matpca <- prcomp(mat, center=TRUE, scale. = TRUE)
pcaDF <- as.data.frame(matpca[[2]]) %>% dplyr::select(PC1, PC2) %>%
         rownames_to_column() %>%
         dplyr::rename(SampleID = rowname) %>%
         mutate(cellSubset = sampleSheet2$type[match(SampleID,
                                                     table = sampleSheet2$SampleID)],
                timePoint = sampleSheet2$timePoint[match(SampleID,
                                                        table = sampleSheet2$SampleID)],
                Batch = sampleSheet2$Batch[match(SampleID,table = sampleSheet2$SampleID)])
Scols2 <- Scols[!names(Scols) %in% "WATERCONTROL"]
plotScat <- ggplot(data = pcaDF,
                   mapping = aes(x = PC1,
                                 y = PC2,
                                 label = SampleID,
                                 color = cellSubset)) +
            geom_text(size = 3) +
            scale_color_manual(values=Scols2) +
            geom_point(size = 4) +
            theme_bw() +
            theme(axis.text.x = element_text(size=8, color = "black"),
                  axis.text.y = element_text(size=8, color = "black"),
                  axis.line = element_line(color = "grey"),
                  panel.background = element_blank(),
                  panel.grid = element_blank())
pdf(file = file.path(figDir, "PCA_withoutWATERCONTROL.pdf"), width = 6, height = 6)
print(plotScat)
dev.off()

## remove outlier samples based on visual inspection and replot PCA
removeSamples <- pcaDF %>% filter(PC1 > -0.05) %>% .$SampleID
removeSamples <- c(removeSamples, c("Monocytes_86","Monocytes_78","Monocytes_94",
                                    "CD3T_RAL07_Mo3", "C306L4I16"))
sampleSheet3 <- sampleSheet2 %>% filter(!SampleID %in% removeSamples)
rownames(sampleSheet3) <- sampleSheet3$SampleID
countDF3 <- countDF2[, rownames(sampleSheet3)]
table(colnames(countDF3) == rownames(sampleSheet3))

# Normalize counts
dge <- DGEList(counts = countDF3, remove.zeros = TRUE)
dge <- calcNormFactors(object = dge, method = "TMM")
normCounts <- cpm(dge, normalized.lib.sizes = TRUE)

# PCA
mat <- log2(normCounts + 0.25)
matpca <- prcomp(mat, center=TRUE, scale. = TRUE)
pcaDF <- as.data.frame(matpca[[2]]) %>% dplyr::select(PC1, PC2) %>%
         rownames_to_column() %>%
         dplyr::rename(SampleID = rowname) %>%
         mutate(cellSubset = sampleSheet3$type[match(SampleID,
                                                    table = sampleSheet3$SampleID)],
                timePoint = sampleSheet3$timePoint[match(SampleID,
                                                        table = sampleSheet3$SampleID)],
                Batch = sampleSheet3$Batch[match(SampleID,table = sampleSheet3$SampleID)])

# summary PCA
summaryPCA <- summary(matpca)$importance


plotScat <- ggplot(data = pcaDF,
                    mapping = aes(x = PC1,
                                  y = PC2,
                                  color = cellSubset)) +
            scale_color_manual(values = Scols2) +
            geom_point(size = 3) +
            labs(x = paste("PC1(",
                            round(summaryPCA["Proportion of Variance", "PC1"] * 100),
                            "%)", sep = ""),
                 y = paste("PC2(",
                            round(summaryPCA["Proportion of Variance", "PC2"] * 100),
                            "%)", sep = "")) +
            theme_bw() +
            theme(axis.text.x = element_text(size = 10, color = "black"),
                  axis.text.y = element_text(size = 10, color = "black"),
                  axis.line = element_line(color = "grey"),
                  panel.background = element_blank(),
                  panel.grid = element_blank())
pdf(file = file.path(figDir, "Fig_S5b.pdf"),
    width = 4, height = 4, useDingbats = F)
print(plotScat)
dev.off()

# save raw/cpm count matrix and eset after all outliers removed
write.table(countDF3, file = file.path(dataDir, "rawCounts.txt"), sep = "\t")
write.table(normCounts, file = file.path(dataDir, "cpmCounts.txt"), sep = "\t")

# save as esetRaw and eset
esetRaw <- newSeqExpressionSet(counts = as.matrix(countDF3),
                               featureData = featuresAnnotation,
                               phenoData = sampleSheet3)
save(esetRaw, file = file.path(dataDir, "esetRaw_V2.RData"))

eset <- newSeqExpressionSet(counts = as.matrix(normCounts),
                            featureData = featuresAnnotation[rownames(normCounts), ],
                            phenoData = sampleSheet3)
save(eset, file = file.path(dataDir, "eset.RData"))



############################################################
# FIG S5C
############################################################
### PCA by cell subset (pre-batch correction)
esetTemp <- eset
tCols <- c("D0" = "black",
           "M3" = "red",
           "M12" = "coral",
           "M24" = "royalblue3")
for(CS in unique(esetTemp$cellSubset)) {
        esetTemp2 <- esetTemp[, esetTemp$cellSubset %in% CS]
        mat <- log2(counts(esetTemp2) + 0.25)
        matpca <- prcomp(mat, center = TRUE, scale. = TRUE)
        pcaDF <- as.data.frame(matpca[[2]]) %>% dplyr::select(PC1, PC2) %>%
                 rownames_to_column() %>%
                 dplyr::rename(SampleID = rowname) %>%
                 mutate(timePoint = esetTemp2$timePoint[match(SampleID,
                                                              table = esetTemp2$SampleID)],
                        Batch = esetTemp2$Batch[match(SampleID,table = esetTemp2$SampleID)])

        # summary PCA
        summaryPCA <- summary(matpca)$importance

        outFile <- paste("PCA_", CS, "_preCombat.pdf", sep = "")
        plotScat <- ggplot(data = pcaDF,
                           mapping = aes(x = PC1,
                                         y = PC2,
                                         color = timePoint,
                                         shape = Batch)) +
                    scale_color_manual(values = tCols) +
                    geom_point(size = 3) +
                    labs(x = paste("PC1(",
                                    round(summaryPCA["Proportion of Variance", "PC1"] * 100),
                                    "%)", sep = ""),
                         y = paste("PC2(",
                                    round(summaryPCA["Proportion of Variance", "PC2"] * 100),
                                    "%)", sep = "")) +
                    theme_bw() +
                    theme(axis.text.x = element_text(size = 8, color = "black"),
                          axis.text.y = element_text(size = 8, color = "black"),
                          axis.line = element_line(color = "grey"),
                          panel.background = element_blank(),
                          panel.grid = element_blank())
        pdf(file = file.path(figDir, outFile), width = 4, height = 3, useDingbats = FALSE)
        print(plotScat)
        dev.off()
    }


############################################################
# FIG S5D
############################################################
### PCA post comBat batch correction
############ combat for batch correction (per subset) ############
esetTemp <- eset
for(CS in unique(esetTemp$cellSubset)) {
    esetTemp2 <- esetTemp[, esetTemp$cellSubset %in% CS]
    esetTemp2 <- esetTemp2[, !esetTemp2$SampleID %in% "H5GYL2I21"]
    
    # build design matrix
    group1 <- factor((esetTemp2$"timePoint"))
    designMat <- model.matrix(~group1)
    rownames(designMat) <- sampleNames(esetTemp2)
    colnames(designMat) <- gsub(pattern = "group1",
                                replacement = "",
                                colnames(designMat))

    # transform count data into a normal distribution by Voom, since comBat expects normally distributed data
    v <- voom(counts(esetTemp2), design = designMat, plot = FALSE)

    # combat
    batch <- esetTemp2$Batch
    mod <- model.matrix(~ 0 + timePoint, data = pData(esetTemp2))
    colnames(mod) <- gsub("timePoint", "", colnames(mod))
    
    # remove batch effect using combat
    rmBatch <- ComBat(dat = v$E, batch = esetTemp2$Batch)
    mat <- rmBatch
    matpca <- prcomp(mat, center = TRUE, scale. = TRUE)
    pcaDF <- as.data.frame(matpca[[2]]) %>% dplyr::select(PC1, PC2) %>%
             rownames_to_column() %>%
             dplyr::rename(SampleID = rowname) %>%
             mutate(timePoint = esetTemp2$timePoint[match(SampleID,
                                                          table = esetTemp2$SampleID)],
                    Batch = esetTemp2$Batch[match(SampleID, table = esetTemp2$SampleID)])

    # summary PCA
    summaryPCA <- summary(matpca)$importance

    outFile <- paste("PCA_", CS, "_postBatchCorrection.pdf", sep = "")
    plotScat <- ggplot(data = pcaDF,
                       mapping = aes(x = PC1,
                                     y = PC2,
                                     color = timePoint,
                                     shape = Batch)) +
                scale_color_manual(values = tCols) +
                geom_point(size = 3) +
                labs(x = paste("PC1(",
                                round(summaryPCA["Proportion of Variance", "PC1"] * 100),
                                "%)", sep = ""),
                     y = paste("PC2(",
                                round(summaryPCA["Proportion of Variance", "PC2"] * 100),
                                "%)", sep = "")) +
                theme_bw() +
                theme(axis.text.x = element_text(size = 8, color = "black"),
                      axis.text.y = element_text(size = 8, color = "black"),
                      axis.line = element_line(color = "grey"),
                      panel.background = element_blank(),
                      panel.grid = element_blank())
    pdf(file = file.path(figDir, outFile), width = 4, height = 3)
    print(plotScat)
    dev.off()
 }


############################################################
# FIG S5E
############################################################
### Which subsets are making these cytokines at baseline
markers <- c(clusters[2], clusters[4], clusters[3]) %>% unlist() %>% unname()
markers <- c(markers, "TGFB1", "CXCL10")

esetTemp <- eset
pData(esetTemp)$BatchTimePoint <- interaction(esetTemp$Batch,
                                              esetTemp$timePoint, sep = "_")
esetTemp <- esetTemp[, !esetTemp$BatchTimePoint %in% "B3_D0"]
esetTemp <- esetTemp[, esetTemp$timePoint %in% "D0"]

# replace cytokines with HGNC symbols
markers[markers == "IL1b"] <- "IL1B"
markers[markers == "MIP1b"] <- "CCL4"
markers[markers == "GCSF"] <- "CSF3"
markers[markers == "MCP1"] <- "CCL2"
markers[markers == "TNFa"] <- "TNF"
markers[markers == "IFNg"] <- "IFNG"
markers[markers == "IL12"] <- "IL12A"


## for markers
markers <- rownames(counts(esetTemp))[grep("^CCL|^CXCL", rownames(counts(esetTemp)))]

# plot heatmap of markers genes
filter <- apply(counts(esetTemp), 1, function(x) mean(x)>0)
esetTemp <- esetTemp[filter, ]
mat <- counts(esetTemp)
comm <- intersect(markers, rownames(mat))
mat <- mat[comm, ]


# rename HGNC symbols back to what cytokine was names
rownames(mat)[rownames(mat) %in% "CCL4"] <- "MIP1B"
rownames(mat)[rownames(mat) %in% "CSF3"] <- "GCSF"
rownames(mat)[rownames(mat) %in% "CCL2"] <- "MCP1"
rownames(mat)[rownames(mat) %in% "CXCL10"] <- "IP10"

# scale and set breaks
mat <- t(scale(t(mat)))
limit <- range(mat)
limit <- c(ceiling(limit[1]), floor(limit[2]))
limit <- min(abs(limit))

# phenotype annotation
matAnnot <- pData(esetTemp) %>%
            dplyr::select(cellSubset, timePoint) %>%
            as.data.frame()
matAnnot <- matAnnot[order(factor(matAnnot$cellSubset, levels = c("Tcells",
                                                                  "mDC",
                                                                  "Monocytes",
                                                                  "NK"))), ]
ann_colors = list(cellSubset = c(Tcells = "royalblue3",
                                 mDC = "springgreen",
                                 Monocytes = "orchid",
                                 NK = "gold"),
                  timePoint = c(D0 = "black"))

# plot
outFile <- "Markers_GeneExpression_Heatmap_D0_V2.pdf"
print(outFile)
pdf(file = file.path(figDir, outFile), width  = 8, height = 5)

colorPalette <- c("blue", "white", "red")
colorPalette <- colorRampPalette(colors = colorPalette)(100)

pheatmap(mat = mat[, rownames(matAnnot)],
         color = colorPalette,
         breaks = c(min(mat),
                    seq(from = -1 * limit,
                        to = limit,
                        length.out = 99),
                    max(mat)),
         cellwidth = 6,
         cellheight = 6,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         treeheight_row = 0,
         annotation = matAnnot,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         border_color = NA,
         fontsize_row = 6,
         fontsize = 7)
dev.off()

############################################################
# FIG S5E
############################################################
######## Pathway expression at baseline across all subsets
# do SLEA (Sample Level Enrichment Analysis) function(z-score per sample)
doSLEA <- function(expressionSet, geneSet) {
           # scale expression
           exprsMat  <- expressionSet
           exprsMat  <- t(scale(t(exprsMat)))
           # extract expression of leGenes of each geneset
           comm <- intersect(geneSet, rownames(expressionSet))
           gsDF <- exprsMat[comm, ]
           # calculate mean expression per sample
           gsM <- colMeans(gsDF)
           # extract random genes of size of the geneSet from full probeset and calculate mean
           # and perform this for 'n' permutations
           nperm <- lapply(1:1000, function(j) {
                       # set seed for every permutation
                       set.seed(j)
                       rGSDF <- exprsMat[sample.int(nrow(exprsMat), length(comm)), ]
                       rGSM <- colMeans(rGSDF)
                       return(value = rGSM)
                   })
           permDF <- do.call(rbind, nperm)
           zscore <- (gsM - colMeans(permDF)) / apply(permDF,2,sd)
           sleaDF <- zscore %>% as.data.frame()
           return(value = sleaDF)
       }

# define geneSet background
gmx <- "h.all.v6.1.symbols.gmt"
dirName <- "gsea"
gmxFile <- file.path(dirName, gmx)
colNames <- max(count.fields(file = gmxFile, sep = "\t"))
colNames <- seq(from = 1, to = colNames)
colNames <- as.character(colNames)
gmx <- read.table(file = gmxFile,
                  sep = "\t",
                  quote = "\"",
                  fill = TRUE,
                  col.names = colNames,
                  row.names = 1)
gmx <- gmx[, -1]
gmx <- apply(gmx, MARGIN = 1, FUN = function(x) {
                return(value = setdiff(unname(x), ""))
           })
names(gmx) <- toupper(names(gmx))

## expression baseline
esetTemp <- eset
pData(esetTemp)$BatchTimePoint <- interaction(esetTemp$Batch,
                                              esetTemp$timePoint, sep = "_")
esetTemp <- esetTemp[, !esetTemp$BatchTimePoint %in% "B3_D0"]
pData(esetTemp)$DonorID_TimePoint <- interaction(esetTemp$subjectID,
                                                 esetTemp$timePoint, sep = "_")
esetTemp <- esetTemp[, esetTemp$timePoint %in% "D0"]
filter <- apply(counts(esetTemp), 1, function(x) mean(x)>0)
esetTemp <- esetTemp[filter, ]

# call SLEA
sleaLS <- lapply(1:length(gmx), function(l) {
               expressionSet = counts(esetTemp)
               geneSet <- gmx[[l]] %>% strsplit(",") %>% unlist(.)
               sDF <- doSLEA(expressionSet = expressionSet, geneSet = geneSet)
               names(sDF) <- names(gmx[l])
               return(value = sDF)
           })
sleaDF <- do.call(cbind, sleaLS)
colnames(sleaDF) <- gsub("HALLMARK_", "", colnames(sleaDF))
sleaDF <- sleaDF %>% t() %>% as.data.frame()

# annotation
matAnnot <- data.frame(SampleID = colnames(sleaDF)) %>%
            mutate(cellSubset = esetTemp$cellSubset[match(SampleID,
                                                         table = esetTemp$SampleID)],
                  TimePoint = "D0") %>%
            column_to_rownames(var = "SampleID")

# annotation colors
matAnnot <- matAnnot[order(factor(matAnnot$cellSubset, levels = c("Tcells",
                                                                 "mDC",
                                                                 "Monocytes",
                                                                 "NK"))), ]
ann_colors = list(cellSubset = c(Tcells = "royalblue3",
                                 mDC = "springgreen",
                                 Monocytes = "orchid",
                                 NK = "gold"),
                  TimePoint = c(D0 = "black"))

# scale and set breaks
mat <- t(scale(t(sleaDF)))
limit <- range(mat)
limit <- c(ceiling(limit[1]), floor(limit[2]))
limit <- min(abs(limit))
mat <- mat[, rownames(matAnnot)]

colorPalette <- c("blue", "white", "red")
colorPalette <- colorRampPalette(colors = colorPalette)(100)

outFile <- "Pathways at baseline across subsets.pdf"
pdf(file = file.path(figDir,outFile), width = 12, height = 12)
pheatmap(mat = mat,
         color = colorPalette,
         breaks = c(min(mat),
                    seq(from = -1 * limit,
                        to = limit,
                        length.out = 99),
                    max(mat)),
         cellwidth = 6,
         cellheight = 6,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         treeheight_row = 0,
         annotation = matAnnot,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         border_color = NA,
         fontsize_row = 6,
         fontsize = 7)
dev.off()


############################################################
# Identify Differential expressed genes between timepoints
############################################################
fits <- list()
SUBSETS <- c("mDC", "Monocytes", "NK", "Tcells")
for(C in SUBSETS) {
  
    esetTemp <- eset
    esetTemp <- esetTemp[, esetTemp$cellSubset %in% C]
    pData(esetTemp)$BatchTimePoint <- interaction(esetTemp$Batch,
                                                  esetTemp$timePoint, sep = "_")
    esetTemp <- esetTemp[, !esetTemp$BatchTimePoint %in% "B3_D0"]
    pData(esetTemp)$DonorID_TimePoint <- interaction(esetTemp$subjectID,
                                                     esetTemp$timePoint, sep = "_")

     # remove rows with mean expression 0
    filter <- apply(counts(esetTemp), 1, function(x) mean(x)>0)
    esetTemp <- esetTemp[filter, ]
    
    # build design matrix
    group <- factor(esetTemp$"timePoint")
    designMat <- model.matrix(~ 0 + group)
    rownames(designMat) <- sampleNames(esetTemp)
    colnames(designMat) <- gsub(pattern = "group",
                                replacement = "",
                                colnames(designMat))
    attr(designMat, "assign") <- attr(designMat, "contrasts") <- NULL
    
    # transform count data into log2 CPM - a normal distribution by Voom
    v <- voom(counts(esetTemp), design = designMat, plot = FALSE)

    # lmFit
    fit <- lmFit(v, design = designMat)
    
    # define contrast
    contrastLS <- c(paste(setdiff(group, "D0"), "-D0", sep = ""),
                    "M12-M3", "M24-M12", "M24-M3")
    contrastMat <- makeContrasts(contrasts = contrastLS, levels = designMat)
    
    fit2 <- contrasts.fit(fit, contrastMat)
    fit2 <- eBayes(fit2)
    fit2$genes <- fData(esetTemp)[rownames(fit$coef), ]
    modelName <- paste(C,"_FTest",sep="")
    fits[[modelName]] <- list(fit = fit, fit2 = fit2)
    
    # print number of genes differently expressed and make heatmaps
    fitsTemp <- fits
    fit2 <- fitsTemp[[modelName]][["fit2"]]
    coefName <- colnames(fit2)
    
    # save sigTags
    sigTags <- topTable(fit = fit2, coef = coefName, n = Inf)
    sigTags <- sigTags[order(sigTags$adj.P.Val, decreasing = FALSE), ]
    
    # write results to degDir
    outputFile <- paste(modelName, ".FTest.topTable.txt", sep="")
    write.table(sigTags,
                file = file.path(degDir, outputFile),
                quote = FALSE,
                sep = "\t",
                row.names = FALSE)
    
    # print # DEGs up and dn
    print(paste("NomPval = ", dim(sigTags %>% filter(`P.Value` < 0.05))[1], sep = ""))
    print(paste("FDR = ", dim(sigTags %>% filter(`adj.P.Val` < 0.05))[1], sep = ""))

    # make heatmaps of top 50 DEGs
    sigTags <- sigTags[1:50, ]
    mat <- v$E %>% as.data.frame()
    mat <- mat[rownames(sigTags), ]
    
    # scale expression and define limits for color gradient on heatmaps
    mat <- t(scale(t(mat)))
    limit <- range(mat)
    limit <- c(ceiling(limit[1]), floor(limit[2]))
    limit <- min(abs(limit))
    
    # phenotype annotation
    matAnnot <- pData(esetTemp)[, c("cellSubset","timePoint","subjectID")]
    ordercol <- matAnnot[order(match(matAnnot$timePoint, table = c("D0", "M3", "M12", "M24"))), ]
    ordercol <- rownames(ordercol)

# colors
    ann_colors = list(cellSubset = c(mDC = "springgreen", Monocytes = "orchid",
                                     NK = "gold", Tcells = "royalblue3"),
                      timePoint = c(D0 = "black",
                                    M3 = "red",
                                    M12 = "coral",
                                    M24 = "royalblue3"))

    # plot heatmap
    fileName <- paste(C,"_FTest_heatmap.pdf",sep = "")
    pdf(file = file.path(degDir, fileName), width = 10, height = 10)
    
    colorPalette <- c("blue", "white", "red")
    colorPalette <- colorRampPalette(colors = colorPalette)(100)
    
    pheatmap(mat = mat[, ordercol],
             color = colorPalette,
             breaks = c(min(mat),
                        seq(from = -1 * limit,
                            to = limit,
                            length.out = 99),
                        max(mat)),
             cellwidth = 5,
             cellheight = 5,
             cluster_cols = FALSE,
             clustering_distance_cols = "euclidean",
             treeheight_row = 0,
             annotation = matAnnot,
             annotation_colors = ann_colors,
             show_rownames = TRUE,
             border_color = NA,
             fontsize = 5)
    dev.off()
}

############################################################
###### Perform Fisher pathway enrichment among FTest genes
############################################################
gmx <- "TH17.gmt"
dirName <- "gsea"
gmxFile <- file.path(dirName, gmx)
colNames <- max(count.fields(file = gmxFile, sep = "\t"))
colNames <- seq(from = 1, to = colNames)
colNames <- as.character(colNames)
gmx <- read.table(file = gmxFile,
                  sep = "\t",
                  quote = "\"",
                  fill = TRUE,
                  col.names = colNames,
                  row.names = 1)
gmx <- gmx[, -1]
gmx <- apply(gmx, MARGIN = 1, FUN = function(x) {
            return(value = setdiff(unname(x), ""))
        })
names(gmx) <- toupper(names(gmx))


### making a GMT of the monocyte subset signatures
### Monocyte subsets signatures
#gmx <- lapply(c(2:4), function(i) {
#           signatures <- read_excel("data/TableS3_gse25913.xlsx", sheet = i) %>%
#                         as.data.frame()
#           colnames(signatures) <- signatures[1,]
#           signatures <- signatures[-1,] %>% .$SYMBOL %>% unique(.)
#           return(value = signatures)
#       })
#names(gmx) <- c("Classical", "Intermediate", "NonClassical")


# read FTest genes
## 1. Genes : M3 and M12 > D0 and M24 < M3 and M12 (Monocytes and Tcells)
fileDF <- read.delim("deg_plot/Tcells_FTest.FTest.topTable.txt") # output from identifying DEGs
gs <- fileDF %>%
      filter(`adj.P.Val` < 0.05) %>%
      filter(`M3.D0` > 0 & M12.D0 > 0 & M24.M3 < 0, M24.M12 < 0) %>%
     .$gene_name %>% unique(.)

## 2. Genes : M3 and M12 < D0 and M24 > M3 and M12 (Monocytes and Tcells)
gs <- fileDF %>%
      filter(`adj.P.Val` < 0.05) %>%
      filter(`M3.D0` < 0 & M12.D0 < 0 & M24.M3 > 0, M24.M12 > 0) %>%
     .$gene_name %>% unique(.)


### do Fisher
# obtain background
bg <- unique(unlist(gmx))
fisherIndex <- NULL
output <- mclapply(gmx, function(x) {
            tab <- table(factor(bg %in% gs, levels = c(TRUE, FALSE)),
                         factor(bg %in% x, levels = c(TRUE, FALSE)))
            fit <- fisher.test(tab, alternative = "greater")
            interSection <- intersect(gs, x)
            interSection <- paste(interSection, collapse = ",")
            return(value = c(RATIO = as.numeric((diag(tab) / colSums(tab))[1]),
                             NOM_pval = fit$p.value,
                             INTERSECT = interSection))
                    })
output <- do.call(what = rbind, args = output)
output <- cbind(output, NAME = names(gmx))
pAdj <- p.adjust(as.numeric(output[, "NOM_pval"]), method = "BH")
output <- cbind(output, ADJ_pval = pAdj)
output <- output[, c("NAME",
                     "RATIO",
                     "NOM_pval",
                     "ADJ_pval",
                     "INTERSECT")]
output2 <- output[order(as.numeric(output[, "NOM_pval"])), ] %>% as.data.frame() %>%
          mutate(NOM_pval = as.numeric(NOM_pval),
                 ADJ_pval = as.numeric(ADJ_pval)) %>%
          filter(NOM_pval < 0.05)

fileName= "Tcells_FTest_Pathways_TH17_M3_M12_higher.txt"
write.table(output2,
            file = file.path("deg_plot", fileName),
            row.names=FALSE,
            sep="\t",
            quote = FALSE)


#### MAKE HEATMAPS of the significant fisher enrichment pathways
###### Th17 signatures heatmap
fileDF <- read.delim("deg_plot/Tcells_FTest_Pathways_TH17_M3_M12_higher.txt")
genes <- fileDF %>%
         filter(ADJ_pval < 0.05) %>%
        .$INTERSECT %>%
         strsplit(",") %>%
         unlist() %>%
         unique()

# subset on eset
esetTemp <- eset
esetTemp <- esetTemp[, esetTemp$cellSubset %in% "Tcells"]
pData(esetTemp)$DonorID_TimePoint <- interaction(esetTemp$subjectID,
                                                 esetTemp$timePoint, sep = "_")
pData(esetTemp)$BatchTimePoint <- interaction(esetTemp$Batch,
                                              esetTemp$timePoint, sep = "_")
esetTemp <- esetTemp[, !esetTemp$BatchTimePoint %in% "B3_D0"]

# remove rows with mean expression 0
filter <- apply(counts(esetTemp), 1, function(x) mean(x)>0)
esetTemp <- esetTemp[filter, ]
mat <- counts(esetTemp)[genes, , drop = F]
mat <- log2(mat + 0.25)
mat <- t(scale(t(mat)))
limit <- range(mat)
limit <- c(ceiling(limit[1]), floor(limit[2]))
limit <- min(abs(limit))

#### mat annotation
matAnnot <- pData(esetTemp)[, c("cellSubset", "timePoint", "DonorID_TimePoint"), drop = FALSE]

## add serratia ratio annotation
matAnnot$Serratia <- pathogenratioDF$serratia_otherbacteria[match(matAnnot$DonorID_TimePoint,
                                                table = pathogenratioDF$DonorID_TimePoint)]
cOrder <- matAnnot[order(match(matAnnot$timePoint, table=c("D0", "M3", "M12", "M24"))),] %>%
          rownames_to_column() %>%
          .$rowname

# annotation colors
ann_colors = list(cellSubset = c(Tcells = "royalblue3"),
                  timePoint = c(D0 = "black",
                                M3 = "red",
                                M12 = "coral",
                                M24 = "royalblue3"),
                  Serratia = c("blue", "white", "red"))
# plot heatmap
fileName <- "TH17signatures.pdf"
pdf(file = file.path(figDir, fileName), width = 10, height = 20)

colorPalette <- c("blue", "white", "red")
colorPalette <- colorRampPalette(colors = colorPalette)(100)
pheatmap(mat = mat[, cOrder],
         color = colorPalette,
         scale = "none",
         breaks = c(min(mat),
                    seq(from = -1 * limit,
                        to = limit,
                        length.out = 99),
                    max(mat)),
         cellwidth = 12,
         cellheight = 12,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation = matAnnot,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         treeheight_row = 0,
         border_color = NA,
         fontsize_row = 3)
dev.off()


## cytokine cluster annotation over heatmap
datDF <- read_excel("PlasmaBiomarkers.xlsx", sheet = 1) %>% as.data.frame()
markers <- clusters[c(2,4,3)] %>% unlist(.) %>% unname()
cytoDF1 <- datDF %>%
           dplyr::filter(TimePoint %in% c("d0", "M3", "M12", "M24")) %>%
           dplyr::select_(.dots = c("SampleID", "TimePoint", unname(markers))) %>%
           gather(Cytokine, value, -SampleID, -TimePoint) %>%
           mutate(Cluster = ifelse(Cytokine %in% clusters[[2]], "C1",
                            ifelse(Cytokine %in% clusters[[3]], "C3", "C2"))) %>%
           spread(TimePoint, value) %>%
           mutate(D0_FC = log2(d0/d0),
                  M3_FC = log2(M3/d0),
                  M12_FC = log2(M12/M3),
                  M24_FC = log2(M24/M12)) %>%
          dplyr::select(-d0, -M3, -M12, -M24) %>%
          gather(TimePoint, value, -SampleID, -Cytokine, -Cluster) %>%
          dplyr::group_by(SampleID, Cluster, TimePoint) %>%
          dplyr::summarize(mFC = mean(value, na.rm = TRUE)) %>%
          as.data.frame() %>%
          mutate(TimePoint = gsub("_FC", "", TimePoint),
                 DonorID_TimePoint = interaction(SampleID, TimePoint, sep = "_"))
cytokineClusterDF <- cytoDF1 %>%
                     dplyr::select(-SampleID, -TimePoint) %>%
                     spread(Cluster, mFC)


# cytokine bar plot for annotation above heatmap
plotDF <- cytokineClusterDF %>%
          mutate(SampleID = esetTemp$SampleID[match(DonorID_TimePoint,
                                                    table = esetTemp$DonorID_TimePoint)],
                 TimePoint = esetTemp$timePoint[match(DonorID_TimePoint,
                                                      table = esetTemp$DonorID_TimePoint)]) %>%
          filter(SampleID %in% cOrder) %>%
          dplyr::select(SampleID, DonorID_TimePoint, C3, TimePoint) %>%
          column_to_rownames(var = "SampleID")
plotDF$SampleID = rownames(plotDF)
plotDF <- plotDF[cOrder, ]
table(rownames(plotDF) == cOrder)

cols <- c("D0" = "black", "M3" = "red", M12 = "coral", M24 = "royalblue3")
outFile <- "Tcells_Heatmap_barplotannotations_C3.pdf"
pdf(file = file.path(figDir, outFile), width = 7, height = 8)
barpl <- ggplot(data = plotDF,
                mapping = aes(x = SampleID, y = C3, color = TimePoint)) +
         geom_bar(stat = "identity", width = 0.1) +
         scale_color_manual(values = cols) +
         scale_x_discrete(limits = rownames(plotDF)) +
         labs(x = NULL, y = "Cluster 3: Mean - log2(Fold Change)") +
         theme_bw() +
         theme(panel.grid = element_blank(),
               axis.text.x = element_blank(),
               axis.text.y = element_text(color = "black", size = 10),
               axis.ticks.x = element_blank())
print(barpl)
dev.off()

#### Same procedure is followed for identifying DEGs and pathways enriched in other cell subsets and plotting figures



###########################################################################
# FIG. S6
###########################################################################
############### In-silico deconvolution of Monocytes  ###############

#### Get monocyte subset frequencies (classical, non-c, int) from Monocyte RNA-Seq subset data using CiberSort
# create reference expression levels of each monocyte subset
load("data/eset.gse25913.RData")
pData(eset) <- pData(eset) %>%
               rownames_to_column() %>%
               mutate(Subset2 = ifelse(Subset %in% "CD14++CD16-", "Classical",
                                ifelse(Subset %in% "CD14++CD16+", "Intermediate",
                                ifelse(Subset %in% "CD14+CD16+", "NonClassical", NA)))) %>%
               column_to_rownames(var = "rowname")

# fetch monocyte subset signatures
monoSigs <- lapply(c(2:4), function(i) {
                signatures <- read_excel("data/TableS3_gse25913.xlsx", sheet = i) %>%
                              as.data.frame()
                colnames(signatures) <- signatures[1,]
                signatures <- signatures[-1,] %>% .$SYMBOL %>% unique(.)
                return(value = signatures)
            })
names(monoSigs) <- c("Classical", "Intermediate", "NonClassical")
sapply(monoSigs, function(x) length(x))
genes <- unique(unlist(monoSigs) %>% as.character())

# get mean gene expression from each subset using the 'GSE25913' gene expression data
mat <- exprs(eset)
rownames(mat) <- make.unique(fData(eset)$Symbol[match(rownames(mat), table = fData(eset)$ID)])
comm <- intersect(rownames(mat), genes)
ref.m <- mat[comm, , drop = FALSE] %>%
         as.data.frame() %>%
         rownames_to_column(var = "SYMBOL") %>%
         gather(SampleID, value, -SYMBOL) %>%
         mutate(cellSubset = eset$Subset2[match(SampleID, table = rownames(pData(eset)))]) %>%
         dplyr::group_by(cellSubset,SYMBOL) %>%
         dplyr::summarize(mValue = mean(value)) %>%
         as.data.frame() %>%
         spread(cellSubset, mValue) %>%
         column_to_rownames(var = "SYMBOL")


##### Get monocyte subset frequencies across time (cibersort)
cellsubset = "Monocytes"
esetTemp <- eset
pData(esetTemp)$BatchTimePoint <- interaction(esetTemp$Batch,
                                              esetTemp$timePoint, sep = "_")
esetTemp <- esetTemp[, !esetTemp$BatchTimePoint %in% "B3_D0"]
pData(esetTemp)$typeTimePoint <- interaction(esetTemp$type,
                                             esetTemp$timePoint, sep = "_")
pData(esetTemp)$DonorID_TimePoint <- interaction(esetTemp$subjectID,
                                                 esetTemp$timePoint, sep = "_")
esetTemp <- esetTemp[, !esetTemp$typeTimePoint %in% "mDC_Mo4"]
esetTemp <- esetTemp[, esetTemp$cellSubset %in% cellsubset]
mat <- counts(esetTemp)

# cibersort
out.l <- epidish(mat, as.matrix(ref.m), method = "CBS")
freqDF <- as.data.frame(out.l$estF) %>%
          rownames_to_column() %>%
          mutate(DonorID = esetTemp$subjectID[match(rowname,table=esetTemp$SampleID)],
                 TimePoint = esetTemp$timePoint[match(rowname,table=esetTemp$SampleID)]) %>%
          dplyr::select(-rowname) %>%
          gather(Monocyte, freq, -DonorID, -TimePoint) %>%
          mutate(freq = freq *100,
                 DonorID_TimePoint = interaction(DonorID, TimePoint, sep = "_"))

# median of each type
xorder <- c("D0", "M3", "M12", "M24")

# plot
cols <- c("Intermediate" = "darkgreen", "Classical" = "red", "NonClassical" = "blue")
pdf(file = file.path(figDir, "Monocyte subset frequencies across time.pdf"), height = 4, width = 5)
plotJit <- ggplot(data = freqDF,
                  mapping = aes(x = TimePoint, y = freq)) +
           geom_boxplot(aes(fill = Monocyte), outlier.colour = NA) +
           scale_fill_manual(values = cols) +
           scale_x_discrete(limits = xorder) +
           labs(x = NULL, y = "Monocyte subset frequency") +
           theme_bw() +
           theme(panel.grid = element_blank(),
                 axis.text.x = element_text(color = "black", size = 10),
                 axis.text.y = element_text(color = "black", size = 10))
print(plotJit)
dev.off()





######## Correlating gene expression signatures with pathogen ratio ########
## gene expression set
esetTemp <- eset
pData(esetTemp)$BatchTimePoint <- interaction(esetTemp$Batch,
                                              esetTemp$timePoint, sep = "_")
esetTemp <- esetTemp[, !esetTemp$BatchTimePoint %in% "B3_D0"]
pData(esetTemp)$DonorID_TimePoint <- interaction(esetTemp$subjectID,
                                                esetTemp$timePoint, sep = "_")
#esetTemp <- esetTemp[, !esetTemp$timePoint %in% c("D0")]
filter <- apply(counts(esetTemp), 1, function(x) mean(x)>0)
esetTemp <- esetTemp[filter, ]

### Th17 signatures
# extract fisher genes
esetTemp2 <- esetTemp[, esetTemp$cellSubset %in% "Tcells"]
filter <- apply(counts(esetTemp2), 1, function(x) mean(x)>0)
esetTemp2 <- esetTemp2[filter, ]
sigs <- read.delim("deg_plot/Tcells_FTest_Pathways_TH17_M3_M12_higher.txt") %>%
        .$INTERSECT %>%
        strsplit(",") %>%
        unlist() %>% unique()
# call SLEA
expressionSet = counts(esetTemp2)
geneSet <- sigs
sDF1 <- doSLEA(expressionSet = expressionSet, geneSet = geneSet)
names(sDF1) <- "TH17"
sDF1 <- sDF1 %>%
        rownames_to_column() %>%
        mutate(DonorID = esetTemp2$subjectID[match(rowname, table = esetTemp2$SampleID)],
               timePoint = esetTemp2$timePoint[match(rowname, table = esetTemp2$SampleID)],
               SampleID = interaction(DonorID, timePoint, sep = "_")) %>%
        dplyr::select(-rowname, -DonorID, -timePoint)


#Th1 + Th2
fileDF1 <- read.delim("deg_plot/g_th1.txt") # refer to supplementary table for Th1 signatures
fileDF2 <- read.delim("deg_plot/g_th2.txt") # refer to supplementary table for Th1 signatures
sigs <- c(fileDF1$gene_name, fileDF2$gene_name) %>% unique()
# call SLEA
expressionSet = counts(esetTemp2)
geneSet <- sigs
sDF2 <- doSLEA(expressionSet = expressionSet, geneSet = geneSet)
names(sDF2) <- "Th1_Th2"
sDF2 <- sDF2 %>%
       rownames_to_column() %>%
       mutate(DonorID = esetTemp2$subjectID[match(rowname, table = esetTemp2$SampleID)],
          timePoint = esetTemp2$timePoint[match(rowname, table = esetTemp2$SampleID)],
          SampleID = interaction(DonorID, timePoint, sep = "_")) %>%
              dplyr::select(-rowname, -DonorID, -timePoint)


### GATA3 + TCF7 taregts signatures
# extract fisher genes
sigs <- read.delim("deg_plot/Tcells_FTest_Pathways_CHEA_M24_higher.txt") %>%
        filter(NAME %in% c("GATA3", "TCF7")) %>%
        .$INTERSECT %>%
        strsplit(",") %>%
        unlist() %>% unique()
# call SLEA
expressionSet = counts(esetTemp2)
geneSet <- sigs
sDF3 <- doSLEA(expressionSet = expressionSet, geneSet = geneSet)
names(sDF3) <- "GATA3_TCF7"
sDF3 <- sDF3 %>%
        rownames_to_column() %>%
        mutate(DonorID = esetTemp2$subjectID[match(rowname, table = esetTemp2$SampleID)],
               timePoint = esetTemp2$timePoint[match(rowname, table = esetTemp2$SampleID)],
               SampleID = interaction(DonorID, timePoint, sep = "_")) %>%
        dplyr::select(-rowname, -DonorID, -timePoint)


## monocyte inflamatory sigs
# extract fisher genes
esetTemp2 <- esetTemp[, esetTemp$cellSubset %in% "Monocytes"]
filter <- apply(counts(esetTemp2), 1, function(x) mean(x)>0)
esetTemp2 <- esetTemp2[filter, ]
sigs <- read.delim("deg_plot/Monocyte_Pathways/Monocytes_FTest_Pathways_Hallmark_M3_M12_higher.txt") %>%
    .$INTERSECT %>%
    strsplit(",") %>%
    unlist() %>% unique()
# call SLEA
expressionSet = counts(esetTemp2)
geneSet <- sigs
sDF4 <- doSLEA(expressionSet = expressionSet, geneSet = geneSet)
names(sDF4) <- "Monocytes_Inf"
sDF4 <- sDF4 %>%
         rownames_to_column() %>%
         mutate(DonorID = esetTemp2$subjectID[match(rowname, table = esetTemp2$SampleID)],
         timePoint = esetTemp2$timePoint[match(rowname, table = esetTemp2$SampleID)],
         SampleID = interaction(DonorID, timePoint, sep = "_")) %>%
         dplyr::select(-rowname, -DonorID, -timePoint)

## Non-Classical monocytes
sigs <- read.delim("deg_plot/Monocyte_Pathways/Monocytes_FTest_Monocyte Subsets_M3 and M12_higher.txt") %>%
       filter(NAME %in% "NonClassical") %>%
      .$INTERSECT %>%
       strsplit(",") %>%
       unlist() %>% unique()
# call SLEA
expressionSet = counts(esetTemp2)
geneSet <- sigs
sDF5 <- doSLEA(expressionSet = expressionSet, geneSet = geneSet)
names(sDF5) <- "Monocytes_NC"
sDF5 <- sDF5 %>%
        rownames_to_column() %>%
        mutate(DonorID = esetTemp2$subjectID[match(rowname, table = esetTemp2$SampleID)],
               timePoint = esetTemp2$timePoint[match(rowname, table = esetTemp2$SampleID)],
               SampleID = interaction(DonorID, timePoint, sep = "_")) %>%
        dplyr::select(-rowname, -DonorID, -timePoint)

## serratia ratio
serratia <- pathogenratioDF[, c("DonorID_TimePoint", "serratia_otherbacteria")] %>%
            dplyr::rename(SampleID = DonorID_TimePoint)


# merge all tables
colnames(cytokineClusterDF)[1] <- "SampleID"
corDF <- (merge(merge(merge(merge(sDF1, sDF2, by = "SampleID"),
                           sDF3, by = "SampleID"),
                      sDF4, by = "SampleID"),
                    sDF5, by = "SampleID")
                serratia, by = "SampleID") %>%
        mutate(timepoint = gsub(".+_(.+)", "\\1", SampleID))
corDF <- merge(corDF, cytokineClusterDF, by = "SampleID")

## CD4/CD8 data frame
CD4 <- read_excel("RALTcellcounts.xlsx", sheet = 1) %>%
       as.data.frame() %>%
       dplyr::rename(CD4CD8Ratio = `CD4:CD8 Ratio`) %>%
       dplyr::select(DonorID, TimePoint, CD4CD8Ratio) %>%
       mutate(CD4CD8Ratio = as.numeric(CD4CD8Ratio)) %>%
       spread(TimePoint, CD4CD8Ratio) %>%
       gather(TimePoint, CD4CD8, -DonorID) %>%
       mutate(TimePoint = gsub("d", "D", TimePoint),
              SampleID = interaction(DonorID, TimePoint, sep = "_")) %>%
       dplyr::select(-DonorID, -TimePoint)

# correlation between intervals
corDF2 <- merge(corDF, CD4, by = "SampleID", all.x = TRUE) %>%
          mutate(TimePoint = gsub(".+_(.+)", "\\1", SampleID)) %>%
          filter(TimePoint %in% c("D0", "M12"))
cor.test(corDF2$CD4CD8, corDF2$C2, method = "spearman")

