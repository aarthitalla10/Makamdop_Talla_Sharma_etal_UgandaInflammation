# load required packages
suppressPackageStartupMessages(library(package = "readxl"))
suppressPackageStartupMessages(library(package = "dplyr"))
suppressPackageStartupMessages(library(package = "tidyr"))
suppressPackageStartupMessages(library(package = "tibble"))
suppressPackageStartupMessages(library(package = "MetaboAnalystR"))
#suppressPackageStartupMessages(library(package = "mzR"))
#suppressPackageStartupMessages(library(package = "xcms"))

######################################
# Define the global options
######################################
options(stringsAsFactors = FALSE)
options(useFancyQuotes = FALSE)
figDir = "Manuscript/Figures_20190223"

#### list of all sheets
fileLS <- c("Metabolon/Thaw1 vs Thaw5/CASE-02-17MD+CDT (HUMAN PLASMA SAMPLES)_V2.XLSX",
            "Metabolon/ProjectFiles-CASE-01-18MD+/CASE-01-18MD+ CDT.XLSX",
            "Metabolon/FINAL/CASE-04-18MD+ CDT 20181130.XLSX")

#### merge sample annotation of all batches
sLS <- lapply(1:length(fileLS), function(i) {
            matDF <- read_excel(fileLS[i], sheet = 2) %>% as.data.frame()
            # sample annotation
            s <- matDF[, c(grep("SAMPLE_NAME", colnames(matDF)):ncol(matDF))]
            s <- s[c(1:grep("Group   HMDB_ID", s[, 1])-1), ] %>% t() %>%
                 as.data.frame() %>%
                 rownames_to_column()
            colnames(s) <- s[1, ]
            s <- s[-1, ] %>%
                 dplyr::select(SAMPLE_NAME, SUBJECT_ID, TIME_POINT) %>%
                 mutate(Batch = paste("B",i,sep=""))
            rownames(s) <- NULL
            return(value = s)
        })
sampleAnnotation = do.call(rbind, sLS) %>%
                   mutate(TIME_POINT = gsub("Day |Day_","D",TIME_POINT)) %>%
                   mutate(TIME_POINT = gsub("Mo|Month_","M",TIME_POINT))
# Annotate EDTA batch
edta <- read_excel("Metabolon/EDTASamples.xlsx", sheet = 1) %>% as.data.frame
sampleAnnotation <- sampleAnnotation %>%
                    mutate(Batch = ifelse(SAMPLE_NAME %in% edta$SAMPLE_NAME,
                                          edta$Batch, Batch)) %>%
                    column_to_rownames(var = "SAMPLE_NAME")
    
    
#### merge features annotation of all batches
fLS <- lapply(1:length(fileLS), function(i) {
            f <- read_excel(fileLS[i], sheet = 2) %>% as.data.frame()
            f <- f[c(grep("PATHWAY_SORTORDER", f[, 1]):nrow(f)), ]
            f <- f[, -grep("^CASE-", colnames(f))]
            colnames(f) <- f[1, ]
            f <- f[-1,]
            rownames(f) <- NULL
            f <- f %>%
                 dplyr::select(BIOCHEMICAL, SUPER_PATHWAY, SUB_PATHWAY, `Group   HMDB_ID`)
            return(value = f)
        })
featuresAnnotation <- Reduce(function(...) merge(..., by = "BIOCHEMICAL", all = T), fLS)
featuresAnnotation$HMDB_ID <- gsub("_NA|NA_",
                                   "",
                                   apply(featuresAnnotation[, c("Group   HMDB_ID.x",
                                                                "Group   HMDB_ID.y",
                                                                "Group   HMDB_ID")],
                                         1, function(x) paste(unique(x), collapse="_")))
featuresAnnotation$PATHWAY_SUPER <- gsub("_NA|NA_",
                                         "",
                                         apply(featuresAnnotation[, c("SUPER_PATHWAY.x",
                                                                      "SUPER_PATHWAY.y",
                                                                      "SUPER_PATHWAY")],
                                               1, function(x) paste(unique(x), collapse="_")))
featuresAnnotation$PATHWAY_SUB <- gsub("_NA|NA_",
                                       "",
                                       apply(featuresAnnotation[, c("SUB_PATHWAY.x",
                                                                    "SUB_PATHWAY.y",
                                                                    "SUB_PATHWAY")],
                                             1, function(x) paste(unique(x), collapse="_")))
featuresAnnotation <- featuresAnnotation %>%
                      dplyr::select(BIOCHEMICAL, PATHWAY_SUPER, PATHWAY_SUB, HMDB_ID) %>%
                      column_to_rownames(var = "BIOCHEMICAL")
    

#### merge expression of all batches
eLS <- lapply(1:length(fileLS), function(i) {
            d <- read_excel(fileLS[i], sheet = 2) %>% as.data.frame()
            d <- d[c(grep("PATHWAY_SORTORDER", d[, 1]):nrow(d)),
                     grep("^...2$|^CASE-", colnames(d))]
            d <- d[-1, ]
            d <- d %>%
                 dplyr::rename("BIOCHEMICAL" = "...2")
            return(value = d)
        })
metaboDF <- Reduce(function(...) merge(..., by = "BIOCHEMICAL", all = T), eLS)
metaboDF[2:ncol(metaboDF)] <- apply(metaboDF[2:ncol(metaboDF)], 2, as.numeric)
metaboDF <- metaboDF%>%
            column_to_rownames(var = "BIOCHEMICAL")
metaboDF[is.na(metaboDF)] <- 0
metaboDF <- log2(metaboDF + 0.25)

## remove unknown compounds
metaboDF <- metaboDF[-grep("^X - ", rownames(metaboDF)), , drop = FALSE]
featuresAnnotation <- featuresAnnotation[rownames(metaboDF), , drop = FALSE]

### check data integrity
table(rownames(featuresAnnotation) == rownames(metaboDF))
table(rownames(sampleAnnotation) == colnames(metaboDF))


## PCA
matpca <- prcomp(metaboDF, center = TRUE, scale. = TRUE)
pcaDF <- as.data.frame(matpca[[2]]) %>% dplyr::select(PC1, PC2)
pcaDF$TIME_POINT <- sampleAnnotation$TIME_POINT
pcaDF$Batch <- sampleAnnotation$Batch
plotScat <- ggplot(data = pcaDF,
                   mapping = aes(x = PC1, y = PC2, color = Batch, shape = Batch)) +
            geom_point(size = 4) +
            theme_bw()
pdf(file = file.path(figDir, "PCA_preBatchcorrection.pdf"), width = 5, height = 4)
print(plotScat)
dev.off()

##### Batch correction
mod <- model.matrix(~SUBJECT_ID+TIME_POINT, data = sampleAnnotation)
cb <- removeBatchEffect(x = as.matrix(metaboDF),
                        batch = sampleAnnotation$Batch,
                        design = mod)

## PCA
matpca <- prcomp(cb, center = TRUE, scale. = TRUE)
pcaDF <- as.data.frame(matpca[[2]]) %>% dplyr::select(PC1, PC2)
pcaDF$TIME_POINT <- sampleAnnotation$TIME_POINT
pcaDF$Batch <- sampleAnnotation$Batch
plotScat <- ggplot(data = pcaDF,
                    mapping = aes(x = PC1, y = PC2, color = TIME_POINT, shape = Batch)) +
            geom_point(size = 4) +
            theme_bw()
pdf(file = file.path(figDir, "PCA_postRmBacthEffects.pdf"), width = 5, height = 4)
print(plotScat)
dev.off()


### avg expression (on batch corrected matrix) across batches for each timepoint
mat <- cb %>% as.data.frame() %>%
       rownames_to_column() %>%
       gather(SampleID, value, -rowname) %>%
       mutate(TIME_POINT = sampleAnnotation$TIME_POINT[match(SampleID,
                                                             table = rownames(sampleAnnotation))],
              Batch = sampleAnnotation$Batch[match(SampleID, table = rownames(sampleAnnotation))],
              SUBJECT_ID = sampleAnnotation$SUBJECT_ID[match(SampleID,
                                                             table = rownames(sampleAnnotation))]) %>%
     # remove M6 and M9 since only few common samples added (M9 = 5, M6 = 1)
       filter(!TIME_POINT %in% c("M6", "M9")) %>%
       dplyr::group_by(SUBJECT_ID, TIME_POINT, rowname) %>%
       dplyr::summarize(mValue = mean(value, na.rm = TRUE)) %>%
       mutate(SampleID = interaction(SUBJECT_ID, TIME_POINT, sep = "_")) %>%
       as.data.frame() %>%
       dplyr::select(SampleID, rowname, mValue) %>%
       spread(SampleID, mValue) %>%
       column_to_rownames(var = "rowname")
mat <- mat^2

####  corrected metabolon matrix
cbMetabolon <- mat
### Keep only the most varying metabolites for down stream analysis
varMet <- data.frame(v = apply(cbMetabolon, 1, var)) %>%
          rownames_to_column() %>%
          filter(v >= median(v)) %>%
          .$rowname
cbMetabolon2 <- cbMetabolon[varMet, , drop = FALSE]

# keep common samples between metabolon and pathseq (RA) of full breakdown
comm <- intersect(pathogenratioDF$DonorID_TimePoint, colnames(cbMetabolon2))
pathogenratioDF2 <- pathogenratioDF %>%
                    filter(DonorID_TimePoint %in% comm) %>%
                    column_to_rownames(var = "DonorID_TimePoint")
pathogenratioDF2 <- pathogenratioDF2[comm, , drop = F]
cbMetabolon2 <- cbMetabolon2[,comm, drop = FALSE]
table(rownames(pathogenratioDF2) == colnames(cbMetabolon2))



######################################################################
############## Correlation between Metabolon and PathSeq ##############
######################################################################
## univariate spearman correlation between pathogen and metabolites
X <- t(cbMetabolon2)
cLS <- lapply(1:ncol(X), function(i) {
            lapply(1:ncol(pathogenratioDF2), function(j) {
            cT <- cor.test(as.numeric(X[, i]),
                           as.numeric(pathogenratioDF2[, j]))
            df <- data.frame(met = colnames(X)[i],
                             pathogen = colnames(pathogenratioDF2)[j],
                             pVal = cT$p.value,
                             rho = cT$estimate)
            })
    })
cDF <- do.call(rbind, lapply(cLS, function(x) do.call(rbind, x))) %>%
       group_by(pathogen) %>%
       mutate(pAdjust = p.adjust(pVal, method = "BH"),
              PATHWAY_SUB = featuresAnnotation$PATHWAY_SUB[match(met,
                                                                 table = rownames(featuresAnnotation))]) %>%
       as.data.frame() %>%
       filter(pVal < 0.05)


# carnitines
carnitines <- unique(cDF$met[grep("Fatty Acid Metabolism", cDF$PATHWAY_SUB)])
# glycolysis
glycolysis <- unique(cDF$met[grep("Glycolysis", cDF$PATHWAY_SUB)])
cDF_Sig <- cDF %>%
           filter(met %in% c(carnitines,glycolysis) & rho > 0) %>%
           filter(pathogen %in% "serratia_otherbacteria")

# taking the mean of carnitines
sigCarnitines <- cDF_Sig$met[grep("carnitine", cDF_Sig$met)]
mcarnitineDF <- cbMetabolon2[sigCarnitines, , drop = FALSE]
mCarnitines <- colMeans(as.matrix(mcarnitineDF))

# data frame of carnitines and glycoltic markers
sigMetDF <- t(cbMetabolon2[c(sigCarnitines, glycolysis), , drop = FALSE]) %>%
            as.data.frame()
sigMetDF$mCarnitines <- mCarnitines
sigMetDF$serratia_otherbacteria <- pathogenratioDF2$serratia_otherbacteria[match(rownames(sigMetDF), table = rownames(pathogenratioDF2))]


### correlation between Th17 signatures and carnitines
fileDF <- read.delim("deg_plot/Tcells_FTest_Pathways_TH17_M3_M12_higher.txt")

    genes <- fileDF %>% .$INTERSECT %>% strsplit(",") %>% unlist() %>% unique()
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
    exprsMat <- counts(esetTemp)

    # define geneSet background
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
    gmx <- gmx[fileDF$NAME]

    # call SLEA
    sleaLS <- lapply(1:length(gmx), function(l) {
                    expressionSet = exprsMat
                    geneSet <- gmx[[l]] %>% strsplit(",") %>% unlist(.)
                    sDF <- doSLEA(expressionSet = expressionSet, geneSet = geneSet)
                    names(sDF) <- names(gmx[l])
                    return(value = sDF)
                })
    sleaDF <- do.call(cbind, sleaLS)
    mTH17 <- rowMeans(as.matrix(sleaDF))
    sleaDF$mTH17 <- mTH17
    sleaDF <- sleaDF %>%
              rownames_to_column() %>%
              mutate(DonorID = esetTemp$subjectID[match(rowname, table = esetTemp$SampleID)],
                     timePoint = esetTemp$timePoint[match(rowname, table = esetTemp$SampleID)],
                     DonorID_TimePoint = interaction(DonorID, timePoint, sep = "_")) %>%
              dplyr::select(-rowname, -DonorID, -timePoint)

# bind to sleaDF
carnitineDF <- cbMetabolon2[unique(cDF_C$met), , drop = FALSE]
mCarnitines <- colMeans(as.matrix(carnitineDF))
carnitineDF <- data.frame(mCarnitines) %>%
               mutate(DonorID_TimePoint = colnames(carnitineDF))
corDF <- merge(carnitineDF, sleaDF, by = "DonorID_TimePoint") %>%
         mutate(Timepoint = gsub(".+_(.+)", "\\1", DonorID_TimePoint))

# cor.test
cor.test(corDF$mCarnitines, corDF$mTH17, method = "spearman")

# correlation plot
plotTheme <- theme(panel.grid = element_blank(),
                    axis.text.x  = element_text(size=10,color="black"),
                    axis.text.y  = element_text(size=10,color="black"),
                    axis.title.x = element_text(size = 10),
                    axis.title.y = element_text(size = 10),
                    panel.background  = element_blank(),
                    axis.line.x = element_line(color = "black"),
                    axis.line.y = element_line(color = "black"))
cols <- c("D0" = "black", "M3" = "red", "M12" = "coral", "M24" = "darkblue")
cPlot <- ggplot(data = corDF, aes(x = mTH17, y = log2(mCarnitines), color = Timepoint)) +
         scale_color_manual(values = cols) +
         geom_point(size = 2) +
         geom_smooth(method = "lm", color = "black") +
         labs(x = "TH17 (slea z-score)", y = "log2(mean carnitine concentration)") +
         plotTheme
outFile <- "carnitines and Th17 correlation.pdf"
pdf(file  = file.path(figDir, outFile), width = 4, height = 3, useDingbats = F)
print(cPlot)
dev.off()
