# load required packages
suppressPackageStartupMessages(library(package = "vegan"))
suppressPackageStartupMessages(library(package = "flux"))
suppressPackageStartupMessages(library(package = "dplyr"))
suppressPackageStartupMessages(library(package = "tidyr"))
suppressPackageStartupMessages(library(package = "tibble"))
suppressPackageStartupMessages(library(package = "readxl"))
suppressPackageStartupMessages(library(package = "ggplot2"))
suppressPackageStartupMessages(library(package = "mixOmics"))
suppressPackageStartupMessages(library(package = "igraph"))
suppressPackageStartupMessages(library(package = "lars"))
suppressPackageStartupMessages(library(package = "sva"))
suppressPackageStartupMessages(library(package = "phyloseq"))
suppressPackageStartupMessages(library(package = "ape"))
suppressPackageStartupMessages(library(package = "ggdendro"))
suppressPackageStartupMessages(library(package = "ggtree"))

######################################
# Define the global options 
######################################
options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)
figDir = "Manuscript/Figures_20190826"

#################### PATHSEQ ############################
# combine batch1 and batch2 data frames into one data frame
b1 <- read_excel("Microbiome_B1.xlsx", sheet = 1) %>%
      as.data.frame()
b2 <- read_excel("Microbiome_B2.xlsx", sheet = 1) %>%
      as.data.frame()
datDF <- merge(b1, b2, by = "Phylum", all.x = TRUE, all.y = TRUE) %>%
         mutate(Classification = ifelse(`Classification.x` %in% NA,
                                        `Classification.y`,
                                        `Classification.x`),
                Superkingdom = ifelse(`Superkingdom.x` %in% NA,
                                      `Superkingdom.y`,
                                      `Superkingdom.x`),
                Group = ifelse(`Group.x` %in% NA, `Group.y`, `Group.x`)) %>%
         dplyr::select(-`Classification.x`,
                       -`Classification.y`,
                       -`Superkingdom.x`,
                       -`Superkingdom.y`,
                       -`Group.x`,
                       -`Group.y`)

# make annotation data frame
# column annotation / sample annotation
CannotDF <- datDF %>%
            dplyr::select(-Superkingdom, -Phylum, -Classification, -Group)
CannotDF <- data.frame(SampleID = colnames(CannotDF)) %>%					
			mutate(DonorID = gsub(pattern = "(.+)-.+_.+",
                                  replacement = "\\1",
								  SampleID),
                   timePoint = gsub(pattern = ".+-(.+)_.+",
						    		replacement = "\\1",
						    		SampleID),
                   batch = gsub(pattern = ".+_(.+)",
                                replacement = "\\1",
                                SampleID)) %>%
            column_to_rownames(var = "SampleID")

# row annotation / pathogen annotation
RannotDF <- datDF %>%
			dplyr::select(Superkingdom, Phylum, Classification, Group) %>%
            column_to_rownames(var = "Phylum")

## FPKM data frame
table(datDF$Phylum == rownames(RannotDF))
table(datDF$Classification == RannotDF$Classification)

matFPKM <- datDF %>% dplyr::select(-Superkingdom, -Phylum, -Classification, -Group)
matFPKM[is.na(matFPKM)] <- 0
rownames(matFPKM) <- rownames(RannotDF)
table(colnames(matFPKM) == rownames(CannotDF))

## relative abundance data frame = abundance/total abundance (per sample as columns)
matRA <- t(t(matFPKM) / colSums(matFPKM))
table(colnames(matFPKM) == colnames(matRA))
table(rownames(matFPKM) == rownames(matRA))


################## CHECKING FOR BATCH EFFECTS ##################
# FPKM (abundance) and RA (relative abundance) data frame
# set metric as "FPKM" (or) "RA"
metric = "FPKM"
    if(metric == "FPKM") {
        mat <- matFPKM
    } else {
        mat <- matRA
    }

#######  Beta-diversity (Bray-Curtis distance) : B1 + B2
# plot DENDROGRAM based on bray-curtis distance
# calculate bray-curtis distance
distMat <- vegdist(as.matrix(t(mat)), method = "bray")

# calculate average BC distance between timepoints (across donors)
# eg: avg(Distance) : D0 and M3 for RAL03, RAL04 ...etc
#      avg(Distance) : D0 and M4 for RAL03, RAL04 ...etc
meanDF <- distMat %>%
          as.matrix() %>%
          as.data.frame() %>%
          rownames_to_column() %>%
          gather(colname, value, -rowname) %>%
          mutate(r.donor = gsub(pattern = "([^\\-]+)\\-.+_.+", replacement = "\\1", rowname),
                 c.donor = gsub(pattern = "([^\\-]+)\\-.+_.+", replacement = "\\1", colname),
                 rowname = gsub(pattern = "[^\\-]+\\-(.+_.+)", replacement = "\\1", rowname),
                 colname = gsub(pattern = "[^\\-]+\\-(.+_.+)", replacement = "\\1", colname)) %>%
         filter(r.donor == c.donor) %>%
         group_by(rowname, colname) %>%
         summarize(mu = mean(value)) %>%
         spread(colname, mu) %>%
         as.data.frame()
rownames(meanDF) <- meanDF$rowname
meanDF$rowname <- NULL

# plot dendrogram
dMat <- as.dist(meanDF)
hC <- hclust(dMat, method = "average")
outFile <- "Bray-Curtis Distance_Dendrogram_B1 and B2.pdf"
pdf(file = file.path(figDir, outFile), width = 4, height = 5)
plot(hC, main = "PathSeq pre-batch correction",
     ylab = "Bray-Curtis Dissimilarity", xlab = NULL)
dev.off()



############ using "combat" to perfrom batch correction #########
metric = "FPKM"
    if(metric == "FPKM") {
        mat <- matFPKM
    } else {
        mat <- matRA
    }

# adjusting for covariates (timepoint and donor)
mod <- model.matrix(~DonorID+timePoint, data = CannotDF)
cb <- ComBat(dat = as.matrix(mat), batch = CannotDF$batch, mod = mod)

#  Beta-diversity (Bray-Curtis distance)
distMat <- vegdist(as.matrix(t(cb)), method = "bray")
meanDF <- distMat %>%
          as.matrix() %>%
          as.data.frame() %>%
          rownames_to_column() %>%
          gather(colname, value, -rowname) %>%
          mutate(r.donor = gsub(pattern = "([^\\-]+)\\-.+_.+", replacement = "\\1", rowname),
                 c.donor = gsub(pattern = "([^\\-]+)\\-.+_.+", replacement = "\\1", colname),
                 rowname = gsub(pattern = "[^\\-]+\\-(.+_.+)", replacement = "\\1", rowname),
                 colname = gsub(pattern = "[^\\-]+\\-(.+_.+)", replacement = "\\1", colname)) %>%
          filter(r.donor == c.donor) %>%
          group_by(rowname, colname) %>%
          summarize(mu = mean(value)) %>%
          spread(colname, mu) %>%
          as.data.frame()
rownames(meanDF) <- meanDF$rowname
meanDF$rowname <- NULL

# plot dendrogram
dMat <- as.dist(meanDF)
hC <- hclust(dMat, method = "average")
outFile <- "Bray-Curtis Distance_Dendrogram_combat.pdf"
pdf(file = file.path(figDir, outFile), width = 4, height = 5)
plot(hC, main = "PathSeq post-batch correction",
     ylab = "Bray-Curtis Dissimilarity", xlab = NULL)
dev.off()

### assess if M9 batch1 and M9 batch2, post batch correction are correlated to each other
mat <- matRA
cb <- ComBat(dat = as.matrix(mat), batch = CannotDF$batch, mod = mod)
b1 <- cb[, grep("M9_B1", colnames(cb)), drop = FALSE]
b2 <- cb[, grep("M9_B2", colnames(cb)), drop = FALSE]
table(rownames(b1) == rownames(b2))
diag(cor(b1, b2))
corDF <- cor(b1, b2)
corDF[upper.tri(corDF)] <- NA
corDF[lower.tri(corDF)] <- NA

# plot correlation heatmap
outFile <- "Correlation between M9 B1 and B2.pdf"
pdf(file = file.path(figDir,outFile), width  = 6, height = 6)
colorPalette <- c("orange", "brown", "red")
colorPalette <- colorRampPalette(colors = colorPalette)(100)
pheatmap(mat  =  corDF,
         color  =  colorPalette,
         cellwidth  =  15,
         cellheight =    15,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_rownames =  TRUE,
         show_colnames = TRUE,
         border_color = "grey",
         fontsize  =  12)
dev.off()

### average M9 Batch1 and Batch2 samples
mat <- matFPKM
cb <- ComBat(dat = as.matrix(mat), batch = CannotDF$batch, mod = mod)
cbFPKM <- cb %>% as.data.frame() %>%
          rownames_to_column() %>%
          gather(SampleID, value, -rowname) %>%
          mutate(TIME_POINT = gsub(".+-(.+)_.+", "\\1", SampleID),
                 Batch = gsub(".+-.+_(.+)", "\\1", SampleID),
                 SubjectID = gsub("(.+)-.+_.+", "\\1", SampleID)) %>%
          group_by(SubjectID, TIME_POINT, rowname) %>%
          summarize(mValue = mean(value, na.rm = TRUE)) %>%
          mutate(SampleID = interaction(SubjectID, TIME_POINT, sep = "_")) %>%
          as.data.frame() %>%
          dplyr::select(SampleID, rowname, mValue) %>%
          spread(SampleID, mValue) %>%
          column_to_rownames(var = "rowname")
cbFPKM[cbFPKM < 0] <- 0

# relative abundance of batch corrected data
cbRA <- t(t(cbFPKM) / colSums(cbFPKM))

##################################################
### FIG 2A
##################################################
#######  Beta-diversity (Bray-Curtis distance)
mat <- cbFPKM
distMat <- vegdist(as.matrix(t(mat)), method = "bray")
meanDF <- distMat %>%
          as.matrix() %>%
          as.data.frame() %>%
          rownames_to_column() %>%
          gather(colname, value, -rowname) %>%
          mutate(r.donor = gsub(pattern = "(.+)_.+", replacement = "\\1", rowname),
                 c.donor = gsub(pattern = "(.+)_.+", replacement = "\\1", colname),
                 rowname = gsub(pattern = ".+_(.+)", replacement = "\\1", rowname),
                 colname = gsub(pattern = ".+_(.+)", replacement = "\\1", colname)) %>%
         filter(r.donor == c.donor) %>%
         group_by(rowname, colname) %>%
         summarize(mu = mean(value)) %>%
         spread(colname, mu) %>%
         as.data.frame()
rownames(meanDF) <- meanDF$rowname
meanDF$rowname <- NULL

# plot dendrogram
dMat <- as.dist(meanDF)
hC <- hclust(dMat, method = "average")
outFile <- "Bray-Curtis Distance_Dendrogram Final FPKM.pdf"
pdf(file = file.path(figDir, outFile), width = 4, height = 5)
plot(hC, main = "Circulating plasma microbiome on ART",
     ylab = "Bray-Curtis Dissimilarity", xlab = NULL)
dev.off()


# dataframe of relative abundance per group
#datDF <- cbFPKM %>% as.data.frame() %>%
#        rownames_to_column() %>%
#        mutate(Group = RannotDF$Group[match(rowname,table = rownames(RannotDF))])
#groups <- c("Bacteria", "Viruses", "Worms","Fungi") ## No parasites bec just 1 present
#mbLS <- lapply(groups, function(GP) {
#           m1 <- datDF %>% filter(Group %in% GP) %>%
#                 dplyr::select(-Group) %>%
#                 column_to_rownames(var = "rowname")
#           m1RA <- t(t(m1) / colSums(m1, na.rm = TRUE))
#        })
#cbRA_Group <- do.call(rbind, mbLS)
#cbRA_Group[cbRA_Group %in% "NaN"] <- 0
#cbRA_Group <- rbind(cbRA_Group, cbRA["Apicomplexa", , drop = FALSE])


##################################################
### FIG 2B
##################################################
## Alpha-diversity (Shannon Diversity)
mat <- cbFPKM
plotDF <- as.matrix(mat)
# calculate shannon index
shannonDF <- vegan::diversity(plotDF,
                              index = "shannon",
                              MARGIN = 2,
                              base = exp(1))
shannonDF <- data.frame(ShannonDiversity = shannonDF) %>%
            rownames_to_column() %>%
            mutate(DonorID = gsub("(.+)_.+", "\\1", rowname),
                   TimePoint = gsub(".+_(.+)", "\\1", rowname)) %>%
            dplyr::select(-rowname) %>%
            mutate(DonorID = c(RAL03 = "2213",
                               RAL04 = "2214",
                               RAL07 = "2215",
                               RAL08 = "2216",
                               RAL10 = "2217",
                               RAL11 = "2218",
                               RAL12 = "2219",
                               RAL13 = "2220",
                               RAL14 = "2221",
                               RAL15 = "2222",
                               RAL16 = "2223")[DonorID])
# define donor colors
dCols <- c("2213" = "#E41A1C",
           "2214" = "#377EB8",
           "2215" = "#4DAF4A",
           "2216" = "#984EA3",
           "2217" = "#FF7F00",
           "2218" = "#000000",
           "2219" = "#A65628",
           "2220" = "#F781BF",
           "2221" = "#999999",
           "2222" = "#FF00FF",
           "2223" = "#0000FF")
xorder <- c("D0","M3","M4","M9","M12","M24")

#  plot of alpha diversity
# median DF
xorder <- c("D0", "M3", "M4", "M9", "M12", "M24")
mDF <- shannonDF %>%
       group_by(TimePoint) %>%
       summarize(m = median(ShannonDiversity)) %>%
       as.data.frame()
mDF <- mDF[order(match(mDF$TimePoint, table = xorder)), ]

# plot
outFile <- "ShannonIndex.pdf"
pdf(file = file.path(figDir, outFile), height = 2.5, width = 3.5, useDingbats = FALSE)
pP <- ggplot(data = shannonDF,
             mapping = aes(x = TimePoint, y = ShannonDiversity)) +
      geom_point(aes(color = DonorID), size = 2) +
      scale_color_manual(values = dCols, name = "DonorID") +
      labs(x = "Duration on ART", y = "Shannon diversity index") +
      geom_path(aes(x = TimePoint, y = m), data = mDF, linetype = 2, group = 1) +
      scale_x_discrete(limits = xorder) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(size = 9, color = "black", angle = 45),
            axis.title.x = element_text(size = 10, color = "black"),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title.y = element_text(size = 10, color = "black"),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 10))
print(pP)
dev.off()

### correlation between plasma cytokine clusters and alpha diversity
## refer to "clusterDF" data frame in "trace_PlasmaCytokines.R:
for(c in unique(clusterDF$Cluster)) {
    print(c)
    cyt <- clusterDF %>%
           filter(Cluster %in% c) %>%
           mutate(DonorID_TimePoint = interaction(SampleID, TimePoint, sep = "_")) %>%
           select(-SampleID, -TimePoint)
    div <- shannonDF %>% mutate(DonorID_TimePoint = interaction(DonorID, TimePoint, sep = "_")) %>%
           select(-DonorID, -TimePoint)
    cytdiv <- merge(cyt, div, by = "DonorID_TimePoint")
    print(cor.test(cytdiv$meanFC, cytdiv$ShannonDiversity))
}



##################################################
### FIG 2C
##################################################
### per timepoint, plot each phylum abundance
datDF <- cbFPKM
mbLS <- lapply(unique(CannotDF$DonorID), function(DONOR) {
         m1 <- datDF %>% as.data.frame() %>%
               rownames_to_column() %>%
               mutate(Classification = RannotDF$Classification[match(rowname,
                                                                   table = rownames(RannotDF))])
         m1 <- m1[, c("Classification",
                      colnames(m1)[grep(DONOR, colnames(m1))]), drop = FALSE] %>%
               gather(SampleID, value, -Classification) %>%
               mutate(TimePoint = gsub(pattern = ".+\\_(.+)", "\\1", SampleID)) %>%
               group_by(Classification, TimePoint) %>%
               mutate(total = sum(value, na.rm = TRUE),
                      DonorID = DONOR) %>%
               as.data.frame() %>%
               dplyr::select(Classification, TimePoint, total, DonorID) %>%
               unique(.)
           })
plotDF <- do.call(rbind, mbLS) %>%
          mutate(DonorID = c(RAL03 = "2213",
                             RAL04 = "2214",
                             RAL07 = "2215",
                             RAL08 = "2216",
                             RAL10 = "2217",
                             RAL11 = "2218",
                             RAL12 = "2219",
                             RAL13 = "2220",
                             RAL14 = "2221",
                             RAL15 = "2222",
                             RAL16 = "2223")[DonorID])
plotDF$TimePoint <- factor(plotDF$TimePoint, levels = c("D0","M3","M4","M9","M12","M24"))
plotDF$Classification <- factor(plotDF$Classification, levels = c("Actinobacteria",
                                                                  "Bacteroidetes",
                                                                  "Firmicutes",
                                                                  "Proteobacteria",
                                                                  "Other bacteria",
                                                                  "Fungi",
                                                                  "Parasites",
                                                                  "Worms - Flat",
                                                                  "Worms - Ringed",
                                                                  "Worms - Round",
                                                                  "GBV-c",
                                                                  "HIV",
                                                                  "Other viruses"))

# Jitter plot
pdf(file = file.path(figDir, "Pathogens per timepoint.pdf"),
     height = 10, width = 7, useDingbats = F)
plotJit <- ggplot(data = plotDF,
                  mapping = aes(x = Classification, y = log2(total))) +
          geom_boxplot(fill = "grey", outlier.colour = NA, size = 0.3) +
          geom_jitter(mapping = aes(color = DonorID), size = 1, width = 0.1) +
          scale_color_manual(values = dCols) +
          scale_y_log10() +
          labs(x = NULL, y = "TPM") +
          facet_wrap(~TimePoint, ncol = 1) +
          theme_bw() +
          theme(panel.grid = element_blank(),
                axis.text.x = element_text(color = "black", size = 8, angle = 90, hjust = 1),
                axis.text.y = element_text(color = "black", size = 5))
print(plotJit)
dev.off()


#### Pie chart of Relative Abundance
datDF <- cbRA %>% as.data.frame() %>%
         rownames_to_column() %>%
         mutate(Classification = RannotDF$Classification[match(rowname,
                                                               table = rownames(RannotDF))])
mbLS <- lapply(unique(CannotDF$DonorID), function(DONOR) {
            m1 <- datDF[, c("Classification",
                            colnames(datDF)[grep(DONOR, colnames(datDF))]), drop = FALSE] %>%
                  gather(SampleID, value, -Classification) %>%
                  mutate(TimePoint = gsub(pattern = ".+_(.+)", "\\1", SampleID)) %>%
                  group_by(Classification, TimePoint) %>%
                  mutate(total = sum(value, na.rm = TRUE),
                         DonorID = DONOR) %>%
                  as.data.frame() %>%
                  dplyr::select(Classification, TimePoint, total, DonorID) %>%
                  unique(.)
            })
plotDF <- do.call(rbind, mbLS) %>%
          group_by(TimePoint, Classification) %>%
          summarize(mPath = mean(total)) %>% as.data.frame()

cols <- c("Proteobacteria" = "red",
          "Actinobacteria" = "coral",
            "Bacteroidetes" = "coral2",
            "Firmicutes" = "coral3",
            "Other bacteria" =  "saddlebrown",

            "GBV-c" = "dodgerblue",
            "HIV" = "navyblue",
            "Other viruses" = "cyan",

            "Worms - Flat" = "orange",
            "Worms - Ringed" = "#fff1aa",
            "Worms - Round" = "#ffd27f",

            "Fungi" = "darkorchid",

            "Parasites" = "darkolivegreen3")
plotDF$Classification <- factor(plotDF$Classification, levels = names(cols))
plotDF$TimePoint <- factor(plotDF$TimePoint, levels = c("D0","M3","M9","M12","M4","M24"))
pP <- ggplot(data = plotDF,
             mapping = aes(x = "", y = mPath, fill = Classification)) +
      geom_bar(width = 1, stat = "identity") +
      scale_fill_manual(values = cols) +
      coord_polar("y") +
      labs(x = NULL, y = NULL) +
      facet_wrap(~TimePoint, ncol = 4) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(color = "black", size= 8),
            axis.text.y = element_text(color = "black", size = 8))
outFile <- "Microbiome Pie Chart Relative Abundance.pdf"
pdf(file = file.path(figDir, outFile), width = 11, height = 4)
print(pP)
dev.off()


##################################################
### FIG S3B
##################################################
## Stacked bar plot of the total microbiome per donor (Relative abundance)
datDF <- cbRA %>% as.data.frame() %>%
         rownames_to_column() %>%
         mutate(Classification = RannotDF$Classification[match(rowname,
                                                                table = rownames(RannotDF))])
mbLS <- lapply(unique(CannotDF$DonorID), function(DONOR) {
            m1 <- datDF[, c("Classification",
                            colnames(datDF)[grep(DONOR, colnames(datDF))]), drop = FALSE] %>%
                  gather(SampleID, value, -Classification) %>%
                  mutate(TimePoint = gsub(pattern = ".+_(.+)", "\\1", SampleID)) %>%
                  group_by(Classification, TimePoint) %>%
                  mutate(total = sum(value, na.rm = TRUE),
                         DonorID = DONOR) %>%
                  as.data.frame() %>%
                  dplyr::select(Classification, TimePoint, total, DonorID) %>%
                  unique(.)
            })
plotDF <- do.call(rbind, mbLS) %>%
          mutate(DonorID = c(RAL03 = "2213",
                             RAL04 = "2214",
                             RAL07 = "2215",
                             RAL08 = "2216",
                             RAL10 = "2217",
                             RAL11 = "2218",
                             RAL12 = "2219",
                             RAL13 = "2220",
                             RAL14 = "2221",
                             RAL15 = "2222",
                             RAL16 = "2223")[DonorID])
plotDF$Classification <- factor(plotDF$Classification, levels = names(cols))
xorder <- c("D0", "M3", "M4", "M9", "M12", "M24")

# plot
lP <- ggplot(data = plotDF,
             mapping = aes(x = factor(TimePoint),
                            y = total,
                            by = Classification,
                            fill = Classification)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = cols) +
      scale_x_discrete(limits = xorder) +
      labs(x = "Duration on ART", y = "Frequency") +
      facet_wrap(~ DonorID,  ncol = 6, scales = "free") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(color = "black", size= 8),
            axis.text.y = element_text(color = "black", size = 8))
outFile <- "Microbiome Stacked Bar plot Relative Abundance.pdf"
pdf(file = file.path(figDir, outFile), width = 11, height = 4)
print(lP)
dev.off()

####################################################################
###### Identifying the phyla associated with alpha diversity
####################################################################
### Keep only the most varying pathogens for down stream analysis
varPt <- data.frame(v = apply(cbFPKM, 1, var)) %>%
         rownames_to_column() %>%
         filter(v >= median(v) & !rowname %in% "Retroviridae") %>%
         .$rowname
cbFPKM2 <- cbFPKM[varPt, , drop = FALSE]
cbRA_2 <- cbRA[varPt, , drop = FALSE]
datDF <- cbRA_2
df1 <- datDF %>% as.data.frame() %>%
       rownames_to_column() %>%
       gather(SampleID, value, -rowname)

# shannonDF
mat <- cbFPKM
plotDF <- as.matrix(mat)

# calculate shannon index
shannonDF <- vegan::diversity(plotDF, index = "shannon", MARGIN = 2, base = exp(1))
shannonDF <- data.frame(ShannonDiversity = shannonDF) %>%
             rownames_to_column() %>%
             mutate(DonorID = gsub("(.+)_.+", "\\1", rowname),
                    TimePoint = gsub(".+_(.+)", "\\1", rowname)) %>%
             dplyr::select(-rowname)
df2 <- shannonDF
rownames(df2) <- interaction(df2$DonorID, df2$TimePoint, sep = "_")
cLS <- lapply(unique(df1$rowname), function(PATH) {
            print(PATH)
            df3 <- df1 %>% filter(rowname %in% PATH) %>%
                   column_to_rownames(var = "SampleID")
            df2 <- df2[rownames(df3), ,drop = FALSE]
            table(rownames(df2)==rownames(df3))
            cT <- cor.test(df2$ShannonDiversity, df3$value, method = "spearman")
            cDF <- data.frame(Pathogen = PATH,
                              pVal = cT$p.value,
                              rho = cT$estimate)
            return(value = cDF)
        })
corDF <- do.call(rbind, cLS) %>%
         mutate(FDR = p.adjust(pVal, method = "BH")) %>%
         filter(FDR < 0.05 & abs(rho) > 0.5) %>%
         arrange(rho)


##################################################
### FIG 2D
##################################################
### calculate correlation between pathogens
datDF <- cbRA_2
# variable 'hOrder' taken from corDF$Pathogen
hOrder <- c("Proteobacteria", "Actinobacteria", "Bacteroidetes", "Firmicutes")
corMat <- datDF[hOrder, , drop = FALSE]

# cor.tests
cor.test(corMat["Proteobacteria",], corMat["Actinobacteria", ], method = "spearman")
cor.test(corMat["Proteobacteria",], corMat["Bacteroidetes", ], method = "spearman")
cor.test(corMat["Proteobacteria",], corMat["Firmicutes", ], method = "spearman")

## correlation matrix
corMat2 <- cor(t(corMat), method = "spearman")
corMat2[lower.tri(corMat2)] <- NA

# plot correlation heatmap
outFile <- "Correlation between pathogens.pdf"
pdf(file = file.path(figDir,outFile), width  = 6, height = 6)
colorPalette <- c("blue", "white", "red")
colorPalette <- colorRampPalette(colors = colorPalette)(100)
breaks <- c(-1,
            seq(from = min(corMat2, na.rm = T),
                to = abs(min(corMat2, na.rm = T)),
                length.out = 99),
            1)
pheatmap(mat = corMat2,
         breaks = breaks,
         color = colorPalette,
         cellwidth = 15,
         cellheight = 15,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         border_color = "grey",
         fontsize = 12)
dev.off()



##################################################
### FIG S3C
##################################################
# Jitter plot comapring the pathogens associated with diversity across time
for(PATH in corDF$Pathogen) {
d1 <- df2
d2 <- df1 %>% filter(rowname %in% PATH) %>% column_to_rownames(var = "SampleID")
d1 <- d1[rownames(d2), ,drop = FALSE]
table(rownames(d1)==rownames(d2))
splotDF <- data.frame(AD = d1$ShannonDiversity,
                      Pathogen = d2$value,
                      TP = d1$TimePoint)
splotDF$TP <- factor(splotDF$TP, levels = c("D0","M3","M4","M9","M12","M24"))

# Jitter plot across time
outFile2 <- paste("Comparing ", PATH, " across timepoints.pdf", sep = "")
pdf(file = file.path(figDir, outFile2),height = 4, width = 5, useDingbats = F)
plotJit <- ggplot(data = splotDF,
                  mapping = aes(x = TP, y = Pathogen)) +
           geom_boxplot(fill = "grey", outlier.colour = NA, size = 0.1) +
           geom_jitter(size = 1) +
           labs(x = "TimePoint",
                y = paste(PATH, "Relative abundance", sep = " ")) +
           theme_bw() +
           theme(panel.grid = element_blank(),
                 axis.text.x = element_text(color = "black", size = 10),
                 axis.text.y = element_text(color = "black", size = 10))
print(plotJit)
dev.off()
}

############################################################
#################### Breakdown analysis ############################
############################################################
# combine batch1 and batch2 data frames into one data frame
# batch 1
b1 <- read_excel("Microbiome_fullbreakdown_B1.xlsx", sheet = 1) %>% as.data.frame()
colnames(b1)[grep("^RAL", colnames(b1))] <- as.character(interaction(colnames(b1)[grep("^RAL",
                                                                                    colnames(b1))],
                                                                    "B1", sep = "_"))
rownames(b1) <- apply(b1[, -grep("^RAL", colnames(b1))], 1, paste, collapse="__")
b1 <- b1[, grep("^RAL", colnames(b1)), drop = F]

# batch 2
b2 <- read_excel("Microbiome_fullbreakdown_B2.xlsx", sheet = 1) %>% as.data.frame()
colnames(b2)[grep("^RAL", colnames(b2))] <- as.character(interaction(colnames(b2)[grep("^RAL",
                                                                                    colnames(b2))],
                                                                     "B2", sep = "_"))
rownames(b2) <- apply(b2[, -grep("^RAL", colnames(b2))], 1, paste, collapse="__")
b2 <- b2[, grep("^RAL", colnames(b2)), drop = F]

# merge batches
datDF <- merge(b1, b2, by = "row.names", all.x = TRUE, all.y = TRUE) %>%
         column_to_rownames(var = "Row.names")

# replace NA's with zeros
matFPKM <- datDF
matFPKM[is.na(matFPKM)] <- 0
filter <- apply(matFPKM, 1, function(x) mean(x) > 0)
matFPKM <- matFPKM[filter, ]
colnames(matFPKM) <- gsub("Mo", "M", colnames(matFPKM))
# make annotation data frame
# column annotation / sample annotation
CannotDF <- data.frame(SampleID = colnames(matFPKM)) %>%
            mutate(DonorID = gsub(pattern = "(.+)-.+_.+",
                                  replacement = "\\1",
                                  SampleID),
                  timePoint = gsub(pattern = ".+-(.+)_.+",
                                   replacement = "\\1",
                                   SampleID),
                  batch = gsub(pattern = ".+-.+_(.+)",
                               replacement = "\\1",
                               SampleID)) %>%
            column_to_rownames(var = "SampleID")

# row annotation / pathogen annotation
RannotDF <- data.frame(Rowname = rownames(matFPKM)) %>%
            mutate(Rowname2 = Rowname) %>%
            separate(Rowname, c("Superkingdom",
                                "Phylum",
                                "Class",
                                "Order",
                                "Family",
                                "Genus",
                                "Species"), "__") %>%
            column_to_rownames(var = "Rowname2")
table(rownames(matFPKM) == rownames(RannotDF))

# Batch correction
# adjusting for covariates (timepoint and donor)
mod <- model.matrix(~DonorID+timePoint, data = CannotDF)
cb <- ComBat(dat = as.matrix(matFPKM), batch = CannotDF$batch, mod = mod)
#  Beta-diversity (Bray-Curtis distance)
distMat <- vegdist(as.matrix(t(cb)), method = "bray")
meanDF <- distMat %>%
          as.matrix() %>%
          as.data.frame() %>%
          rownames_to_column() %>%
          gather(colname, value, -rowname) %>%
          mutate(r.donor = gsub(pattern = "([^\\-]+)\\-.+_.+", replacement = "\\1", rowname),
                 c.donor = gsub(pattern = "([^\\-]+)\\-.+_.+", replacement = "\\1", colname),
                 rowname = gsub(pattern = "[^\\-]+\\-(.+_.+)", replacement = "\\1", rowname),
                 colname = gsub(pattern = "[^\\-]+\\-(.+_.+)", replacement = "\\1", colname)) %>%
         filter(r.donor == c.donor) %>%
         group_by(rowname, colname) %>%
         summarize(mu = mean(value)) %>%
         spread(colname, mu) %>%
         as.data.frame()
rownames(meanDF) <- meanDF$rowname
meanDF$rowname <- NULL

### average M9 Batch1 and Batch2 samples
cbFPKM <- cb %>% as.data.frame() %>%
          rownames_to_column() %>%
          gather(SampleID, value, -rowname) %>%
          mutate(TIME_POINT = gsub(".+-(.+)_.+", "\\1", SampleID),
                 Batch = gsub(".+-.+_(.+)", "\\1", SampleID),
                 SubjectID = gsub("(.+)-.+_.+", "\\1", SampleID)) %>%
         group_by(SubjectID, TIME_POINT, rowname) %>%
         summarize(mValue = mean(value, na.rm = TRUE)) %>%
         mutate(SampleID = interaction(SubjectID, TIME_POINT, sep = "_")) %>%
         as.data.frame() %>%
         dplyr::select(SampleID, rowname, mValue) %>%
         spread(SampleID, mValue) %>%
         column_to_rownames(var = "rowname")
cbFPKM[cbFPKM < 0] <- 0

# relative abundance
# relative abundance of batch corrected data
cbRA <- t(t(cbFPKM) / colSums(cbFPKM)) %>% as.data.frame()


##################################################
### FIG 3A
##################################################
## Per phylum, order by the most abundant Genus : all in one plot
phyla <- c("Proteobacteria", "Actinobacteria", "Bacteroidetes", "Firmicutes")

    # top 10 abundant per genus
    gOrder <- cbFPKM %>%
              rownames_to_column() %>%
              gather(SampleID, FPKM, -rowname) %>%
              mutate(Phylum = RannotDF$Phylum[match(rowname, table = rownames(RannotDF))],
                     Genus = RannotDF$Genus[match(rowname, table = rownames(RannotDF))],
                     DonorID = gsub("(.+)_.+", "\\1", SampleID),
                     TimePoint = gsub(".+_(.+)", "\\1", SampleID)) %>%
              filter(Phylum %in% phyla) %>%
              dplyr::select(-rowname, -SampleID) %>%
              group_by(Phylum, Genus, DonorID, TimePoint) %>%
              summarize(s = sum(FPKM)) %>%
              as.data.frame() %>%
             # calculate median abundance per Genus across all donors and timepoints
              group_by(Phylum, Genus) %>%
              summarize(m = median(log10(s))) %>%
              as.data.frame() %>%
              filter(!m < 1) %>%
              group_by(Phylum) %>%
              top_n(10, m) %>%
              as.data.frame() %>%
              arrange(Phylum, desc(m))
gOrder <- gOrder[order(match(gOrder$Phylum, table = c("Proteobacteria",
                                                      "Actinobacteria",
                                                      "Bacteroidetes",
                                                      "Firmicutes"))), ] %>% .$Genus
    
    # bubble plot per donor across time
    plotDF <- cbFPKM %>%
              rownames_to_column() %>%
              gather(SampleID, FPKM, -rowname) %>%
              mutate(Phylum = RannotDF$Phylum[match(rowname, table = rownames(RannotDF))],
                     Genus = RannotDF$Genus[match(rowname, table = rownames(RannotDF))],
                     DonorID = gsub("(.+)_.+", "\\1", SampleID),
                     TimePoint = gsub(".+_(.+)", "\\1", SampleID)) %>%
              filter(Phylum %in% phyla) %>%
              dplyr::select(-rowname, -SampleID) %>%
              group_by(Phylum, Genus, DonorID, TimePoint) %>%
              summarize(s = sum(FPKM)) %>%
              as.data.frame() %>%
            # calculate median abundance per Genus across timepoints per donor
             group_by(Genus, DonorID) %>%
             summarize(m = median(log10(s))) %>%
             as.data.frame() %>%
            # filter on the top 10 most abundant
             filter(Genus %in% gOrder) %>%
             mutate(DonorID = c(RAL03 = "2213",
                                RAL04 = "2214",
                                RAL07 = "2215",
                                RAL08 = "2216",
                                RAL10 = "2217",
                                RAL11 = "2218",
                                RAL12 = "2219",
                                RAL13 = "2220",
                                RAL14 = "2221",
                                RAL15 = "2222",
                                RAL16 = "2223")[DonorID])

    outFile <- paste("Fig_3a.pdf", sep = "")
    
    bP <- ggplot(plotDF, aes(x = DonorID,
                             y = Genus,
                             size = m)) +
           geom_point(color = "black") +
           scale_y_discrete(limits = rev(gOrder)) +
           scale_size(range = c(0, 6)) +
           labs(y = "Genus") +
           theme_bw() +
           theme(panel.grid = element_blank(),
                 legend.title = element_text(size = 10),
                 axis.text.x = element_text(color = "black", size= 8, angle = 90),
                 axis.text.y = element_text(color = "black", size = 8))
    pdf(file = file.path(figDir,outFile), width = 8, height = 9, useDingbats = FALSE)
    print(bP)
    dev.off()

## save top 10 Genus per phylum
topGenus <- cbFPKM %>%
            rownames_to_column() %>%
            gather(SampleID, FPKM, -rowname) %>%
            mutate(Phylum = RannotDF$Phylum[match(rowname, table = rownames(RannotDF))],
                   Genus = RannotDF$Genus[match(rowname, table = rownames(RannotDF))],
                   DonorID = gsub("(.+)_.+", "\\1", SampleID),
                   TimePoint = gsub(".+_(.+)", "\\1", SampleID)) %>%
            filter(Phylum %in% phyla) %>%
            dplyr::select(-rowname, -SampleID) %>%
            group_by(Phylum, Genus, DonorID, TimePoint) %>%
            summarize(s = sum(FPKM)) %>%
            as.data.frame() %>%
            # calculate mean abundance per Genus across all donors and timepoints
            group_by(Phylum, Genus) %>%
            summarize(m = mean(s)) %>%
            as.data.frame() %>%
            group_by(Phylum) %>%
            top_n(10, m) %>%
            as.data.frame() %>%
            .$Genus


##################################################
### FIG 3B
##################################################
#### MIXOMICS : Bacterial genus and cytokines
# Genus data frame
genDF <- cbRA %>%
         rownames_to_column() %>%
         gather(SampleID, value, -rowname) %>%
         mutate(Genus = RannotDF$Genus[match(rowname, table = rownames(RannotDF))],
                Phylum = RannotDF$Phylum[match(rowname, table = rownames(RannotDF))],
                DonorID = gsub("(.+)_.+", "\\1", SampleID),
                TimePoint = gsub(".+_(.+)", "\\1", SampleID)) %>%
        filter(Phylum %in% phyla & Genus %in% topGenus) %>%
        dplyr::select(-rowname, -SampleID, -Phylum) %>%
        group_by(Genus, DonorID, TimePoint) %>%
        summarize(s = sum(value)) %>%
        as.data.frame() %>%
        mutate(SampleID = interaction(DonorID, TimePoint, sep = "_")) %>%
        dplyr::select(-DonorID, -TimePoint) %>%
        spread(SampleID, s) %>%
        column_to_rownames(var = "Genus")

# add ratio of Serratia/other Genus
otherGen <- genDF[setdiff(rownames(genDF), "Serratia"), ]
otherMean <- colMeans(otherGen)
serratia <- genDF["Serratia", , drop = F]
serratia_others <- log2(serratia/otherMean)
others_serratia <- log2(otherMean/serratia)
genDF1 <- rbind(genDF,
                serratiaRatio = serratia_others,
                serratiaRatio_rev = others_serratia)


# cytokine data frame
cytDF <- read_excel("PlasmaBiomarkers.xlsx", sheet = 1) %>% as.data.frame()
markers <- c("IL1b","IL6","IL8","MIP1b",
             "IL17A","TNFa",
             "IL2","IL4","IL5","IL7","IL10","IL12","IL13","GCSF","IFNg","MCP1")
Y <- cytDF %>%
     filter(TimePoint %in% c("d0", "M3","M4","M9","M12","M24")) %>%
     dplyr::select_(.dots = c("SampleID", "TimePoint",  markers)) %>%
     mutate(TimePoint = gsub("d", "D", TimePoint),
            DonorID_TimePoint = interaction(SampleID, TimePoint, sep = "_")) %>%
     dplyr::select(-SampleID, -TimePoint) %>%
     column_to_rownames(var = "DonorID_TimePoint")
Y <- Y[complete.cases(Y), ]

# mixomics
### Define X and Y
## X = Predictor matrix ; Y = response vector matrix
X <- genDF1
X <- t(X)
comm <- intersect(rownames(X), rownames(Y))
X <- X[comm, , drop = FALSE]
Y <- Y[comm, , drop = FALSE]
table(rownames(X) == rownames(Y))

# Fit a sparse regression model (spls object)
splsFull <- spls(X,
                 Y,
                 mode = "regression",
                 scale = TRUE)
# perform leave-one-out cross validation of the fitted sparse regression model, to identify the best principal components to project on
set.seed(seed = 1)
splsCV <- perf(splsFull,
               validation = "loo")
nComp <- which.max(splsCV$Q2.total)

# keep variables whose absolute sum of principal components > 0
keep.X <- apply(abs(splsFull$loadings$X), 1, sum) > 0
keep.Y <- apply(abs(splsFull$loadings$Y), 1, sum) > 0

# project variables on the best principal components
cord.X <- cor(splsFull$X[, keep.X], splsFull$variates$X[, 1:nComp], use = "pairwise")
cord.Y <- cor(splsFull$Y[, keep.Y], splsFull$variates$X[, 1:nComp], use = "pairwise")

# now calculate correlation between the two projected omics
simFullMat <- cord.X %*% t(cord.Y)
simFullMat <- round(simFullMat, 1)

# put a cut-off on correlation coefficient thats in the 95th quantile
rThreshold <- quantile(as.numeric(simFullMat), probs = 0.95)

# keep interactions >= rThreshold
simFullMat1 <- simFullMat %>%
               as.data.frame() %>%
               rownames_to_column() %>%
               dplyr::rename(Genus = rowname) %>%
               gather(cytokine, r, -Genus) %>%
               filter(abs(r) >= rThreshold) %>%
               mutate(Phylum = RannotDF$Phylum[match(Genus, table = RannotDF$Genus)])

### make network of the significant associations from above
graphDF <- simFullMat1 %>%
           mutate(dir = ifelse(r > 0, "pos", "neg"),
                  cGroup = ifelse(cytokine %in% c("IL1b","IL6","IL8","MIP1b"),
                                  "C1",
                           ifelse(cytokine %in% c("IL17A", "TNFa"), "C2","C3"))) %>%
           group_by(Genus,dir,cGroup) %>%
           summarize(m = mean(r)) %>%
           as.data.frame() %>%
           dplyr::select(-dir) %>%
           dplyr::rename(from = Genus, to = cGroup, r = m)

## make graph data frame
g <- graph_from_data_frame(graphDF, directed = FALSE)

# define edge colors based on correlation coefficient for graph
minVal <- min(graphDF$r)
maxVal <- max(graphDF$r)

# define color of edge based on correlation coefficient
colorPalette <- breaks <- NULL
if (0 > minVal & 0 < maxVal) {
    breakVal <- min(abs(range(graphDF$r)))
    breaks <- c(floor(minVal),
                seq(from = -breakVal,
                    to = breakVal,
                    length.out = 99),
                 ceiling(maxVal))
    colorPalette <- c("#D73027", "#FC8D59", "#FEE090", "#FFFFFF",
                     "#E0F3F8", "#91BFDB", "#4575B4")
    colorPalette <- colorRampPalette(colors = rev(colorPalette))(100)
        } else {
            if (0 <= minVal) {
                breaks <- seq(from = 0, to = maxVal, length.out = 101)
                colorPalette <- c("#D73027", "#FC8D59", "#FEE090","#FFFFFF")
                colorPalette <- colorRampPalette(colors = rev(colorPalette))(100)
            } else {
                breaks <- seq(from = minVal, to = 0, length.out = 101)
                colorPalette <- c("#FFFFFF", "#E0F3F8", "#91BFDB","#4575B4")
                colorPalette <- colorRampPalette(colors = rev(colorPalette))(100)
            }
        }
keyColor <- findInterval(graphDF$r, vec = breaks, all.inside = TRUE)
keyColor <- colorPalette[keyColor]
minLeg <- which(graphDF$r==minVal)[1]
maxLeg <- which(graphDF$r==maxVal)[1]

# save graph
outFile <- "Mixomics_Genus_Cytokines.pdf"
pdf(file = file.path(figDir, outFile), width = 6, height = 6, useDingbats = F)
plot(g,
     edge.color = keyColor,
     edge.width = 0.5,
     vertex.size = 8,
     vertex.label.color = "black",
     vertex.frame.color = NA,
     vertex.label.cex = 0.8)
legend(x = 0.5,
       y = -1.1,
       lwd = range(graphDF$r) * 5,
       col = c(keyColor[minLeg], keyColor[maxLeg]),
       legend = signif(range(graphDF$r), digits = 2),
       bty = "n",
       title = "correlation coefficient")
dev.off()


##################################################
### FIG 3C
##################################################
### Network of serrtia/others ~ cytokines and others/serratia ~ cytokines
graphDF <- simFullMat1 %>%
           filter(Genus %in% c("serratiaRatio")) %>%
           dplyr::select(-Phylum) %>%
           dplyr::rename(from = Genus, to = cytokine)

## make graph data frame
g <- graph_from_data_frame(graphDF, directed = FALSE)

# define edge colors based on correlation coefficient for graph
minVal <- min(graphDF$r)
maxVal <- max(graphDF$r)

# define color of edge based on correlation coefficient
colorPalette <- breaks <- NULL
if (0 > minVal & 0 < maxVal) {
    breakVal <- min(abs(range(graphDF$r)))
    breaks <- c(floor(minVal),
                seq(from = -breakVal,
                    to = breakVal,
                    length.out = 99),
                 ceiling(maxVal))
    colorPalette <- c("#D73027", "#FC8D59", "#FEE090", "#FFFFFF",
                     "#E0F3F8", "#91BFDB", "#4575B4")
    colorPalette <- colorRampPalette(colors = rev(colorPalette))(100)
        } else {
            if (0 <= minVal) {
                breaks <- seq(from = 0, to = maxVal, length.out = 101)
                colorPalette <- c("#D73027", "#FC8D59", "#FEE090","#FFFFFF")
                colorPalette <- colorRampPalette(colors = rev(colorPalette))(100)
            } else {
                breaks <- seq(from = minVal, to = 0, length.out = 101)
                colorPalette <- c("#FFFFFF", "#E0F3F8", "#91BFDB","#4575B4")
                colorPalette <- colorRampPalette(colors = rev(colorPalette))(100)
            }
        }
keyColor <- findInterval(graphDF$r, vec = breaks, all.inside = TRUE)
keyColor <- colorPalette[keyColor]
minLeg <- which(graphDF$r==minVal)[1]
maxLeg <- which(graphDF$r==maxVal)[1]

# save graph
outFile <- "Mixomics_GenusRatio_Cytokines.pdf"
pdf(file = file.path(figDir, outFile), width = 6, height = 6, useDingbats = F)
plot(g,
     edge.color = keyColor,
     edge.width = 0.5,
     vertex.size = 8,
     vertex.label.color = "black",
     vertex.frame.color = NA,
     vertex.label.cex = 0.8,
     layout = posMat)
legend(x = 0.5,
       y = -1.1,
       lwd = range(graphDF$r) * 5,
       col = c(keyColor[minLeg], keyColor[maxLeg]),
       legend = signif(range(graphDF$r), digits = 2),
       bty = "n",
       title = "correlation coefficient")
dev.off()


### Table of serratia/other bacteria genus and serratia/other flavivirus ratio values
pathogenratioDF <- cbind(t(genDF1[c("Serratia","serratiaRatio"), ]), t(genDF2["serratiaRatio", ]))
colnames(pathogenratioDF) <- c("Serratia", "serratia_otherbacteria")
pathogenratioDF <- pathogenratioDF %>%
                   as.data.frame() %>%
                   rownames_to_column(var = "DonorID_TimePoint")


##################################################
### FIG S4
##################################################
#### line plots of serratia, Pseudomonas and other genues separately
# Genus data frame
genDF <- cbRA %>%
         rownames_to_column() %>%
         gather(SampleID, value, -rowname) %>%
         mutate(Genus = RannotDF$Genus[match(rowname, table = rownames(RannotDF))],
                Phylum = RannotDF$Phylum[match(rowname, table = rownames(RannotDF))],
                DonorID = gsub("(.+)_.+", "\\1", SampleID),
                TimePoint = gsub(".+_(.+)", "\\1", SampleID)) %>%
        filter(Phylum %in% phyla & Genus %in% topGenus) %>%
        dplyr::select(-rowname, -SampleID, -Phylum) %>%
        group_by(Genus, DonorID, TimePoint) %>%
        summarize(s = sum(value)) %>%
        as.data.frame() %>%
        mutate(SampleID = interaction(DonorID, TimePoint, sep = "_")) %>%
        dplyr::select(-DonorID, -TimePoint) %>%
        spread(SampleID, s) %>%
        column_to_rownames(var = "Genus")

# serratia
serratia <- as.data.frame(t(genDF["Serratia", ]))

# Pseudomonas
pseudomonas <- as.data.frame(t(genDF["Pseudomonas", ]))

# other Genus (data frame - mean acorss all other genues per donor per timepoint)
otherGen <- genDF[setdiff(rownames(genDF), c("Serratia", "Pseudomonas")), ]
otherMean <- data.frame(OtherGenuses = colMeans(otherGen)) %>% rownames_to_column(var = "Row.names")
plotDF <- merge(merge(serratia, pseudomonas, by = "row.names"),
                      otherMean, by = "Row.names") %>%
          mutate(DonorID = gsub("(.+)_.+", "\\1", `Row.names`),
                 TimePoint = gsub(".+_(.+)", "\\1", `Row.names`)) %>%
          gather(pathogen, value, -DonorID, -TimePoint, -`Row.names`)
plotDF$TimePoint <- factor(plotDF$TimePoint,
                         levels = c("D0","M3","M4","M9","M12","M24"))
plotDF$pathogen <- factor(plotDF$pathogen,
                    levels = c("Serratia", "Pseudomonas", "OtherGenuses"))
cols <- c("Serratia" = "red",
          "Pseudomonas" = "blue",
          "OtherGenuses" = "green")
outFile <- "Serratia and other genuses over time_RA.pdf"
pdf(file = file.path(figDir, outFile), width = 10, height = 2, useDingbats = FALSE)
pDF <- ggplot(data = plotDF, aes(TimePoint,
                                 value*100,
                                 group = DonorID,
                                 color = pathogen)) +
       scale_color_manual(values = cols) +
       geom_point(size = 0.5) +
       geom_line(size = 0.2) +
       facet_wrap(~pathogen, scales = "free") +
       labs(x = "Duration on ART", y = "Relative Abundance") +
       theme_bw() +
       theme(panel.grid = element_blank(),
             axis.text.x = element_text(size = 8, color = "black"),
             axis.text.y = element_text(size = 8, color = "black"))
print(pDF)
dev.off()
