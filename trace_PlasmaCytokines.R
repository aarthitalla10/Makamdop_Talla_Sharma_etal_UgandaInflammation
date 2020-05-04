# load required packages
suppressPackageStartupMessages(library(package = "readxl"))
suppressPackageStartupMessages(library(package = "dplyr"))
suppressPackageStartupMessages(library(package = "tidyr"))
suppressPackageStartupMessages(library(package = "tibble"))
suppressPackageStartupMessages(library(package = "reshape2"))
suppressPackageStartupMessages(library(package = "ggplot2"))
suppressPackageStartupMessages(library(package = "edgeR"))
suppressPackageStartupMessages(library(package = "EDASeq"))
suppressPackageStartupMessages(library(package = "cluster"))
suppressPackageStartupMessages(library(package = "pheatmap"))
suppressPackageStartupMessages(library(package = "quantmod"))
#suppressPackageStartupMessages(library(package = "robustbase"))

######################################
# Define the global options 
######################################
options(stringsAsFactors = FALSE, useFancyQuotes = FALSE)

######################################
# Initializing directories
######################################
figDir = "Manuscript/Figures_20191129"


######################################
### FIG S1A
######################################
### Plasma VL
plotDF <- read_excel("RALTcellcounts.xlsx", sheet = 1) %>%
          as.data.frame() %>%
          select(DonorID, TimePoint, VL) %>%
          filter(!VL %in% NA) %>%
          mutate(TimePoint = gsub("d", "D", TimePoint)) %>%
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
order <- unique(plotDF$TimePoint)

# define donor colors
dCols <- c("2213" = "#E41A1C",
			"2214" = "#377EB8",
            "2215" = "#4DAF4A",
            "2216" = "#984EA3",
            "2217" = "#FF7F00",
            "2218" = "#000000",
            "2219" = "#A65628",
            "2220" = "#F781BF",
            "2221"  = "#999999",
            "2222" = "#FF00FF",
            "2223"  = "#0000FF")

#ggplot
outFile <- paste("ViralLoad.pdf")
pdf(file = file.path(figDir, outFile), width = 5, height = 4)
pDF <- ggplot(data = plotDF, aes(TimePoint, 
								 log10(VL),
								 group = DonorID,
                                 color = DonorID)) +
       scale_color_manual(values = dCols) +
       geom_point(size=4) +
       geom_line() +
       scale_x_discrete(limits = order) +
       labs(x = "Duration on ART", y = "log10(HIV RNA copies/ml)") +
       theme_bw() +
       theme(panel.grid = element_blank(),
			axis.text.x = element_text(size = 12, color = "black"),
			axis.title.x = element_text(size = 12, color = "black"),
			axis.text.y = element_text(size = 12, color = "black"),
			axis.title.y = element_text(size = 12, color = "black"))
print(pDF)
dev.off()

# lm
plotDF1 <- plotDF %>%
           mutate(TimePoint = gsub("\\D", "", TimePoint),
				  TimePoint = as.numeric(TimePoint))
fit <- lm(VL ~ TimePoint, data = plotDF1)



######################################
### FIG S1B
######################################
### CD4/CD8 Ratio
plotDF <- read_excel("RALTcellcounts.xlsx", sheet = 1) %>%
          as.data.frame() %>%
          select(DonorID, TimePoint, `CD4:CD8 Ratio`) %>%
          dplyr::rename(CD4CD8Ratio = `CD4:CD8 Ratio`) %>%
          mutate(CD4CD8Ratio = as.numeric(CD4CD8Ratio)) %>%
		  filter(!CD4CD8Ratio %in% NA) %>%
          mutate(TimePoint = gsub("d", "D", TimePoint)) %>%
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
order <- unique(plotDF$TimePoint)

#ggplot
outFile <- paste("CD4CD8 Ratio.pdf")
pdf(file = file.path(figDir, outFile), width = 5, height = 4)
pDF <- ggplot(data = plotDF, aes(TimePoint, 
								 CD4CD8Ratio,
								 group = DonorID,
								 color = DonorID)) +
       scale_color_manual(values = dCols) +
       geom_point(size=4) +
       geom_line() +
       scale_x_discrete(limits = order) +
       labs(x = "Duration on ART", y = "CD4:CD8 Ratio") +
       theme_bw() +
       theme(panel.grid = element_blank(),
             axis.text.x = element_text(size = 12, color = "black"),
             axis.title.x = element_text(size = 12, color = "black"),
             axis.text.y = element_text(size = 12, color = "black"),
             axis.title.y = element_text(size = 12, color = "black"))
print(pDF)
dev.off()

# lm
plotDF1 <- plotDF %>%
           mutate(TimePoint = gsub("\\D", "", TimePoint),
                  TimePoint = as.numeric(TimePoint))
fit <- lm(CD4CD8Ratio ~ TimePoint, data = plotDF1)
summary(fit)



########################## Plasma Cytokines ##########################
######################################
### FIG 1A
######################################
datDF1 <- read_excel("PlasmaBiomarkers.xlsx", sheet = 1) %>% as.data.frame() %>%
          filter(!TimePoint %in% "LT") %>%
          dplyr::rename(DonorID = SampleID) %>%
          dplyr::select(-sCD14, -iFABP, -LPS) %>%
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
			    			RAL16 = "2223")[DonorID],
		 DonorID_TimePoint = interaction(DonorID, TimePoint, sep="_")) %>%
		 dplyr::select(-VL, -CD4, -CD4_CD8_RATIO, -duration, -Ddimer)
rownames(datDF1) <- as.character(datDF1$DonorID_TimePoint)

# PCA on plasma cytokines
mat <- datDF1
mat <- mat %>% dplyr::select(-DonorID, -TimePoint, -DonorID_TimePoint)
mat2 <- as.data.frame(t(mat))

## imupte missing values for PCA visualization of data
# since TGF-b since is NA for M24
set.seed(seed = 1)
mat2_imp <- impute.knn(data = as.matrix(mat2),
                        k = 10,
                        rowmax = 0.5,
                        colmax = 0.8)
mat2_imp <- log2(mat2_imp$data)
exprsMat <- t(scale(t(mat2_imp))) %>% as.matrix
matpca <- prcomp(exprsMat, center = TRUE, scale = TRUE)
pcaDF <- as.data.frame(matpca[[2]]) %>% dplyr::select(PC1, PC2) %>%
         rownames_to_column() %>%
         dplyr::rename(SampleID = rowname) %>%
         mutate(TimePoint = gsub(".+_(.+)", "\\1", SampleID))
pcaDF$TimePoint <- factor(pcaDF$TimePoint, levels = unique(pcaDF$TimePoint))

# get rotation coordinates
pcaDF2 <- matpca$x %>%
          as.data.frame() %>%
          dplyr::select(PC1, PC2) %>%
          rownames_to_column(var = "cytokine")

# define timpoint colors
tCols <- c("d0" = "black",
           "d7" = "grey",
           "W2" = "#898989",
           "M3" = "red",
           "M4" = "lightblue",
           "M6" = "orange",
           "M9" = "yellow",
           "M12" = "coral",
           "M24" = "royalblue3")

# summary PCA
summaryPCA <- summary(matpca)$importance

### PCA biplot
x <- pcaDF %>%
     dplyr::select(SampleID, PC1, PC2) %>%
     column_to_rownames(var = "SampleID")
x <- scale(x)

y <- pcaDF2 %>% column_to_rownames(var = "cytokine")
y <- scale(y)

plotDF <- rbind(x, y) %>%
          as.data.frame() %>%
          rownames_to_column() %>%
          mutate(lab = ifelse(rowname %in% rownames(y), rowname, ""))

# segment data frame for arrows
segmentDF <- data.frame(x1 = rep(0, nrow(y)),
                        x2 = y[, "PC1"],
                        y1 = rep(0, nrow(y)),
                        y2 = y[, "PC2"])

# biplot
plotPCA <- ggplot(data = plotDF,
                  mapping = aes(x = PC1,
                                y = PC2)) +
           geom_text(mapping = aes(label = lab),
                     size = 3,
                     hjust = 0,
                     nudge_x = 0.05,
                     color = "black") +
           geom_point(size = 2, mapping = aes(color = c(as.character(pcaDF$TimePoint),
                                                        rep(NA, length(pcaDF2$cytokine))))) +
           scale_color_manual(values = tCols) +
           labs(color = "Duration on ART") +
           labs(x = paste("PC1(",
                           round(summaryPCA["Proportion of Variance", "PC1"] * 100),
                           "%)", sep = ""),
                y = paste("PC2(",
                           round(summaryPCA["Proportion of Variance", "PC2"] * 100),
                           "%)", sep = "")) +
            geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
            geom_vline(xintercept = 0, linetype = "dashed", size = 0.2) +
            geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
                         data = segmentDF,
                         arrow = arrow(length = unit(0.2, "cm")),
                         size = 0.3) +
            theme_bw() +
            theme(axis.text.x = element_text(size=8, color = "black"),
                  axis.text.y = element_text(size=8, color = "black"),
                  axis.line = element_line(color = "grey"),
                  panel.background = element_blank(),
                  panel.grid = element_blank())
pdf(file = file.path(figDir, "Fig_1a.pdf"), width = 6, height = 5, useDingbats = F)
print(plotPCA)
dev.off()


######################################
### FIG S2B
######################################
# detecting clusters using the gapStatistic
# include TGF-b
mat <- datDF1[complete.cases(datDF1), ]
mat2 <- mat %>% dplyr::select(-DonorID, -TimePoint, -DonorID_TimePoint)
mat2 <- as.data.frame(t(mat2))  ## cyt * samples
exprsMat <- t(scale(t(mat2))) %>% as.matrix

# defining the function for Gap statistic
set.seed(seed = 2)
hclustFun <- function(x, k) {
				return(value = list(cluster = cutree(hclust(d = as.dist(1 - cor(t(x),
																        method = "spearman"))),
                                                            k)))}
fit <- clusGap(exprsMat, hclustFun, K.max = 10, B = 100)
nC <- maxSE(fit$Tab[, "gap"], fit$Tab[, "SE.sim"], method = "globalmax")

# plot gap
pdf(file = file.path(figDir, "GapStatistic_PlasmaCytokines.pdf"), width = 5, height = 5)
plot(fit, main = "Number of optimal plasma cytokine clusters")
abline(v = nC, col = "blue")
dev.off()

# identify clusters
clLS <- hclustFun(exprsMat, nC) %>%
		as.data.frame %>%
		rownames_to_column() %>% arrange(cluster)
markers <- clLS$rowname
clusters <- unstack(clLS)

# Test significance of expression of the plasma cytokines (loess + piecewise regression)
datDF <- read_excel("PlasmaBiomarkers.xlsx", sheet =1) %>% as.data.frame()

# Reference lines for uninfected UG and USA donors
uiLS <- lapply(c(1:2), function(i) {
				datDF1 <- read_excel("Uninfected_UG_USA.xlsx", sheet = i) %>% 
                          as.data.frame()
				sName <- excel_sheets("Uninfected_UG_USA.xlsx")[i]
				colnames(datDF1) <- gsub(pattern = " pg\\/ml| ng\\/ml| pg\\/mL|Hu ",
                                         replacement = "",
										 colnames(datDF1))
				datDF2 <- datDF1 %>% select_(.dots = c("SampleID", markers))
				datDF2[, c(2:ncol(datDF2))] <- apply(datDF2[, c(2:ncol(datDF2))], 
													 MARGIN = 2,
													 FUN=as.numeric) %>%
												as.data.frame()
				uimDF <- colMedians(as.matrix(datDF2[, c(2:ncol(datDF2))]), 
                                    na.rm = TRUE) %>%
						 as.data.frame()
				rownames(uimDF) <- markers
				colnames(uimDF) = sName				
				return(value = uimDF)
			 })
uiDF <- do.call(cbind, uiLS) %>%
		rownames_to_column() %>%
		dplyr::rename(markers = rowname) %>%
        mutate(UI_UG = log10(UI_UG), UI_USA = log10(UI_USA))


######################################
### FIG 1B/C/D
######################################

### Plot cytokine expression
datDF <- read_excel("PlasmaBiomarkers.xlsx", sheet =1) %>% as.data.frame()
clusters2 <- clusters[c(2,3,4)]
datDF[, as.character(unlist(clusters))] <- log10(datDF[, as.character(unlist(clusters))])

for(CLUSTER in clusters2) {
    cluster <- CLUSTER
    # print(cluster)
    for(i in 1:length(cluster)) {
        mrk <- cluster[i]
        datDF1 <- datDF %>%
                  filter(!TimePoint %in% "LT") %>%
                  select_(.dots = c("SampleID", "TimePoint", "duration", mrk))
        colnames(datDF1)[colnames(datDF1) %in% mrk] <- "cytokine"
        uiDFm <- uiDF %>% filter(markers %in% mrk)
        xorder <- unique(datDF1$TimePoint)

       # median of cytokine
       mDF <- datDF1 %>%
              group_by(TimePoint) %>%
              summarize(medianHIV_UG = median(cytokine, na.rm = TRUE)) %>%
              as.data.frame() %>%
              dplyr::select(TimePoint, medianHIV_UG)
       rownames(mDF) <- mDF$"TimePoint"
       mDF <- mDF[xorder, ]
      
      # plot
       outFile <- paste(mrk, ".pdf", sep = "")
       pdf(file = file.path(figDir, outFile), height = 4, width = 4, useDingbats = FALSE)
       p <- ggplot(data = datDF1,
                   mapping = aes(x = TimePoint, y = cytokine)) +
	        geom_point(size = 0.8) +
            labs(x = "Duration on ART",
                 y = paste("log10(",mrk,")", " pg/ml", sep="")) +
            geom_hline(yintercept = uiDFm[, "UI_UG"], color = "blue") +
            geom_path(data = mDF, mapping = aes(x = TimePoint, y = medianHIV_UG),
                      color = "red", group = 1) +
            scale_x_discrete(limits = xorder) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.text.x = element_text(size = 8, color = "black"),
                  axis.title.x = element_text(size = 8, color = "black"),
                  axis.title.y = element_text(size = 8, color = "black"),
                  axis.text.y = element_text(size = 8, color = "black"),
                  axis.line = element_line(color = "black"))
             print(p)
        dev.off()
    }
}

#  piecewise regression for clusters 2,3 and 4
datDF <- read_excel("PlasmaBiomarkers.xlsx", sheet =1) %>% as.data.frame()
clusters2 <- clusters[c(2,3)]
datDF_log <- datDF
datDF_log[, as.character(unlist(clusters))] <- log10(datDF_log[, as.character(unlist(clusters))])

for(CLUSTER in clusters2) {
	cluster <- CLUSTER
	# print(cluster)
	   for(i in 1:length(cluster)) {
			 mrk <- cluster[i]
			 print(mrk)
    # loess regression to detect peaks (until M12)
		     datDF1 <- datDF %>%
                       filter(!TimePoint %in% c("LT", "M24")) %>%
                       select_(.dots = c("SampleID", "TimePoint", "duration", mrk))
             datDF1 <- datDF1[complete.cases(datDF1), ]
             markSel <- datDF1[, mrk] %>% as.numeric()
             reg1 <- loess(markSel ~ duration, data = datDF1)
             
             # predict on new data
             dNew <- seq(min(datDF1$duration), max(datDF1$duration), length.out = 500)
             yhat <- predict(reg1, newdata = data.frame(duration = dNew))
             
             # find peaks and valleys
             pks <- findPeaks(yhat)
             pks <- dNew[pks]
             val <- findValleys(yhat)
             val <- dNew[val]
             K <- c(pks, val) %>% sort()

            # add M12 to K
             K <- c(K, 365)
  
            # fit piecewise regression with within determined knots range
            datDF2 <- datDF_log %>%
                      select_(.dots = c("SampleID", "TimePoint", "duration", mrk)) %>%
                      mutate(duration = ifelse(TimePoint == "M24", 730, duration))
            datDF2 <- datDF2[complete.cases(datDF2), ]
            markSel <- datDF2[, mrk] %>% as.numeric()
            
            if(length(K) == 2) {
                lm1a <- lm(markSel ~ ifelse(duration < K[1], duration-K[1], 0) +
                                     ifelse(duration >= K[1], duration-K[1], 0) +
                                     ifelse(duration >= K[2], duration-K[2], 0) ,
                           data = datDF2)
                    } else if(length(K) == 3) {
                lm1a <- lm(markSel ~ ifelse(duration < K[1], duration-K[1], 0) +
                                     ifelse(duration >= K[1], duration-K[1], 0) +
                                     ifelse(duration >= K[2], duration-K[2], 0) +
                                     ifelse(duration >= K[3], duration-K[3], 0),
                           data = datDF2)} else {
                lm1a <- lm(markSel ~ ifelse(duration < K[1], duration-K[1], 0) +
                                     ifelse(duration >= K[1], duration-K[1], 0) +
                                     ifelse(duration >= K[2], duration-K[2], 0) +
                                     ifelse(duration >= K[3], duration-K[3], 0) +
                                     ifelse(duration >= K[4], duration-K[4], 0),
                          data = datDF2)
                        }

         # wilcox tests : D0-M3; M3-M4; M4-M12; M12-M24
         print(paste("D0-M3", wilcox.test(datDF2[,mrk][datDF2$TimePoint %in% "d0"],
                                          datDF2[,mrk][datDF2$TimePoint %in% "M3"])$p.value,
                      sep = " : "))
         print(paste("M3-M4", wilcox.test(datDF2[,mrk][datDF2$TimePoint %in% "M3"],
                                          datDF2[,mrk][datDF2$TimePoint %in% "M4"])$p.value,
                      sep = " : "))
         print(paste("M4-M12", wilcox.test(datDF2[,mrk][datDF2$TimePoint %in% "M4"],
                                           datDF2[,mrk][datDF2$TimePoint %in% "M12"])$p.value,
                      sep = " : "))
         print(paste("M12-M24", wilcox.test(datDF2[,mrk][datDF2$TimePoint %in% "M12"],
                                            datDF2[,mrk][datDF2$TimePoint %in% "M24"])$p.value,
                      sep = " : "))

         # predict on new data and calculate confidence interval
         lineDF <- data.frame(duration = 1:datDF2$duration[length(datDF2$duration)])
         lineDF$yhat <- predict(lm1a,newdata = lineDF)
         yhat1 <- predict(lm1a,
                          newdata = data.frame(duration = 1:datDF2$duration[length(datDF2$duration)]),
                          se.fit = T)
         ciDF <- data.frame(up = yhat1$fit + 2 * yhat1$se.fit,
                            dn = yhat1$fit - 2 * yhat1$se.fit,
                            duration = 1:datDF2$duration[length(datDF2$duration)])

        # fetch UG and USA un-infected plasma cytokine levels
        uiDFm <- uiDF %>% filter(markers %in% mrk)

        # plot
        outFile2 <- paste(mrk, "_loess_piecewise_withUninfected.pdf", sep = "")
        pdf(file = file.path(figDir, outFile2),height = 4, width = 5, useDingbats = F)
        p <- ggplot(data = datDF2,
                    mapping = aes(x = duration, y = datDF2[, mrk])) +
             geom_point(size = 0.8) +
             labs(x = "Duration on ART",
                  y = paste("log10(", mrk, ")", " pg/ml", sep = "")) +
             geom_hline(yintercept = uiDFm[, "UI_UG"], color = "blue", size = 0.5) +
             geom_hline(yintercept = uiDFm[, "UI_USA"], color = "green", size = 0.5) +
             geom_line(data = lineDF,
                      mapping = aes(y = yhat), color = "red", size = 0.5) +
             scale_x_continuous(breaks = unique(datDF2$duration),
                                labels = unique(datDF2$TimePoint)) +
             geom_line(data = ciDF, mapping = aes(y = up), linetype = 2, color = "blue") +
             geom_line(data = ciDF, mapping = aes(y = dn), linetype = 2, color = "blue") +
             theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   axis.text.x = element_text(color = "black", size = 8, angle = 45),
                   axis.text.y = element_text(color = "black", size = 10),
                   axis.line = element_line(color = "black"))
        print(p)
        dev.off()
       }
	}

######################################
### FIG S2C/D
######################################
### For TGF-b and IP10
datDF <- read_excel("PlasmaBiomarkers.xlsx", sheet =1) %>% as.data.frame()
mrk = "IP10"
        datDF1 <- datDF %>%
                  filter(!TimePoint %in% "LT") %>%
                  select_(.dots = c("SampleID", "TimePoint", "duration", mrk))
        colnames(datDF1)[colnames(datDF1) %in% mrk] <- "cytokine"
        datDF1$cytokine <- log10(datDF1$cytokine)
        uiDFm <- uiDF %>% filter(markers %in% mrk)
        xorder <- unique(datDF1$TimePoint)

       # median of cytokine
       mDF <- datDF1 %>%
              group_by(TimePoint) %>%
              summarize(medianHIV_UG = median(cytokine, na.rm = TRUE)) %>%
              as.data.frame() %>%
              dplyr::select(TimePoint, medianHIV_UG)
       rownames(mDF) <- mDF$"TimePoint"
       mDF <- mDF[xorder, ]
      
       outFile <- paste(mrk, ".pdf", sep = "")
       pdf(file = file.path(figDir, outFile), height = 4, width = 4, useDingbats = FALSE)
       p <- ggplot(data = datDF1,
                   mapping = aes(x = TimePoint, y = cytokine)) +
            geom_point(size = 0.8) +
            labs(x = "Duration on ART",
                 y = paste("log10(",mrk,")", " pg/ml", sep="")) +
            geom_hline(yintercept = uiDFm[, "UI_UG"], color = "blue") +
            geom_path(data = mDF, mapping = aes(x = TimePoint, y = medianHIV_UG),
                      color = "red", group = 1) +
            scale_x_discrete(limits = xorder) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.text.x = element_text(size = 8, color = "black"),
                  axis.title.x = element_text(size = 8, color = "black"),
                  axis.title.y = element_text(size = 8, color = "black"),
                  axis.text.y = element_text(size = 8, color = "black"),
                  axis.line = element_line(color = "black"))
        print(p)
        dev.off()



######################################
### FIG S2E
######################################
#### Correlation between plasma cytokine expression clusters
datDF <- read_excel("PlasmaBiomarkers.xlsx", sheet = 1) %>% as.data.frame()
markers <- clusters[c(2,4,3)] %>% unlist(.) %>% unname()
datDF1 <- datDF %>%
          filter(TimePoint %in% c("d0", "M3", "M4", "M12", "M24")) %>%
          select(markers)

# correlation data frame
corDF <- cor(datDF1, use = "complete.obs", method = "spearman")
corDF[lower.tri(corDF)] <- NA

# plot correlation heatmap
outFile <- "Correlation between cytokine clusters Heatmap.pdf"
pdf(file = file.path(figDir,outFile), width  = 6, height = 6)
colorPalette <- c("blue", "white", "red")
colorPalette <- colorRampPalette(colors = colorPalette)(100)
breaks <- c(-1,
 	  		seq(from = min(corDF, na.rm = T),
 	  			to = abs(min(corDF, na.rm = T)),
 	  			length.out = 99),
 	  		1)
pheatmap(mat = corDF,
         color = colorPalette,
         breaks = breaks,
         cellwidth = 15,
         cellheight = 15,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         gaps_col = c(4,6),
         gaps_row = c(4,6),
         border_color = "grey",
         fontsize = 12)
dev.off()


######################################
### FIG 1E
######################################
#### Summarizing fold change across cytokines per cluster
# 1: log2FC for each cytokine and donor - Cytokine1:Donor1 - D0/D0, M3/D0, M4/M3, M12/M4, M24/M12
# 2: mean fold change across cytokines per donor and timepoint
datDF <- read_excel("PlasmaBiomarkers.xlsx", sheet = 1) %>% as.data.frame()
markers <- clusters[c(2,4,3)] %>% unlist(.) %>% unname()
plotDF <- datDF %>%
          dplyr::filter(TimePoint %in% c("d0", "M3", "M4", "M9", "M12", "M24")) %>%
          dplyr::select_(.dots = c("SampleID", "TimePoint", unname(markers))) %>%
          gather(Cytokine, value, -SampleID, -TimePoint) %>%
          mutate(Cluster = ifelse(Cytokine %in% clusters[[2]], "C1",
                           ifelse(Cytokine %in% clusters[[3]], "C3", "C2"))) %>%
          spread(TimePoint, value) %>%
          mutate(D0_FC = log2(d0/d0),
                 M3_FC = log2(M3/d0),
                 M4_FC = log2(M4/M3),
                 M9_FC = log2(M9/M4),
                 M12_FC = log2(M12/M9),
                 M24_FC = log2(M24/M12)) %>%
          dplyr::select(-d0, -M3, -M4, -M9, -M12, -M24) %>%
          gather(TimePoint, value, -SampleID, -Cytokine, -Cluster) %>%
          group_by(SampleID, Cluster, TimePoint) %>%
          mutate(meanFC = mean(value, na.rm = TRUE)) %>%
          dplyr::select(-value, -Cytokine) %>%
          unique() %>%
          as.data.frame() %>%
          mutate(TimePoint = gsub("_FC", "", TimePoint),
                 SampleID = c(RAL03 = "2213",
                              RAL04 = "2214",
                              RAL07 = "2215",
                              RAL08 = "2216",
                              RAL10 = "2217",
                              RAL11 = "2218",
                              RAL12 = "2219",
                              RAL13 = "2220",
                              RAL14 = "2221",
                              RAL15 = "2222",
                              RAL16 = "2223")[SampleID])
clusterDF <- plotDF

# median DF
xorder <- c("D0", "M3", "M4", "M9", "M12", "M24")
mDF <- plotDF %>%
       group_by(TimePoint, Cluster) %>%
       mutate(medFC = median(meanFC)) %>%
       dplyr::select(-meanFC, -SampleID) %>%
       unique(.) %>%
       as.data.frame()
mDF <- mDF[order(match(mDF$TimePoint, table = xorder)), ]
mDF1 <- mDF %>% filter(Cluster %in% "C1")
mDF2 <- mDF %>% filter(Cluster %in% "C2")
mDF3 <- mDF %>% filter(Cluster %in% "C3")

# plot summarized expression (point plot)
outFile <- "Summarized Cytokine Expression_FoldChange.pdf"
pdf(file = file.path(figDir, outFile), height = 3, width = 8, useDingbats = FALSE)
sumP <- ggplot(data = plotDF,
               mapping = aes(x = TimePoint, y = meanFC)) +
        geom_point(aes(color = SampleID), size = 2) +
        scale_color_manual(values = dCols, name = "DonorID") +
        labs(x = "Duration on ART", y = "Mean - log2(Fold Change)") +
        geom_hline(yintercept = 0) +
        facet_wrap(~Cluster, scales = "free") +
        geom_path(aes(x = TimePoint, y = medFC), data = mDF1, linetype = 2, group = 1) +
        geom_path(aes(x = TimePoint, y = medFC), data = mDF2, linetype = 2, group = 1) +
        geom_path(aes(x = TimePoint, y = medFC), data = mDF3, linetype = 2, group = 1) +
        scale_x_discrete(limits = xorder) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text.x = element_text(size = 9, color = "black", angle = 45),
              axis.title.x = element_text(size = 10, color = "black"),
              axis.text.y = element_text(size = 10, color = "black"),
              axis.title.y = element_text(size = 10, color = "black"),
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10))
print(sumP)
dev.off()

###### calculate correlation between clusters
# between C1 and C2
cor.test(plotDF$meanFC[plotDF$Cluster %in% "C1"],
         plotDF$meanFC[plotDF$Cluster %in% "C2"], method = "spearman")



######################################
### FIG 1F
######################################
###### Correlation between Plasma cytokines and CD4 counts
##### correlation between summarized logFC plasma cytokines (per cluster) and delta cd4/cd8 ratio
# taking timepoints D0, M4,
datDF <- read_excel("PlasmaBiomarkers.xlsx", sheet = 1) %>% as.data.frame()
markers <- clusters[c(2,4,3)] %>% unlist(.) %>% unname()
plotDF <- datDF %>%
          dplyr::filter(TimePoint %in% c("d0", "M12", "M24")) %>%
          dplyr::select_(.dots = c("SampleID", "TimePoint", unname(markers))) %>%
          gather(Cytokine, value, -SampleID, -TimePoint) %>%
          mutate(Cluster = ifelse(Cytokine %in% clusters[[2]], "C1",
                           ifelse(Cytokine %in% clusters[[3]], "C3", "C2"))) %>%
          spread(TimePoint, value) %>%
          mutate(D0_FC = log2(d0/d0),
                 M12_FC = log2(M12/d0),
                 M24_FC = log2(M24/M12)) %>%
          dplyr::select(-d0, -M12, -M24) %>%
          gather(TimePoint, value, -SampleID, -Cytokine, -Cluster) %>%
          group_by(SampleID, Cluster, TimePoint) %>%
          mutate(meanFC = mean(value, na.rm = TRUE)) %>%
          dplyr::select(-value, -Cytokine) %>%
          unique() %>%
          as.data.frame() %>%
          mutate(TimePoint = gsub("_FC", "", TimePoint),
                 SampleID = c(RAL03 = "2213",
                              RAL04 = "2214",
                              RAL07 = "2215",
                              RAL08 = "2216",
                              RAL10 = "2217",
                              RAL11 = "2218",
                              RAL12 = "2219",
                              RAL13 = "2220",
                              RAL14 = "2221",
                              RAL15 = "2222",
                              RAL16 = "2223")[SampleID],
                 SampleID_TimePoint = interaction(SampleID, TimePoint, sep = "_"))
# median DF
xorder <- c("D0", "M12", "M24")
mDF <- plotDF %>%
       group_by(TimePoint, Cluster) %>%
       mutate(medFC = median(meanFC)) %>%
       dplyr::select(-meanFC, -SampleID) %>%
       unique(.) %>%
       as.data.frame()

mDF <- mDF[order(match(mDF$TimePoint, table = xorder)), ]
mDF1 <- mDF %>% filter(Cluster %in% "C1")
mDF2 <- mDF %>% filter(Cluster %in% "C2")
mDF3 <- mDF %>% filter(Cluster %in% "C3")

# plot summarized expression (point plot)
outFile <- "Summarized Cytokine Expression_FoldChange_forCorwithCD4CD8.pdf"
pdf(file = file.path(figDir, outFile), height = 4, width = 10)
sumP <- ggplot(data = plotDF,
               mapping = aes(x = TimePoint, y = meanFC)) +
        geom_point(aes(color = SampleID), size = 4) +
        scale_color_manual(values = dCols, name = "DonorID") +
        labs(x = "Duration on ART", y = "log2 change between timepoints") +
        geom_hline(yintercept = 0) +
        facet_wrap(~Cluster, scales = "free") +
        geom_path(aes(x = TimePoint, y = medFC), data = mDF1, linetype = 2, group = 1) +
        geom_path(aes(x = TimePoint, y = medFC), data = mDF2, linetype = 2, group = 1) +
        geom_path(aes(x = TimePoint, y = medFC), data = mDF3, linetype = 2, group = 1) +
        scale_x_discrete(limits = xorder) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text.x = element_text(size = 9, color = "black", angle = 45),
              axis.title.x = element_text(size = 10, color = "black"),
              axis.text.y = element_text(size = 10, color = "black"),
              axis.title.y = element_text(size = 10, color = "black"),
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10))
print(sumP)
dev.off()


# CD4 data frame : (D0, M12, M24)
plotDF2 <- read_excel("RALTcellcounts.xlsx", sheet = 1) %>%
         as.data.frame() %>%
         dplyr::rename(CD4CD8Ratio = `CD4:CD8 Ratio`) %>%
         dplyr::select(DonorID, TimePoint, CD4CD8Ratio) %>%
         filter(TimePoint %in% c("d0", "M12", "M24")) %>%
         mutate(CD4CD8Ratio = as.numeric(CD4CD8Ratio)) %>%
         spread(TimePoint, CD4CD8Ratio) %>%
         mutate(D0_FC = d0-d0,
                M12_FC = M12-d0,
                M24_FC = M24-M12) %>%
         dplyr::select(-d0, -M12, -M24) %>%
         gather(TimePoint, CD4CD8, -DonorID) %>%
         mutate(TimePoint = gsub("_FC", "", TimePoint),
                SampleID = c(RAL03 = "2213",
                             RAL04 = "2214",
                             RAL07 = "2215",
                             RAL08 = "2216",
                             RAL10 = "2217",
                             RAL11 = "2218",
                             RAL12 = "2219",
                             RAL13 = "2220",
                             RAL14 = "2221",
                             RAL15 = "2222",
                             RAL16 = "2223")[DonorID],
               SampleID_TimePoint = interaction(SampleID, TimePoint, sep = "_")) %>%
               dplyr::select(-DonorID)
# median DF
xorder <- c("D0", "M12", "M24")
mDF <- plotDF2 %>%
       group_by(TimePoint) %>%
       summarize(medFC = median(CD4CD8, na.rm = T)) %>%
       as.data.frame()

# plot summarized expression (point plot)
outFile <- "CD4CD8 Delta.pdf"
pdf(file = file.path(figDir, outFile), height = 4, width = 4)
sumP <- ggplot(data = plotDF2,
               mapping = aes(x = TimePoint, y = CD4CD8)) +
        geom_point(aes(color = SampleID), size = 4) +
        scale_color_manual(values = dCols, name = "DonorID") +
        labs(x = "Duration on ART", y = "change between timepoints") +
        geom_path(aes(x = TimePoint, y = medFC), data = mDF, linetype = 2, group = 1) +
        scale_x_discrete(limits = xorder) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text.x = element_text(size = 9, color = "black", angle = 45),
              axis.title.x = element_text(size = 10, color = "black"),
              axis.text.y = element_text(size = 10, color = "black"),
              axis.title.y = element_text(size = 10, color = "black"),
              legend.title = element_text(size = 8),
              legend.text = element_text(size = 10))
print(sumP)
dev.off()
