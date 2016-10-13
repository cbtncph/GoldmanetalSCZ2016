
# NormalizationAndDEAnalysis.R
# For publication by Goldman et al. 2016 on schizophrenia-derived glia
# Author and contact: Mikhail Osipovitch (mxo@sund.ku.dk)
#
# The following script performs RUVSeq-based normalization and differential expression analysis
# of RNA-Seq data produced from schizophrenia- and control-derived human glial cells. The procedure
# is accomplished in three steps:
#
# 1 - first-pass differential expression analysis for determination of in silico negative control
#     genes (FDR-corrected P Value > 0.75) that are not affected by the condition of interest
# 2 - calculation of variance normalization factors by RUVs function of RUVSeq
# 3 - second-pass differential expression analysis (5% FDR and log2 fold change > 1) for
#     determination of disease-dysregulated genes using the original counts but with adjusting for
#     RUVs-calculated variance factors by multi-factor GLM models implemented in edgeR and DESeq2
#
# NOTE: Steps 2 and 3 are accomplished in a loop that runs as many iterations as there are samples
# in the dataset. At each iteration, the loop increases the number of additional covariates by 1
# (the variable k in RUVs function), produces diagnostic plots and performs the second-pass
# differential expression analysis with adjustment for these covariates. Output is produced
# corresponding to each iteration in the form of PCA plots, RLE plots, P Value histogram plots,
# and differential expression tables containing calculated fold changes and P Values.

### load required libraries:
library(gplots)
library(ape)
library(scatterplot3d)
library(RColorBrewer)
library(DESeq2)
library(edgeR)
library(RUVSeq)
library(EDASeq)


##########################################################################################################################
### DATASET SPECIFICATION ################################################################################################
##########################################################################################################################

experiment_title = "SCZ RNA-Seq - Goldman et al. 2016"

### specify colors for cell lines to be used in plots and figures:
col_22  <- rep("#636363", 3)
col_37  <- rep("#636363", 4)
col_205 <- rep("#636363", 7)
col_08  <- rep("#528881", 4)
col_29  <- rep("#5e3c58", 3)
col_51  <- rep("#886451", 7)
col_164 <- rep("#2a334f", 8)

### read count data and sample sheet:
countData_allSamples  <- read.table("countData_allSamples.txt", sep = "\t")
sampleData_allSamples <- read.table("sampleData_allSamples.txt", sep = "\t")

### The lines below specify differential comparison to be made. The variable
### "which_comparison" serves as a switch between 4 differential comparisons.
### The following values are to be assigned to "which_comparison":
### 1 : 4 pooled SCZ lines vs. 3 pooled CTR lines
### 2 : SCZ line 8 vs. 3 pooled CTR lines
### 3 : SCZ line 29 vs. 3 pooled CTR lines
### 4 : SCZ line 51 vs. 3 pooled CTR lines
### 5 : SCS line 164 vs. 3 pooled CTR lines

which_comparison <- 5 ### accepts values between 1 and 5

if (which_comparison == 1) { ### 1 : 4 pooled SCZ lines vs. 3 pooled CTR lines
  
  ### suffix is used in filenames to distinguish output for different comparisons:
  suffix <- " - SCZ vs CTR"
  
  ### subset countData:
  countData <- countData_allSamples
  
  ### subset sampleData:
  sampleData <- sampleData_allSamples
  
  ### combined vector of colors to be used in plots and figures:
  cell_line_colors <- c(col_22, col_37, col_205, col_08, col_29, col_51, col_164)
  
  ### factors for edgeR differential expression analysis (1 = CTR, 2 = SCZ):
  x <- factor(c(rep(1, 14), rep(2, 22)))
  
  ### groups for RUVg function of RUVSeq package for normalization (first vector is CTR, second in SCZ):
  groups <- rbind(c(1:14), c(15:36))
  
} else if (which_comparison == 2) { ### 2 : SCZ line 8 vs. 3 pooled CTR lines
  
  suffix <- " - SCZ.8 vs CTR"
  countData <- countData_allSamples[,c(c(1:14), c(15:18))]
  sampleData <- sampleData_allSamples[c(c(1:14), c(15:18)),]
  cell_line_colors <- c(col_22, col_37, col_205, col_08)
  x <- factor(c(rep(1, 14), rep(2, 4)))
  groups <- rbind(c(1:14), c(15:18, rep(-1, 10)))
  
} else if (which_comparison == 3) { ### 3 : SCZ line 29 vs. 3 pooled CTR lines
  
  suffix <- " - SCZ.29 vs CTR"
  countData <- countData_allSamples[,c(c(1:14), c(19:21))]
  sampleData <- sampleData_allSamples[c(c(1:14), c(19:21)),]
  cell_line_colors <- c(col_22, col_37, col_205, col_29)
  x <- factor(c(rep(1, 14), rep(2, 3)))
  groups <- rbind(c(1:14), c(15:17, rep(-1, 11)))
  
} else if (which_comparison == 4) { ### 4 : SCZ line 51 vs. 3 pooled CTR lines
  
  suffix <- " - SCZ.51 vs CTR"
  countData <- countData_allSamples[,c(c(1:14), c(22:28))]
  sampleData <- sampleData_allSamples[c(c(1:14), c(22:28)),]
  cell_line_colors <- c(col_22, col_37, col_205, col_51)
  x <- factor(c(rep(1, 14), rep(2, 7)))
  groups <- rbind(c(1:14), c(15:21, rep(-1, 7)))
  
} else if (which_comparison == 5) { ### 5 : SCZ line 164 vs.3 pooled CTR lines
  
  suffix <- " - SCZ.164 vs CTR"
  countData <- countData_allSamples[,c(c(1:14), c(29:36))]
  sampleData <- sampleData_allSamples[c(c(1:14), c(29:36)),]
  cell_line_colors <- c(col_22, col_37, col_205, col_164)
  x <- factor(c(rep(1, 14), rep(2, 8)))
  groups <- rbind(c(1:14), c(15:22, rep(-1, 6))) 

} else {
  
  print("Please provide a legal value for which_comparison to make :/")

} # end else


##########################################################################################################################
### START of ANALYSIS: FILTER DATA and MAKE PRELIMINARY PLOTS  ###########################################################
##########################################################################################################################

### make necessary folders to store output:

output_folder <- paste("Output", suffix, "/", sep = "")
dir.create(output_folder)

output_rle <- paste(output_folder, "Output - RLE Plots/", sep = "")
output_pca <- paste(output_folder, "Output - PCA Plots/", sep = "")
output_matrices <- paste(output_folder, "Output - Normalized Matrices/", sep = "")
dir.create(output_rle)
dir.create(output_pca)
dir.create(output_matrices)

### filter out lowly expressed transcripts:
countData <- countData[apply(countData, 1, function(x) length(x[x >= 5]) >= 5),]

### write out filtered count matrix (NN = No Normalization):
write.table(countData, paste(output_matrices,"countData - NN", suffix, ".txt", sep = ""), sep = "\t", quote = FALSE)

### write out filtered, shifted, log2-transformed count matrix (NN = No Normalization):
write.table(log2(countData + 1), paste(output_matrices,"countData - NN_Log2", suffix, ".txt", sep = ""), sep = "\t", quote = FALSE)

### UQ Normalization to account for library sizes:
uq <- as.data.frame(betweenLaneNormalization(as.matrix(countData), which = "upper"))

### write out UQ count matrix (UQ = Upper Quartile):
write.table(uq, paste(output_matrices, "countData - UQ", suffix, ".txt", sep = ""), sep = "\t", quote = FALSE)

### write out UQ, shifted, log2-transformed count matrix, *recommended for visualization* (UQ = Upper Quartile):
write.table(log2(uq + 1), paste(output_matrices, "countData - UQ_log2", suffix, ".txt", sep = ""), sep = "\t", quote = FALSE)

### RLE plot (NN = No Normalization):
png(paste(output_rle, "RLE - NN", suffix, ".png", sep = ""), width = 660, height = 440)
par(cex = 1.2)
plotRLE(as.matrix(countData), outline = FALSE, ylim = c(-2, 2.5), col = cell_line_colors, las = 2, ylab = "Relative Log Expression")
graphics.off()

### 2D PCA Plot (NN = No Normalization):
png(paste(output_pca, "PCA 2D - NN", suffix, ".png", sep = ""), width = 660, height = 440)
par(cex = 1.3)
plotPCA(as.matrix(countData), col = cell_line_colors, xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
graphics.off()

### RLE plot (UQ = Upper Quartile):
png(paste(output_rle, "RLE - UQ", suffix, ".png", sep = ""), width = 660, height = 440)
par(cex = 1.2)
plotRLE(as.matrix(uq), outline = FALSE, ylim = c(-2, 2.5), col = cell_line_colors, las = 2, ylab = "Relative Log Expression")
graphics.off() 

### 2D PCA Plot (UQ = Upper Quartile):
png(paste(output_pca, "PCA 2D - UQ", suffix, ".png", sep = ""), width = 660, height = 440)
par(cex = 1.3)
plotPCA(as.matrix(uq), col = cell_line_colors, xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
graphics.off()


##########################################################################################################################
### STEP 1: FIRST PASS DIFFERENTIAL EXPRESSION ANALYSIS ##################################################################
##########################################################################################################################

### make necessary folders to store output:
dir.create(paste(output_folder, "Output - DE First Pass/", sep  = ""))

output_path_int <- paste(output_folder, "Output - DE First Pass/DEGs_INTs/", sep = "")
output_path_deseq2 <- paste(output_folder, "Output - DE First Pass/DEGs_DESeq2/", sep = "")
output_path_edger <- paste(output_folder, "Output - DE First Pass/DEGs_edgeR/", sep = "")
output_path_allde <- paste(output_folder, "Output - DE First Pass/All DE Results Tables/", sep = "")
output_path_pvalhist <- paste(output_folder, "Output - P Value Histograms/", sep = "")

dir.create(output_path_int)
dir.create(output_path_deseq2)
dir.create(output_path_edger)
dir.create(output_path_allde)
dir.create(output_path_pvalhist)

##### DESeq2 Differential Expresstion:

### Create a DESeq2 dataset - no normalization
labels <- unique(sampleData$condition)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = ~ condition)
dds$condition <- factor(dds$condition, levels = labels)

### DESEq2 single factor: ~ condition
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
res_matrix <- as.matrix(res)

### write out complete differential expression results from DESeq2:
write.table(res_matrix, paste(output_path_allde, "All DE Results - DESeq2 - First Pass", suffix, ".txt"), sep = "\t", quote = FALSE)

degs_deseq2_001 <- as.data.frame(res[which(res$padj < 0.01), ])
degs_deseq2_005 <- as.data.frame(res[which(res$padj < 0.05), ])
degs_deseq2_010 <- as.data.frame(res[which(res$padj < 0.10), ])

### P Value histogram for DESeq2 results:
h1 <- hist(res_matrix[,6], breaks=0:50/50, plot=FALSE)
png(paste(output_path_pvalhist ,"PVal Hist DESeq2 - NN", suffix, ".png", sep = ""),
    width = 660, height = 440)
par(cex = 1.3)
barplot(height = h1$counts, beside = FALSE, ylim = c(0, 5000),
        col = "#d2d4dc", space = 0, ylab="Frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
graphics.off()

##### edgeR Differential Expression:

### edgeR single factor: ~ condition (~ x)
design <- model.matrix(~ x, data = sampleData)
y <- DGEList(counts = uq, group = x)
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
all_tags <- topTags(lrt, n = nrow(uq))
all_tags <- all_tags$table
all_tags <- all_tags[order(all_tags$FDR),]

### write out complete differential expression results from edgeR:
write.table(all_tags, paste(output_path_allde, "All DE Results - edgeR - First Pass", suffix, ".txt"), sep = "\t", quote = FALSE)

degs_edger_001 <- all_tags[which(all_tags$FDR < 0.01), ]
degs_edger_005 <- all_tags[which(all_tags$FDR < 0.05), ]
degs_edger_010 <- all_tags[which(all_tags$FDR < 0.10), ]

### P Value histogram for edgeR results:
h1 <- hist(all_tags$FDR, breaks=0:50/50, plot=FALSE)
png(paste(output_path_pvalhist, "PVal Hist edgeR - NN", suffix, ".png", sep = ""),
    width = 660, height = 440)
par(cex = 1.3)
barplot(height = h1$counts, beside = FALSE, ylim = c(0, 5000),
        col = "#d2d4dc", space = 0, ylab="Frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
graphics.off()

### write out individual and intersection tables from DESeq2 and edgeR:  
write.table(as.data.frame(degs_deseq2_001), paste(output_path_deseq2, "degs_deseq2_001", suffix,".txt", sep = ""), sep = "\t", quote = FALSE)
write.table(as.data.frame(degs_deseq2_005), paste(output_path_deseq2, "degs_deseq2_005", suffix,".txt", sep = ""), sep = "\t", quote = FALSE)
write.table(as.data.frame(degs_deseq2_010), paste(output_path_deseq2, "degs_deseq2_010", suffix,".txt", sep = ""), sep = "\t", quote = FALSE)

write.table(degs_edger_001, paste(output_path_edger, "degs_edgeR_001", suffix,".txt", sep = ""), sep = "\t", quote = FALSE)
write.table(degs_edger_005, paste(output_path_edger, "degs_edgeR_005", suffix,".txt", sep = ""), sep = "\t", quote = FALSE)
write.table(degs_edger_010, paste(output_path_edger, "degs_edgeR_010", suffix,".txt", sep = ""), sep = "\t", quote = FALSE)

int_001 <- intersect(rownames(degs_deseq2_001), rownames(degs_edger_001))
int_005 <- intersect(rownames(degs_deseq2_005), rownames(degs_edger_005))
int_010 <- intersect(rownames(degs_deseq2_010), rownames(degs_edger_010))

write.table(cbind(degs_deseq2_001[int_001,], degs_edger_001[int_001,]), paste(output_path_int, "DEGs_001", suffix, ".txt", sep = ""),
            quote = FALSE, sep = "\t")
write.table(cbind(degs_deseq2_005[int_005,], degs_edger_005[int_005,]), paste(output_path_int, "DEGs_005", suffix, ".txt", sep = ""),
            quote = FALSE, sep = "\t")
write.table(cbind(degs_deseq2_010[int_010,], degs_edger_010[int_010,]), paste(output_path_int, "DEGs_010", suffix, ".txt", sep = ""),
            quote = FALSE, sep = "\t")

### initialize matrix to hold summary of DE results:

deResults <- matrix(nrow = ncol(uq) + 1, ncol = 10)
colnames(deResults) <- c("K","DESeq2_01","DESeq2_05","DESeq2_10","edgeR_01","edgeR_05","edgeR_10","INT_01","INT_05","INT_10")

### store the results of the initial first-pass comparison:
deResults[1, 1] <- 0
deResults[1, 2] <- nrow(degs_deseq2_001)
deResults[1, 3] <- nrow(degs_deseq2_005)
deResults[1, 4] <- nrow(degs_deseq2_010)
deResults[1, 5] <- nrow(degs_edger_001)
deResults[1, 6] <- nrow(degs_edger_005)
deResults[1, 7] <- nrow(degs_edger_010)
deResults[1, 8] <- length(int_001)
deResults[1, 9] <- length(int_005)
deResults[1,10] <- length(int_010)

##### select Negative Controls:

PVal_T = 0.75

cutoff_deseq2 <- length(which(res_matrix[,6] < PVal_T)) 
cutoff_edger <- length(which(all_tags$FDR < PVal_T))

### sorted by increasing p value, take the bottom rows:
emp_deseq2 <- row.names(res_matrix[c(cutoff_deseq2:nrow(res_matrix)),])
emp_edger <- row.names(all_tags[c(cutoff_edger:nrow(all_tags)),])

### final list of control genes, intersection of DESeq2 and edgeR:
empirical <- intersect(emp_deseq2, emp_edger)

### write out the empirical control genes with DESeq2 (res_matrix) and edgeR (all_tags) DE data:
write.table(cbind(res_matrix[empirical,], all_tags[empirical,]),
            paste(output_folder, "Output - DE First Pass/Negative Empiricals", suffix, ".txt", sep = ""),
            quote = FALSE, sep = "\t")


##########################################################################################################################
### STEP 2: RUVs NORMALIZATION BY SAMPLES AND NEGATIVE CONTROLS ##########################################################
##########################################################################################################################

### make necessary folders to store output:
dir.create(paste(output_folder, "Output - DE Second Pass/", sep = ""))

output_path_int <- paste(output_folder, "Output - DE Second Pass/DEGs_INTs/", sep = "")
output_path_deseq2 <- paste(output_folder, "Output - DE Second Pass/DEGs_DESeq2/", sep = "")
output_path_edger <- paste(output_folder, "Output - DE Second Pass/DEGs_edgeR/", sep = "")
output_path_allde <- paste(output_folder, "Output - DE Second Pass/All DE Results Tables/", sep = "")

dir.create(output_path_int)
dir.create(output_path_deseq2)
dir.create(output_path_edger)
dir.create(output_path_allde)

for (i in 1:ncol(uq)) {
  
  ### call to RUVs function for normalization:
  s <- RUVs(as.matrix(uq), empirical, k = i, groups)
  
  ### 2D PCA plot:
  png(paste(output_pca,"PCA 2D - UQ_RUVs_", i, suffix, ".png", sep = ""), width = 660, height = 440)
  par(cex = 1.4)
  plotPCA(s$normalizedCounts, col = cell_line_colors, xlim = c(-0.5, 0.7), ylim = c(-0.5, 0.5))
  graphics.off()
  
  ### RLE Plot:
  png(paste(output_rle, "RLE - UQ_RUVs_", i, suffix, ".png", sep = ""), width = 660, height = 440)
  par(cex = 1.2)
  plotRLE(s$normalizedCounts, outline = FALSE, ylim = c(-2, 2.5), col = cell_line_colors, las = 2, ylab = "Relative Log Expression")
  graphics.off()
  
  ### write out RUVs-normalized matrix:
  write.table(s$normalizedCounts, paste(output_matrices, "countData - UQ_RUVs_", i, suffix, ".txt", sep = ""), sep = "\t", quote = FALSE)
  
  ##########################################################################################################################
  ### STEP 3: SECOND PASS DIFFERENTIAL EXPRESSION ##########################################################################
  ##########################################################################################################################
  
  if (TRUE) { # change to FALSE to only generate diagnostic plots
    
    ### make a copy of sampleData with appended normalization factores from RUVs:
    sampleDataW <- cbind(sampleData, s$W)
    
    ##### DESeq2 Differential Expresstion:
    
    ### create a DESeq2 dataset:
    ddsW <- DESeqDataSetFromMatrix(countData = countData, colData = sampleDataW, design = ~ condition)
    ddsW$condition <- factor(ddsW$condition, levels = labels)
    
    ### dinamically compile formula for DESeq2 and edgeR GLMs:
    form_string <- ""
    for (w in 1:i) {
      if (w == i){
        form_string <- paste(form_string, "W_", w, sep = "")
      } else {
        form_string <- paste(form_string, "W_", w, " + ", sep = "")
      } # end else
    } # end for
    
    design(ddsW) <- formula(paste("~", form_string, " + condition", sep = ""))
    
    ### DESEq2 single factor: ~ W_k + condition
    ddsW <- DESeq(ddsW)
    res <- results(ddsW)
    res <- res[order(res$padj),]
    res_matrix <- as.matrix(res)
    
    ### write out complete differential expression results from DESeq2:
    write.table(res_matrix, paste(output_path_allde, "All DE Results - DESeq2 - UQ_RUVs_", i, ".txt", sep = ""), sep = "\t", quote = FALSE)
    
    degs_deseq2_001 <- as.data.frame(res[which(res$padj < 0.01), ])
    degs_deseq2_005 <- as.data.frame(res[which(res$padj < 0.05), ])
    degs_deseq2_010 <- as.data.frame(res[which(res$padj < 0.10), ])
    
    ### P Value histogram for DESeq2 results:
    h1 <- hist(res_matrix[,6], breaks=0:50/50, plot=FALSE)
    png(paste(output_path_pvalhist, "PVal Hist DESeq2 - UQ_RUVs_", i, ".png", sep = ""),
        width = 660, height = 440)
    par(cex = 1.3)
    barplot(height = h1$counts, beside = FALSE, ylim = c(0, 5000),
            col = "#d2d4dc", space = 0, ylab="Frequency")
    text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
         adj = c(0.5,1.7), xpd=NA)
    graphics.off()
    
    ##### edgeR Differential Expression:
    
    ### edgeR single factor: ~ condition (~ x) + W_k
    design <- model.matrix(formula(paste("~ x +", form_string, sep = "")), data = sampleDataW)
    y <- DGEList(counts = uq, group = x)
    y <- calcNormFactors(y)
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef = 2)
    all_tags <- topTags(lrt, n = nrow(uq))
    all_tags <- all_tags$table
    all_tags <- all_tags[order(all_tags$FDR),]
    
    ### write out complete differential expression results from edgeR:
    write.table(all_tags, paste(output_path_allde, "All DE Results - edgeR - UQ_RUVs_", i, ".txt", sep = ""), sep = "\t", quote = FALSE)
    
    degs_edger_001 <- all_tags[which(all_tags$FDR < 0.01), ]
    degs_edger_005 <- all_tags[which(all_tags$FDR < 0.05), ]
    degs_edger_010 <- all_tags[which(all_tags$FDR < 0.10), ]
    
    ### p value histogram for edgeR results:
    h1 <- hist(all_tags$FDR, breaks=0:50/50, plot=FALSE)
    png(paste(output_path_pvalhist, "PVal Hist edgeR - UQ_RUVs_", i, ".png", sep = ""),
        width = 660, height = 440)
    par(cex = 1.3)
    barplot(height = h1$counts, beside = FALSE, ylim = c(0, 5000),
            col = "#d2d4dc", space = 0, ylab="Frequency")
    text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
         adj = c(0.5,1.7), xpd=NA)
    graphics.off()
    
    ### write out individual and intersection tables from DESeq2 and edgeR:
    write.table(as.data.frame(degs_deseq2_001), paste(output_path_deseq2, "degs_deseq2_001 - W_", i, ".txt", sep = ""), sep = "\t", quote = FALSE)
    write.table(as.data.frame(degs_deseq2_005), paste(output_path_deseq2, "degs_deseq2_005 - W_", i, ".txt", sep = ""), sep = "\t", quote = FALSE)
    write.table(as.data.frame(degs_deseq2_010), paste(output_path_deseq2, "degs_deseq2_010 - W_", i, ".txt", sep = ""), sep = "\t", quote = FALSE)
    
    write.table(degs_edger_001, paste(output_path_edger, "degs_edgeR_001 - W_", i, ".txt", sep = ""), sep = "\t", quote = FALSE)
    write.table(degs_edger_005, paste(output_path_edger, "degs_edgeR_005 - W_", i, ".txt", sep = ""), sep = "\t", quote = FALSE)
    write.table(degs_edger_010, paste(output_path_edger, "degs_edgeR_010 - W_", i, ".txt", sep = ""), sep = "\t", quote = FALSE)
    
    ### obtain intersects:
    int_001 <- intersect(rownames(degs_deseq2_001), rownames(degs_edger_001))
    int_005 <- intersect(rownames(degs_deseq2_005), rownames(degs_edger_005))
    int_010 <- intersect(rownames(degs_deseq2_010), rownames(degs_edger_010))
    
    write.table(cbind(degs_deseq2_001[int_001,], degs_edger_001[int_001,]), paste(output_path_int, "DEGs_001 - W_", i, ".txt", sep = ""),
                quote = FALSE, sep = "\t")
    write.table(cbind(degs_deseq2_005[int_005,], degs_edger_005[int_005,]), paste(output_path_int, "DEGs_005 - W_", i, ".txt", sep = ""),
                quote = FALSE, sep = "\t")
    write.table(cbind(degs_deseq2_010[int_010,], degs_edger_010[int_010,]), paste(output_path_int, "DEGs_010 - W_", i, ".txt", sep = ""),
                quote = FALSE, sep = "\t")
    
    ### store results in DE summary matrix:
    deResults[i+1, 1] <- i
    deResults[i+1, 2] <- nrow(degs_deseq2_001)
    deResults[i+1, 3] <- nrow(degs_deseq2_005)
    deResults[i+1, 4] <- nrow(degs_deseq2_010)
    deResults[i+1, 5] <- nrow(degs_edger_001)
    deResults[i+1, 6] <- nrow(degs_edger_005)
    deResults[i+1, 7] <- nrow(degs_edger_010)
    deResults[i+1, 8] <- length(int_001)
    deResults[i+1, 9] <- length(int_005)
    deResults[i+1,10] <- length(int_010)
    
    ### write out the DE summary matrix:
    write.table(deResults, paste(output_folder, "Output - DE Second Pass/deResuls_Summary", suffix, ".txt", sep = ""),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    ### write out the sampleData with appended normalization factors from RUVs:
    write.table(sampleDataW, paste(output_folder, "Output - DE Second Pass/sampleDataW", suffix, ".txt", sep = ""),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
  } # end if
  
} # end for
