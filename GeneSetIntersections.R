
# GeneSetIntersections.R
# For publication by Goldman et al. 2016 on schizophrenia-derived glia
# Author and contact: Mikhail Osipovitch (mxo@sund.ku.dk)
#
# The following script opens the specified lists of differentially expressed genes generated
# by NormalizationAndDEAnalysis.R script, performs set intersections, generates Venn diagrams
# with shared up- and down-regulated genes, and outputs the resulting intersection gene list
# with DE information (fold changed and P Values) pulled from individual SCZ line comparisons
# to pooled controls and the comparison of pooled SCZ lines to pooled controls.

### load required libraries:
library(VennDiagram)

### read the lists of differentially expressed genes obtained from "NormalizationAndDEAnalysis.R" script:
degs_005_08  <- read.table("Output - SCZ.8 vs CTR/Output - DE Second Pass/DEGs_INTs/DEGs_005 - W_3.txt", sep = "\t")
degs_005_29  <- read.table("Output - SCZ.29 vs CTR/Output - DE Second Pass/DEGs_INTs/DEGs_005 - W_1.txt", sep = "\t")
degs_005_51  <- read.table("Output - SCZ.51 vs CTR/Output - DE Second Pass/DEGs_INTs/DEGs_005 - W_7.txt", sep = "\t")
degs_005_164 <- read.table("Output - SCZ.164 vs CTR/Output - DE Second Pass/DEGs_INTs/DEGs_005 - W_3.txt", sep = "\t")

### filter genes by fold change (Log2 FC > 1) and split into up- and down-regulated:
degs_005_fc2_08_up    <- degs_005_08[which(degs_005_08$logFC >  1),]
degs_005_fc2_08_down  <- degs_005_08[which(degs_005_08$logFC < -1),]
degs_005_fc2_29_up    <- degs_005_29[which(degs_005_29$logFC >  1),]
degs_005_fc2_29_down  <- degs_005_29[which(degs_005_29$logFC < -1),]
degs_005_fc2_51_up    <- degs_005_51[which(degs_005_51$logFC >  1),]
degs_005_fc2_51_down  <- degs_005_51[which(degs_005_51$logFC < -1),]
degs_005_fc2_164_up   <- degs_005_164[which(degs_005_164$logFC >  1),]
degs_005_fc2_164_down <- degs_005_164[which(degs_005_164$logFC < -1),]

output_folder <- "Output - DE Genes Intersections/"
dir.create(output_folder)

labels <- c("SCZ08","SCZ29","SCZ51","SCZ164")

###### Gene set intersection and Venn diagram for up-regulated genes:

l1 <- rownames(degs_005_fc2_08_up)
l2 <- rownames(degs_005_fc2_29_up)
l3 <- rownames(degs_005_fc2_51_up)
l4 <- rownames(degs_005_fc2_164_up)

png(paste(output_folder, "Venn 005_FC2 - UP.png", sep = ""), width = 800, height = 700, units = "px")  
vd <- venn.diagram(list(v1 = l1, v2 = l2, v3 = l3, v4 = l4),
                   main = "Up-regulated genes in each of 4 SCZ hGPC lines\ncompared to 3 pooled control hGPC lines\n5% FDR, L2FC > 1",
                   main.cex = 2.5, main.fontfamily = "Verdana",
                   category.names = labels, fill = c("#528881","#5e3c58","#886459","#2a334f"),
                   fontfamily = "Verdana", cat.fontfamily = "Verdana", cex = 3, cat.cex = 2.7,
                   margin = c(0.1,0.1), filename = NULL)
grid.draw(vd)
graphics.off()

### compile the intersection gene list with fold changes and P Values calculated by edgeR:
degs_int_up <- cbind(degs_005_fc2_08_up[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3, v4 = l4)),c(7,11)],
                     degs_005_fc2_29_up[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3, v4 = l4)),c(7,11)],
                     degs_005_fc2_51_up[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3, v4 = l4)),c(7,11)],
                     degs_005_fc2_164_up[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3, v4 = l4)),c(7,11)])

###### Gene set intersection and Venn diagram for down-regulated genes:

l1 <- rownames(degs_005_fc2_08_down)
l2 <- rownames(degs_005_fc2_29_down)
l3 <- rownames(degs_005_fc2_51_down)
l4 <- rownames(degs_005_fc2_164_down)

png(paste(output_folder, "Venn 005_FC2 - DOWN.png", sep = ""), width = 800, height = 700, units = "px")  
vd <- venn.diagram(list(v1 = l1, v2 = l2, v3 = l3, v4 = l4),
                   main = "Down-regulated genes in each of 4 SCZ hGPC lines\ncompared to 3 pooled control hGPC lines\n5% FDR, L2FC < -1",
                   main.cex = 2.5, main.fontfamily = "Verdana",
                   category.names = labels, fill = c("#528881","#5e3c58","#886459","#2a334f"),
                   fontfamily = "Verdana", cat.fontfamily = "Verdana", cex = 3, cat.cex = 2.7,
                   margin = c(0.1,0.1), filename = NULL)
grid.draw(vd)
graphics.off()

### compile the intersection gene list with fold changes and P Values calculated by edgeR:
degs_int_down <- cbind(degs_005_fc2_08_down[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3, v4 = l4)),c(7,11)],
                       degs_005_fc2_29_down[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3, v4 = l4)),c(7,11)],
                       degs_005_fc2_51_down[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3, v4 = l4)),c(7,11)],
                       degs_005_fc2_164_down[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3, v4 = l4)),c(7,11)])

###### Read data from DE comparison of pooled SCZ to pooled CTR:

degs_005_pooled <- read.table("Output - SCZ vs CTR/Output - DE Second Pass/DEGs_INTs/DEGs_005 - W_9.txt", sep = "\t")

### bind together the intersection of up- and down-regulated genes: 
degs_int <- rbind(degs_int_up, degs_int_down)

colnames(degs_int) <- c("logFC_SCZ08","FDR_SCZ08","logFC_SCZ29","FDR_SCZ29","logFC_SCZ51","FDR_SCZ51","logFC_SCZ164","FDR_SCZ164")

### Description of column names:
# logFC_SCZ08:  edgeR-calculated log2 fold change from DE comparison
#               of SCZ line 8 vs. pooled CTR lines (22, 37, and 205)
# FDR_SCZ08:    edgeR-calculated FDR-corrected P Value from DE comparison
#               of SCZ line 8 vs. pooled CTR lines (22, 37, and 205)
# logFC_SCZ29:  edgeR-calculated log2 fold change from DE comparison
#               of SCZ line 29 vs. pooled CTR lines (22, 37, and 205)
# FDR_SCZ29:    edgeR-calculated FDR-corrected P Value from DE comparison
#               of SCZ line 29 vs. pooled CTR lines (22, 37, and 205)
# logFC_SCZ51:  edgeR-calculated log2 fold change from DE comparison
#               of SCZ line 51 vs. pooled CTR lines (22, 37, and 205)
# FDR_SCZ51:    edgeR-calculated FDR-corrected P Value from DE comparison
#               of SCZ line 51 vs. pooled CTR lines (22, 37, and 205)
# logFC_SCZ164: edgeR-calculated log2 fold change from DE comparison
#               of SCZ line 164 vs. pooled CTR lines (22, 37, and 205)
# FDR_SCZ164:   edgeR-calculated FDR-corrected P Value from DE comparison
#               of SCZ line 164 vs. pooled CTR lines (22, 37, and 205)
# logFC_pooled: edgeR-calculated log2 fold change from DE comparison
#               of pooled SCZ lines (8, 29, and 51) vs. pooled CTR lines (22, 37, and 205)
# FDR_pooled:   edgeR-calculated FDR-corrected P Value from DE comparison
#               of pooled SCZ lines (8, 29, and 51) vs. pooled CTR lines (22, 37, and 205)

### append data from pooled comparison:
degs_int$logFC_pooled <- degs_005_pooled[rownames(degs_int), "logFC"]
degs_int$FDR_pooled <- degs_005_pooled[rownames(degs_int), "FDR"]

### order genes by FDR-corrected P Value from pooled comparison:
degs_int <- degs_int[order(degs_int$FDR_pooled),]

### write out the compiled intersection DE gene list:
write.table(degs_int, paste(output_folder, "DEGs Intersection - FDR005_FC2.txt", sep = ""), sep = "\t", quote = FALSE)
