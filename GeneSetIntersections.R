
# GeneSetIntersections.R
# For publication by Goldman et al. 2016 on schizophrenia-derived glia
# Author and contact: Mikhail Osipovitch (mxo@sund.ku.dk)
#
# The following script opens the specified lists of differentially expressed gened generated
# by NormalizationAndDEAnalysis.R script, performs set intersections, and generates the resulting
# Venn diagrams.

### load required libraries:
library(VennDiagram)

### read the lists of differentially expressed genes obtained from "Normalization 
### and DE Analysis.R" script:
degs_005_08 <- read.table("Output - SCZ.8 vs CTR/Output - DE Second Pass/DEGs_INTs/DEGs_005 - W_3.txt", sep = "\t")
degs_005_29 <- read.table("Output - SCZ.29 vs CTR/Output - DE Second Pass/DEGs_INTs/DEGs_005 - W_1.txt", sep = "\t")
degs_005_51 <- read.table("Output - SCZ.51 vs CTR/Output - DE Second Pass/DEGs_INTs/DEGs_005 - W_7.txt", sep = "\t")

### filter genes by fold change (Log2 FC > 1) and split into up- and down-regulated:
degs_005_fc2_08_up   <- degs_005_08[which(degs_005_08$logFC >  1),]
degs_005_fc2_08_down <- degs_005_08[which(degs_005_08$logFC < -1),]
degs_005_fc2_29_up   <- degs_005_29[which(degs_005_29$logFC >  1),]
degs_005_fc2_29_down <- degs_005_29[which(degs_005_29$logFC < -1),]
degs_005_fc2_51_up   <- degs_005_51[which(degs_005_51$logFC >  1),]
degs_005_fc2_51_down <- degs_005_51[which(degs_005_51$logFC < -1),]

output_folder <- "Output - DE Genes Intersections/"
dir.create(output_folder)

labels <- c("SCZ.08","SCZ.29","SCZ.51")

###### Gene set intersection and Venn diagram for up-regulated genes:

l1 <- rownames(degs_005_fc2_08_up)
l2 <- rownames(degs_005_fc2_29_up)
l3 <- rownames(degs_005_fc2_51_up)

a1 <- length(unique(l1))
a2 <- length(unique(l2))
a3 <- length(unique(l3))

n12 <- length(unique(intersect(l1, l2)))
n23 <- length(unique(intersect(l2, l3)))
n13 <- length(unique(intersect(l1, l3)))

n123 <- length(Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3)))

png(paste(output_folder, "Venn 005_FC2 - UP.png", sep = ""), width = 800, height = 600, units = "px")

draw.triple.venn(a1, a2, a3, n12, n23, n13, n123,
                 category = labels, fill = c("#528881","#5e3c58","#886451"), euler.d = TRUE, scaled = FALSE,
                 cat.pos = c(340, 20, 180), cat.dist = c(0.07, 0.07, 0.05), cex = 2.9, cat.cex = 2.4,
                 fontfamily = c(rep("Verdana", 7)), cat.fontfamily = (c(rep("Verdana", 3))))
dev.off()
graphics.off()

degs_int_up <- cbind(degs_005_fc2_08_up[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3)),c(7,11)],
                     degs_005_fc2_29_up[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3)),c(7,11)],
                     degs_005_fc2_51_up[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3)),c(7,11)])

colnames(degs_int_up) <- c("logFC_SCZ08","FDR_SCZ08","logFC_SCZ29","FDR_SCZ29","logFC_SCZ51","FDR_SCZ51")

write.table(degs_int_up, paste(output_folder, "DEGs Intersection - FDR005_FC2 - UP.txt", sep = ""), sep = "\t", quote = FALSE)

rm(l1, l2, l3, a1, a2, a3, n12, n23, n13, n123)

###### Gene set intersection and Venn diagram for down-regulated genes:

l1 <- rownames(degs_005_fc2_08_down)
l2 <- rownames(degs_005_fc2_29_down)
l3 <- rownames(degs_005_fc2_51_down)

a1 <- length(unique(l1))
a2 <- length(unique(l2))
a3 <- length(unique(l3))

n12 <- length(unique(intersect(l1, l2)))
n23 <- length(unique(intersect(l2, l3)))
n13 <- length(unique(intersect(l1, l3)))

n123 <- length(Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3)))

png(paste(output_folder, "Venn 005_FC2 - DOWN.png", sep = ""), width = 800, height = 600, units = "px")

draw.triple.venn(a1, a2, a3, n12, n23, n13, n123,
                 category = labels, fill = c("#528881","#5e3c58","#886451"), euler.d = TRUE, scaled = FALSE,
                 cat.pos = c(340, 20, 180), cat.dist = c(0.07, 0.07, 0.05), cex = 2.9, cat.cex = 2.4,
                 fontfamily = c(rep("Verdana", 7)), cat.fontfamily = (c(rep("Verdana", 3))))
dev.off()
graphics.off()

degs_int_down <- cbind(degs_005_fc2_08_down[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3)),c(7,11)],
                       degs_005_fc2_29_down[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3)),c(7,11)],
                       degs_005_fc2_51_down[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3)),c(7,11)])

colnames(degs_int_up) <- c("logFC_SCZ08","FDR_SCZ08","logFC_SCZ29","FDR_SCZ29","logFC_SCZ51","FDR_SCZ51")

write.table(degs_int_up, paste(output_folder, "DEGs Intersection - FDR005_FC2 - DOWN.txt", sep = ""), sep = "\t", quote = FALSE)

rm(l1, l2, l3, a1, a2, a3, n12, n23, n13, n123)
