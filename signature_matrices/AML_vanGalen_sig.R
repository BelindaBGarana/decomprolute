# prep van Galen binary signature matrix
library(readxl); library(dplyr); library(reshape2)
setwd("~/OneDrive - PNNL/Documents/GitHub/decomprolute/signature_matrices")

#### 1. Load signatures ####
vg <- readxl::read_excel("aml_vanGalen_cellTypes.xlsx")

# just extract "tumor-derived, per cell-type"
vg <- vg[1:nrow(vg), 8:ncol(vg)]
colnames(vg) <- vg[1,]
vg <- vg[2:nrow(vg), -3]
vg <- vg[1:50,1:5] # only keep rows with gene names (na.omit didn't give desired result)

#### 2. Reformat for decomprolute ####
cellType <- colnames(vg)
Gene <- unique(unlist(vg))

sig <- data.frame(row.names = Gene)
sig[,cellType] <- 0

for (i in cellType) {
  cell.markers <- unique(unlist(vg[,i]))
  sig[cell.markers,i] <- 1
}

sig$NAME <- rownames(sig)
sig <- sig[,c("NAME", colnames(sig)[1:ncol(sig)-1])]
write.table(sig, "AML_vanGalen.txt", sep="\t", row.names = FALSE)
