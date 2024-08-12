# combine sorted AML proteomics signatures for deconvolution purposes
library(synapser); library(dplyr); library(reshape2)
setwd("~/OneDrive - PNNL/Documents/GitHub/decomprolute/signature_matrices")

#### 1. Load signatures ####
synapser::synLogin()
syn.id <- "syn59448202"
sig <- read.csv(synapser::synGet(syn.id)$path)

#### 2. Reformat for decomprolute ####
### extract sigs for each cell type vs. all others
no.filter <- sig[is.na(sig$Filter),]

## global
global <- no.filter[no.filter$Feature_type == "Gene",]
global.DIA <- dplyr::distinct(global[global$method == "DIA",])
global.TMT <- dplyr::distinct(global[global$method == "TMT",])
global.DIA <- reshape2::dcast(global.DIA, Feature ~ Contrast, value.var = "Log2FC")
global.TMT <- reshape2::dcast(global.TMT, Feature ~ Contrast, value.var = "Log2FC")

# for now, merge DIA & TMT sigs because CD14 & CD34 are only in TMT and MSC are only in DIA
global <- merge(global.TMT, global.DIA, all=TRUE)
colnames(global) <- c("NAME", "Monocyte", "Progenitor", "MSC")
write.table(global, "globalSortedAML.txt", sep="\t", row.names = FALSE)

## phospho - only have TMT
phospho <- no.filter[no.filter$Feature_type == "SUB_SITE",]
#phospho.DIA <- dplyr::distinct(phospho[phospho$method == "DIA",])
phospho.TMT <- dplyr::distinct(phospho[phospho$method == "TMT",])
#phospho.DIA <- reshape2::dcast(phospho.DIA, Feature ~ Contrast, value.var = "Log2FC")
phospho.TMT <- reshape2::dcast(phospho.TMT, Feature ~ Contrast, value.var = "Log2FC")

# for now, only have DIA which only has CD14 & CD34
#phospho <- merge(phospho.TMT, phospho.DIA, all=TRUE)
phospho <- phospho.TMT
#colnames(phospho) <- c("NAME", "Monocyte", "Progenitor", "MSC")
colnames(phospho) <- c("NAME", "Monocyte", "Progenitor")
write.table(phospho, "phosphoSortedAML.txt", sep="\t", row.names = FALSE)

#### 3. Filter for drug sensitivity ####
filters <- na.omit(unique(sig$Filter))
drugs <- c("Aza", "Ven")
drug.filters <- c()
for (i in drugs) {
  drug.filters <- unique(c(drug.filters, filters[grepl(i, filters)])) 
}
# name.map <- list("NAME" = "Feature", 
#                 "Monocyte" = "CD14_Pos_vs_Neg", 
#                  "Progenitor" = "CD34_Pos_vs_Neg",
#                 "MSC" = "MSC_MSC_vs_Non_MSC")
name.map <- list("Feature" = "NAME", 
                 "CD14_Pos_vs_Neg" = "Monocyte", 
                 "CD34_Pos_vs_Neg" = "Progenitor",
                 "MSC_MSC_vs_Non_MSC" = "MSC")
data.types <- list("global" = "Gene",
                   "phospho" = "SUB_SITE")
for (i in drug.filters) {
  ### extract sigs for each cell type vs. all others
  no.filter <- sig[sig$Filter == i,]
  
  ## for each data type
  for (j in 1:length(data.types)) {
    global <- no.filter[no.filter$Feature_type == data.types[[j]],]
    global.DIA <- dplyr::distinct(global[global$method == "DIA",])
    global.TMT <- dplyr::distinct(global[global$method == "TMT",])
    if (nrow(global.DIA) > 0) {global.DIA <- reshape2::dcast(global.DIA, Feature ~ Contrast, value.var = "Log2FC")}
    if (nrow(global.TMT) > 0) {global.TMT <- reshape2::dcast(global.TMT, Feature ~ Contrast, value.var = "Log2FC")}
    
    # for now, merge DIA & TMT sigs because CD14 & CD34 are only in TMT and MSC are only in DIA
    if (nrow(global.DIA) > 0 & nrow(global.TMT) > 0) {
      global <- merge(global.TMT, global.DIA, all=TRUE)
    } else if (nrow(global.DIA) > 0) {
      global <- global.DIA
    } else {
      global <- global.TMT
    }
    
    # fix colnames
    colnames(global) <- unlist(name.map[colnames(phospho.TMT)])
    write.table(global, paste(names(data.types)[j], i, "sortedAML.txt", sep = "_"), sep="\t", row.names = FALSE)
  }
}

#### 4. Select top gene symbols ####
n.genes <- c(1, 5, 10, 25, 50, 100)
filters <- na.omit(unique(sig$Filter))
drugs <- c("Aza", "Ven")
drug.filters <- c()
for (i in drugs) {
  drug.filters <- unique(c(drug.filters, filters[grepl(i, filters)])) 
}
drug.filters <- c(drug.filters, "")
name.map <- list("Feature" = "NAME", 
                 "CD14_Pos_vs_Neg" = "Monocyte", 
                 "CD34_Pos_vs_Neg" = "Progenitor",
                 "MSC_MSC_vs_Non_MSC" = "MSC")
data.types <- list("global" = "Gene",
                   "phospho" = "SUB_SITE")
for (i in drug.filters) {
  message("drug filter ", i)
  ### extract sigs for each cell type vs. all others
  if (i == "") {
    no.filter <- sig[is.na(sig$Filter),]
    temp.fname <- "sortedAML.txt"
  } else {
    no.filter <- sig[sig$Filter == i,] 
    temp.fname <- paste(i, "sortedAML.txt", sep = "_")
  }
  
  ## for each data type
  for (j in 1:length(data.types)) {
    message("data type ", names(data.types)[j])
    global <- na.omit(no.filter[no.filter$Feature_type == data.types[[j]] & 
                                  no.filter$Contrast %in% names(name.map)[2:4],])
    global.DIA <- na.omit(dplyr::distinct(global[global$method == "DIA",]))
    global.TMT <- na.omit(dplyr::distinct(global[global$method == "TMT",]))
    
    for (k in 1:length(n.genes)) {
      if (nrow(global.DIA) > 0) {
        temp.sigs.DIA <- list()
        for (m in 2:length(name.map)) {
          temp.global.DIA <- global.DIA[global.DIA$Contrast == names(name.map)[m],] 
          if (nrow(temp.global.DIA) > 0) {
            temp.sigs.DIA[[names(name.map)[m]]] <- temp.global.DIA %>% 
              slice_max(order_by = abs(Log2FC), n = n.genes[k])
          }
          temp.global.DIA.sig <- data.table::rbindlist(temp.sigs.DIA, use.names = FALSE)
        }
        temp.global.DIA <- reshape2::dcast(temp.global.DIA.sig, Feature ~ Contrast, value.var = "Log2FC")
        message("compiled DIA sigs")
      }
      if (nrow(global.TMT) > 0) {
        temp.sigs.TMT <- list()
        for (m in 2:length(name.map)) {
          temp.global.TMT <- global.TMT[global.TMT$Contrast == names(name.map)[m],] 
          if (nrow(temp.global.TMT) > 0) {
            temp.sigs.TMT[[names(name.map)[m]]] <- temp.global.TMT %>% 
              slice_max(order_by = abs(Log2FC), n = n.genes[k])
          }
          temp.global.TMT.sig <- data.table::rbindlist(temp.sigs.TMT, use.names = FALSE)
        }
        temp.global.TMT <- reshape2::dcast(temp.global.TMT.sig, Feature ~ Contrast, value.var = "Log2FC")
        message("compiled TMT sigs")
      }
      
      # for now, merge DIA & TMT sigs because CD14 & CD34 are only in TMT and MSC are only in DIA
      if (nrow(global.DIA) > 0 & nrow(global.TMT) > 0) {
        temp.global <- merge(temp.global.TMT, temp.global.DIA, all=TRUE)
      } else if (nrow(global.DIA) > 0) {
        temp.global <- temp.global.DIA
      } else {
        temp.global <- temp.global.TMT
      }
      
      # fix colnames
      colnames(temp.global) <- unlist(name.map[colnames(temp.global)])
      write.table(temp.global, paste0("top_", n.genes[k], "_", 
                                      names(data.types)[j], "_", temp.fname), 
                  sep="\t", row.names = FALSE) 
    }
  }
}
