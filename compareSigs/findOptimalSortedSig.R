# how should we define our sorted proteomics signature for Monocyte?
# assuming trends with WV deconvolution will hold true with other deconvolution methods
# 1- what is the best filter for accurately scoring monocytes (based on t-test, correct category per sample)
# 2- what is the best filter for predicting drug sensitivity (based on Ven, Aza+Ven AUC)

library(synapser)
library(DMEA)
library(tidyr)
library(tibble)
library(reshape2)
library(plyr)
library(dplyr)
synapser::synLogin()

evalMonoSig <- function(temp.sig, temp.name) {
  ## sorted proteomics
  # perform weighted voting on sorted proteomics
  sorted.wv <- DMEA::WV(global.df100, temp.sig, weight.values = "Log2FC")
  temp.wv <- sorted.wv$scores
  temp.wv$filter <- temp.name
  
  # evaluate accuracy based on t-test
  temp.test <- stats::t.test(temp.wv[grepl("CD14", temp.wv$Sample),]$WV,
                             temp.wv[!grepl("CD14", temp.wv$Sample),]$WV, 
                             "greater")
  temp.test.df <- as.data.frame(unlist(temp.test))
  colnames(temp.test.df)[1] <- "value"
  temp.test.df$variable <- rownames(temp.test.df)
  temp.test.df$filter <- temp.name
  test.df <- rbind(test.df, temp.test.df)
  
  # evaluate accuracy based on categorization: are monocyte scores highest for CD14 samples?
  
  
  
  
  return(list(wv = temp.wv, test = temp.test.df, corr = corr.df))
}

# load sorted proteomics
synapser::synLogin()
global.df <- read.csv(synapser::synGet("syn58895933")$path) # DIA
global.df <- global.df[ , which(colMeans(!is.na(global.df)) >= 0.75)] # 36 out of 48 samples are kept
# syn.test <- synapser::synGet("syn58914135") # for some reason, can't access pre-filtered version ???
rownames(global.df) <- global.df$Gene
global.df$Gene <- NULL
global.df <- as.data.frame(t(global.df))
global.df$Sample <- rownames(global.df)
global.df <- global.df[,c("Sample", colnames(global.df)[1:(ncol(global.df)-1)])]
global.df100 <- global.df[,colSums(is.na(global.df)) == 0]

# load sorted proteomics signature
cd14.sig <- read.csv(synapser::synGet("syn59424916")$path) # DIA
cd14.sig$sig <- FALSE
cd14.sig[cd14.sig$adj.P.Val <= 0.05,]$sig <- TRUE
cd14.sig$score <- -log(cd14.sig$adj.P.Val, 10) * cd14.sig$Log2FC
sig.filter <- cd14.sig[cd14.sig$sig,]
n.top <- c(1,5,10,25,50,100,200)
input.sigs <- list("Monocyte" = cd14.sig, "Progenitor" = cd34.sig, "MSC" = msc.sig)
filters <- list("noFilter" = NA, "maxFDR0.05" = NA,
                "absLog2FC" = n.top, "maxFDR0.05-absLog2FC" = n.top,
                "Log2FC" = n.top, "maxFDR0.05-Log2FC" = n.top, "FDR" = n.top, 
                "FDR-Log2FC" = n.top, "FDR-absLog2FC" = n.top)
# keep track of which filter produces the best result
test.df <- data.frame()
all.corr.df <- data.frame()
wv.df <- data.frame()
gmt <- readRDS("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/data/gmt_BeatAML_drug_MOA_2024-02-22.rds")
for (i in 1:length(filters)) {
  # filter signature appropriately
  if (grepl("absLog2FC", names(filters)[i])) {
    sigs <- list()
    for (k in 1:length(n.top)) {
      temp.name <- paste0("Top_", n.top[k],"_", names(filters[i]))
      if (grepl("maxFDR0.05", names(filters)[i])) {
        sigs[[temp.name]] <- sig.filter %>% slice_max(abs(Log2FC), n = n.top[k])
      } else if (grepl("FDR", names(filters)[i])) {
        sigs[[temp.name]] <- sig.filter %>% slice_max(abs(score), n = n.top[k])
      } else {
        sigs[[temp.name]] <- cd14.sig %>% slice_max(abs(Log2FC), n = n.top[k]) 
      }
    }
  } else if (grepl("Log2FC", names(filters)[i])) {
    sigs <- list()
    for (k in 1:length(n.top)) {
      temp.name <- paste0("Top_", n.top[k],"_", names(filters)[i])
      if (grepl("maxFDR0.05", names(filters)[i])) {
        sigs[[temp.name]] <- sig.filter %>% slice_max(Log2FC, n = n.top[k])
      } else if (grepl("FDR", names(filters)[i])) {
        sigs[[temp.name]] <- sig.filter %>% slice_max(score, n = n.top[k])
      } else {
        sigs[[temp.name]] <- cd14.sig %>% slice_max(Log2FC, n = n.top[k]) 
      }
    }
  } else if (is.data.frame(filters[[i]])) { # or use pre-filtered signature
    sigs <- list(filters[[i]])
    names(sigs) <- names(filters)[i]
  }
  
  for (j in 1:length(sigs)) {
    # evaluate signature
    evalResults <- evalMonoSig(sigs[[j]], names(sigs)[j]) 
    
    # store results
    wv.df <- rbind(wv.df, evalResults$wv)
    test.df <- rbind(test.df, evalResults$test)
    all.corr.df <- rbind(all.corr.df, evalResults$corr)
  }
}
write.csv(wv.df, "CD14_DIA_filters_WV.csv", row.names = FALSE)
write.csv(test.df, "Ttest_CD14_DIA_filters_WV.csv", row.names = FALSE)
write.csv(all.corr.df, "Correlation_BeatAML_drug_AUC_CD14_DIA_filters_WV.csv", row.names = FALSE)
write.csv(test.df[test.df$variable=="p.value",], "pValue_CD14_DIA_filters_WV.csv", row.names = FALSE)
write.csv(all.corr.df[all.corr.df$Drug == "Venetoclax",], "Venetoclax_BeatAML_drug_AUC_CD14_DIA_filters_WV.csv", row.names = FALSE)
write.csv(all.corr.df[all.corr.df$Drug == "Azacytidine - Venetoclax",], "AzaVen_BeatAML_drug_AUC_CD14_DIA_filters_WV.csv", row.names = FALSE)
write.csv(all.corr.df[all.corr.df$Drug == "Azacytidine",], "Azacytidine_BeatAML_drug_AUC_CD14_DIA_filters_WV.csv", row.names = FALSE)

## plot overview of results
# p-value
test.df <- test.df|>
  tidyr::separate(filter,into=c('Top','N_proteins','filterType'),sep='_',remove=FALSE)
test.df$N_proteins <- as.numeric(test.df$N_proteins)
test.df$`Filtered for FDR <= 0.05` <- FALSE
test.df[grepl("maxFDR0.05",test.df$filter),]$`Filtered for FDR <= 0.05` <- TRUE
test.df$`Sorting Method` <- test.df$filterType
test.df[test.df$`Filtered for FDR <= 0.05`, ]$`Sorting Method` <- stringr::str_split_fixed(test.df[test.df$`Filtered for FDR <= 0.05`, ]$filterType, "-", 2)[,2]
# any remaining "sorting method" with a dash is using the score column
test.df[which(test.df$`filterType` == "FDR-Log2FC"),]$`Sorting Method` <- "-LogFDR*Log2FC"
test.df[which(test.df$`filterType` == "FDR-absLog2FC"),]$`Sorting Method` <- "-LogFDR*|Log2FC|"
test.df[which(test.df$`Sorting Method` == "absLog2FC"),]$`Sorting Method` <- "|Log2FC|"
library(ggplot2)
ggplot(test.df[test.df$variable=="p.value" & !is.na(test.df$N_proteins),], 
       aes(x=N_proteins, y=-log(as.numeric(value), 10), color = `Sorting Method`, shape = `Filtered for FDR <= 0.05`)) + 
  scale_x_continuous(trans="log10") +
  geom_point() + theme_minimal() + xlab("Number of Proteins") + ylab("-Log(P-value)") + ggtitle("T-test: Monocyte scores are higher in CD14+ samples")
ggsave("pValue_CD14_DIA_filters_WV.pdf", width = 5, height = 5)

# Venetoclax
all.corr.df <- all.corr.df|>
  tidyr::separate(filter,into=c('Top','N_proteins','filterType'),sep='_',remove=FALSE)
all.corr.df$N_proteins <- as.numeric(all.corr.df$N_proteins)
all.corr.df$`Filtered for FDR <= 0.05` <- FALSE
all.corr.df[grepl("maxFDR0.05",all.corr.df$filter),]$`Filtered for FDR <= 0.05` <- TRUE
all.corr.df$`Sorting Method` <- all.corr.df$filterType
all.corr.df[all.corr.df$`Filtered for FDR <= 0.05`, ]$`Sorting Method` <- stringr::str_split_fixed(all.corr.df[all.corr.df$`Filtered for FDR <= 0.05`, ]$filterType, "-", 2)[,2]
# any remaining "sorting method" with a dash is using the score column
all.corr.df[which(all.corr.df$`filterType` == "FDR-Log2FC"),]$`Sorting Method` <- "-LogFDR*Log2FC"
all.corr.df[which(all.corr.df$`filterType` == "FDR-absLog2FC"),]$`Sorting Method` <- "-LogFDR*|Log2FC|"
all.corr.df[which(all.corr.df$`Sorting Method` == "absLog2FC"),]$`Sorting Method` <- "|Log2FC|"
all.corr.df$`-LogFDR` <- -log(all.corr.df$Pearson.q,10)
library(ggplot2)
ggplot(all.corr.df[all.corr.df$Drug == "Venetoclax" & !is.na(all.corr.df$N_proteins),], 
       aes(x=N_proteins, y=Pearson.est, color = `Sorting Method`, shape = `Filtered for FDR <= 0.05`)) + 
  scale_x_continuous(trans="log10") +
  geom_point() + theme_minimal() + xlab("Number of Proteins") + ylab("Pearson Correlation Estimate") + ggtitle("Venetoclax Sensitivity is Predicted by CD14+ Signature")
ggsave("Venetoclax_BeatAML_drug_AUC_CD14_DIA_filters_WV.pdf", width = 5, height = 5)

# AV
ggplot(all.corr.df[all.corr.df$Drug == "Azacytidine - Venetoclax" & !is.na(all.corr.df$N_proteins),], 
       aes(x=N_proteins, y=Pearson.est, color = `Sorting Method`, shape = `Filtered for FDR <= 0.05`)) + 
  scale_x_continuous(trans="log10") +
  geom_point() + theme_minimal() + xlab("Number of Proteins") + ylab("Pearson Correlation Estimate") + ggtitle("Aza+Ven Sensitivity is Predicted by CD14+ Signature")
ggsave("AzaVen_BeatAML_drug_AUC_CD14_DIA_filters_WV.pdf", width = 5, height = 5)

# Aza
ggplot(all.corr.df[all.corr.df$Drug == "Azacytidine" & !is.na(all.corr.df$N_proteins),], 
       aes(x=N_proteins, y=Pearson.est, color = `Sorting Method`, shape = `Filtered for FDR <= 0.05`)) + 
  scale_x_continuous(trans="log10") +
  geom_point() + theme_minimal() + xlab("Number of Proteins") + ylab("Pearson Correlation Estimate") + ggtitle("Azacytidine Sensitivity is Predicted by CD14+ Signature")
ggsave("Azacytidine_BeatAML_drug_AUC_CD14_DIA_filters_WV.pdf", width = 5, height = 5)

# combo
ggplot(na.omit(all.corr.df[all.corr.df$Drug == "Azacytidine" | 
                     all.corr.df$Drug == "Azacytidine - Venetoclax" | 
                     all.corr.df$Drug == "Venetoclax" & !is.na(all.corr.df$N_proteins),]), 
       aes(x=N_proteins, y=Pearson.est, color = `Sorting Method`, shape = `Filtered for FDR <= 0.05`)) + 
  scale_x_continuous(trans="log10") + facet_grid(rows = vars(Drug)) +
  geom_point() + xlab("Number of Proteins") + ylab("Pearson Correlation Estimate") + ggtitle("Drug Sensitivity is Predicted by CD14+ Signature")
ggsave("BeatAML_drug_AUC_CD14_DIA_filters_WV.pdf")

ggplot(na.omit(all.corr.df[all.corr.df$Drug == "Azacytidine" | 
                             all.corr.df$Drug == "Azacytidine - Venetoclax" | 
                             all.corr.df$Drug == "Venetoclax" & !is.na(all.corr.df$N_proteins),]), 
       aes(x=N_proteins, y=Pearson.est, color = `Sorting Method`, shape = `Filtered for FDR <= 0.05`)) + 
  scale_x_continuous(trans="log10") + facet_grid(cols = vars(Drug)) +
  geom_point() + xlab("Number of Proteins") + ylab("Pearson Correlation Estimate") + ggtitle("Drug Sensitivity is Predicted by CD14+ Signature")
ggsave("BeatAML_drug_AUC_CD14_DIA_filters_WV_v2.pdf")

## also try incorporating significance with size
ggplot(all.corr.df[all.corr.df$Drug == "Venetoclax" & !is.na(all.corr.df$N_proteins),], 
       aes(x=N_proteins, y=Pearson.est, color = `Sorting Method`, shape = `Filtered for FDR <= 0.05`, size = `-LogFDR`)) + 
  scale_x_continuous(trans="log10") +
  geom_point() + theme_minimal() + xlab("Number of Proteins") + ylab("Pearson Correlation Estimate") + ggtitle("Venetoclax Sensitivity is Predicted by CD14+ Signature")
ggsave("Venetoclax_BeatAML_drug_AUC_CD14_DIA_filters_WV_v2.pdf", width = 5, height = 5)

# AV
ggplot(all.corr.df[all.corr.df$Drug == "Azacytidine - Venetoclax" & !is.na(all.corr.df$N_proteins),], 
       aes(x=N_proteins, y=Pearson.est, color = `Sorting Method`, shape = `Filtered for FDR <= 0.05`, size = `-LogFDR`)) + 
  scale_x_continuous(trans="log10") +
  geom_point() + theme_minimal() + xlab("Number of Proteins") + ylab("Pearson Correlation Estimate") + ggtitle("Aza+Ven Sensitivity is Predicted by CD14+ Signature")
ggsave("AzaVen_BeatAML_drug_AUC_CD14_DIA_filters_WV_v2.pdf", width = 5, height = 5)

# Aza
ggplot(all.corr.df[all.corr.df$Drug == "Azacytidine" & !is.na(all.corr.df$N_proteins),], 
       aes(x=N_proteins, y=Pearson.est, color = `Sorting Method`, shape = `Filtered for FDR <= 0.05`, size = `-LogFDR`)) + 
  scale_x_continuous(trans="log10") +
  geom_point() + theme_minimal() + xlab("Number of Proteins") + ylab("Pearson Correlation Estimate") + ggtitle("Azacytidine Sensitivity is Predicted by CD14+ Signature")
ggsave("Azacytidine_BeatAML_drug_AUC_CD14_DIA_filters_WV_v2.pdf", width = 5, height = 5)

write.csv(test.df, "Ttest_CD14_DIA_filters_WV_v2.csv", row.names = FALSE)
write.csv(all.corr.df, "Correlation_BeatAML_drug_AUC_CD14_DIA_filters_WV_v2.csv", row.names = FALSE)
write.csv(test.df[test.df$variable=="p.value",], "pValue_CD14_DIA_filters_WV_v2.csv", row.names = FALSE)
write.csv(all.corr.df[all.corr.df$Drug == "Venetoclax",], "Venetoclax_BeatAML_drug_AUC_CD14_DIA_filters_WV_v2.csv", row.names = FALSE)
write.csv(all.corr.df[all.corr.df$Drug == "Azacytidine - Venetoclax",], "AzaVen_BeatAML_drug_AUC_CD14_DIA_filters_WV_v2.csv", row.names = FALSE)
write.csv(all.corr.df[all.corr.df$Drug == "Azacytidine",], "Azacytidine_BeatAML_drug_AUC_CD14_DIA_filters_WV_v2.csv", row.names = FALSE)