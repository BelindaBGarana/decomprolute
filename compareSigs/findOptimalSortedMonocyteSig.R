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
  
  # evaluate accuracy based on t-test
  temp.test <- stats::t.test(temp.wv[grepl("CD14", temp.wv$Sample),]$WV,
                             temp.wv[!grepl("CD14", temp.wv$Sample),]$WV, 
                             "greater")
  temp.test.df <- as.data.frame(unlist(temp.test))
  colnames(temp.test.df)[1] <- "value"
  temp.test.df$variable <- rownames(temp.test.df)
  temp.test.df$filter <- temp.name
  test.df <- rbind(test.df, temp.test.df)
  
  # evaluate accuracy based on categorization: are monocyte scores highest for CD14 samples? - actually, would have to do CD34 & MSC scores also
  
  ## Beat AML
  # perform DMEA on Beat AML global proteomics
  DMEA.result <- DMEA::DMEA(BeatAML$drug, gmt, BeatAML$global, temp.sig)
  corr.df <- DMEA.result$corr.result
  corr.df$filter <- temp.name
  all.corr.df <- rbind(all.corr.df, corr.df)
  return(list(wv = temp.wv, test = temp.test.df, corr = corr.df))
}

load_BeatAML_for_DMEA <- function() {
  # load Beat AML drug AUC
  drug.BeatAML <- read.csv(synapser::synGet("syn51674470")$path)
  drug.BeatAML$X <- NULL
  colnames(drug.BeatAML)[1] <- "sample"
  drug.BeatAML.wide <- reshape2::dcast(drug.BeatAML, sample ~ inhibitor, value.var = "auc")
  colnames(drug.BeatAML)[1] <- "sample" # to match sample column name in camilo
  
  # other questions
  # does ratio of mono/prog in Beat AML correlate with sorted cell CD14/CD34 signature?
  globalFileDIA <- synapser::synGet("syn59429685") # filtered for max FDR of 0.05
  global.sig.DIA <- read.csv(globalFileDIA$path) # 1842 proteins
  global.sig.25 <- global.sig.DIA %>% slice_max(abs(Log2FC), n=25)
  
  # load Beat AML global data
  meta.BeatAML <- read.table(synapser::synGet("syn25807733")$path, 
                             sep = "\t", header = TRUE)
  global.BeatAML <- read.table(synapser::synGet("syn25714248")$path,
                               sep = "\t", header = TRUE)
  
  # change global.BeatAML column names from SampleID.abbrev to 
  # Barcode.ID to match drug.BeatAML
  global.ids <- names(global.BeatAML)
  
  # remove X and any 0's from start of each column name and then
  # replace SampleID.abbrev with Barcode.ID to match drug.BeatAML
  for(i in seq_len(length(global.ids))){
    global.ids[i] <- substr(global.ids[i], 2, nchar(global.ids[i]))
    
    if(substring(global.ids[i], 1, 1) == 0){
      global.ids[i] <- substr(global.ids[i], 2, nchar(global.ids[i]))
    }
    
    if(global.ids[i] %in% meta.BeatAML$SampleID.abbrev){
      global.ids[i] <- meta.BeatAML[meta.BeatAML$SampleID.abbrev == global.ids[i], ]$Barcode.ID
    }
  }
  
  # replace global.BeatAML column names 
  names(global.BeatAML) <- global.ids
  
  # transpose global.BeatAML so that first column is Barcode.ID and 
  # rest of columns are gene symbols
  global.BeatAML <- as.data.frame(t(global.BeatAML))
  
  # make first column sample
  global.BeatAML[, "sample"] <- rownames(global.BeatAML)
  global.BeatAML <- 
    global.BeatAML[ , c("sample", 
                        names(global.BeatAML[ , 1:(ncol(global.BeatAML)-1)]))]
  
  # filter for genes with 100% coverage in Beat AML database
  global100 <- global.BeatAML[,colSums(is.na(global.BeatAML)) == 0]
  return(list(drug = drug.BeatAML.wide, global = global100))
}
BeatAML <- load_BeatAML_for_DMEA()

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
filters <- list("noFilter" = cd14.sig, "maxFDR0.05" = sig.filter,
                "absLog2FC" = n.top, "maxFDR0.05-absLog2FC" = n.top,
                "Log2FC" = n.top, "maxFDR0.05-Log2FC" = n.top, "FDR" = n.top, 
                "FDR-Log2FC" = n.top, "FDR-absLog2FC" = n.top)
# keep track of which filter produces the best result
test.df <- data.frame()
all.corr.df <- data.frame()
gmt <- readRDS("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/data/gmt_BeatAML_drug_MOA_2024-02-22.rds")
for (i in 1:length(filters)) {
  if (grepl("-", names(filters)[i])) {
    
  } else if (grepl("absLog2FC", names(filters)[i])) {
    sigs <- list()
    for (k in 1:length(n.top)) {
      temp.name <- paste0("Top_", n.top[k],"_absLog2FC")
      sigs[[temp.name]] <- cd14.sig %<% slice_max(abs(Log2FC), n = n.top[k])
    }
  } else if (grepl("Log2FC", names(filters)[i])) {
    sigs <- list()
    for (k in 1:length(n.top)) {
      
    }
  } else if (is.data.frame(filters[[i]])) {
    sigs <- list(filters[[i]])
    names(sigs) <- names(filters)[i]
  }
  for (j in 1:length(sigs)) {
    evalResults <- evalMonoSig(sigs[[j]], names(sigs)[j]) 
  }
}
write.csv(test.df, "Ttest_CD14_DIA_filters_WV.csv", row.names = FALSE)
write.csv(all.corr.df, "Correlation_BeatAML_drug_AUC_CD14_DIA_filters_WV.csv", row.names = FALSE)
write.csv(test.df[test.df$variable=="p.value",], "pValue_CD14_DIA_filters_WV.csv", row.names = FALSE)
write.csv(all.corr.df[all.corr.df$Drug == "Venetoclax",], "Venetoclax_BeatAML_drug_AUC_CD14_DIA_filters_WV.csv", row.names = FALSE)
write.csv(all.corr.df[all.corr.df$Drug == "Azacytidine - Venetoclax",], "AzaVen_BeatAML_drug_AUC_CD14_DIA_filters_WV.csv", row.names = FALSE)
write.csv(all.corr.df[all.corr.df$Drug == "Azacytidine",], "Azacytidine_BeatAML_drug_AUC_CD14_DIA_filters_WV.csv", row.names = FALSE)