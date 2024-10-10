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

evalOneSig <- function(global.df100, frac.df, temp.sig, temp.marker, temp.name) {
  ## sorted proteomics
  # perform weighted voting on sorted proteomics
  temp.sig$Gene <- rownames(temp.sig)
  temp.sig <- temp.sig[,c("Gene", colnames(temp.sig)[1])]
  sorted.wv <- DMEA::WV(global.df100, temp.sig, weight.values = "Log2FC")
  temp.wv <- sorted.wv$scores
  temp.wv$filter <- temp.name
  
  # evaluate accuracy based on t-test
  temp.test <- stats::t.test(temp.wv[grepl(temp.marker, temp.wv$Sample),]$WV,
                             temp.wv[!grepl(temp.marker, temp.wv$Sample),]$WV, 
                             "greater")
  temp.test.df <- as.data.frame(unlist(temp.test))
  colnames(temp.test.df)[1] <- "value"
  temp.test.df$variable <- rownames(temp.test.df)
  temp.test.df$filter <- temp.name
  
  # compare to known cell fractions
  temp.wv.frac <- merge(frac.df, temp.wv, by="Sample")
  frac.corr <- DMEA::rank_corr(temp.wv.frac)
  frac.corr.result <- frac.corr$result
  frac.corr.result$filter <- temp.name

  ## Beat AML
  # perform DMEA on Beat AML global proteomics
  DMEA.result <- DMEA::DMEA(BeatAML$drug, gmt, BeatAML$global, temp.sig)
  corr.df <- DMEA.result$corr.result
  corr.df$filter <- temp.name
  return(list(wv = temp.wv, test = temp.test.df, drug.corr = corr.df, 
              frac.corr = frac.corr.result))
}

evalSig <- function(global.df100, frac.df, sig.matrix,
                    markers = colnames(sig.matrix),
                    temp.name = "Top_25_FDR-absLog2FC") {
  # evaluate each signature
  wv.list <- list()
  test.list <- list()
  drug.corr.list <- list()
  frac.corr.list <- list()
  for (i in 1:length(markers)) {
    temp.result <- evalOneSig(global.df100, frac.df, 
                              as.data.frame(sig.matrix[,markers[i]]), 
                              markers[i], temp.name = temp.name)
    wv.list[[markers[i]]] <- temp.result$wv
    test.list[[markers[i]]] <- temp.result$test
    drug.corr.list[[markers[i]]] <- temp.result$drug.corr
    frac.corr.list[[markers[i]]] <- temp.result$frac.corr
  }
  wv.df <- data.table::rbindlist(wv.list, use.names = TRUE, idcol = "Cell Type")
  test.df <- data.table::rbindlist(test.list, use.names = TRUE, idcol = "Cell Type")
  drug.corr.df <- data.table::rbindlist(drug.corr.list, use.names = TRUE, 
                                        idcol = "Cell Type")
  frac.corr.df <- data.table::rbindlist(frac.corr.list, use.names = TRUE, 
                                        idcol = "Cell Type")
  p.df <- test.df[test.df$variable == "p.value",]
  
  # determine classification accuracy - is CD14 highest for CD14 samples?
  class.df <- reshape2::dcast(wv.df, Sample ~ `Cell Type`, value.var = "WV")
  class.df$Classification <- NA
  for (i in 1:nrow(class.df)) {
    temp.data <- class.df[i,]
    marker.names <- colnames(class.df)[2:(ncol(class.df)-1)]
    class.df$Classification[i] <- marker.names[which(max(temp.data))]
  }
  class.df$Accurate <- FALSE
  for (i in markers) {
    class.df[grepl(i, class.df$Sample) & 
               grepl(i, class.df$Classification),]$Accurate <- TRUE
  }
  
  # compile accuracy for whole signature
  agg.tab1 <- plyr::ddply(class.df, .(`Cell Type`), summarize,
                          N_correct = sum(Accurate, na.rm = TRUE),
                          N_total = length(na.omit(Accurate)))
  agg.tab1$Accuracy <- 100*agg.tab1$N_correct/agg.tab1$N_total
  
  # count % of samples correctly guessed for each signature, marker
  agg.tab2 <- plyr::ddply(class.df, .(`Cell Type`, Classification), summarize,
                          N_correct = sum(Accurate, na.rm = TRUE),
                          N_total = length(na.omit(Accurate)))
  agg.tab2$Accuracy <- 100*agg.tab2$N_correct/agg.tab2$N_total

  return(list(wv = wv.df, p = p.df, 
              drug.corr = drug.corr.df, frac.corr = frac.corr.df,
              accuracy = agg.tab1, marker.accuracy = agg.tab2))
}

testNprots <- function(global.df100, frac.df, sig.list, select.by = "Log2FC", 
                       value.var = select.by, n.prots = 25) {
  # filter each signature for top N proteins
  sig.names <- names(sig.list)
  filtered.sigs <- list()
  for (i in 1:length(sig.list)) {
    temp.sig <- sig.list[[i]]
    temp.sig$criteria <- temp.sig[,select.by]
    filtered.sigs[[sig.names[i]]] <- temp.sig %>% slice_max(criteria, n=n.prots)
  }
  filtered.sigs.df <- data.table::rbindlist(filtered.sigs, use.names = TRUE, idcol = "Signature")
  sig.matrix <- reshape2::dcast(filtered.sigs.df, Gene ~ Signature, mean,
                                value.var = value.var)
  rownames(sig.matrix) <- sig.matrix$Gene
  sig.matrix$Gene <- NULL
  
  # test signature matrix
  temp.name <- paste0("Top_", n.prots, "_", select.by)
  sigResults <- evalSig(global.df100, frac.df, sig.matrix,
                        markers = colnames(sig.matrix),
                        temp.name = temp.name)
  return(sigResults)
}

compareSigs <- function(global.df100, frac.df, all.sig.list, 
                        select.by = "Log2FC", value.var = select.by, 
                        all.n.prots = c(1,5,10,15,20,25,50,100),
                        Filter = rep(NA, length(all.sig.list))) {
  # for each signature type, test various # of proteins
  wv.list <- list()
  p.list <- list()
  drug.corr.list <- list()
  frac.corr.list <- list()
  acc.list <- list()
  marker.acc.list <- list()
  for (i in 1:length(all.n.prots)) {
    temp.wv.list <- list()
    temp.p.list <- list()
    temp.drug.corr.list <- list()
    temp.frac.corr.list <- list()
    temp.acc.list <- list()
    temp.marker.acc.list <- list()
    for (j in 1:length(all.sig.list)) {
      sig.name <- names(all.sig.list)[j]
      temp.result <- testNprots(global.df100, frac.df, all.sig.list[[j]],
                                select.by = select.by, value.var = value.var,
                                n.prots = all.n.prots[i])
      temp.wv.list[[sig.name]] <- temp.result$wv
      temp.p.list[[sig.name]] <- temp.result$p
      temp.drug.corr.list[[sig.name]] <- temp.result$drug.corr
      temp.frac.corr.list[[sig.name]] <- temp.result$frac.corr
      temp.acc.list[[sig.name]] <- temp.result$accuracy
      temp.marker.acc.list[[sig.name]] <- temp.result$marker.accuracy
    }
    wv.list[[as.character(all.n.prots[i])]] <- 
      data.table::rbindlist(temp.wv.list, use.names = TRUE, idcol = "Signature")
    p.list[[as.character(all.n.prots[i])]] <- 
      data.table::rbindlist(temp.p.list, use.names = TRUE, idcol = "Signature")
    drug.corr.list[[as.character(all.n.prots[i])]] <- 
      data.table::rbindlist(temp.drug.corr.list, use.names = TRUE, 
                            idcol = "Signature")
    frac.corr.list[[as.character(all.n.prots[i])]] <- 
      data.table::rbindlist(temp.frac.corr.list, use.names = TRUE, 
                                          idcol = "Signature")
    acc.list[[as.character(all.n.prots[i])]] <- 
      data.table::rbindlist(temp.acc.list, use.names = TRUE, 
                                               idcol = "Signature")
    marker.acc.list[[as.character(all.n.prots[i])]] <- 
      data.table::rbindlist(temp.marker.acc.list, use.names = TRUE, 
                                         idcol = "Signature")
  }
  wv.df <- data.table::rbindlist(wv.list, use.names = TRUE, 
                                 idcol = "Number of Proteins")
  p.df <- data.table::rbindlist(p.list, use.names = TRUE, 
                                idcol = "Number of Proteins")
  drug.corr.df <- data.table::rbindlist(drug.corr.list, use.names = TRUE, 
                                        idcol = "Number of Proteins")
  frac.corr.df <- data.table::rbindlist(frac.corr.list, use.names = TRUE, 
                                        idcol = "Number of Proteins")
  acc.df <- data.table::rbindlist(acc.list, use.names = TRUE, 
                                        idcol = "Number of Proteins")
  marker.acc.df <- data.table::rbindlist(marker.acc.list, use.names = TRUE, 
                                        idcol = "Number of Proteins")
  
  # add filter info
  Signature <- names(all.sig.list)
  sig.df <- data.frame(Signature, Filter)
  wv.df <- merge(sig.df, wv.df, by="Signature", all = TRUE)
  p.df <- merge(sig.df, p.df, by="Signature", all = TRUE)
  drug.corr.df <- merge(sig.df, drug.corr.df, by="Signature", all = TRUE)
  frac.corr.df <- merge(sig.df, frac.corr.df, by="Signature", all = TRUE)
  acc.df <- merge(sig.df, acc.df, by="Signature", all = TRUE)
  marker.acc.df <- merge(sig.df, marker.acc.df, by="Signature", all = TRUE)
  
  # plot results
  ggplot(p.df, aes(x=`Number of Proteins`, y=-log(p.value, 10), 
                   color = Signature, shape = Filter)) + 
    geom_point() + theme_minimal() + ylab("-Log(P-value)") + 
    ggtitle("T-test: Cell type scores are higher in their respective samples")
  ggsave("pValue_DIA_WV.pdf", width = 5, height = 5)
  
  doi <- c("Azacytidine", "Venetoclax", "Azacytidine - Venetoclax")
  drug.corr.df$`Drug Treatment` <- NA
  drug.corr.df[drug.corr.df$Drug == "Azacytidine",]$`Drug Treatment` <- "Aza"
  drug.corr.df[drug.corr.df$Drug == "Azacytidine - Venetoclax",]$`Drug Treatment` <- "Aza + Ven"
  drug.corr.df[drug.corr.df$Drug == "Venetoclax",]$`Drug Treatment` <- "Ven"
  ggplot(drug.corr.df[drug.corr.df$Drug %in% doi], aes(x=`Number of Proteins`, y=Pearson.est, 
                           color = Signature, shape = Filter)) + 
    geom_point() + theme_minimal() + ylab("Pearson Correlation Estimate") + 
    facet_wrap(`Cell Type` ~ `Drug Treatment`) +
    ggtitle("Cell-specific signatures predict drug response")
  ggsave("drugCorr_DIA_WV.pdf", width = 5, height = 5)
  
  ggplot(drug.corr.df[drug.corr.df$Drug == "Venetoclax",], 
         aes(x=`Number of Proteins`, y=Pearson.est, 
             color = Signature, shape = Filter)) + 
    geom_point() + theme_minimal() + ylab("Pearson Correlation Estimate") + 
    facet_wrap(~ `Cell Type`) +
    ggtitle("Cell-specific signatures predict Venetoclax response")
  ggsave("VenCorr_DIA_WV.pdf", width = 5, height = 5)
  
  ggplot(drug.corr.df[drug.corr.df$Drug == "Azacytidine - Venetoclax",], 
         aes(x=`Number of Proteins`, y=Pearson.est, 
             color = Signature, shape = Filter)) + 
    geom_point() + theme_minimal() + ylab("Pearson Correlation Estimate") + 
    facet_wrap(~ `Cell Type`) +
    ggtitle("Cell-specific signatures predict Aza + Ven response")
  ggsave("AzaVenCorr_DIA_WV.pdf", width = 5, height = 5)
  
  ggplot(drug.corr.df[drug.corr.df$Drug == "Azacytidine",], 
         aes(x=`Number of Proteins`, y=Pearson.est, 
             color = Signature, shape = Filter)) + 
    geom_point() + theme_minimal() + ylab("Pearson Correlation Estimate") + 
    facet_wrap(~ `Cell Type`) +
    ggtitle("Cell-specific signatures predict Azacytidine response")
  ggsave("AzaCorr_DIA_WV.pdf", width = 5, height = 5)
  
  ggplot(frac.corr.df, aes(x=`Number of Proteins`, y=Pearson.est,
                           color = Signature, shape = Filter)) + 
    geom_point() + theme_minimal() + ylab("Pearson Correlation Estimate") + 
    facet_wrap(`Measured Cell Type` ~ `Cell Type`) +
    ggtitle("Cell-specific signatures predict cell fraction")
  ggsave("fracCorr_DIA_WV.pdf", width = 5, height = 5)
  
  ggplot(acc.df, aes(x=`Number of Proteins`, y=Accuracy, 
             color = Signature, shape = Filter)) + 
    geom_point() + theme_minimal() + ylab("% Correctly Classified") + 
    facet_wrap(~ `Cell Type`) +
    ggtitle("Cell-specific signatures classify cell type")
  ggsave("Accuracy_DIA_WV.pdf", width = 5, height = 5)
  
  ggplot(marker.acc.df, aes(x=`Number of Proteins`, y=Accuracy,
                           color = Signature, shape = Filter)) + 
    geom_point() + theme_minimal() + ylab("% Correctly Classified") + 
    facet_wrap(`Measured Cell Type` ~ `Cell Type`) +
    ggtitle("Cell-specific signatures classify cell type")
  ggsave("AccuracyPerCellType_DIA_WV.pdf", width = 5, height = 5)
  
  return(list(wv = wv.df, p = p.df, drug.corr = drug.corr.df, 
              frac.corr = frac.corr.df, 
              accuracy = acc.df, marker.accuracy = marker.acc.df))
}

load_not_norm_BeatAML_for_DMEA3 <- function(BeatAML.path = "BeatAML_DMEA_inputs_not_normalized",
                                            exclude.samples = c()) {
  message("Loading Beat AML data for DMEA")
  BeatAML_synapse_id <- list("drug_response.csv" = "syn51674470", 
                             "Ex10_metadata.txt" = "syn25807733",
                             "ptrc_ex10_crosstab_global_gene_original.txt" = "syn25714254",
                             "ptrc_ex10_crosstab_phospho_siteID_original.txt" = "syn25714936")
  
  ### download files if any not already downloaded
  if (!file.exists(BeatAML.path)) {
    lapply(BeatAML_synapse_id, synapser::synGet, downloadLocation = BeatAML.path)
  } else if (!any(FALSE %in% lapply(names(BeatAML_synapse_id), file.exists))) {
    lapply(BeatAML_synapse_id, synapser::synGet, downloadLocation = BeatAML.path)
  }
  
  ### load files
  drug.BeatAML <- read.csv(file.path(BeatAML.path, names(BeatAML_synapse_id)[1]))
  meta.BeatAML <- read.table(file.path(BeatAML.path, names(BeatAML_synapse_id)[2]), 
                             sep = "\t", header = TRUE)
  global.BeatAML <- read.table(file.path(BeatAML.path, names(BeatAML_synapse_id)[3]),
                               sep = "\t", header = TRUE)
  phospho.BeatAML <- read.table(file.path(BeatAML.path, names(BeatAML_synapse_id)[4]),
                                sep = "\t", header = TRUE)
  rna.BeatAML <- synapser::synTableQuery("select * from syn26545877")$asDataFrame()
  
  ### format BeatAML data for DMEA
  sample.names <- "Barcode.ID"
  
  ## format drug sensitivity data frame
  # format drug.BeatAML wide (samples in first column, drug names for rest of columns)
  drug.BeatAML <- reshape2::dcast(drug.BeatAML, sample_id ~ inhibitor, 
                                  value.var = "auc", fill = NA)
  
  # change sample column name to match expression data
  names(drug.BeatAML)[1] <- sample.names
  
  ## format global proteomics data frame
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
  
  # subtract sample medians
  sample.names <- colnames(dplyr::select_if(global.BeatAML, is.numeric))
  #global.BeatAML[,sample.names] <- log(global.BeatAML[,sample.names], 2)
  global_sample_coef <- apply(global.BeatAML[,sample.names], 2, median, na.rm = T)
  global.BeatAML[,sample.names] <- sweep(global.BeatAML[,sample.names], 2, global_sample_coef, FUN = '-')
  
  # transpose global.BeatAML so that first column is Barcode.ID and 
  # rest of columns are gene symbols
  global.BeatAML <- as.data.frame(t(global.BeatAML))
  
  # make first column Barcode.ID
  global.BeatAML[,"Barcode.ID"] <- rownames(global.BeatAML)
  global.BeatAML <- 
    global.BeatAML[ , c("Barcode.ID", 
                        names(global.BeatAML[ , 1:(ncol(global.BeatAML)-1)]))]
  
  ## format phospho-proteomics data frame
  # change global.BeatAML column names from SampleID.abbrev to Barcode.ID to match drug.BeatAML
  phospho.ids <- names(phospho.BeatAML)
  
  # remove X and any 0's from start of each column name and then
  # replace SampleID.abbrev with Barcode.ID to match drug.BeatAML
  for(i in seq_len(length(phospho.ids))){
    phospho.ids[i] <- substr(phospho.ids[i], 2, nchar(phospho.ids[i]))
    
    if(substring(phospho.ids[i], 1, 1) == 0){
      phospho.ids[i] <- substr(phospho.ids[i], 2, nchar(phospho.ids[i]))
    }
    
    if(phospho.ids[i] %in% meta.BeatAML$SampleID.abbrev){
      phospho.ids[i] <- meta.BeatAML[
        meta.BeatAML$SampleID.abbrev == phospho.ids[i], ]$Barcode.ID
    }
  }
  
  # replace phospho.BeatAML column names
  names(phospho.BeatAML) <- phospho.ids
  
  # subtract sample medians
  sample.names <- colnames(dplyr::select_if(phospho.BeatAML, is.numeric))
  #phospho.BeatAML[,sample.names] <- log(phospho.BeatAML[,sample.names], 2)
  phospho_sample_coef <- apply(phospho.BeatAML[,sample.names], 2, median, na.rm = T)
  phospho.BeatAML[,sample.names] <- sweep(phospho.BeatAML[,sample.names], 2, phospho_sample_coef, FUN = '-')
  
  # transpose phospho.BeatAML so that first column is Barcode.ID and rest of columns are gene symbols
  phospho.BeatAML <- as.data.frame(t(phospho.BeatAML))
  
  # make first column Barcode.ID
  phospho.BeatAML[, "Barcode.ID"] <- rownames(phospho.BeatAML)
  phospho.BeatAML <- phospho.BeatAML[ , c("Barcode.ID", names(phospho.BeatAML[ , 1:(ncol(phospho.BeatAML)-1)]))]
  
  ## format rnaSeq data frame
  rna.BeatAML$Barcode.ID <- rna.BeatAML$labId
  rna.BeatAML <- reshape2::dcast(rna.BeatAML, Barcode.ID ~ display_label, mean,
                                 value.var = "RNA counts")
  rownames(rna.BeatAML) <- rna.BeatAML$Barcode.ID
  
  return(list(meta = meta.BeatAML[!(meta.BeatAML$Barcode.ID %in% exclude.samples),], 
              drug = drug.BeatAML[!(drug.BeatAML$Barcode.ID %in% exclude.samples),], 
              rna = rna.BeatAML[!(rna.BeatAML$Barcode.ID %in% exclude.samples),],
              global = global.BeatAML[!(global.BeatAML$Barcode.ID %in% exclude.samples),],
              phospho = phospho.BeatAML[!(phospho.BeatAML$Barcode.ID %in% exclude.samples),]))
}
sorted.patients <- c("18-00105", "21-00839", "22-00571", "22-00117", "16-01184",
                     "19-00074", "18-00103", "21-00432", "17-01060", "22-00251")
BeatAML <- load_not_norm_BeatAML_for_DMEA3(exclude.samples = sorted.patients)
gmt.drug <- readRDS("gmt_BeatAML_drug_MOA_2024-02-22.rds")

# load sorted proteomics
synapser::synLogin()
global.df <- read.csv(synapser::synGet("syn58895933")$path) # DIA
global.df <- global.df[ , which(colMeans(!is.na(global.df)) >= 0.75)] # 36 out of 48 samples are kept
outliers <- c("X00839_CD34plusFlow", "X00117_CD34plus", 
              "X00432_CD14plus", "X00251_CD14plus", "X00105_CD14plusFlow")
# syn.test <- synapser::synGet("syn58914135") # for some reason, can't access pre-filtered version ???
rownames(global.df) <- global.df$Gene
global.df$Gene <- NULL
global.df <- as.data.frame(t(global.df))
global.df$Sample <- rownames(global.df)
global.df <- global.df[,c("Sample", colnames(global.df)[1:(ncol(global.df)-1)])]
global.df <- global.df[!(global.df$Sample %in% outliers),]
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
    evalResults <- evalMonoSig(global.df100, sigs[[j]], names(sigs)[j]) 
    
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