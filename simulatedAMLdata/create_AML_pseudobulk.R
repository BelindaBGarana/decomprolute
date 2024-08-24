# combine sorted AML proteomics signatures for deconvolution purposes
library(synapser); library(dplyr); library(reshape2); library(readxl)
setwd("~/OneDrive - PNNL/Documents/GitHub/decomprolute/signature_matrices")

#### 1. Load sorted proteomics ####
# load raw data
synapser::synLogin()
sorted.dia <- read.csv(synapser::synGet("syn58895933")$path)
sorted.tmt <- read.csv(synapser::synGet("syn53493077")$path)

# require proteins to be quantified in at least 75% of samples
sorted.dia <- sorted.dia[rowSums(is.na(sorted.dia)) < 0.25*ncol(sorted.dia),]
sorted.tmt <- sorted.tmt[rowSums(is.na(sorted.tmt)) < 0.25*ncol(sorted.tmt),]

# load metadata
meta <- readxl::read_excel(synapser::synGet("syn53461075"))
meta$id <- paste0("X", rownames(meta))
cd14.ids <- meta[grep("CD14", meta$SampleType), "id"]
cd34.ids <- meta[grep("CD34", meta$SampleType), "id"]
msc.ids <- meta[grep("MSC", meta$SampleType), "id"]

# #### 2. Create representative distributions for each cell type ####
# ### identify mean, stdev for each cell type
# 
# ### plot real vs. simulated distributions for each cell type

#### 2. Generate pseudobulk data ####
n.cells <- 1000
frac.dom <- 0.8 # fraction of cells which are the dominant cell type
cell.types <- c("Monocyte", "Progenitor", "MSC")
methods <- list("DIA" = sorted.dia, "TMT" = sorted.tmt)
for (j in names(methods)) {
  temp.df <- methods[j]
  pseudo.df <- data.frame(row.names=cell.types, col.names = rownames(temp.df))
  for (i in cell.types) {
    # separate out desired dominant cell type
    # cell.df <- t()
    # other.df <- t()
    
    ### sample 1000 times with replacement, 80% of cells being i, 20% being the rest
    # use dplyr::slice_sample with replacement option
    sampled.cell.df <- cell.df %>% slice_sample(n=frac.dom*n.cells, replace = TRUE)
    sampled.other.df <- other.df %>% slice_sample(n=(1-frac.dom)*n.cells, replace = TRUE)
    sampled.df <- rbind(sampled.cell.df, sampled.other.df)
    
    ### sum across the samples to create 1 pseudobulk profile
    #final.df <- data.frame(col.names = colnames(sampled.df), row.names = i)
    pseudo.df[i,] <- colSums(sampled.df)
  } 
  #pseudo.df <- data.table::rbindlist(all.data, use.names = TRUE, fill = NA, idcol = "dominant_cell_type")
  pseudo.df <- t(pseudo.df)
  
  ### median-normalize the data ?
  ## across samples
  global_row_medians <- apply(pseudo.df, 1, median, na.rm = T)
  pseudo.df <- sweep(pseudo.df, 1, global_row_medians, FUN = '-')
  
  ## across proteins
  global_sample_coef <- apply(pseudo.df, 2, median, na.rm = T)
  pseudo.df <- sweep(pseudo.df, 2, global_sample_coef, FUN = '-')
  
  ### save pseudobulk
  write.csv(pseudo.df, paste0(j, "_pseudobulk_global_proteomics.csv"), row.names = FALSE)
}
