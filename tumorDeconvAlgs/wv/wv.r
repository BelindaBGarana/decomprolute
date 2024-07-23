#!/usr/local/bin/env Rscript --vanilla
args <- commandArgs(TRUE)
if (is.null(args[1])) {
  stop("No expression matrix provided!")
}

if (is.null(args[2])) {
  stop("No signature matrix provided!")
}

# wv R script (last updated 2024-07-22)
# Description: Weighted Voting
# Author: Belinda B. Garana, Pacific Northwest National Laboratory (belinda.garana@pnnl.gov)
# Requirements:
#       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#       install.packages('DMEA')
# Windows users using the R GUI may need to Run as Administrator to install or update packages.
#
# Usage:
#       Navigate to directory containing R script
#
#   In R:
#       source('wv.R')
#       results <- wv('sig_matrix_file.txt','mixture_file.txt')
#
# Input: signature matrix and mixture file, formatted as specified at https://github.com/BelindaBGarana/DMEA/blob/main/man/WV.Rd
# Output: dataframe containing all results
# License: CC0


#dependencies
library(DMEA)

#main function
wv <- function(sig_matrix, mixture_file){
  #read in data
  #X <- read.table(sig_matrix,header=T,sep="\t",row.names=1,check.names=F)
  X <- read.csv(sig_matrix, sep = "\t", row.names = 1,check.names=F)
  #Y <- read.table(mixture_file, header=T, sep="\t",check.names=F)
  Y <- read.csv(mixture_file, sep = "\t",check.names=F)
  
  #to prevent crashing on duplicated gene symbols, add unique numbers to identical names
  if (!is.numeric(Y[,1])) { # if first column contains gene symbols
    dups <- dim(Y)[1] - length(unique(Y[,1]))
    if(dups > 0) {
      warning(paste(dups," duplicated gene symbol(s) found in mixture file!",sep=""))
      rownames(Y) <- make.names(Y[,1], unique=TRUE)
    }else {rownames(Y) <- Y[,1]}  
    Y <- Y[,2:ncol(Y)] # remove column of gene symbols after setting rownames
  }
  
  # reformat Y so that gene symbols are along columns and samples along rows
  Y <- as.data.frame(t(Y))
  Y$Sample <- rownames(Y)
  Y <- Y[,c("Sample", colnames(Y)[1:(ncol(Y)-1)])]
  
  # make sure first column of sig_matrix contains gene symbols
  if (ncol(X) == ncol(dplyr::select_if(X, is.numeric))) {
    X$Gene <- rownames(X)
    X <- X[,c("Gene", colnames(X)[1:(ncol(X)-1)])]
  }
  
  # prepare result dataframe
  cell.types <- colnames(X)[2:ncol(X)]
  obj <- as.data.frame(cell.types)
  obj[,Y$Sample] <- NA
  
  for (i in seq_len(length(cell.types))) {
    # run WV for each cell type signature
    temp.sig <- X[,c(1,(i+1))]
    temp.WV <- DMEA::WV(Y, temp.sig)$scores
    
    # organize results to store in obj
    rownames(temp.WV) <- temp.WV[,1]
    temp.WV <- temp.WV[Y$Sample,]
    obj[i,2:ncol(obj)] <- temp.WV$WV
  }
  
  # format output
  rownames(obj) <- obj[,1]
  obj <- obj[,2:ncol(obj)]
  return(obj)
}

# ## test wv
# # load test data
# mixture_file <- "~/OneDrive - PNNL/Documents/GitHub/decomprolute/toy_data/ov-all-mrna-reduced.tsv"
# 
# # load example signature matrix
# sig_matrix <- "~/OneDrive - PNNL/Documents/GitHub/decomprolute/signature_matrices/AML.txt"
# 
# # run wv
# test.output <- wv(sig_matrix, mixture_file)

tryCatch(
    expr = {
        wv.result <- wv(args[2], args[1])
        write.table(wv.result, file="deconvoluted.tsv", quote = FALSE, col.names = NA, sep = "\t")
    },
    error = function(e){
      # (Optional)
      # Do this if an error is caught...
      print(e)
      X <- read.csv(args[2], sep = "\t", row.names = 1)
      Y <- read.csv(args[1], sep = "\t")
      wv <- matrix(0, nrow = length(colnames(Y)) - 1, ncol = length(colnames(X)), dimnames = list(colnames(Y)[2:length(colnames(Y))], colnames(X)))
      write.table(t(wv), file="deconvoluted.tsv", quote = FALSE, col.names = NA, sep = "\t")
    }
)