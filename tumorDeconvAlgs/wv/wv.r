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

# format input data
# input X: signature dataframe
# input Y: expression dataframe
# output X: signature matrix where rownames are gene names and all colnames are cell types
# output Y: expression matrix where rownames are sample names and all colnames are gene names in X
prep_for_wv <- function(X, Y) {
  if (nrow(Y) > ncol(Y)) { # checks if genes are along rows, assuming # of genes > # of samples
    # assumes gene symbols are in first column or rownames
    # to prevent crashing on duplicated gene symbols, add unique numbers to identical names
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
  }
  
  # make sure gene symbols are rownames for signature matrix
  if (!is.numeric(X[,1])) {
    rownames(X) <- X[,1]
    X <- X[,2:ncol(X)]
  }
  
  # only keep genes in signatures
  Y <- Y[,rownames(X)[rownames(X) %in% colnames(Y)]]
  
  # make sure genes are quantified in all samples even if all were in signature
  Y <- t(na.omit(t(Y)))
  
  while (nrow(X) != ncol(Y)) { # X & Y must have same # of genes
    # filter signature matrix for genes also in filtered expression matrix
    X <- X[rownames(X) %in% colnames(Y),]
    
    # filter expression matrix for genes also in filtered signature matrix
    Y <- Y[,rownames(X)]
  }
  
  return(list(X = as.matrix(X), Y = as.matrix(Y)))
}

#main function
wv <- function(sig_matrix, mixture_file){
  #read in data
  X <- read.csv(sig_matrix, sep = "\t", row.names = 1,check.names=F)
  Y <- read.csv(mixture_file, sep = "\t",check.names=F)
  
  new.inputs <- prep_for_wv(X, Y)
  X <- new.inputs$X
  Y <- new.inputs$Y
  
  # prepare result dataframe
  cell.types <- colnames(X)
  obj <- matrix(NA, nrow = ncol(X), ncol = nrow(Y), 
                dimnames = list(cell.types, rownames(Y)))
  
  if (nrow(X) > 0) {
    # run WV for each cell type signature
    for (i in cell.types) {
      message("Running signature for ", i, " cells")
      temp.sig <- na.omit(as.matrix(X[,i]))
      temp.expr <- as.matrix(Y[,rownames(temp.sig)])
      obj[i,] <- temp.expr %*% temp.sig
    } 
  } else {
    message("no genes in weights match in expression")
  }

  return(obj)
}

## test wv
# input type: either "prot" to test proteomics or "mrna" to test mRNA input
test_wv <- function(type="prot") {
  base.path <- "~/OneDrive - PNNL/Documents/GitHub/decomprolute"
  # load test data
  mixture_file <- file.path(base.path, "toy_data", 
                            paste0("ov-all-",type,"-reduced.tsv"))
  
  # load example signature matrix
  sig_matrix <- file.path(base.path,"signature_matrices/AML.txt")
  
  # run wv
  test.output <- wv(sig_matrix, mixture_file) 
  return(test.output)
}

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
      wv.result <- matrix(NA, nrow = length(colnames(Y)) - 1, ncol = length(colnames(X)), dimnames = list(colnames(Y)[2:length(colnames(Y))], colnames(X)))
      write.table(t(wv.result), file="deconvoluted.tsv", quote = FALSE, col.names = NA, sep = "\t")
    }
)
