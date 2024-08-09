#' merge results across matrices, algorithms and cancers
#' currently we have two assessments:
#' 1 - assess correlations for each patient/disease
#' 2- assess correlations for each celltype/matrix
#'
#' We generally only use correlation right now
#' to run:
#' Rscript  combine_results.R --metricType sample --metric correlation --repNumber 0 */*corr.tsv
#' Rscript  combine_results.R --metricType cellType --metric correlation  --repNumber 0  */*corrXcelltypes.tsv
#' Rscript  combine_results.R --metricType sample --metric meanCorrelation  --repNumber 0 */*corr.tsv
#' Rscript  combine_results.R --metricType cellType --metric meanCorrelation  --repNumber 0  */*corrXcelltypes.tsv
#' Rscript  combine_results.R --metricType js --metric distance  --repNumber 0  */*dist.tsv
library(plyr)
library(dplyr)
library(argparser)
library(ggplot2)
library(nationalparkcolors)
library(reshape2)
library(plotly)
library(htmlwidgets)
#library(webshot)
#webshot::install_phantomjs()
# pal<-c(park_palette('GeneralGrant'), park_palette('Redwoods'))

##here is the color scheme
pal <- unlist(park_palettes[c(7, 15, 25, 22, 19, 16, 13, 11, 12, 3, 1, 2, 4, 5, 6, 8, 9, 10, 14, 17, 18, 20, 21, 23, 24)], use.names = F)


###combineResultsFiles
## this is an important function that is used by all figure files
## combines all ersults into a single file
combineResultsFiles<-function(file.list,colnamevals=c('patient','correlation')){
  
  message(paste0('Combining ',length(file.list),' files'))
  
  ##first we read in each file and make int a single table
  full.tab<-do.call(rbind,lapply(file.list,function(file){
    vars <- unlist(strsplit(basename(file),split='-')) #split into pieces
    message("vars are ", paste0(vars, sep=" "))
    tissue=vars[1]
    disease=vars[2]
    signature1=vars[3]
    signature2=vars[5]
    #matrix=vars[6]
    #sample=vars[7]
    #rep = vars[9]
    tab<-read.table(file,fill=TRUE,sep = '\t',check.names=FALSE)
    
    ##here we have to do some tomfoolery
    if (ncol(tab) > 1) {
      colnames(tab)<-colnamevals
    } else {
      tab$values <- NaN
      colnames(tab)<-colnamevals
    }
    return(data.frame(tab,tissue,disease,signature1,signature2))
  }))
  full.tab$algorithm <- "wv"
  return(full.tab)
}

combineValFiles<-function(file.list,colnamevals=c('cellType')){
  
  message(paste0('Combining ',length(file.list),' files'))
  
  ##first we read in each file and make int a single table
  full.tab<-do.call(rbind,lapply(file.list,function(file){
    # extract metadata from filename
    vars <- unlist(strsplit(basename(file),split='-')) #split into pieces
    message("vars are ", paste0(vars, sep=" "))
    tissue=vars[2]
    disease=vars[1]
    signature=vars[6]
    
    # import data
    tab<-read.table(file,fill=TRUE,sep = '\t',check.names=FALSE)
    colnames(tab)[1]<-colnamevals[1]
    
    # convert data from wide to long format
    tab <- reshape2::melt(tab, id.vars=colnamevals[1], variable.name="sample",
                          value.name="wv")
    
    # extract more metadata
    disease.info <- strsplit(disease, split="_")[[1]]
    tab$cancerType <- disease.info[1] # "AML"
    tab$knownCellType <- disease.info[2] # e.g., "Monocyte"
    tab$method <- disease.info[3] # e.g., "DIA" (MS method)
    
    return(data.frame(tab,tissue,disease,signature))
  }))
  full.tab[full.tab$signature == "AML_sorted_100.tsv",]$signature <- "sorted proteomics"
  full.tab[full.tab$signature == "AML_vanGalen_100.tsv",]$signature <- "single-cell transcriptomics" # van Galen et al, 2019
  return(full.tab)
}

#' combine list of files by cell type values
combineCellTypeVals<-function(file.list){
  # import & format data
  full.tab<-combineValFiles(file.list,colnamevals=c('cellType'))
  print(head(full.tab))
  
  # plot data
  fc<-ggplot(full.tab,aes(x=cellType,y=wv,fill=knownCellType))+geom_bar(stat='identity')+scale_fill_manual(values=pal)+
    facet_grid(rows=vars(signature),cols=vars(method))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(paste0('cellTypeValueAllSamples.pdf'),fc,width=10)
  
  # calculate p-values (e.g., if monocyte score is higher than others for known monocytes)
  cellType <- sort(unique(full.tab$cellType))
  knownCellTypes <- sort(unique(full.tab$knownCellType))
  p.tab <- data.frame(cellType)
  p.types <- c("overall","DIA","TMT","van_Galen","sorted")
  p.tab[,p.types] <- NA
  for (i in 1:length(knownCellTypes)) {
    overall.tab <- full.tab[full.tab$knownCellType == knownCellTypes[i],]
    p.overall <- t.test(overall.tab[overall.tab$cellType == cellType[i],], 
                        overall.tab[overall.tab$cellType != cellType[i],])$p.value
    
    DIA.tab <- overall.tab[overall.tab$method == "DIA",]
    p.DIA <- t.test(DIA.tab[DIA.tab$cellType == cellType[i],], 
                        DIA.tab[DIA.tab$cellType != cellType[i],])$p.value
    
    TMT.tab <- overall.tab[overall.tab$method == "TMT",]
    p.TMT <- t.test(TMT.tab[TMT.tab$cellType == cellType[i],], 
                        TMT.tab[TMT.tab$cellType != cellType[i],])$p.value
    
    vg.tab <- overall.tab[overall.tab$signature == "single-cell transcriptomics",]
    p.vg <- t.test(vg.tab[vg.tab$cellType == cellType[i],], 
                    vg.tab[vg.tab$cellType != cellType[i],])$p.value
    
    sorted.tab <- overall.tab[overall.tab$signature == "sorted proteomics",]
    p.sorted <- t.test(sorted.tab[sorted.tab$cellType == cellType[i],], 
                       sorted.tab[sorted.tab$cellType != cellType[i],])$p.value
    
    p.tab[p.tab$cellType == cellType[i],p.types] <- c(p.overall, p.DIA, p.TMT, p.vg, p.sorted)
  }
  write.table(p.tab,'pValues_knownCellType.tsv',row.names=F,col.names=T)
  
  # identify samples correctly classified by each method
  methods <- unique(full.tab$method)
  signatures <- unique(full.tab$signature)
  full.tab$predictedCellType <- NA
  full.tab$correctPrediction <- NA
  for (i in methods) {
    method.tab <- full.tab[full.tab$method == i,]
    samples <- unique(method.tab$sample)
    for (j in samples) {
      sample.tab <- method.tab[method.tab$sample == j,]
      for (k in signatures) {
        sig.tab <- sample.tab[sample.tab$signature == k,]
        
        # identify prediction
        prediction <- sig.tab[sig.tab$wv == max(sig.tab$wv),]$cellType
        full.tab[full.tab$method == i & 
                   full.tab$sample == j & 
                   full.tab$signature == k,]$predictedCellType <- prediction
        
        # determine if it was accurate
        acc <- which(knownCellTypes == prediction) == which(cellType == prediction) # assuming order of cell types match
        full.tab[full.tab$method == i & 
                     full.tab$sample == j & 
                     full.tab$signature == k,]$correctPrediction <- acc
      }
    }
  }
  write.table(p.tab,'accuracy_knownCellType.tsv',row.names=F,col.names=T)
  
  # create ternary plot (corners are Monocyte, Progenitor, MSC)
  tern.df <- reshape2::dcast(full.tab, sample + method + knownCellType + signature ~ cellType, value.var="wv")
  axis <- function(title) {
    list(
      title = title,
      titlefont = list(
        size = 20
      ),
      tickfont = list(
        size = 15
      ),
      tickcolor = 'rgba(0,0,0,0)',
      ticklen = 5
    )
  }
  
  tern.fig <- tern.df %>% plot_ly() %>% add_trace(
    type = 'scatterternary', mode = 'markers',
    a = ~Monocytle-like,
    b = ~Progenitor-like,
    c = ~MSC-like,
    marker = list(
      #colorscale = pal,
      color = ~knownCellType,
      shape = ~signature,
      symbol = 100,
      size = 14,
      line = list('width' = 2),
      showscale = TRUE
    )
  ) %>% layout(
    title = "",
    ternary = list(
      sum = 100,
      aaxis = axis('Monocyte-like'),
      baxis = axis('Progenitor-like'),
      caxis = axis('MSC-like')
    )
  )
  temp.fname <- "cellTypeTernary"
  saveWidget(tern.fig, paste0(temp.fname,".html"))
  #webshot(paste0(temp.fname,".html"), paste0(temp.fname,".pdf"))
  
  for (i in methods) {
    filtered.df <- tern.df[tern.df$method==i,]
    tern.fig <- filtered.df %>% plot_ly() %>% add_trace(
        type = 'scatterternary', mode = 'markers',
        a = ~Monocytle-like,
        b = ~Progenitor-like,
        c = ~MSC-like,
        marker = list(
          #colorscale = pal,
          color = ~knownCellType,
          shape = ~method,
          symbol = 100,
          size = 14,
          line = list('width' = 2),
          showscale = TRUE
        )
      ) %>% layout(
        title = i,
        ternary = list(
          sum = 100,
          aaxis = axis('Monocyte-like'),
          baxis = axis('Progenitor-like'),
          caxis = axis('MSC-like')
        )
      )
    temp.fname <- paste0("cellTypeTernary_",i)
    saveWidget(tern.fig, paste0(temp.fname,".html"))
    #webshot(paste0(temp.fname,".html"), paste0(temp.fname,".pdf")) 
  }
  
  # count % of samples correctly guessed for each signature
  agg.tab1 <- plyr::ddply(full.tab, .(signature), summarize,
                         N_correct = sum(correctPrediction, na.rm = TRUE),
                         N_total = length(na.omit(correctPrediction)))
  agg.tab1$percent_correct <- 100*agg.tab1$N_correct/agg.tab1$N_total
  
  perc.agg1 <- ggplot2::ggplot(agg.tab1, aes(x=signature, y=percent_correct))+
    geom_bar()+scale_fill_manual(values=pal)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave('cellTypeAccuracy_by_signature.pdf',perc.agg1,width=10)
  
  # count % of samples correctly guessed for each signature, method
  agg.tab2 <- plyr::ddply(full.tab, .(signature, method), summarize,
                          N_correct = sum(correctPrediction, na.rm = TRUE),
                          N_total = length(na.omit(correctPrediction)))
  agg.tab2$percent_correct <- 100*agg.tab2$N_correct/agg.tab2$N_total
  
  perc.agg2 <- ggplot2::ggplot(agg.tab2, aes(x=signature, y=percent_correct))+
    geom_bar()+scale_fill_manual(values=pal)+
    facet_grid(rows=vars(method))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave('cellTypeAccuracy_by_signature-method.pdf',perc.agg2,width=10)
  
  # also count % of correctly categorized samples for each cell type
  agg.tab3 <- plyr::ddply(full.tab, .(signature, method, knownCellType), summarize,
                          N_correct = sum(correctPrediction, na.rm = TRUE),
                          N_total = length(na.omit(correctPrediction)))
  agg.tab3$percent_correct <- 100*agg.tab3$N_correct/agg.tab3$N_total
  
  perc.agg3 <- ggplot2::ggplot(agg.tab3, aes(x=signature, y=percent_correct))+
    geom_bar()+scale_fill_manual(values=pal)+
    facet_grid(rows=vars(method),cols=vars(knownCellType))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave('cellTypeAccuracy_by_signature-method-knownCellType.pdf',perc.agg3,width=10)
  
  return(full.tab)
}

#' combine list of files of patient correlations
#' patient correlations represent how similar the cell type distribution is on a per-patient value, by cancer
combinePatientCors<-function(file.list,metric='correlation'){

  full.tab<-combineResultsFiles(file.list,colnamevals=c('patient',metric))|>
	dplyr::rename(value=metric)
  mats<-unique(full.tab$matrix)

  lapply(mats,function(mat){
    p<-full.tab%>%
      subset(matrix==mat)%>%
      ggplot()+
      geom_violin(aes(x=tissue,y=value,fill=disease))+
      facet_grid(rows=vars(signature1),cols=vars(signature2))+
      scale_fill_manual(values=pal)+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    ggsave(paste0(mat,'patient',metric,'s.pdf'),p,width=12,height=12)
  })

  p2<-full.tab%>%
    ggplot(aes(x=matrix,y=value,fill=disease))+geom_violin()+
    facet_grid(rows=vars(signature1),cols=vars(signature2))+scale_fill_manual(values=pal)
  ggsave(paste0('allSigsPatient',metric,'.pdf'),p2,width=12,height=12)

  mean.tab<-full.tab%>%
    group_by(tissue,disease,signature1,signature2)%>%
    summarize(meanVal=mean(value,na.rm=T))

  p3<-ggplot(mean.tab,aes(x=matrix,shape=tissue,y=meanVal,col=disease))+geom_jitter()+
    facet_grid(rows=vars(mrna.algorithm),cols=vars(prot.algorithm))+scale_color_manual(values=pal)
  ggsave(paste0('patient',metric,'averages.pdf'),p2,width=12,height=12)

  return(full.tab)
}

#' combine list of files by cell type correlations
combineCellTypeCors<-function(file.list,metric='correlation'){
  
  full.tab<-combineResultsFiles(file.list,colnamevals=c('cellType',metric))|>
    dplyr::rename(value=metric)
  print(head(full.tab))

  fc<-ggplot(full.tab,aes(x=cellType,y=value,fill=disease))+geom_boxplot()+scale_fill_manual(values=pal)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggsave(paste0('cellType',metric,'AllSamples.pdf'),fc,width=10)
  
  return(full.tab)
}

#' combine list of files of patient mean correlations
combineCorsMean<-function(file.list,metric='meanCorrelation', metricType='patient'){
  
  full.tab<-combineResultsFiles(file.list,colnamevals = c(metricType,'correlation'))
  
  full.tab<-full.tab|>
    dplyr::rename(newMet=metric)|>
    group_by(tissue,disease,mrna.algorithm,prot.algorithm,matrix,sample)|>
    summarize(meanCorr=mean(newMet))
  
  # message(paste0('Combining ',length(file.list),' files'))
  # 
  # full.tab<-do.call(rbind,lapply(file.list,function(file){
  #   vars <- unlist(strsplit(basename(file),split='-')) #split into pieces
  #   tissue=vars[1]
  #   disease=vars[2]
  #   mrna.algorithm=vars[3]
  #   prot.algorithm=vars[5]
  #   matrix=vars[6]
  #   sample=vars[7]
  #   tab<-read.table(file,fill=TRUE, sep = '\t',check.names=FALSE)
  #   if (ncol(tab) > 1) {
  #     colnames(tab)<-(c(metricType,metric))
  #     meanCorr <- mean(tab[[metric]])
  #   } else {
  #     meanCorr <- NaN
  #   }
  #   return(data.frame(meanCorr,tissue,disease,mrna.algorithm,prot.algorithm,matrix,sample))
  # }))
  # full.tab<-full.tab%>%
  #   mutate(algorithm=paste(mrna.algorithm,prot.algorithm,sep='-'))
  # 
 
  p2<-full.tab%>%
    ggplot(aes(x = mrna.algorithm, y = prot.algorithm, fill = meanCorr)) + geom_tile() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    # ggplot(aes(x=matrix,y=value,fill=disease))+geom_violin()+
    facet_grid(rows=vars(matrix),cols=vars(disease))+scale_fill_viridis_c()
  ggsave(paste0('heatmaps-', metricType, '-', metric,'.pdf'),p2)

#  p3<-full.tab%>%
#    ggplot(aes(x=matrix,shape=tissue,y=meanCorr,col=disease))+geom_jitter()+
#    facet_grid(rows=vars(mrna.algorithm),cols=vars(prot.algorithm))+scale_color_viridis_d() + ylab(metric)
#  ggsave(paste0('scatters-', metricType, '-', metric,'.pdf'),p3)

  p4<-full.tab%>%
    ggplot(aes(x=disease,y=meanCorr,fill=matrix))+geom_bar(stat='identity',position='dodge')+
    facet_grid(cols=vars(prot.algorithm),rows=vars(mrna.algorithm))+scale_fill_viridis_d()
  ggsave(paste0('barplot-matrix-', metricType, '-', metric,'.pdf'),p4,width=20)

  p5<-full.tab%>%
    ggplot(aes(x=matrix,y=meanCorr,fill=disease))+geom_bar(stat='identity',position='dodge')+
    facet_grid(cols=vars(prot.algorithm),rows=vars(mrna.algorithm))+scale_fill_viridis_d()
  ggsave(paste0('barplot-disease-', metricType, '-', metric,'.pdf'),p5,width=20)

  return(full.tab)
}

#' combine list of files of patient distances
combineDists<-function(file.list,metric='distance', metricType='js'){
 
  full.tab<-combineResultsFiles(file.list,colnamevals = c(metricType,'correlation'))
  full.tab<-full.tab|>
    dplyr::rename(newMet=metric)|>
    group_by(tissue,disease,mrna.algorithm,prot.algorithm,matrix,sample)|>
    summarize(distance=mean(newMet))
  
  # message(paste0('Combining ',length(file.list),' files'))
  # 
  # full.tab<-do.call(rbind,lapply(file.list,function(file){
  #   vars <- unlist(strsplit(basename(file),split='-')) #split into pieces
  #   tissue=vars[1]
  #   disease=vars[2]
  #   mrna.algorithm=vars[3]
  #   prot.algorithm=vars[5]
  #   matrix=vars[6]
  #   sample=vars[7]
  #   tab <- NULL
  # 
  #   try(tab<-read.table(file,fill=TRUE,sep = '\t', check.names=FALSE))
  #   if(is.null(tab))
  #     return(NULL)
  #   if (ncol(tab) > 1) {
  #     colnames(tab)<-(c('patient',metric))
  #     distance <- mean(tab[[metric]])
  #   } else {
  #     distance <- NaN
  #   }
  #   return(data.frame(distance,tissue,disease,mrna.algorithm,prot.algorithm,matrix,sample))
  # }))
  # full.tab<-full.tab%>%
  #   mutate(algorithm=paste(mrna.algorithm,prot.algorithm,sep='-'))

  p2<-full.tab%>%
    ggplot(aes(x = mrna.algorithm, y = prot.algorithm, fill = distance)) + geom_tile(height=1,width=1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    # ggplot(aes(x=matrix,y=value,fill=disease))+geo_violin()+
    facet_grid(rows=vars(matrix),cols=vars(disease))+scale_fill_gradient(low=pal[1],high=pal[3])
  ggsave(paste0('heatmaps-', metricType, '-', metric,'.pdf'),p2)

  p4<-full.tab%>%
    ggplot(aes(x=disease,y=distance,fill=matrix))+geom_bar(stat='identity',position='dodge')+
    facet_grid(cols=vars(prot.algorithm),rows=vars(mrna.algorithm))+scale_fill_manual(values=pal)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#+scale_fill_viridis_d()
  ggsave(paste0('barplot-matrix-', metricType, '-', metric,'.pdf'),p4)


  p5<-full.tab%>%
    ggplot(aes(x=matrix,y=distance,fill=as.factor(disease)))+geom_bar(stat='identity',position='dodge')+
    facet_grid(cols=vars(prot.algorithm),rows=vars(mrna.algorithm))+scale_fill_manual(values=pal)#+scale_fill_viridis_d()
  ggsave(paste0('barplot-disease-', metricType, '-', metric,'.pdf'),p5)


  p5<-full.tab%>%
    ggplot(aes(x=matrix,y=distance,fill=mrna.algorithm))+geom_boxplot()+
    facet_grid(cols=vars(prot.algorithm))+scale_fill_manual(values=pal)#+scale_fill_viridis_d()
  ggsave(paste0('boxplot-disease-', metricType, '-', metric,'mrna.pdf'),p5)

  p5<-full.tab%>%
    ggplot(aes(x=matrix,y=distance,fill=prot.algorithm))+geom_boxplot()+
    facet_grid(cols=vars(mrna.algorithm))+scale_fill_manual(values=pal)#+scale_fill_viridis_d()
  ggsave(paste0('boxplot-disease-', metricType, '-', metric,'prot.pdf'),p5)

  return(full.tab)
}




main<-function(){

  ##todo: store in synapse
  argv <- commandArgs(trailingOnly = TRUE)
  file.list<-argv[4:length(argv)]
  metric=argv[2]
  metricType = argv[1]
  repNumber = argv[3]
  if(metric=='correlation' && metricType=='sample'){
    tab<-combinePatientCors(file.list,metric)
    print(dim(tab))
    write.table(tab,paste0('combined-', metricType, '-', metric,'-',repNumber,'.tsv'),row.names=F,col.names=T)
  }else if(metric=='correlation' && metricType=='cellType'){
    tab<-combineCellTypeCors(file.list,metric)
    print(dim(tab))
    write.table(tab,paste0('combined-', metricType, '-', metric,'-',repNumber,'.tsv'),row.names=F,col.names=T)
  } else if (metric=='meanCorrelation'){
    tab<-combineCorsMean(file.list,metric, metricType)
    print(dim(tab))
    write.table(tab,paste0('combined-', metricType, '-', metric,'-',repNumber,'.tsv'),row.names=F,col.names=T)
  }else if(metric=='distance'){
    tab<-combineDists(file.list,metric,metricType)
    print(dim(tab))
    write.table(tab,paste0('combined-', metricType, '-', metric,'-',repNumber,'.tsv'),row.names=F,col.names=T)
  }else if(metric =='value'){
    tab <- combineCellTypeVals(file.list)
    print(dim(tab))
    write.table(tab,paste0('combined-', metricType, '-', metric,'-',repNumber,'.tsv'),row.names=F,col.names=T)
  }else{
    print("First argument must be metricType and second must be metric name")

  }

}

main()
