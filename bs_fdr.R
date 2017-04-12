bs_fdr<-function(community.matrix,init.threshold=0.5,iters=1000,fdr=1,risk=0.05,cor.method=c('spearman','pearson','kendall','manhattan','euclidean',
  'canberra','bray','kulczynski','jaccard','gower','altGower','morisita','horn','mountford','raup','binomial','chao',
  'cao','mahalanobis')) 
{
  options(warn=-1,scipen=12)
  if(is.null(community.matrix)==TRUE) stop(c('Must provide community abundance matrix.'))
  y<-community.matrix
  if(is.null(cor.method)==FALSE) cor.method<-match.arg(cor.method)
  else cor.method<-'spearman'
  base.methods<-c('spearman','pearson','kendall')
  if(cor.method %in% base.methods) x<-as.dist(cor(y,method=cor.method),diag=F,upper=F)
  else {
    require(vegan)
    x<-vegdist(t(y),method=cor.method,diag=F,upper=F,binary=F)
  }
  #
  temp<-t(combn(colnames(y),2))
  #
  #run code in 10,000-row blocks
  block<-floor(length(x)/10000)  
  edge_data<-data.frame('to'=character(),'from'=character(),'cor'=double(),stringsAsFactors=F)
  write.csv(edge_data,'edge_data.csv',row.names=F,quote=F)
  start<-0
  #
  #This builds the data frame in 10,000 row blocks, and represents the first selection of significant values
  message('Initial edge selection. Processing ',length(x),' edges...',appendLF=F)
  for(n in 1:block) {
    end<-n*10000
    flush.console()
    window<-seq(from=(start+1),to=end)
    otu_names<-temp[window,]
    cor<-x[window]
    cor<-c(cor)
    keep<-which(abs(cor)>init.threshold & is.na(cor)==F)
    edge_data<-data.frame('to'=otu_names[keep,1],'from'=otu_names[keep,2],'cor'=cor[keep],stringsAsFactors=F)
    write.csv(edge_data,'edge_data.csv',append=T,col.names=F,row.names=F,quote=F)
    start<-start+10000
  }
  #
  #And get the remaining edges
  remaining<-abs((block*10000)-length(x))
  last.window<-seq(from=(length(x)-remaining),to=length(x))
  otu_names<-temp[window,]
  cor<-x[last.window]
  cor<-c(cor)
  keep<-which(abs(cor)>init.threshold & is.na(cor)==F)
  edge_data<-data.frame('to'=otu_names[keep,1],'from'=otu_names[keep,2],'cor'=cor[keep],stringsAsFactors=F)
  write.csv(edge_data,'edge_data.csv',append=T,col.names=F,row.names=F,quote=F)
  message('complete')
  #
  edge_data<-read.csv('edge_data.csv',header=T)
  message('Creating bootstrapped association scores (',iters,' iterations)...',appendLF=F)
  #first subset community matrix following first selection
  edges_to<-as.character(unique(edge_data$to))
  edges_from<-as.character(unique(edge_data$from))
  unique_nodes<-unique(c(edges_to,edges_from))
  kept<-nrow(edge_data)
  rm(edges_to,edges_from)
  y<-y[,unique_nodes]
  if(cor.method %in% base.methods) x<-as.dist(cor(y,method=cor.method),diag=F,upper=F)
  else x<-vegdist(t(y),method=cor.method,diag=F,upper=F,binary=F)
  x<-abs(c(x))
  x[x<0.5]<-NA
  adjust<-double()
  #resampling to get adjustment distribution
  for(i in 1:iters) {
    y2<-y[sample(1:nrow(y),replace=T),]
    if(cor.method %in% base.methods) x2<-as.dist(cor(y2,method=cor.method),diag=F,upper=F)
    else x2<-vegdist(t(y2),method=cor.method,diag=F,upper=F,binary=F)
    x2<-abs(c(x2))
    x.temp<-x[which(x2<1 & x<1)]
    x2<-x2[which(x2<1 & x<1)]
    diff<-x2-x.temp
    false.pos<-which(x2>max(x.temp)) #may need to be which(x2>x)
    diff2<-diff[false.pos]
    diff2<-diff2[which(is.na(diff2)==F)]
    #cutoff<-quantile(diff2,(1-fdr/length(diff2)),names=FALSE,na.rm=TRUE)
    cutoff<-sort(diff2,decreasing=T)[fdr+1]
    adj<-min(diff[diff>=cutoff])
    adjust<-c(adjust,adj)
  }
  message('complete')
  rm(y2,x2,x.temp,diff,false.pos,diff2,unique_nodes)
  final.adj<-quantile(adjust,(1-risk),names=F,na.rm=T)
  final.threshold<-init.threshold+final.adj
  keep<-which(abs(edge_data$cor)>final.threshold)
  edge_data<-edge_data[keep,]
  message('Kept ',length(keep),' significant edges out of ',kept,' from first selection')
  message('Adjusted threshold: ',round(final.threshold,4),' at risk ',risk,' (adjustment: ',round(final.adj,4),')\n')
  file.remove('./edge_data.csv')
  return(edge_data)
}
