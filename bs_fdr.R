bs_fdr<-function(community.matrix,multiple.adjust=c('BY','BSFD'),output.filename='edge_data',
  init.threshold=0.5,iters=1000,fdr=1,risk=0.05,cor.method=c('spearman','pearson','kendall','manhattan','euclidean',
  'canberra','bray','kulczynski','jaccard','gower','altGower','morisita','horn','mountford','raup','binomial','chao',
  'cao','mahalanobis')) 
{
  options(warn=-1,scipen=12)
  if(is.null(community.matrix)==TRUE) stop(c('Must provide community abundance matrix.'))
  y<-community.matrix
  if(is.null(cor.method)==FALSE) cor.method<-match.arg(cor.method)
  else cor.method<-'spearman'
  base.methods<-c('spearman','pearson','kendall')
  if(cor.method %in% base.methods) x<-as.dist(cor(y,method=cor.method),diag=FALSE,upper=FALSE)
  else {
    require(vegan)
    x<-vegdist(t(y),method=cor.method,diag=FALSE,upper=FALSE,binary=FALSE)
  }
  if(is.null(multiple.adjust)==TRUE) stop(c('Must provide an adjustment method for multiple tests, either "BY" or "BSFD".'))
  multiple.adjust<-match.arg(multiple.adjust)
  #
  temp<-t(combn(colnames(y),2))
  #
  #run code in 10,000-row blocks
  block<-floor(length(x)/10000)
  edge_data<-data.frame('to'=character(),'from'=character(),'cor'=double(),stringsAsFactors=FALSE)
   write.table(edge_data,output.filename,sep='\t',row.names=FALSE,quote=FALSE)
  all_p_values<-double()
  start<-0
  #
  #This builds the data frame in 10,000 row blocks, and represents the first selection
  #of significant values
  for(n in 1:block) {
    end<-n*10000
    cat(paste('\r','Initial edge selection. Processing edges',end,'of',length(x),sep=' '))
    flush.console()
    window<-seq(from=(start+1),to=end)
    otu_names<-temp[window,]
    cor<-x[window]
    if(multiple.adjust=='BY') {
      t<-abs(cor*sqrt((nrow(y)-2)/(1-cor^2))) 
      p_value<-2*pt(t,(nrow(y)-1),lower.tail=FALSE)
      keep<-which(p_value<0.1)
      all_p_values<-c(all_p_values,p_value[keep])
     }
    else if(multiple.adjust=='BSFD') {
      cor<-c(cor)
      keep<-which(abs(cor)>init.threshold & is.na(cor)==FALSE)
    }
    edge_data<-data.frame('to'=otu_names[keep,1],'from'=otu_names[keep,2],'cor'=cor[keep],stringsAsFactors=FALSE)
    write.table(edge_data,output.filename,append=TRUE,sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
    start<-start+10000
  }
  #
   #And get the remaining edges
  cat(paste('\r','Initial edge selection. Processing edges',length(x),'of',length(x),'\r',sep=' '))
  flush.console()
  remaining<-abs((block*10000)-length(x))
  last.window<-seq(from=(length(x)-remaining),to=length(x))
  #otu_names<-unlist(strsplit(temp[last.window],' '))
  #otu_names<-matrix(otu_names,ncol=2,byrow=TRUE)
  otu_names<-temp[window,]
  cor<-x[last.window]
  if(multiple.adjust=='BY') {
    t<-abs(cor*sqrt((nrow(y)-2)/(1-cor^2)))
    p_value<-2*pt(t,(nrow(y)-1),lower.tail=FALSE)
    keep<-which(p_value<0.1)
    all_p_values<-c(all_p_values,p_value[keep])
  }
  else if(multiple.adjust=='BSFD') {
    cor<-c(cor)
    keep<-which(abs(cor)>init.threshold & is.na(cor)==FALSE)
  }
  edge_data<-data.frame('to'=otu_names[keep,1],'from'=otu_names[keep,2],'cor'=cor[keep],stringsAsFactors=FALSE)
  write.table(edge_data,output.filename,append=TRUE,sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
  #
  #If using BY method, adjust significance of p-values based on all rows from the first selection
  #This will vastly cut down on the network size
  #
  edge_data<-read.delim(output.filename,header=TRUE)
  #
  #Benjamini-Yekutieli method
  if(multiple.adjust=='BY') {
    cat(paste('\n','Adjustment of p-values using the Benjamini & Yekutieli method...\r',sep=' '))
    flush.console()
    all_p_values<-p.adjust(all_p_values,method=multiple.adjust)
    keep<-which(p_value<=0.05)
    edge_data<-edge_data[keep,]
    cat(paste('\n','Kept',length(keep),'significant edges out of',length(all_p_values),'from first selection.\r\n',sep=' '))
  }
  #Bootstrap False Discovery method	
  else if(multiple.adjust=='BSFD') {
    cat(paste('\n','Adjustment of association scores using the bootstrap false discovery method...\r\n',sep=' '))
    flush.console()
    #first subset community matrix following first selection
    edges_to<-as.character(unique(edge_data$to))
    edges_from<-as.character(unique(edge_data$from))
    unique_nodes<-unique(c(edges_to,edges_from))
    kept<-nrow(edge_data)
    rm(edges_to,edges_from)
    y<-y[,unique_nodes]
    if(cor.method %in% base.methods) x<-as.dist(cor(y,method=cor.method),diag=FALSE,upper=FALSE)
    else x<-vegdist(t(y),method=cor.method,diag=FALSE,upper=FALSE,binary=FALSE)
    x<-abs(c(x))
    x[x<0.5]<-NA
    adjust<-double()
    #resampling to get adjustment distribution
    for(i in 1:iters) {
      cat(paste('\r','Iteration',i,'of',iters,sep=' '))
      flush.console()
      y2<-y[sample(1:nrow(y),replace=TRUE),]
      if(cor.method %in% base.methods) x2<-as.dist(cor(y2,method=cor.method),diag=FALSE,upper=FALSE)
      else x2<-vegdist(t(y2),method=cor.method,diag=FALSE,upper=FALSE,binary=FALSE)
      x2<-abs(c(x2))
      x.temp<-x[which(x2<1 & x<1)]
      x2<-x2[which(x2<1 & x<1)]
      diff<-x2-x.temp
      false.pos<-which(x2>max(x.temp)) #may need to be which(x2>x)
      diff2<-diff[false.pos]
      diff2<-diff2[which(is.na(diff2)==FALSE)]
      #cutoff<-quantile(diff2,(1-fdr/length(diff2)),names=FALSE,na.rm=TRUE)
      cutoff<-sort(diff2,decreasing=TRUE)[fdr+1]
      adj<-min(diff[diff>=cutoff])
      adjust<-c(adjust,adj)
    }
    rm(y2,x2,x.temp,diff,false.pos,diff2,unique_nodes)
    final.adj<-quantile(adjust,(1-risk),names=FALSE,na.rm=TRUE)
    final.threshold<-init.threshold+final.adj
    keep<-which(abs(edge_data$cor)>final.threshold)
    edge_data<-edge_data[keep,]
    cat(paste('\n','Kept',length(keep),'significant edges out of',kept,'from first selection\n\r',sep=' '))
    cat(paste('\r','Adjusted threshold:',round(final.threshold,4),'at risk',risk,'( adjustment:',round(final.adj,4),')\n\r',sep=' '))
    cat(paste0('\r'))
    flush.console()
  }
  write.table(edge_data,output.filename,sep='\t',row.names=FALSE,quote=FALSE)
}
