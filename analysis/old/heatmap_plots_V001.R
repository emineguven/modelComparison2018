results.sub<-data.frame(cbind(sderr),cbind(R),cbind(G),cbind(Delta_LL))


R.els = unlist( unique(results.sub$R))
colnum = length(R.els)

tmp = unlist( unique(results.sub$sderr))
noise.els = tmp[order(tmp)]
rownum = length(noise.els)

mat = matrix( data=NA, nrow= rownum, ncol=colnum) #noise as row, alpha as columns
rownames(mat) = noise.els
colnames(mat) = R.els
for (k in c(0.05,0.08, 0.1,0.15,0.17, 0.2, 0.25)){
  data = results.sub[results.sub[,3]==k, 4]
  
  
  
  
  data<-unlist(data)
  
  heat_mat<-matrix(data,ncol=colnum,nrow=rownum)
  
  #rownames(heat_mat, do.NULL = TRUE, prefix = "row")
  rownames(heat_mat) <- c("0","0.5","1","2","3","4","5")
  
  colnames(heat_mat) <- R.els
  library(gplots)
  hM <- format(round(heat_mat, 2))
  data_mat<-scale(heat_mat,scale=TRUE,center=FALSE)
  
  
  #paste(file = "~/github/model.comparison/plots/heatplot_zero_noise_G",k,".jpeg",sep="") 
  jpeg(paste("plots/",k, ".batch.jpg", sep=''))
  
  #paste(“myplot_”, i, “.jpeg”, sep=””)
  
  heatmap.2(data_mat, cellnote=hM,col = cm.colors(256), scale="none", notecol="black",  margins=c(5,10),
            dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
            xlab     = "R parameters",
            ylab     = "noise", main = bquote(paste("R vs.noise of dLL fixed " ~ G==.(k))),srtCol=315, adjCol = c(0,1),cexRow=0.8,cexCol=0.8)
  
  
  dev.off()
  
}

