library(MASS)
# y_n is the phenotype matrix, each row is a subject, each column represents a phenotype
HCDC<-function(y_n){
  sumcol=ncol(y_n)
  colnames(y_n)=1:sumcol
  lis=list()
  h=rep(NA,sumcol-1)
  comb=rep(NA,sumcol-1)
  comb_NULL=rep(NA,sumcol-1)
  for (i in 1:(sumcol-1)) {
    if(i==1) {
      lis[i]=list(matrix(NA,choose(sumcol,2),3))
      lis[[i]][,1:2]=t(combn(sumcol,2))
      lis[[i]][,3]=sapply(1:choose(sumcol,2), function(x) {
        (1-cor(y_n[,lis[[i]][x,1]],y_n[,lis[[i]][x,2]]))/(length(lis[[i]][x,1])*length(lis[[i]][x,2]))
      })
      lis_one=lis[[i]][which(lis[[i]][,3]==min(lis[[i]][,3])),]
      h[i]=lis_one[3]
      comb_NULL[i]=comb[i]=paste(lis_one[1],lis_one[2],sep='_')
    }
    else {
      lis[i]=list(matrix(NA,choose((sumcol-i+1),2),3))
      lis[[i]][,1:2]=t(combn(c(unique(comb_NULL[1:(i-1)]),(1:sumcol)[
        -as.numeric(unlist(strsplit(unique(comb_NULL[1:(i-1)]),split='_')))]),2))
      lis[[i]][,3]=sapply(1:choose(sumcol-i+1,2), function(x) {
        s1=unlist(strsplit(lis[[i]][x,1],split='_'))
        s2=unlist(strsplit(lis[[i]][x,2],split='_'))
        1-(cancor(y_n[,s1],y_n[,s2])$cor[1])
      })
      lis_two=lis[[i]][which(lis[[i]][,3]==min(lis[[i]][,3])),]
      h[i]=lis_two[3]
      comb_NULL[i]=comb[i]=paste(lis_two[1],lis_two[2],sep='_')
      index=sapply(1:(i-1),function(x) {
        (unlist(strsplit(comb[x],split = '_')) %in% unlist(strsplit(comb[i],split = '_')))[1]})
      comb_NULL[which(index)]=comb[i]
    }
  }
  h=as.numeric(h)
  D_value=do.call('c',sapply(1:(sumcol-1), function(x) {
    if (x>1) {h[x]-h[x-1]}
  }))
  lis_o=lis[[which(D_value==max(D_value))+1]]
  cluster_o=unique(rbind(matrix(lis_o[,1]),matrix(lis_o[,2])))
  clus=list()
  for (i in 1:nrow(cluster_o)) {
    if (length(grep("_",cluster_o[i,1])==1)) {clus[[i]]=as.numeric(unlist(strsplit(cluster_o[i,1],split = '_')))}
    else {clus[[i]]=as.numeric(cluster_o[i,1])}
  }
  return(list(n=nrow(cluster_o),cluster=clus))
}