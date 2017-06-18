insulationLevel<-function(bed,beta=10000,lambda=1e6,ins.steps,g.start,g.end){
  vec.point<-seq(from=g.start+lambda,to=g.end-lambda,by = ins.steps)
  new.coord<-cbind(rowMin(as.matrix(bed[,c(2,5)])),rowMax(as.matrix(bed[,c(2,5)])))
  vec.ins<-vector()
  for(i in  1:length(vec.point)){
    sel.square<-new.coord[new.coord[,1]<=(vec.point[i]-beta),]
    sel.square<-sel.square[sel.square[,2]>=(vec.point[i]+beta),]
    sel.square<-sel.square[sel.square[,1]>=(vec.point[i]-beta-lambda),]
    sel.square<-sel.square[sel.square[,2]<=(vec.point[i]+beta+lambda),]
    vec.ins[i]<-nrow(sel.square)
  }
  return(cbind(vec.point,vec.ins))
}
