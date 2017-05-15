library(regioneR)
library(lattice)
library(colorRamps)

clustMatrix<-function(mat,...){
  d<-dist(mat)
  hc<-hclust(d,...)
  ind<-hc$order

  d1<-dist(t(mat))
  hc1<-hclust(d1,...)
  ind1<-hc1$order

  mat1<-mat[ind,ind1]
  return(mat1)
}
clustMatrix.double<-function(mat,...){

  ordx<-order(colnames(mat))
  ordy<-order(rownames(mat))
  mat<-mat[ordy,ordx]
  d<-dist(mat)
  hc<-hclust(d,...)
  ind<-hc$order

  mat1<-mat[ind,ind]
  return(mat1)
}
listOverlapPermTest.true<-function(Alist,Blist,ntimes,...){
  list.tabs<-list()
  list.pt<-list()
  for ( i in 1:length(Alist)){
    A<-Alist[[i]]
    new.names<-names(Blist)
    func.list <- createFunctionsList(FUN=numOverlaps, param.name="B", values=Blist)

    ptm <- proc.time()
    pt <- permTest(A=A, ntimes=ntimes, evaluate.function=func.list,
                   randomize.function=randomizeRegions,count.once=T,per.chromosome=T)
    time<-proc.time() - ptm
    time<-time[3]/60
    print(paste0(" run in ",time,"  minute"))
    tab<-vector()
    for (j in 1:length(pt)){
      if(pt[[j]]$zscore>0){
        zscore.norm<-pt[[j]]$zscore/sqrt(length(A))
        max.val<-min(length(A),length(Blist[[j]]))
        max.z<-(max.val-mean(pt[[j]]$permuted))/sd(pt[[j]]$permuted)
        max.z.norm<-max.z/sqrt(length(A))
        zscore.std<-zscore.norm/max.z.norm
      }
      if(pt[[j]]$zscore<=0){
        zscore.norm<-pt[[j]]$zscore/sqrt(length(A))
        max.val<-0
        max.z<-(max.val-mean(pt[[j]]$permuted))/sd(pt[[j]]$permuted)
        max.z.norm<-max.z/sqrt(length(A))
        zscore.std<-(-(zscore.norm/max.z.norm))
      }
      vec<-c(j,new.names[j],length(Blist[[j]]),pt[[j]]$zscore,zscore.norm,zscore.std,pt[[j]]$pval,pt[[j]]$observed,
             mean(pt[[j]]$permuted), sd(pt[[j]]$permuted))
      tab<-rbind(vec,tab)
    }
    colnames(tab)<-c("order.id","name","nº region","z-score","norm.z-score","standard.z-score","p-value","n.overlaps","mean.perm","sd.perm")
    print(tab)
    list.tabs[[i]]<-tab
    #list.pt[[i]]<-pt
    names(list.tabs)[i]<-names(Alist)[i]
    #names(list.pt)[i]<-names(Alist)[i]

  }
  return(list.tabs)
}
matListOver<-function(listOverlapPermTest.obj,zs.type="std"){
  A.obj<-listOverlapPermTest.obj
  mat<-vector()
  for (i in 1:length(A.obj)){
    if (zs.type=="std"){
      mat<-cbind(mat,as.numeric(A.obj[[i]][,6]))
    }
    if (zs.type=="norm"){
      mat<-cbind(mat,as.numeric(A.obj[[i]][,5]))
    }
  }
  colnames(mat)<-names(A.obj)
  rownames(mat)<-A.obj[[1]][,2]
  mat<-as.matrix(mat)
  return(mat)
}

listOverlapPermTest.true.circ<-function(Alist,Blist,ntimes,...){
  list.tabs<-list()
  list.pt<-list()
  for ( i in 1:length(Alist)){
    A<-Alist[[i]]
    A<-filterChromosomes(A)
    new.names<-names(Blist)
    func.list <- createFunctionsList(FUN=numOverlaps, param.name="B", values=Blist)

    ptm <- proc.time()
    pt <- permTest(A=A, ntimes=ntimes, evaluate.function=func.list,
                   randomize.function=circularRandomizeRegions,...)
    time<-proc.time() - ptm
    time<-time[3]/60
    print(paste0(" run in ",time,"  minute"))
    tab<-vector()
    for (j in 1:length(pt)){
      if(pt[[j]]$zscore>0){
        zscore.norm<-pt[[j]]$zscore/sqrt(length(A))
        max.val<-min(length(A),length(Blist[[j]]))
        max.z<-(max.val-mean(pt[[j]]$permuted))/sd(pt[[j]]$permuted)
        max.z.norm<-max.z/sqrt(length(A))
        zscore.std<-zscore.norm/max.z.norm
      }
      if(pt[[j]]$zscore<=0){
        zscore.norm<-pt[[j]]$zscore/sqrt(length(A))
        max.val<-0
        max.z<-(max.val-mean(pt[[j]]$permuted))/sd(pt[[j]]$permuted)
        max.z.norm<-max.z/sqrt(length(A))
        zscore.std<-(-(zscore.norm/max.z.norm))
      }
      vec<-c(j,new.names[j],length(Blist[[j]]),pt[[j]]$zscore,zscore.norm,zscore.std,pt[[j]]$pval,pt[[j]]$observed,
             mean(pt[[j]]$permuted), sd(pt[[j]]$permuted))
      tab<-rbind(vec,tab)
    }
    colnames(tab)<-c("order.id","name","nº region","z-score","norm.z-score","standard.z-score","p-value","n.overlaps","mean.perm","sd.perm")
    print(tab)
    list.tabs[[i]]<-tab
    #list.pt[[i]]<-pt
    names(list.tabs)[i]<-names(Alist)[i]
    #names(list.pt)[i]<-names(Alist)[i]

  }
  return(list.tabs)
}
ReduxAndCompareMatrix<-function(mat.A,mat.B,nameA="",nameB="",...){

  colnames(mat.A)<-names.A<-sub(".*-", "", colnames(mat.A))
  colnames(mat.B)<-names.B<-sub(".*-", "", colnames(mat.B))
  rownames(mat.A)<-sub(".*-", "",rownames(mat.A))
  rownames(mat.B)<-sub(".*-", "",rownames(mat.B))

  mat1<-mat.A[which(rownames(mat.A)%in%rownames(mat.B)),which(colnames(mat.A)%in%colnames(mat.B))]
  mat1<-mat1[!duplicated(rownames(mat1)),!duplicated(colnames(mat1))]
  mat2<-mat.B[which(rownames(mat.B)%in%rownames(mat.A)),which(colnames(mat.B)%in%colnames(mat.A))]
  mat2<-mat2[!duplicated(rownames(mat2)),!duplicated(colnames(mat2))]
  mat1<-mat1[order(rownames(mat1)),order(colnames(mat1))]
  mat2<-mat2[order(rownames(mat2)),order(colnames(mat2))]

  mat1<-clustMatrix(mat1)
  mat2<-mat2[match(rownames(mat1),rownames(mat2)),match(colnames(mat1),colnames(mat2))]
  mat3<-mat1-mat2

  mat.tot<-list(mat1,mat2,mat3)
  names(mat.tot)<-c(nameA,nameB,paste0(nameA,"-",nameB))
  return(mat.tot)

}
TabToListHMM<-function(tab,vect.test){
  lis<-list()
  names<-as.vector(unique(vect.test))
  for(i in 1:length(names)){
    lis[[i]]<-toGRanges(tab[which(vect.test==names[i]),])
  }
  names(lis)<-names

  for(i in 1:length(lis)){
    if (length(lis[[i]])>50000){
      lis[[i]]<-lis[[i]][sample(length(lis[[i]]),50000)]
    }
  }

  for(i in 1:length(lis)){
    names(lis)[i]<-paste0(i,"-",names(lis)[i])
  }
  return(lis)
}
listOverlapPermTest.resampling<-function(Alist,Blist,ntimes,universe="auto",mc.cores=4,...){
  list.tabs<-list()
  list.pt<-list()
  if (universe=="auto"){
    DFF<-data.frame()
    for(i in 1:length(Alist)){
      DF<-toDataframe(Alist[[i]])[,1:3]
      DFF<-rbind(DFF,DF)
    }
    universe<-DFF
  }
  for ( i in 1:length(Alist)){
    A<-Alist[[i]]
    new.names<-names(Blist)
    func.list <- createFunctionsList(FUN=numOverlaps, param.name="B", values=Blist)

    ptm <- proc.time()
    pt <- permTest(A=A, ntimes=ntimes, evaluate.function=func.list, count.once=T,
                   randomize.function=resampleRegions,universe=universe,per.chromosome=F,mc.cores=mc.cores)
    time<-proc.time() - ptm
    time<-time[3]/60
    print(paste0(" run in ",time,"  minute"))
    tab<-vector()
    for (j in 1:length(pt)){
      if(pt[[j]]$zscore>0){
        zscore.norm<-pt[[j]]$zscore/sqrt(length(A))
        max.val<-min(length(A),length(Blist[[j]]))
        max.z<-(max.val-mean(pt[[j]]$permuted))/sd(pt[[j]]$permuted)
        max.z.norm<-max.z/sqrt(length(A))
        zscore.std<-zscore.norm/max.z.norm
      }
      if(pt[[j]]$zscore<=0){
        zscore.norm<-pt[[j]]$zscore/sqrt(length(A))
        max.val<-0
        max.z<-(max.val-mean(pt[[j]]$permuted))/sd(pt[[j]]$permuted)
        max.z.norm<-max.z/sqrt(length(A))
        zscore.std<-(-(zscore.norm/max.z.norm))
      }
      vec<-c(j,new.names[j],length(Blist[[j]]),pt[[j]]$zscore,zscore.norm,zscore.std,pt[[j]]$pval,pt[[j]]$observed,
             mean(pt[[j]]$permuted), sd(pt[[j]]$permuted))
      tab<-rbind(vec,tab)
    }


    tab<-cbind(tab,p.adjust(as.numeric(tab[,7]),method="BH"))
    colnames(tab)<-c("order.id","name","nº region","z-score","norm.z-score","standard.z-score","p-value","n.overlaps","mean.perm","sd.perm","adj.p-value")
    print(tab)
    list.tabs[[i]]<-tab
    #list.pt[[i]]<-pt
    names(list.tabs)[i]<-names(Alist)[i]
    #names(list.pt)[i]<-names(Alist)[i]

  }
  return(list.tabs)
}
plotMatrix<-function(mat,palette="bwr",type="normal",name="",limit=TRUE,gr.range=c(-1.2,1.2)){

  atx<-seq(min(mat),max(mat),length.out=100)

  if(limit==TRUE){
    atx<-seq(gr.range[1],gr.range[2],length.out=100)
  } else {

  }
  if (palette=="bwr"){rgb.palette <- colorRampPalette(c("darkblue","blue","white","red", "darkred"),space = "rgb")}
  if (palette=="gbr"){rgb.palette <- colorRampPalette(c("green","black","red"),space = "rgb")}
  if (palette=="bwy"){rgb.palette <- colorRampPalette(c("black","white","yellow"),space = "rgb")
  }
  trellis.par.set(regions=list(col=rgb.palette(100)))
  if (type=="normal"){
    pl<-levelplot(mat,at=atx, scale=list(x=list(rot=45)),main=name, aspect="fill")
  }
  if (type=="compare"){
    pl<-levelplot(mat, scale=list(x=list(rot=45)),main=name, aspect="fill")
  }
  return(pl)
}
clustMatrix.double.2<-function(mat1,method="complete",...){
  mat<-mat1
  mat<-mat^1
  ordx<-order(colnames(mat))
  ordy<-order(rownames(mat))
  mat<-mat[ordy,ordx]
  mat1<-mat1[ordy,ordx]
  d<-dist(mat)
  hc<-hclust(d,method=method)
  ind<-hc$order

  mat1<-mat1[ind,ind]
  return(mat1)
}
selectCalls <- function(calls){
  for (i in 1:length(calls)){
    if (length(calls[[i]]) > 5000){
      ind <- sample(length(calls[[i]]), 5000)
      calls[[i]] <- calls[[i]][ind]
    }
  }
  return (calls)
}
globalAssociations <- function(call, GRanges, count.once = TRUE, lz = TRUE, ntimes = 100, mc.cores=8, window = 10000, step = 100, genome = "hg19", evaluate.func = "numOverlaps"){
  if (evaluate.func == "numOverlaps"){
    function.list <- createFunctionsList(FUN = numOverlaps, param.name = "B", values = GRanges)
    m_rmsk.pt <- permTest(A = call,
                          ntimes = ntimes,
                          evaluate.function = function.list,
                          randomize.function = randomizeRegions,
                          genome = genome,
                          per.chromosome = TRUE,
                          count.once = count.once,
                          mc.cores = mc.cores)
  }
  if (evaluate.func == "ATGC.eval"){
    m_rmsk.pt <- permTest(A = call,
                          ntimes = ntimes,
                          evaluate.function = ATGC.eval,
                          randomize.function = randomizeRegions,
                          genome = genome,
                          per.chromosome = TRUE,
                          count.once = count.once,
                          mc.cores = mc.cores)
  }
  if (lz == TRUE){
    lz_m_rmsk <- lzglobalAs(call, m_rmsk.pt, window = window, step = step, count.once=count.once)
    return (list(m_rmsk.pt,lz_m_rmsk))
  } else {
    return (m_rmsk.pt)
  }
}
parse.bed<-function(bed,chr,st=0,en=1e20){
  parsed<-bed[bed[,1]==chr,]
  parsed<-parsed[parsed[,2]>=st,]
  parsed<-parsed[parsed[,3]<=en,]
  return(parsed)
}
ATGC.eval<-function(A,type="AT",genome="hg19",...){
  A<-toDataframe(A)
  ATGCper<-data.frame()
  for(i in nrow(A)){
    if (genome=="hg19"){gs<-toString(getSeq(Hsapiens, as.character(A[i,1]), as.integer(A[i,2]), as.integer(A[i,3])))}
    if (genome=="mm9"){gs<-toString(getSeq(Mmusculus, as.character(A[i,1]), as.integer(A[i,2]), as.integer(A[i,3])))}
    ATGCper<-rbind(ATGCper,perBase(gs))
  }
  if (type=="GC"){return(as.numeric(colMeans(ATGCper)[1]))}
  if (type=="AT"){return(as.numeric(colMeans(ATGCper)[2]))}
}
lzglobalAs <- function(call, m_rmsk.pt, window = 10000, step = 100, count.once = TRUE){
  lz_m_rmsk <- list()
  for(i in 1:length(m_rmsk.pt)){
    lz_m_rmsk[[i]] <- localZScore(A = call,  pt = m_rmsk.pt[[i]], window = window, step = step, count.once = count.once)
    names(lz_m_rmsk)[i]<-names(m_rmsk.pt[i])
  }
  return(lz_m_rmsk)
}
matCreate <- function(call, lz_m_rmsk){
  lz_m_rmsk_shifted <- vector()
  lz.tabs <- lz_m_rmsk[2][[1]]
  for(i in 1:length(lz.tabs)){
    lz_m_rmsk_shifted <- rbind(lz_m_rmsk_shifted, lz.tabs[[i]]$shifted.z.scores/sqrt(length(call)))
  }
  #lz_m_rmsk_shifted[(lz_m_rmsk_shifted < cutoff) & (lz_m_rmsk_shifted > -cutoff)] <- 0
  rownames(lz_m_rmsk_shifted) <- names(lz.tabs)
  lz_m_rmsk_shifted <- as.matrix(lz_m_rmsk_shifted)
  return (lz_m_rmsk_shifted)
}
matCreate_matrix <- function(call, lz_m_rmsk){
  # argument is a matrix
  lz_m_rmsk_shifted <- vector()
  #lz.tabs <- lz_m_rmsk[2][[1]]
  for(i in 1:length(lz_m_rmsk)){
    lz_m_rmsk_shifted <- rbind(lz_m_rmsk_shifted, lz_m_rmsk[[i]]$shifted.z.scores/sqrt(length(call)))
  }
  rownames(lz_m_rmsk_shifted) <- names(lz_m_rmsk)
  lz_m_rmsk_shifted <- as.matrix(lz_m_rmsk_shifted)
  return (lz_m_rmsk_shifted)
}
clustMatrix.column<-function(mat,...){
  d<-dist(mat)
  hc<-hclust(d)
  ind<-hc$order

  mat1<-mat[ind,]
  return(list(mat1, ind))
}
getMatrixAssociation <- function(call, lz_m_rmsk, name = "", mat.plot = FALSE, ind = NULL, limit=TRUE, gr.range=c(-0.8,0.8)){
  lz_m_rmsk_shifted <- matCreate(call, lz_m_rmsk)
  if(is.null(ind)){
    list_lz_m_rmsk_shifted.HC_ind <- clustMatrix.column(lz_m_rmsk_shifted)
    lz_m_rmsk_shifted.HC <- list_lz_m_rmsk_shifted.HC_ind[[1]]
  } else {
    lz_m_rmsk_shifted.HC <- lz_m_rmsk_shifted[ind,]
  }
  names <- rownames(lz_m_rmsk_shifted.HC)
  if(mat.plot == TRUE){
    matplot(t(lz_m_rmsk_shifted.HC), type="l", col = primary.colors(16), main = name)
    #legend("topleft",legend=c(names), col = primary.colors(16), lty=1:16, cex=0.8)
  }
  plot.matrix <- plotMatrix(t(lz_m_rmsk_shifted.HC),palette="gbr",type="normal",name = name,limit=TRUE,gr.range=gr.range)
  if(is.null(ind)){
    return (list(list_lz_m_rmsk_shifted.HC_ind, lz_m_rmsk_shifted.HC, plot.matrix))
  } else {
    return (list(lz_m_rmsk_shifted.HC, plot.matrix))
  }
}
getMatrixAssociation_matrix <- function(call, lz_m_rmsk, name = "", mat.plot = FALSE, ind = NULL, limit=TRUE, gr.range=c(-0.8,0.8)){
  lz_m_rmsk_shifted <- matCreate_matrix(call, lz_m_rmsk)
  if(is.null(ind)){
    list_lz_m_rmsk_shifted.HC_ind <- clustMatrix.column(lz_m_rmsk_shifted)
    lz_m_rmsk_shifted.HC <- list_lz_m_rmsk_shifted.HC_ind[[1]]
  } else {
    lz_m_rmsk_shifted.HC <- lz_m_rmsk_shifted[ind,]
  }
  names <- rownames(lz_m_rmsk_shifted.HC)
  if(mat.plot == TRUE){
    matplot(t(lz_m_rmsk_shifted.HC), type="l", col = primary.colors(16), main = name)
    #legend("topleft",legend=c(names), col = primary.colors(16), lty=1:16, cex=0.8)
  }
  plot.matrix <- plotMatrix(t(lz_m_rmsk_shifted.HC),palette="gbr",type="normal",name = name,limit=TRUE,gr.range=gr.range)
  if(is.null(ind)){
    return (list(list_lz_m_rmsk_shifted.HC_ind, lz_m_rmsk_shifted.HC, plot.matrix))
  } else {
    return (list(lz_m_rmsk_shifted.HC, plot.matrix))
  }
}
changenames <- function(file){
  for(i in 1:length(file)){
    names(file)[i]<-strsplit(names(file)[i],"UniPk")[[1]][1]
    names(file)[i]<-paste0(i,"-",strsplit(names(file)[i],"Hepg2")[[1]][2])
  }
  return (file)
}
chromAutosomal <- function(list, organism = "hg"){
  for (i in 1:length(list)){
    list[[i]] <- filterChromosomes(list[[i]], organism = organism, chr.type = "autosomal")
  }
  return (list)
}
filterFunction <- function(listGRanges, vecmin = 500){
  vec <- vector()
  for(i in 1:length(listGRanges)){
    vec[i] <- length(listGRanges[[i]])
  }
  listGRanges <- listGRanges[vec > vecmin]
  return(listGRanges)
}
dirtInputCleaner<-function(test,test.q="0.995",input.call){
  call.test<-gscreen(paste0(test," > ",test.q),iterator=150)
  call.test.input<-extendRegions(input.call,extend.start=150,extend.end=150)
  ind<-(which(overlapRegions(call.test,call.test.input,only.boolean=T)))
  call.test<-toGRanges(call.test[-ind,])
  return(call.test)
}
multiple_merges <- function(Region){
  merges_Region <- Region[[1]]
  for(i in 2:length(Region) - 1){
    test <- mergeRegions(Region[[i]], Region[[i+1]])
    merges_Region <- mergeRegions(merges_Region,test)
  }
  return(merges_Region)
}
parseResults <- function(toGrep = "", call, mat, limit=TRUE, gr.range= c(-1,1), ind = NULL, mat.plot = TRUE){
  if (!is.null(toGrep)){
    pind <- grep(toGrep, names(mat))
    if(is.null(ind)){
      plot.mat <- getMatrixAssociation_matrix(call, mat[pind], name = deparse(substitute(mat)), mat.plot = TRUE, limit = limit, gr.range = gr.range)
    } else {
      plot.mat <- getMatrixAssociation_matrix(call, mat[pind], name = deparse(substitute(mat)), mat.plot = TRUE, ind = ind[[1]][[2]], limit = limit, gr.range= gr.range)
    }
  } else {
    if(is.null(ind)){
      plot.mat <- getMatrixAssociation_matrix(call, mat, name = deparse(substitute(mat)), mat.plot = TRUE, limit = limit, gr.range = gr.range)
    } else {
      plot.mat <- getMatrixAssociation_matrix(call, mat, name = deparse(substitute(mat)), mat.plot = TRUE, ind = ind[[1]][[2]], limit = limit, gr.range= gr.range)
    }
  }
  return (plot.mat)
}
