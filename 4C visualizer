  nred<-adjustcolor( "red", alpha.f = 0.12)
  norange<-adjustcolor( "orange", alpha.f = 0.3)
  ngreen<-adjustcolor( "green", alpha.f = 0.3)
  
  treatment<-"WT"
  len=20000
  spst<-30000
  pos1<-c(45193963,45193963+14000)
  range.win<-c(44e6,46e6)
  vec.wt<-vectorContacts(bed=bed.wt,g.start=g.start,g.end=g.end,pos=pos1,bin=bin,len=len)
  vec.ko<-vectorContacts(bed=bed.ko,g.start=g.start,g.end=g.end,pos=pos1,bin=bin,len=len)

  if(treatment=="WT"){
    name.file<-paste0("/imppc/labs/mblab/share/Cristina/Hi-C/images/",chr,":",pos1[1],"-",pos1[2],"_WT_pow.jpg")
  }
  if(treatment=="KO"){
    name.file<-paste0("/imppc/labs/mblab/share/Cristina/Hi-C/images/",chr,":",pos1[1],"-",pos1[2],"_KO_pow.jpg")
  }
  plotRegionOnProbe(probe.annotation.df,chr=chr,
                    start=range.win[1],end=range.win[2],
                    tab.genes=tab.genes,
                    ChrState=ChrState,adjust=-0.05,PoW=pos1)
  plotVec(vec.wt,cut=6,col="green",g.start=g.start,g.end=g.end)
  plotVec(vec.ko,cut=6,col="red",g.start=g.start,g.end=g.end)
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  segments(x0=pos1[1],x1=pos1[2],y0=0.05,y1=0.05,lwd=2)
  ins.plot<-ins.ko[(ins.ko[,1]>range.win[1] & ins.ko[,1]<range.win[2]),]
  ins.plot.2<-ins.wt[(ins.wt[,1]>range.win[1] & ins.wt[,1]<range.win[2]),]
  
  lines(ins.plot[,1],range01(ins.plot[,2])-1,col="brown")
  lines(ins.plot.2[,1],range01(ins.plot.2[,2])-1,col="green")
  
  norm.vec.wt<-vec.wt/mean(vec.wt)
  norm.vec.ko<-vec.ko/mean(vec.ko)
  diff.vec<-(norm.vec.wt-norm.vec.ko)/(((norm.vec.wt+norm.vec.ko)/2)+1)  

  m<-mean(diff.vec)
  s<-sd(diff.vec)
  col.vec<-rep("white",length(diff.vec))
  col.vec[diff.vec>=m+s*2.5]<-"green"
  col.vec[diff.vec<=m-s*2.5]<-"red"
  points((g.start*bin)+(1:length(diff.vec)*len),rep(1.05,length(diff.vec)),col = col.vec,pch=15,cex=1.2)
  vec.s<-vector()
  vec.e<-vector()
  if(treatment=="WT"){
    vec.s<-(bed.wt[,2]+bed.wt[,5])/sqrt(2)
    vec.e<-(bed.wt[,5]-bed.wt[,2])/sqrt(2)
    colhc<-nred
  }
  if(treatment=="KO"){
    vec.s<-(bed.ko[,2]+bed.ko[,5])/sqrt(2)
    vec.e<-(bed.ko[,5]-bed.ko[,2])/sqrt(2)
    colhc<-nred
  }
  vec.s<-sqrt((vec.s^2)/2)
  vec.e<-sqrt((vec.e^2)/2)
  
  #plot(g.start:g.end,main="Test draw.arc",xlim=c(min(vec.s),max(vec.s)),ylim=c(-max(abs(vec.e)),0),type="n")
  for (i in 1:length(vec.s)){
    
    start<-min(vec.s[i],vec.e[i])
    end<-max(vec.s[i],vec.e[i])
    lenx<-(end-start)/2
    x<-start+lenx
    
    cutter<-max(vec.e)/5
    points(x=vec.s[i],y=((-abs(vec.e[i]))/cutter)-1,col=colhc,pch=18,cex=1.1)
    #points(x=vec.s[i],y=((-abs(vec.e[i]))/cutter)-1,col=nred,pch=18,cex=1.1)
    #draw.arc(x, 1, len, deg1=1, deg2=180, col=nred, lwd=20, lend=1)
    
  }
  dev.off()
