library(misha)
library(devtools)
source("funct.need_assoc.R")
gsetroot("/imppc/labs/mblab/share/misha_roots/hg19_misha")

gvtrack.create("M1_HepG2.vir", "M1_HepG2", "global.percentile")
gvtrack.create("M2_HepG2.vir", "M2_HepG2", "global.percentile")
gvtrack.create("input_HepG2.vir", "input_HepG2", "global.percentile")

M1.Tracks.vir<-gextract("M1_HepG2.vir", gintervals(1:22),iterator=150)
M2.Tracks.vir<-gextract("M2_HepG2.vir", gintervals(1:22),iterator=150)
input.Tracks.vir<-gextract("input_HepG2.vir", gintervals(1:22),iterator=150)

png("/imppc/labs/mblab/share/Cristina/ABRIL/misha.0995.png", width = 800, height = 800)
plot(ecdf(-log2(1-M1.Tracks.vir[,4])),main="macroH2A")
lines(ecdf(-log2(1-M2.Tracks.vir[,4])), col="blue")
lines(ecdf(-log2(1-input.Tracks.vir[,4])), col="red")
abline(h=0.995, lty = 3)
dev.off()
