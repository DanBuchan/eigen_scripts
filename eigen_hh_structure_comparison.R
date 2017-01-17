library(ggplot2)
library(reshape2) 
library(Rmisc)
tm_comparison <- read.csv(file="/mnt/bioinf/archive0/eigen_thread/results/t1_tm_eigen_hh_comparison.csv", header=TRUE, check.names=FALSE, strip.white = TRUE, sep=",",na.strings= c("999", "NA", " ", ""))
scatter_data <- dcast(tm_comparison, formula = pdb~method, value.var="score", fill=0)

ggplot(scatter_data, aes(x=hh, y=eigen))+geom_point(shape=16, size=5)+geom_abline(colour="grey")+scale_y_continuous(limits=c(0,1))+scale_x_continuous(limits=c(0,1))+ylab("EigenThreader")+xlab("HHSearch")+theme(axis.text=element_text(size=20), )