library(ggplot2)
library(reshape2) 
library(Rmisc)
vector_performance <- read.csv(file="/mnt/bioinf/archive0/eigen_thread/results/eigenvectors_performance.csv", header=TRUE, check.names=FALSE, strip.white = TRUE, sep=",",na.strings= c("999", "NA", " ", ""))
distance_performance <- read.csv(file="/mnt/bioinf/archive0/eigen_thread/results/distance_performance.csv", header=TRUE, check.names=FALSE, strip.white = TRUE, sep=",",na.strings= c("999", "NA", " ", ""))

vector_performance$top<-as.factor(vector_performance$top)

distance_performance$top<-as.factor(distance_performance$top)
colnames(distance_performance)<-c("top", "distance", "class", "fold", "superf", "family")
v1<-ggplot(vector_performance, aes(x=vectors, y=class, group=top, color=top))+geom_line(size=0.8)+theme_set(theme_gray(base_size = 25)) + ggtitle("Class")+ylab("TPR")+xlab("No. of Eigenvectors")+theme(plot.title= element_text(hjust=0.5))+scale_y_continuous(limits=c(0.25,1))
v2<-ggplot(vector_performance, aes(x=vectors, y=fold, group=top, color=top))+geom_line(size=0.8)+theme_set(theme_gray(base_size = 25)) + ggtitle("Fold")+ylab("TPR")+xlab("No. of Eigenvectors")+theme(plot.title= element_text(hjust=0.5))+scale_y_continuous(limits=c(0.25,1))
v3<-ggplot(vector_performance, aes(x=vectors, y=superf, group=top, color=top))+geom_line(size=0.8)+theme_set(theme_gray(base_size = 25)) + ggtitle("Superfamily")+ylab("TPR")+xlab("No. of Eigenvectors")+theme(plot.title= element_text(hjust=0.5))+scale_y_continuous(limits=c(0.25,1))
multiplot(v1,v3,v2,cols=2)

d1<-ggplot(distance_performance, aes(x=distance, y=class, group=top, color=top))+geom_line(size=0.8)+theme_set(theme_gray(base_size = 25)) + ggtitle("Class")+ylab("TPR")+xlab("Distance (angstroms)")+theme(plot.title= element_text(hjust=0.5))+scale_y_continuous(limits=c(0,1))
d2<-ggplot(distance_performance, aes(x=distance, y=fold, group=top, color=top))+geom_line(size=0.8)+theme_set(theme_gray(base_size = 25)) + ggtitle("Fold")+ylab("TPR")+xlab("Distance (angstroms)")+theme(plot.title= element_text(hjust=0.5))+scale_y_continuous(limits=c(0,1))
d3<-ggplot(distance_performance, aes(x=distance, y=superf, group=top, color=top))+geom_line(size=0.8)+theme_set(theme_gray(base_size = 25)) + ggtitle("Superfamily")+ylab("TPR")+xlab("Distance (angstroms)")+theme(plot.title= element_text(hjust=0.5))+scale_y_continuous(limits=c(0,1))
multiplot(d1,d3,d2,cols=2)

vect_class<-vector_performance
vect_class$fold<-NULL
vect_class$superf<-NULL
vect_class$family<-NULL
colnames(vect_class)[3]<-"value"
vect_class$label="Class"
head(vect_class)

vect_fold<-vector_performance
vect_fold$class<-NULL
vect_fold$superf<-NULL
vect_fold$family<-NULL
colnames(vect_fold)[3]<-"value"
vect_fold$label="Fold"
head(vect_fold)

vect_superf<-vector_performance
vect_superf$class<-NULL
vect_superf$fold<-NULL
vect_superf$family<-NULL
colnames(vect_superf)[3]<-"value"
vect_superf$label="Superfamily"
head(vect_superf)

whole_performance<-rbind(vect_class,vect_fold,vect_superf)
whole_performance

ggplot(whole_performance, aes(x=vectors, y=value, group=top, color=top))+geom_line(size=0.8)+geom_point(aes(shape=top), size=4)+theme_set(theme_gray(base_size = 15)) + ylab("TPR")+xlab("No. of Eigenvectors")+theme(strip.text.x=element_text(size=15), strip.background = element_rect(colour="white", fill="white"), axis.line = element_line(), panel.background=element_blank(), plot.title= element_text(hjust=0.5), legend.key.size=unit(2, "cm"))+scale_y_continuous(limits=c(0.25,1))+facet_wrap(~label, scales="free",ncol=2)
ggsave("/home/dbuchan/Dropbox/EigenThreader/proofed_figures/figure_2_eigenvector_performance.png", width=10, height=7, dpi=300)

dist_class<-distance_performance
dist_class$fold<-NULL
dist_class$superf<-NULL
dist_class$family<-NULL
colnames(dist_class)[3]<-"value"
dist_class$label="Class"
head(dist_class)

dist_fold<-distance_performance
dist_fold$class<-NULL
dist_fold$superf<-NULL
dist_fold$family<-NULL
colnames(dist_fold)[3]<-"value"
dist_fold$label="Fold"
head(dist_fold)

dist_superf<-distance_performance
dist_superf$class<-NULL
dist_superf$fold<-NULL
dist_superf$family<-NULL
colnames(dist_superf)[3]<-"value"
dist_superf$label="Superfamily"
head(dist_superf)

whole_distance<-rbind(dist_class,dist_fold,dist_superf)
whole_distance
ggplot(whole_distance, aes(x=distance, y=value, group=top, color=top))+geom_line(size=0.8)+geom_point(aes(shape=top), size=4)+theme_set(theme_gray(base_size = 15)) + ylab("TPR")+xlab("No. of Eigenvectors")+theme(strip.text.x=element_text(size=15), strip.background = element_rect(colour="white", fill="white"), axis.line = element_line(), panel.background=element_blank(), plot.title= element_text(hjust=0.5), legend.key.size=unit(2, "cm"))+scale_y_continuous(limits=c(0,1))+facet_wrap(~label, scales="free",ncol=2)
ggsave("/home/dbuchan/Dropbox/EigenThreader/proofed_figures/figure_3_distance_performance.png", width=10, height=7, dpi=300)
