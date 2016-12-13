library(ggplot2)
require(scales)
library(Rmisc)

comparison_data <- read.csv(file="/mnt/bioinf/archive0/eigen_thread/results/comparison_performance.csv", header=TRUE,  strip.white = TRUE, sep=",",na.strings= c("999", "NA", " ", ""))
comparison_data$superf<-NULL
comparison_data$family<-NULL
colnames(comparison_data)<-c("Top","Method","Class","Fold")
comparison_data$Top<-as.factor(comparison_data$Top)
comparison_data$Method<-as.factor(comparison_data$Method)

sf_comparison_data <- read.csv(file="/mnt/bioinf/archive0/eigen_thread/results/comparison_performance_superfamily_allowed.csv", header=TRUE,  strip.white = TRUE, sep=",",na.strings= c("999", "NA", " ", ""))
sf_comparison_data$superf<-NULL
sf_comparison_data$family<-NULL
colnames(sf_comparison_data)<-c("Top","Method","Class","Fold")
sf_comparison_data$Top<-as.factor(sf_comparison_data$Top)
sf_comparison_data$Method<-as.factor(sf_comparison_data$Method)

p1<-ggplot(comparison_data, aes(x=Top, y=Class, group=Method, fill=Method))+geom_bar(stat="identity", position="dodge")+ggtitle("Comparison of Class Recognition Performance")+ylab("TPR")+xlab("Top X")+theme(plot.title= element_text(hjust=0.5),text=element_text(size=25))
p2<-ggplot(comparison_data, aes(x=Top, y=Fold, group=Method, fill=Method))+geom_bar(stat="identity", position="dodge")+ggtitle("Superfamily and family hits excluded")+ylab("TPR")+xlab("Top X")+theme(plot.title= element_text(hjust=0.5),text=element_text(size=25))

p3<-ggplot(sf_comparison_data, aes(x=Top, y=Class, group=Method, fill=Method))+geom_bar(stat="identity", position="dodge")+ggtitle("Comparison of Class Recognition Performance")+ylab("TPR")+xlab("Top X")+theme(plot.title= element_text(hjust=0.5),text=element_text(size=25))
p4<-ggplot(sf_comparison_data, aes(x=Top, y=Fold, group=Method, fill=Method))+geom_bar(stat="identity", position="dodge")+ggtitle("Family hits excluded")+ylab("TPR")+xlab("Top X")+theme(plot.title= element_text(hjust=0.5),text=element_text(size=25))
multiplot(p4,p2,cols=2)