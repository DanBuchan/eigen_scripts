library(ggplot2)
require(scales)

comparison_data <- read.csv(file="/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/results/comparison_performance.csv", header=TRUE,  strip.white = TRUE, sep=",",na.strings= c("999", "NA", " ", ""))
comparison_data$superf<-NULL
comparison_data$family<-NULL
colnames(comparison_data)<-c("Top","Method","Class","Fold")
comparison_data$Top<-as.factor(comparison_data$Top)
comparison_data$Method<-as.factor(comparison_data$Method)

ggplot(comparison_data, aes(x=Top, y=Class, group=Method, fill=Method))+geom_bar(stat="identity", position="dodge")+ggtitle("Comparison of Class Recognition Performance")+ylab("TPR")+xlab("Top X")
ggplot(comparison_data, aes(x=Top, y=Fold, group=Method, fill=Method))+geom_bar(stat="identity", position="dodge")+ggtitle("Comparison of Fold Recognition Performance")+ylab("TPR")+xlab("Top X")