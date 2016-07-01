library(ggplot2)
require(scales)

comparison_data <- read.csv(file="/Users/dbuchan/Code/eigen_scripts/comparison_performance.csv", header=TRUE,  strip.white = TRUE, sep=",",na.strings= c("999", "NA", " ", ""))
comparison_data$superf<-NULL
comparison_data$family<-NULL
colnames(comparison_data)<-c("Top","Method","Class","Fold")
comparison_data$top<-as.factor(comparison_data$top)
comparison_data$Method<-as.factor(comparison_data$Method)

ggplot(comparison_data, aes(x=Top, y=Class, group=Method, colour=Method))+geom_line(size=0.8)+theme_set(theme_gray(base_size=18))+ggtitle("Comparison of Class Recognition Performance")
ggplot(comparison_data, aes(x=Top, y=Fold, group=Method, colour=Method))+geom_line(size=0.8)+theme_set(theme_gray(base_size=18))+ggtitle("Comparison of Fold Recognition Performance")