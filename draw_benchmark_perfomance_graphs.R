library(ggplot2)
library(reshape2) 
vector_performance <- read.csv(file="/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/results/eigenvectors_performance.csv", header=TRUE, check.names=FALSE, strip.white = TRUE, sep=",",na.strings= c("999", "NA", " ", ""))
distance_performance <- read.csv(file="/cs/research/bioinf/home1/green/dbuchan/archive0/eigen_thread/results/distance_performance.csv", header=TRUE, check.names=FALSE, strip.white = TRUE, sep=",",na.strings= c("999", "NA", " ", ""))

vector_performance$top<-as.factor(vector_performance$top)

distance_performance$top<-as.factor(distance_performance$top)
colnames(distance_performance)<-c("top", "distance", "class", "fold", "superf", "family")
ggplot(vector_performance, aes(x=vectors, y=class, group=top, color=top))+geom_line(size=0.8)+theme_set(theme_gray(base_size = 18)) + ggtitle("Class")
ggplot(distance_performance, aes(x=distance, y=class, group=top, color=top))+geom_line(size=0.8)+theme_set(theme_gray(base_size = 18)) + ggtitle("Class")

ggplot(vector_performance, aes(x=vectors, y=fold, group=top, color=top))+geom_line(size=0.8)+theme_set(theme_gray(base_size = 18)) + ggtitle("Fold")
ggplot(distance_performance, aes(x=distance, y=fold, group=top, color=top))+geom_line(size=0.8)+theme_set(theme_gray(base_size = 18)) + ggtitle("Fold")

ggplot(vector_performance, aes(x=vectors, y=superf, group=top, color=top))+geom_line(size=0.8)+theme_set(theme_gray(base_size = 18)) + ggtitle("Superfamily")
ggplot(distance_performance, aes(x=distance, y=superf, group=top, color=top))+geom_line(size=0.8)+theme_set(theme_gray(base_size = 18)) + ggtitle("Superfamily")
