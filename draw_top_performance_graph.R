library(ggplot2)
library(reshape2)
library(Rmisc)
top_performance <- read.csv(file="/Users/dbuchan/Downloads/L_top_performance.csv", header=TRUE, check.names=FALSE, strip.white = TRUE, sep=",",na.strings= c("999", "NA", " ", ""))
random_performance <- read.csv(file="/Users/dbuchan/Downloads/L_random_performance.csv", header=TRUE, check.names=FALSE, strip.white = TRUE, sep=",",na.strings= c("999", "NA", " ", ""))

top_performance$top<-as.factor(top_performance$top)
random_performance$top<-as.factor(random_performance$top)

colnames(random_performance)<-c("top", "L", "class", "fold", "superf", "family")
colnames(top_performance)<-c("top", "L", "class", "fold", "superf", "family")

top_performance$L<-as.factor(top_performance$L)
random_performance$L<-as.factor(random_performance$L)

t1<-ggplot(top_performance, aes(x=L, y=fold, group=top, color=top))+geom_line(size=0.8)+theme_set(theme_gray(base_size = 18)) + ggtitle("Top")+ylab("TPR")+xlab("Contact Proportion (L/x)")
r1<-ggplot(random_performance, aes(x=L, y=fold, group=top, color=top))+geom_line(size=0.8)+theme_set(theme_gray(base_size = 18)) + ggtitle("Random")+ylab("TPR")+xlab("Contact Proportion (L/x)")
multiplot(t1,r1,cols=2)
