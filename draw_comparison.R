library(ggplot2)
require(scales)
library(Rmisc)
# require(cowplot)

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

p1<-ggplot(comparison_data, aes(x=Top, y=Class, group=Method, fill=Method))+geom_bar(stat="identity", position="dodge")+scale_fill_brewer(palette='Blues')+scale_colour_brewer(palette='Blues')+ggtitle("Comparison of Class Recognition Performance")+ylab("TPR")+xlab("Top X")+theme(plot.title= element_text(hjust=0.5),text=element_text(size=25))
p2<-ggplot(comparison_data, aes(x=Top, y=Fold, group=Method, fill=Method))+geom_bar(stat="identity", position="dodge")+scale_fill_brewer(palette='Blues')+scale_colour_brewer(palette='Blues')+ylab("TPR")+xlab("Top X")+theme(plot.title= element_text(hjust=0.5),text=element_text(size=25))

p3<-ggplot(sf_comparison_data, aes(x=Top, y=Class, group=Method, fill=Method))+geom_bar(stat="identity", position="dodge")+scale_fill_brewer(palette='Blues')+scale_colour_brewer(palette='Blues')+ggtitle("Comparison of Class Recognition Performance")+ylab("TPR")+xlab("Top X")+theme(plot.title= element_text(hjust=0.5),text=element_text(size=25))
p4<-ggplot(sf_comparison_data, aes(x=Top, y=Fold, group=Method, fill=Method))+geom_bar(stat="identity", position="dodge", aes(colour = Method))+scale_fill_brewer(palette='Blues')+scale_colour_brewer(palette='Blues')+ylab("TPR")+xlab("Top X")+theme(legend.position="none", plot.title= element_text(hjust=0.5),text=element_text(size=25))
multiplot(p4,p2,cols=2)

# plot_grid(p4,p2,nrow=1,ncol=2, rel_widths= c(1/3,1/2), label_size = 20, labels=c("Family hits excluded","Superfamily and family hits excluded"))
comp_data_with_type <- comparison_data
comp_data_with_type$Label = "Superfamily & family hits excluded"
sf_data_with_type<- sf_comparison_data
sf_data_with_type$Label = "Family hits excluded"
whole_data<-rbind(sf_data_with_type,comp_data_with_type)
ggplot(whole_data, aes(x=Top, y=Fold, group=Method, fill=Method))+geom_bar(stat="identity", position="dodge", aes(colour = Method))+scale_fill_brewer(palette='Blues')+scale_colour_brewer(palette='Blues')+ylab("TPR")+xlab("Top X")+theme(strip.text.x=element_text(size=25), strip.background = element_rect(colour="white", fill="white"), axis.line = element_line(), panel.background=element_blank(), plot.title= element_text(hjust=0.5),text=element_text(size=25))+facet_wrap(~Label, scales="free")+scale_y_continuous(limits=c(0,1))
