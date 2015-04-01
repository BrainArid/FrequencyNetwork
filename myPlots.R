modStats <- read.table(file="../../freqNet0.3ModuleStats.txt",sep="\t",header=TRUE)

library("ggplot2");
library("reshape2")

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

d_formatter <- function(x) {
  print(as.character(x))
  x<-unlist(lapply(X=as.character(x),FUN=function(X){if(X=="0") return("0.00") else if (X=="0.5") return("0.50") else if(X=="1") return("1.00")else return(X)}));
  print(x)
  lab <- paste0('d = ',x);
}

mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="spearmanCut") {
    value[value=="0.1"] <- "c = 0.1"
    value[value=="0.2"] <- "c = 0.2"
    value[value=="0.3"] <- "c = 0.3"
    value[value=="0.35"] <- "c = 0.35"
    value[value=="0.4"] <- "c = 0.4"
  }
  return(value)
}

myTheme<-   theme(panel.margin=unit(x=0.1,units="in"),
                  axis.ticks.margin=unit(x=0.1,units="in"),
                  axis.text=element_text(colour="BLACK"),
                  axis.title=element_text(face="bold"),
                  title=element_text(face="bold"),
                  axis.text.x=element_text(angle=90));
plotBase <-  ggtitle("c (Spearman's Rank correlation threshold)")

  

filter<- modStats$E==2
p5 <- ggplot(modStats[filter,], aes(x=D, y=numModules)) + plotBase +
  facet_grid(. ~ spearmanCut, labeller=mf_labeller)+
  scale_x_continuous(label=d_formatter)+
  geom_point() +
  geom_line(weight=5, colour="BLACK")+
  scale_y_log10(limits=c(1,2000))+
  xlab("d (minimum module density threshold)")+ 
  ylab("Number of resulting modules (LOG scale)")+
  myTheme
p5

modStats2 <- melt(modStats, id.vars=c("spearmanCut", "D", "C","B","S","Q","E","G","numModules"))
filter<- modStats2$E==2
levels(modStats2$variable)<-factor(c("MIN","MAX","AVGERAGE","MEDIAN"))
modStats2$Statistic<-modStats2$variable
modStats2$variable<-NULL
p6 <- ggplot(modStats2[filter,], aes(x=D, y=value, colour=Statistic, group=Statistic)) +plotBase +
  scale_x_continuous(label=d_formatter) +
  scale_y_continuous(limits=c(0,85)) +
  geom_point() +
  geom_line(weight=5) +
  facet_grid(. ~ spearmanCut, labeller=mf_labeller) +
  xlab("d (minimum module density threshold)")+ 
  ylab("Genes count") + 
  theme(legend.position="none")+
  myTheme;
p6


filter<- modStats$E==3

p7 <- ggplot(modStats[filter,], aes(x=D, y=numModules)) +plotBase +
  scale_x_continuous(label=d_formatter) +
  scale_y_log10(limits=c(1,2000)) +
  geom_point() +
  geom_line(weight=5, colour="BLACK") +
  facet_grid(. ~ spearmanCut, labeller=mf_labeller,) +
  xlab("d (minimum module density threshold)")+ 
  ylab("Number of resulting modules (LOG scale)")+
  myTheme;
p7

filter<- modStats2$E==3
p8 <- ggplot(modStats2[filter,], aes(x=D, y=value, colour=Statistic, group=Statistic)) +plotBase +
  scale_x_continuous(label=d_formatter) +
  scale_y_continuous(limits=c(0,85)) +
  geom_point() +
  geom_line(weight=5) +
  facet_grid(. ~ spearmanCut, labeller=mf_labeller) +
  xlab("d (minimum module density threshold)")+
  ylab("Genes count")+
  theme(legend.position="none")+
  myTheme;
p8

grid.arrange(p5,p6,p7,p8, cols=2, main = "Main title")
dataGrid<-c(1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,3,4,3,4,3,4,3,4,3,4,3,4,3,4,3,4,4)
p9<-multiplot(plotlist=list(p5,p6,p7,p8,NULL), layout=matrix(data=c(1,2,3,4),nrow=2,ncol=2,byrow=FALSE),file="my.png")

ggsave(filename="SC_numMods_E=2.png",path="../../Figures/",plot = p5,width = 6,height=4)
ggsave(filename="SC_numGenes_E=2.png",path="../../Figures/",plot = p6,width = 6,height=4)
ggsave(filename="SC_numMods_E=3.png",path="../../Figures/",plot = p7,width = 6,height=4)
ggsave(filename="SC_numGenes_E=3.png",path="../../Figures/",plot = p8,width = 6,height=4)

plots<-list();
modStats2 <- melt(modStats, id.vars=c("spearmanCut", "D", "C","B","S","Q","E","G","numModules"))
levels(modStats2$variable)<-factor(c("MIN","MAX","AVGERAGE","MEDIAN"))
modStats2$Statistic<-modStats2$variable
modStats2$variable<-NULL
d_formatter <- function(x) {
  lab <- paste0('d = ',x);
}

# i<- 1;
# for(sFilter in c(0.1,0.2,0.3,0.35,0.4))
# {
#   for(eFilter in 2:3)
#   {
#   filter<- modStats$E==eFilter & modStats$spearmanCut==sFilter
#   maxy<-2000;
#   if(i!=11&&i!=12){color1="grey20";color2="white"} else {color1="yellow";color2="black"}
#   rect <- data.frame (xmin=-Inf, xmax=Inf, ymin=maxy-maxy/50, ymax=Inf)
#   rect2 <- data.frame (xmin=0.58, xmax=0.62, ymin=0.00001, ymax=Inf)
#   
#   p1 <- ggplot(modStats[filter,], aes(x=D, y=numModules)) +
#     scale_x_continuous(label=d_formatter) +
#     scale_y_log10(limits=c(1, 2000)) +
#     geom_point() +
#     geom_line(weight=5, colour="BLACK") +
#     geom_rect(data=rect,aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=color1,color=color1, inherit.aes = FALSE)+
#     ggtitle(label=paste0("c = ", sFilter))+
#     theme(plot.title = element_text(vjust=-2, color=color2, face="bold"))+
#     xlab("")+ 
#     ylab("")+
#     theme(panel.margin=unit(x=0,units="in"), axis.ticks.y = element_blank(),axis.text.y=element_blank(), axis.text.x=element_text(angle=90, color="BLACK"))
#   if(i==11||i==12){
#     p1<-p1+geom_rect(data=rect2,aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=color1,color=color1, alpha=0.5, inherit.aes = FALSE)
#   }
# 
#   
#   plots[[i]]<-p1;
#   i<-i+1
# 
#   filter<- modStats2$E==eFilter & modStats2$spearmanCut==sFilter
#   maxy<-85;
#   if(i!=11&&i!=12){color1="grey20";color2="white"} else {color1="yellow";color2="black"}
#   rect <- data.frame (xmin=-Inf, xmax=Inf, ymin=maxy-maxy/50, ymax=Inf)
#   rect2 <- data.frame (xmin=0.58, xmax=0.62, ymin=-Inf, ymax=Inf)
#   p2 <- ggplot(modStats2[filter,], aes(x=D, y=value, colour=Statistic, group=Statistic)) +
#     scale_x_continuous(label=d_formatter) +
#     scale_y_continuous(limits=c(0, maxy))+
#     geom_point() +
#     geom_line(weight=5) +
#     geom_rect(data=rect,aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=color1,color=color1, inherit.aes = FALSE)+
#     ggtitle(label=paste0("c = ", sFilter))+
#     theme(plot.title = element_text(vjust=-2, color=color2, face="bold"))+
#     xlab("")+ 
#     ylab("") + 
#     theme(legend.position="none")+
#     theme(panel.margin=unit(x=0,units="in"),axis.ticks.y = element_blank(),axis.text.y=element_blank(), axis.text.x=element_text(angle=90, color="BLACK"))
#   if(i==11||i==12){
#     p2<-p2+geom_rect(data=rect2,aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=color1,color=color1, alpha=0.5, inherit.aes = FALSE)
#   }
#   
#   plots[[i]]<-p2;
#   i<-i+1
#   
#   }
# }
# blank <- grid.rect(gp=gpar(col="white"))
# blank
# library(gridExtra)
# grid.arrange(plots[[1]],plots[[2]])
# 
# grid.arrange(gpLine, gpBar, heights = c(2/3, 1/3), 
#              sub = textGrob("this is my signature", x=1, hjust=1, vjust=0))
# 
# plots
# p9<-multiplot(plotlist=plots, layout=matrix(data=seq(1,20),nrow=2,ncol=10,byrow=FALSE),file="my.png")


Data <- readRDS(file="../../normalization.RDS")