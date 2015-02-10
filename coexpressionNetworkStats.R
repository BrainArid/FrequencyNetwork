filename="Data//coexpressionNetworks//frequencyNetwork_0.2.txt"
data <- as.matrix(read.table(file=filename))
name<-basename(path=filename)
name <- substr(name,start=0,stop=regexpr(".txt",text=name)[1]-1);

#list of edges = 6
di<-which(data==3)
edgeID1<-di%/%dim(data)[2]+1-1;
edgeID2<-di%%dim(data)[2]-1;

#geneMap <- read.table(file="Data//coexpressionNetworks//geneOrder.txt")
#genes<- geneMap[union(edgeID1,edgeID2)+1,1];
#constantNet <- matrix(nrow=dim(geneMap)[1],ncol=dim(geneMap)[1]);
#constantNet[is.na(constantNet)]<-0
#constantNet[edgeID1+1,edgeID2+1]<-1

png(filename=paste0("Figures/",name,"_distribution.png"));
membershipCountDistribution <- hist(data[upper.tri(data)],xlab="Membership Count", main=paste0("Frequency Net Distribution: ",name),breaks=seq(from=-0.5,to=6.5,by=1));
dev.off();

membershipLogCountDistribution <- hist(x=data[upper.tri(data)],xlab="Membership Count", main=paste0("Log Frequency Net Distribution: ",name),breaks=seq(from=-0.5,to=6.5,by=1));
png(filename=paste0("Figures/",name,"_logDistribution.png"));
plot(log(membershipLogCountDistribution$counts),x=membershipLogCountDistribution$mids,xlab="Membership Count", ylab="Log(frequency)", type="h",lwd=50, lend=1, main=paste0("Log Frequency Net Distribution: ",name))
dev.off();

nodeDegreeCounts <- apply(X=data,MARGIN=2,FUN=function(x){return(sum(x!=0))})-1;
png(filename=paste0("Figures/",name,"_degreeDistribution.png"));
nodeDegreeDistribution <- hist(x=nodeDegreeCounts,xlab="node degree",main=paste0("Node Degree Distribution: ",name),breaks=seq(from=min(nodeDegreeCounts)-0.5,to=max(nodeDegreeCounts)+0.5,by=1))
dev.off();

#density
numNodes<-dim(data)[1];
numEdges<-sum(data>0);
density<-(numEdges)/(numNodes*(numNodes-1));
