thresholdNetwork <- function(network, threshold=NULL, targetFrequency=NULL, binaryOutput=FALSE)
{
  if(is.null(x=threshold))
  {
    if(is.null(x=targetFrequency))
      return(NULL);
    threshold=quantile(abs(network[upper.tri(network)]),1-targetFrequency);
  }
  if(binaryOutput)
  {
    network <- 1*(abs(network) >= threshold);
  }else
  {
    network[abs(network) < threshold] <- 0;
  }

  return(network);
}

constructSupportVector <- function(networkAddresses)
{   
  net <- read.table(file=networkAddresses[[1]],sep=",",header=TRUE);
  row.names(net) <- net[,1];
  net<- net[,-1];
  net[is.na(net)]<-0;
  numFeatures <- dim(net)[1];
  numNetworks <- length(networkAddresses)
  features <- names(net);
  supportVector <- matrix(data=0,nrow=(numFeatures*(numFeatures-1))/2,ncol=numNetworks+2);
  k<-1
  for(i in 1:(numFeatures-1))
  {
    for(j in (i+1):numFeatures)
    {  
      supportVector[k,1] <- i-1
      supportVector[k,2] <- j-1
      supportVector[k,3] <- as.integer(net[i,j]*1000)
      k<- k+1;
    }
  }
  
  if(length(networkAddresses) > 1)
  {
    for(n in 2:length(networkAddresses))
    {
      net <- read.table(file=networkAddresses[[n]],sep=",",header=TRUE);
      row.names(net) <- net[,1];
      net<- net[,-1];
      net[is.na(net)]<-0;
      order <- match(features, names(net));
      net <- net[order,order];
      k<-1
      for(i in 1:(numFeatures-1))
      {
        for(j in (i+1):numFeatures)
        {
          supportVector[k,n+2] <- as.integer(net[i,j]*1000);
          k<-k+1;
        }
      }
    }
  }
  rm(net);
  supportVector <- as.data.frame(supportVector)
  return(supportVector);
}
#edit output to fool CODENSE algorithm which sigfaults when only inputting a single network
doubling=TRUE;
outputSupportVector=FALSE;

inFile<-"/home/barand/Single_Cell_RNASeq/Data/coexpressionNetworks/GSE48865_norm_TPM_spearman_int.txt";
if(!doubling)
{
  outFile<-"/home/barand/Single_Cell_RNASeq/Data/coexpressionNetworks/GSE48865_norm_TPM_spearman_filt_0.0002.txt";
  supOutFile <-"/home/barand/Single_Cell_RNASeq/Data/coexpressionNetworks/GSE48865_supportVector.txt";
}else
{
  outFile<-"/home/barand/Single_Cell_RNASeq/Data/coexpressionNetworks/GSE48865_norm_TPM_spearman_filt_0.0002_double.txt";
  supOutFile <-"/home/barand/Single_Cell_RNASeq/Data/coexpressionNetworks/GSE48865_supportVector_double.txt";
}

#messed this up originally by setting targetFrequency to 0.04141572.
#this was frequency before CODENSE filteres out edges by threshold 'e'.
#Should have been 0.00158862374195553 or 0.0002718476  for 'e' of 2 or 3 respectively.

targetFrequency<-0.0002718476;
#threshold<- quantile(network, 1-targetFrequency);
binaryOutput<-TRUE;

print("Threshold Network");
network<- read.table(file=inFile,header=TRUE,sep=",");
row.names(network)<-network[,1];
network <- network[,-1]
network <- thresholdNetwork(network=network, targetFrequency=targetFrequency, binaryOutput=binaryOutput)  


if(!doubling)
{
  write.table(file=outFile, x=network,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE);
}else
{
  #multiply by 2 to get around CODENSE's crash when e=1.
  write.table(file=outFile, x=2*network,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE);
}

print("Output GeneOrder");
write.table(x=row.names(network),file=paste("/home/barand/Single_Cell_RNASeq/Data//coexpressionNetworks/GSE48865_geneOrder.txt"),sep="\t",quote=FALSE);
rm(network);

print("Generate Support Vector");
if(outputSupportVector)
{
  supportVector <- constructSupportVector(networkAddresses=c(inFile));
  if(doubling)
  {
    #add duplicate column to support vector get around CODENSE's crash when e=1.
    supportVector[,4] <- supportVector[,3]
  }
  write.table(x=supportVector,file=supOutFile,sep="\t",quote=FALSE,row.names = FALSE,col.names = FALSE);
}