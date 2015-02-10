constructFrequencyNetwork <- function(networkAddresses, threshold=0.1)
{
  net <- read.table(file=networkAddresses[[1]],sep=",",header=TRUE);
  row.names(net) <- net[,1];
  net<- net[,-1];
  net[is.na(net)]<-0;
  numFeatures <- dim(net)[1];
  features <- names(net);
  freqNet <- matrix(data=0,nrow=numFeatures,ncol=numFeatures);
  freqNet <- freqNet + as.integer(net > threshold);
  
  if(length(networkAddresses > 1))
  {
    for(i in 2:length(networkAddresses))
    {
      net <- read.table(file=networkAddresses[[i]],sep=",",header=TRUE);
      row.names(net) <- net[,1];
      net<- net[,-1];
      net[is.na(net)]<-0;
      order <- match(features, names(net));
      net <- net[order,order];
      freqNet <- freqNet + as.integer(net > threshold);
    }
  }
  rm(net);
  freqNet <- as.data.frame(freqNet)
  names(freqNet) <-features;
  row.names(freqNet) <- features;
  return(freqNet);
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
  
  if(length(networkAddresses > 1))
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

#networkAddresses <- c("Data//coexpressionNetworks/test262.txt",
#                      "Data//coexpressionNetworks/test26.txt",
#                      "Data//coexpressionNetworks/test28.txt",
#                      "Data//coexpressionNetworks/test29.txt",
#                      "Data//coexpressionNetworks/test30.txt",
#                      "Data//coexpressionNetworks/test31.txt");

networkAddresses <- c("Data/coexpressionNetworks/MGH26_spearman_int.txt",
                      "Data/coexpressionNetworks/MGH262_spearman_int.txt",
                      "Data/coexpressionNetworks/MGH28_spearman_int.txt",
                      "Data/coexpressionNetworks/MGH29_spearman_int.txt",
                      "Data/coexpressionNetworks/MGH30_spearman_int.txt",
                      "Data/coexpressionNetworks/MGH31_spearman_int.txt");

threshold<- 0.3;
freqNet <- constructFrequencyNetwork(networkAddresses,threshold);
write.table(x=freqNet,file=paste("Data/coexpressionNetworks/frequencyNetwork_",threshold,".txt"),sep="\t",quote=FALSE,row.names = FALSE,col.names = FALSE);
write.table(x=names(freqNet),file=paste("Data//coexpressionNetworks/geneOrder_",threshold,".txt"),sep="\t",quote=FALSE);
rm(freqNet);

#supportVector <- constructSupportVector(networkAddresses);
#write.table(x=supportVector,file=paste("Data/coexpressionNetworks/supportVector_",threshold,".txt"),sep="\t",quote=FALSE,row.names = FALSE,col.names = FALSE);

