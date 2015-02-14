topGo_get_geneID2GO <- function()
{
  #source("http://bioconductor.org/biocLite.R")
  #biocLite("topGO")
  library("topGO")
  
  #read geneOntologyAnnotation
  # currently this list contains doubled-up genes like: UNQ5830/PRO19650/PRO19816 and MNB/DYRK
  # if you don't have this file get it here: http://geneontology.org/page/download-annotations
  print(getwd());
  temp <- read.table(file = "../..//Gene_GO_Annotation//gene_association.goa_human",blank.lines.skip = TRUE, comment.char = "!",header = FALSE,sep="\t",fill=TRUE)
  temp <- data.frame(Gene=temp[,3],GOTerm=temp[,5])
  allGenes <- unique(temp[,1])
  geneID2GO <- list()
  for(i in 1:length(allGenes))
    geneID2GO[[i]]<- as.character(temp[temp$Gene==allGenes[i],2]);
  names(geneID2GO)<-allGenes;
  rm(temp);
  rm(allGenes)
  return(geneID2GO);
}

topGOEnrich <- function(genesOfInterest, geneID2GO)
{ 
  #construct geneList contains 0 and 1 for not-interesting genes and interesting genes
  geneNames<-names(geneID2GO)
  geneList <- factor(as.integer(geneNames %in% genesOfInterest))
  names(geneList) <- geneNames
  str(geneList)
  rm(geneNames)
  
  #construct three part enrichment objects
  GO_BP <- new("topGOdata", description="Biological Process enrichment", 
               ontology = "BP", allGenes = geneList,
               annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
  GO_MF <- new("topGOdata", description="Molecular Function enrichment", 
               ontology = "MF", allGenes = geneList,
               annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
  GO_CC <- new("topGOdata", description="Cell Cycle enrichment", 
               ontology = "CC", allGenes = geneList,
               annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
  GO_BP.results.cf <- runTest(GO_BP,algorithm="classic",statistic = "fisher");
  GO_MF.results.cf <- runTest(GO_MF,algorithm="classic",statistic = "fisher");
  GO_CC.results.cf <- runTest(GO_CC,algorithm="classic",statistic = "fisher");
  GO_BP.results.ct <- runTest(GO_BP,algorithm="classic",statistic = "t");
  GO_MF.results.ct <- runTest(GO_MF,algorithm="classic",statistic = "t");
  GO_CC.results.ct <- runTest(GO_CC,algorithm="classic",statistic = "t");
  GO_BP.results.w1f <- runTest(GO_BP,algorithm="weight01",statistic = "fisher");
  GO_MF.results.w1f <- runTest(GO_MF,algorithm="weight01",statistic = "fisher");
  GO_CC.results.w1f <- runTest(GO_CC,algorithm="weight01",statistic = "fisher");
  
  topGOsCount<-15;
  GO_BP.allRes <- GenTable(GO_BP, classicFisher=GO_BP.results.cf, classicT=GO_BP.results.ct, weight01Fisher=GO_BP.results.w1f, orderBy="weight01Fisher",topNodes=topGOsCount);
  GO_MF.allRes <- GenTable(GO_MF, classicFisher=GO_MF.results.cf, classicT=GO_MF.results.ct, weight01Fisher=GO_MF.results.w1f, orderBy="weight01Fisher",topNodes=topGOsCount);
  GO_CC.allRes <- GenTable(GO_CC, classicFisher=GO_CC.results.cf, classicT=GO_CC.results.ct, weight01Fisher=GO_CC.results.w1f, orderBy="weight01Fisher",topNodes=topGOsCount);
  
  return(list(BP=GO_BP.allRes,MF=GO_MF.allRes,CC=GO_CC.allRes));
}
setwd("Data/codensedModules/frequencyNetwork_0.3/")
files<- list.files(pattern='*FO$')
geneMap <- read.table(file="../../coexpressionNetworks/geneOrder.txt")

for(file in files)
{
  print(paste0("Working on file: ", file));
  
  maxColumns<- max(count.fields(file=file))
  if(maxColumns==-Inf)
    next;
  
  moduleIndecies <- read.table(file=file,fill=TRUE,header=FALSE,colClasses=as.character(rep("numeric",maxColumns)))
  
  HGNCModules <-list();
  #geneID2GO <- topGo_get_geneID2GO();
  #enrichments <- list();
  numRows <-dim(moduleIndecies)[2]-1;
  numCols <-dim(moduleIndecies)[1]
  HGNCModules.m <- matrix(data=rep(x="",times=numRows*numCols),nrow=numRows,ncol=numCols,
                          dimnames=list(
                            NULL,
                            unlist(lapply(X=seq(1,numCols),FUN=function(x){return(paste0("module",x))}))));#data.frame(row.names=maxColumns-1,stringsAsFactors=FALSE);

  for(i in 1:dim(moduleIndecies)[1])
  {
    genes<-as.character(geneMap[as.vector(as.matrix(moduleIndecies[i,c(FALSE,!is.na(moduleIndecies[i,-1]))]))+1,1]);#the '+1' is to convert from CODENSE's base 0 to R's base 1 gene indexing scheme
    HGNCModules.m[1:length(genes),i] <- genes;
  }  
  text<- apply(X=as.matrix(apply(X=HGNCModules.m,MARGIN=1,FUN=function(x){paste(x,collapse="\t")})),MARGIN=2,FUN=function(x){paste(x,collapse="\n")})
  text<-paste0(apply(X=as.matrix(colnames(HGNCModules.m)),MARGIN=2,FUN=function(x){paste(x,collapse="\t")}),"\n",text);
  fileConn<-file(paste0(file,".hgnc.david.module"));
  writeLines(text=text, fileConn)
  close(fileConn);
}

setwd("../../..");
