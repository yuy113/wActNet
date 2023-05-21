#' @title Generate node scores,edge scores and network
#' @description Generate both node and edge scores for biological correlation network
#' based on P values of clinical outcome and nodes
#' and also P values by multiple permulation tests
#' on all the possible edges in biological correlaton network
#' This function calculates edge scores based on P values in multiple testing
#' assuming P values follows a mixture of two beta distributions-
#'  beta(a,1),beta(1,1)-Uniform(0,1)
#' with those in permutated samples through nonparametric permutation tests.
#' Also this function calculates node scores based on P values in multiple testing
#' assuming P values follows a mixture of two beta distributions-
#'  beta(a,1),beta(1,1)-Uniform(0,1)
#' @param pval.node  P values from multiple hypothesis testing in vector format
#'  hypothesis testing on clinical outcome with each individual node (variable)
#'  The names of the edges (names of two connecting nodes(variables)) must be specified.
#' @param pval.edge  P values from multiple hypothesis testing in vector format
#'  The names of the edges (names of two connecting nodes(variables)) must be specified.
#' @param FDR.node The overall level of False Discovery Rate level for multiple testing
#' on all the nodes(variables) ,typically 0.05,0.01
#' @param FDR.edge The level of False Discovery Rate level for multiple testing
#' on all possible edges,typically 0.01,0.001
#' @param dat  Dataset containing the observations for p covariates (nodes)
#'  The column names of the dat must be specified.
#' @return A list of five elements as decribed below:
#'
#'This function returns an object with class \emph{NetworkScore}. The items in the object are:
#'\item{NodeScore}{The vector of scores of the nodes in network}
#'\item{EdgeScore}{The vector of scores of all possible edges in network}
#'\item{FDR.node}{False Discovery Rate for multiple testing on the nodes in the network}
#'\item{FDR.edge}{False Discovery Rate for multiple testing on all possible edges in the network}
#'\item{Network}{An igraph object-Network with attributes in nodes-node score and name,and edges-edge scores and names}
#'
#' @examples
#' dat1<-matrix(rnorm(4000),ncol=40,nrow=100)
#' colnames(dat1)<-paste("Var", as.character(1:40),sep="")
#'
#' # simulate the p values for all the possible edges in the network
#' ind.pos.pval.edge<-rbinom(40*39/2,1,0.5)
#' pval.edge<-(1-ind.pos.pval.edge)*runif(40*39/2)+ind.pos.pval.edge*rbeta(40*39/2,0.1,1)
#' names(pval.edge)<-unlist(sapply(1:39,function(i){sapply((i+1):40, function(j){paste(paste("Var",
#' as.character(i),sep=""),paste("Var",as.character(j),sep=""),sep="_")})}))
#'
#' # simulate p values for all the nodes in the network
#' ind.pos.pval.node<-rbinom(40,1,0.2)
#' pval.node<-(1-ind.pos.pval.node)*runif(40)+ind.pos.pval.node*rbeta(40,0.1,1)
#' names(pval.node)<-paste("Var", as.character(1:40),sep="")
#'
#' # generate the node score-NodeScore, edge scores-EdgeScore and igraph object-Network
#' network.test<-uniform.beta.node.edge.score(pval.node,pval.edge,0.05,0.05,dat1)
#' @export
#'
#'
uniform.beta.node.edge.score<-function(pval.node,pval.edge,FDR.node,FDR.edge,dat)
{

  fdr1<-FDR.node
  fdr2<-FDR.edge



uniform.beta.edge.score<-function(pval,fdr){
  #fit mixture of uniform-beta model to p-values of the edges
  #input,parameters,pval-the vector of the p-values of the edges in the graph
  #input,fdr-pre-specified FDR for the multiple testing framework
  #a-the proportion of the noise
  #b-the signal distribution-Beta(b,1)
  #pval-the vector of p-values
  if(is.null(names(pval))){
    warning("unnamed P values,please specify the names")
    break
  }
  if(any(is.na(names(pval)))) {
    warning("P values contain unnamed values")
    pval<-pval[!is.na(names(pval))]
  }

  log.like.edge<-function(par){
    a<- par[1]
    b<-par[2]
    a1<-exp(a)/(1+exp(a))
    b1<-exp(b)/(1+exp(b))
    like.edge<-a1+(1-a1)*b1*pval^(b1-1)
    log.like.edge1<-sum(-log(like.edge))
    return(log.like.edge1)}
  logit.inv<-function(x){
    return(exp(x)/(1+exp(x)))
  }
  init<-runif(2,0.1,0.9)
  param.opt<-optim(par=init,log.like.edge)$par
  a.est<-logit.inv(param.opt[2])
  lambda.est<-logit.inv(param.opt[1])
#p(p-1)/2=length(pval)
  p<-sqrt(1+2*length(pval))
  fdr2<-fdr
  pi.est<-lambda.est+(1-lambda.est)*a.est
  tau.fdr<-((pi.est-fdr2*lambda.est)/(fdr2*(1-lambda.est)))^(1/(a.est-1))
  edge.score<-(a.est-1)*(log(pval)-log(tau.fdr))
  names(edge.score)<-names(pval)
  return(edge.score)
}
edge.scores<-uniform.beta.edge.score(pval.edge,fdr2)

#Contruct the igraph object from the names of the covariates in the studies
#incoporate the node scores as the attributes of the nodes in the network
#incoporate the edge scores as the attributes of the edges in the network
#use the function-graph_from_data_frame from the library-igraph
#input:dat-the dataset of the covariates in the study, without any nonmissing names of the covariates
#input:edge.score-the edge scores from the function-edge_score_func
#input:node.weight, the weights or scores for the covariates as the nodes in the network
#input:edge.weight, the weights for the relationships among the covariates in the dataset-dat
#the inputs, dat and edge.score are requried for this function
#other input variables are optional
induced.graph.data.frame<-function(dat,node.score=NA,edge.score,node.weight=NA,edge.weight=NA){
#library(igraph)
  if(is.null(names(node.score)) && (!is.null(node.score))){
    warning("Unnamed node scores")
  }


  if(!is.null(node.score)){
    if(any(is.na(names(node.score)))){
    warning("NA in names of node scores, the node scores will be removed")
    node.score<-node.score[!is.na(node.score)]
  }}

  if(is.null(colnames(dat))){
    warning("Unnamed variables in the dataset")
  }
  if(any(is.na(colnames(dat)))){
    warning("NA in variable names of dataset, the corresponding variables will be removed")
    dat<-dat[,colnames(dat)[!is.na(colnames(dat))]]
    }

  if(!is.null(node.score) && length(node.score)!=length(colnames(dat))){
    warning("The number of nodes not same in data and node scores")
    names.node<-intersect(names(node.score),colnames(dat))
    node.score<-node.score[names.node]
    dat<-dat[,names.node]
  }

  if(is.null(names(edge.score)) && (!is.null(edge.score))){
    warning("Unnamed edge scores")
  }

  if(!is.null(edge.score)){
  if(any(is.na(names(edge.score)))){
    warning("NA in names of edge scores, the edge scores will be removed")
    edge.score<-edge.score[!is.na(edge.score)]
  }}

  ## Compare lengths of colnames(dat),node.score,node.weight
 # if(!is.null(node.score) && !is.null(node.weight) && !is.null(colnames(dat))) {
#  if (length(unique(c(length(node.score), length(node.weight), length(colnames(dat)))))!=1){
 #   warning("node score, node weight and variables in the data must have same length") }
#  }
## The typical case is that these tables are read in from files....
  if(is.null(node.score) && is.null(node.weight)){
node <- data.frame(name=colnames(dat)[order(colnames(dat))])}

  if(!is.null(node.score) && is.null(node.weight)){
    node <- data.frame(name=names(node.score),score=node.score)}

  if(is.null(node.score) && !is.null(node.weight)){
    node <- data.frame(name=colnames(dat)[order(colnames(dat))],weight=node.weight)}

if (!is.null(names(edge.score)) && (!is.null(edge.weight))){

  from.name<-unlist(strsplit(names(edge.score),"_"))[seq(1,2*length(names(edge.score)),by=2)]
  to.name<-unlist(strsplit(names(edge.score),"_"))[seq(2,2*length(names(edge.score)),by=2)]
relations <- data.frame(from=from.name,
                        to=to.name,
                        score=edge.score,
                       weight=edge.weight)
}

  if (!is.null(names(edge.score)) && is.null(edge.weight)){

    from.name<-unlist(strsplit(names(edge.score),"_"))[seq(1,2*length(names(edge.score)),by=2)]
    to.name<-unlist(strsplit(names(edge.score),"_"))[seq(2,2*length(names(edge.score)),by=2)]
    relations <- data.frame(from=from.name,
                            to=to.name,
                            score=edge.score)
  }

  g <- graph_from_data_frame(relations, directed=F, vertices=node)

  g


}
da.igraph<-induced.graph.data.frame(dat,node.score=NULL,edge.score=edge.scores,node.weight=NULL,edge.weight=NULL)
#get node scores from the functions in BioNet library
#input, da.igraph-the graph object with the format-igraph from igraph library
#input:the vector of p-values of the nodes ,
#and the name of the vector must equal to names of nodes in the network
#adjusted node score based on p-values
node.score<-function(da.igraph,pval,fdr){
p<-length(pval)
fdr1<-fdr
#library(BioNet)
fb.bm.node<-fitBumModel(pval,plot=F)
scoreNodes(da.igraph,fb=fb.bm.node,fdr=fdr1)

}
node.scores<-node.score(da.igraph,pval.node,fdr1)

network.dat<-induced.graph.data.frame(dat,node.score=node.scores,edge.score=edge.scores,node.weight=NULL,edge.weight=NULL)
z <- list(NodeScore=node.scores, EdgeScore=edge.scores,Network=network.dat)
class(z) <- "NetworkScore"

return(z)

}



