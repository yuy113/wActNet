#' @title Two heuristic algorithms to detect multiple optimized subnetwork
#' @description Detect multiple possible optimized subnetworks in network object
#' based on either combined node scores and edge scores or node scores only in the network
#' by heuristic algorithm by maximizing the weighted sum of node scores and edge scores in all possible subnetworks
#' or by maximizing sum of node scores in all possible subnetworks
#' @param network  igraph or graphNEL object with nodes and edges
#'  with default node scores and node names, the edges with names and scores
#' @param node.scores  The vector of scores of the corresponding nodes in network
#'  the names of the node.scores are the names of the nodes in the network
#' @param edge.scores The vector of scores of the connecting edges of all the nodes in network
#'  the names of the edge.scores are the names of the edges in the network, which is
#'  the names of two connecting nodes in the network joined with "_"
#' @param weightratio.edge.node The weight ratio of edges scores versus node scores
#' with default value 1, reflecting the weights of the nodes and edges in the network at object function
#' @param method Two possible fixed input values-"NodeOnly",apply Dittrich's heuristic algorithm based on node scores only(Dittrich 2008)
#'  or "NodeEdge" optimizition based on both node scores and edge scores
#' @param ncluster The number of prespecified signaling clusters/optimized networks to detect from the whole network
#' Iteratively running our proposed algorithm to identify optimized subnetwork after removing all the previous detected subnetworks
#' from the input whole network
#' @return The list of multiple possible optimized subnetworks in igraph format with attributes node and edge
#' with both name and scores for both nodes and edges, the length of the list equal to ncluster
#' Please note when the function completes execution and fails to identify any further optimized or active subnetworks, 
#' a warning message "No positive nodes" will be generated.
#'
#' @references 	MT. Dittrich, GW. Klau, A. Rosenwald, T. Dandekar, and T. Muller
#' Identifying functional modules in protein-protein interaction
#' networks: an integrated exact approach. \emph{Bioinformatics}, 24 (13): i223-231, 2008
#'
#'
#'  D.Beisser, GW. Klau, T. Dandekar, T. Muller, and MT. Dittrich. Bionet: an r-package for the functional analysis of biological
#'  networks. \emph{Bioinformatics}, 26 (8): 1129-1130, 2010.
#'
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
#'
#' network.test1<-network.test$Network
#' node.scores<-network.test$NodeScore
#' edge.scores<-network.test$EdgeScore
#'
#' #identify all possible signaling modules in simulated network
#' # using Dittrich's method using node scores only
#' multi.mod.e<-MultiModuleFind(network.test1,node.scores,edge.scores,
#' weightratio.edge.node=1,ncluster=3,method="NodeEdge")
#' # identify all possible signaling modules in simulated network
#' # using our proposed method based on node scores and edge scores
#' multi.mod.n<-MultiModuleFind(network.test1,node.scores,edge.scores,
#' weightratio.edge.node=1,ncluster=3,method="NodeOnly")
#' @export
#'
#'
#'
MultiModuleFind<-function(network,node.scores,edge.scores,weightratio.edge.node=1,ncluster,method){
  #subnetwork to induced subgraph from graph object considering both the vertices and edges,
  #with remove.vertex=T,F; if True, then remove the vertices which is not in the edges, otherwise keep them
  #vid-the names of the vertices in the graph
  #eid-the names of the edges in the graph, with the format-Name(from_vertex)_Name(to_vertex)
  #output:the subgraph containing only the vertices of vid and the edges


.subNetwork0 <- function(nodeList, network)
{
  if(is(network, "igraph"))
  {
    mapping <- seq(1, (length(V(network))))
    if(is.null(V(network)$name))
    {
      V(network)$name <- as.character(V(network))
    }
    names(mapping) <- V(network)$name
    nodeList = mapping[nodeList]
    if(any(is.na(nodeList)))
    {
      nodeList = na.omit(nodeList)
      warning("Not all nodes found in network")
    }
    subgr <- induced.subgraph(network, vids=nodeList) 
  }
  else
  {
    subgr <- subGraph(nodes(network)[nodes(network) %in% nodeList], network)
  }
  return(subgr)
}

subnetwork.e<-function(graph,vid,eid,remove.vertex=F){
  
  
  if(is.null(vid) || all(is.na(vid))){
    warning("No nodes for subnetwork")
    break
  }
  
  if(is.null(E(graph)$name) || all(is.na(E(graph)$name)) ){
    warning("No edge names in the network")
    E(graph)$name<-names(edge.scores)
    
    names(E(graph)$score)
  }
  
  
  if(!(is.null(V(graph)$score)  || all(is.na(E(graph)$name)))        ){
    node.scores<-V(graph)$score
    
    
    names(node.scores)<-V(graph)$name
    node.score.sub<-node.scores[vid]
    
  }
  if(!(is.null(V(graph)$weight)  || all(is.na(V(graph)$weight)))         ){
    node.weight<-V(graph)$weight
    
    names(node.weight)<-V(graph)$name
    node.weight.sub<-node.weight[vid]
    
  }
  
  
  
  
  
  if(! (is.null(E(graph)$weight) || all(is.na(E(graph)$weight)))    ){
    edge.weight<-E(graph)$weight
    names(edge.weight)<-E(graph)$name
    
  }
  

    if( any(is.na(eid)) ){
    eid<-as.vector(na.omit(eid))
    
  }
  
  
  if(!is.null(E(graph)$score)){
    edge.scores<-E(graph)$score
    names(edge.scores)<-E(graph)$name
    edge.score.sub<-edge.scores[eid]
    
  }
  
  
  if(is.null(eid)  || all(is.na(eid))  ){
    warning("No edges for subnetwork")
    g<-make_empty_graph(length(vid),directed=F)
   V(g)$name<-vid
    
  
    
    if(  (is.null(V(graph)$score) | all(is.na(V(graph)$score)))  && !(is.null(V(graph)$weight)| all(is.na(V(graph)$weight))) ){
      
      
      V(g)$weight= node.weight.sub
      
      
      }
    
    if(!(is.null(V(graph)$score) | all(is.na(V(graph)$score))) && (is.null(V(graph)$weight)| all(is.na(V(graph)$weight))) ) {
      
      
      V(g)$score= node.score.sub
      
      }
    
   if(!(is.null(V(graph)$score) | all(is.na(V(graph)$score))) && !(is.null(V(graph)$weight)| all(is.na(V(graph)$weight))) ) {
     
     V(g)$weight= node.weight.sub
     V(g)$score= node.score.sub
     
   }
   
   
   # g <- graph_from_data_frame(eid, directed=F, vertices=pos.nodes.dat)
    return(g)
    break
  }
  

  
  names.edge.score.sub<-names(edge.score.sub)
  
  from.name.sub <- unlist(strsplit(names.edge.score.sub, "_"))[seq(1,
                                                                   2 * length(names.edge.score.sub), by = 2)]
  to.name.sub <- unlist(strsplit(names.edge.score.sub, "_"))[seq(2,
                                                                 2 * length(names.edge.score.sub), by = 2)]
  edge.name.sub<-as.matrix(cbind(from.name.sub,to.name.sub))
  if(dim( edge.name.sub)[1]>=1){
    edge.name.sub.node<-matrix(NA,nrow=dim(edge.name.sub)[1],ncol=dim(edge.name.sub)[2])
    
    
    
    for(i in 1:dim( edge.name.sub)[1]){
      if(length(intersect(vid,edge.name.sub[i,]))==2){
        edge.name.sub.node[i,] <-edge.name.sub[i,]
      }
      else
        edge.name.sub.node[i,] <-c(NA,NA)
      
    }
    edge.name.sub.node<-na.omit(edge.name.sub.node)
    from.name.edge.sub.nodes<-edge.name.sub.node[,1]
    to.name.edge.sub.nodes<-edge.name.sub.node[,2]
  }
  
  if( remove.vertex){
    if(length(unique(c( from.name.edge.sub.nodes, to.name.edge.sub.nodes)))<length(unique(vid))){
      vid<-unique(c( from.name.edge.sub.nodes, to.name.edge.sub.nodes))
      if ( dim(edge.name.sub.node)[1] >= 1){
        edge.name.sub.node2<-matrix(NA,nrow=dim(edge.name.sub.node)[1],ncol=dim(edge.name.sub.node)[2])
        
        for(i in 1:dim( edge.name.sub.node)[1]){
          if(length(intersect(vid,edge.name.sub.node[i,]))==2){
            edge.name.sub.node2[i,] <-edge.name.sub[i,]
          }
          else
            edge.name.sub.node2[i,] <-c(NA,NA)
          
        }
        edge.name.sub.node2<-na.omit(edge.name.sub.node2)
        from.name.edge.sub.nodes<-edge.name.sub.node2[,1]
        to.name.edge.sub.nodes<-edge.name.sub.node2[,2]
      }
    }
  }
  
  
  
  #  remove.vertex=F
  
  
  
  ################################################################################
  ##in case remove.vertex=T, remove some nodes from input parameter-vid#########
  ###############################################################################
  node.score.sub<-node.score.sub[vid]
  if( ! (is.null(V(graph)$weight) || all(is.na(V(graph)$weight)))  ){
    node.weight.sub<-node.weight[vid]
  }
  names.edge.sub.nodes<-paste( from.name.edge.sub.nodes,to.name.edge.sub.nodes,sep="_")
  
  if(! (is.null(E(graph)$score) || all(is.na(E(graph)$score)))            ){
    edge.score.sub.nodes<-edge.scores[names.edge.sub.nodes]
  }
  
  
  if(!(is.null(E(graph)$weight) || all(is.na(E(graph)$weight)))         ){
    edge.weight.sub.nodes<-edge.weight[names.edge.sub.nodes]}
  
  
  if(  (is.null(V(graph)$score) || all(is.na(V(graph)$score))) &&  (is.null(V(graph)$weight) || all(is.na(V(graph)$weight))) ){
    pos.nodes.dat<-data.frame(name=vid)}
  
  if(  (is.null(V(graph)$score) || all(is.na(V(graph)$score))) && !(is.null(V(graph)$weight) || all(is.na(V(graph)$weight))) ){
    pos.nodes.dat<-data.frame(name=vid,weight= node.weight.sub)}
  
  if(  !(is.null(V(graph)$score) || all(is.na(V(graph)$score))) && (is.null(V(graph)$weight) || all(is.na(V(graph)$weight))) ){
    pos.nodes.dat<-data.frame(name=vid,score= node.score.sub)}
  
  if(  !(is.null(V(graph)$score) || all(is.na(V(graph)$score))) && !(is.null(V(graph)$weight) || all(is.na(V(graph)$weight))) ){
    pos.nodes.dat<-data.frame(name=vid,score= node.score.sub,weight=node.weight.sub)}
  
  if(   !(is.null(E(graph)$score) || all(is.na(E(graph)$score))) && !(is.null(E(graph)$weight) || all(is.na(E(graph)$weight)))  ){
    edge.pos.nodes.dat<-data.frame(from= from.name.edge.sub.nodes,to=to.name.edge.sub.nodes,
                                   score=edge.score.sub.nodes,name= names.edge.sub.nodes,weight=edge.weight.sub.nodes)}
  
  if(  !(is.null(E(graph)$score) || all(is.na(E(graph)$score))) && (is.null(E(graph)$weight) || all(is.na(E(graph)$weight)))  ){
    edge.pos.nodes.dat<-data.frame(from= from.name.edge.sub.nodes,to=to.name.edge.sub.nodes,
                                   score=edge.score.sub.nodes,name= names.edge.sub.nodes)}
  
  
  if(  (is.null(E(graph)$score) || all(is.na(E(graph)$score))) && (is.null(E(graph)$weight) || all(is.na(E(graph)$weight)))  ){
    edge.pos.nodes.dat<-data.frame(from= from.name.edge.sub.nodes,to=to.name.edge.sub.nodes,
                                   name= names.edge.sub.nodes)}
  
  
  if((is.null(E(graph)$score) || all(is.na(E(graph)$score))) && !(is.null(E(graph)$weight) || all(is.na(E(graph)$weight)))  ){
    edge.pos.nodes.dat<-data.frame(from= from.name.edge.sub.nodes,to=to.name.edge.sub.nodes,
                                   weight=edge.weight.sub.nodes,name= names.edge.sub.nodes)}
  
  
  
  g <- graph_from_data_frame(edge.pos.nodes.dat, directed=F, vertices=pos.nodes.dat)
  return(g)

}


 num.cluster <- 1
  lst.modules <- list()
  
  
  if (is.null(E(network)$name)) {
    E(network)$name <- names(edge.scores)
  }
  if (method == "NodeOnly") {
    module.1 <- runFastHeinz(network, node.scores)
  }
  if (method == "NodeEdge") {
    module.1 <- runFastHeinz.e(network, node.scores, edge.scores, 
                               weightratio.edge.node = weightratio.edge.node)
  }
  
 # names(node.scores)
#  E(network)$name
  
  lst.modules[[1]]<-module.1

  
  while( !(is.null(V(module.1)$name))  && num.cluster<ncluster ){
    lst.modules[[num.cluster]]<-module.1

    vid.after.cluster1<-setdiff(names(node.scores)
                                ,V(module.1)$name)

    network.after.cluster1<-subnetwork.e(network,vid=vid.after.cluster1,eid=names(edge.scores),remove.vertex=F)

    node.scores2<-node.scores[vid.after.cluster1]
    eid.after.cluster1<-E(network.after.cluster1)$name



    edge.scores2<-edge.scores[eid.after.cluster1]


    if (method=="NodeOnly"){
      module.2<-runFastHeinz(network.after.cluster1,node.scores2)

    }
    if (method=="NodeEdge"){
      module.2<-runFastHeinz.e(network.after.cluster1,node.scores2,edge.scores2,weightratio.edge.node=weightratio.edge.node)

    }



    num.cluster<-num.cluster+1

    if(num.cluster==ncluster && (!is.null(V(module.2)$name)) ){
      lst.modules[[num.cluster]]<-module.2
    }
    network<-network.after.cluster1
    node.scores<-node.scores2

    edge.scores<-edge.scores2
    module.1<-module.2


  }

  lst.modules

}
