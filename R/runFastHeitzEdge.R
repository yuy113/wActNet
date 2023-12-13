#' @title Heuristic algorithm to detect optimized subnetwork in correlation network
#' @description Detect the optimized subnetwork in one network object
#' based on combined node scores and edge scores in the netowrk
#' by heuristic algorithm by maximized the weighted sum of node scores and edge scores in all possible subnetworks
#' @param network  An igraph or graphNEL object with nodes and edges
#'  with default node scores and node names, the edges with names and scores
#' @param node.scores  The vector of scores of the corresponding nodes in network
#'  the names of the node.scores are the names of the nodes in the network
#' @param edge.scores The vector of scores of connecting edges of all the nodes in network
#'  the names of the edge.scores are the names of the edges in the network, which is
#'  the names of two connecting nodes in the network joined with _
#' @param weightratio.edge.node The weight ratio of edges scores versus node scores
#' with default value 1, reflecting the weights of the nodes and edges in the network at object function
#' @return The optimized subnetwork in igraph format with attributes node and edge
#' with both name and scores for nodes and also edges
#' Please note when the function completes execution and fails to identify any further optimized or active subnetworks, 
#' a warning message "No positive nodes" will be generated.
#' @examples
#' dat1<-matrix(rnorm(20000),ncol=200,nrow=100)
#' colnames(dat1)<-paste("Var", as.character(1:200),sep="")
#'
#' # simulate the p values for all the possible edges in the network
#' ind.pos.pval.edge<-rbinom(200*199/2,1,0.5)
#' pval.edge<-(1-ind.pos.pval.edge)*runif(200*199/2)+ind.pos.pval.edge*rbeta(200*199/2,0.1,1)
#' names(pval.edge)<-unlist(sapply(1:199,function(i){sapply((i+1):200, function(j){paste(paste("Var",
#' as.character(i),sep=""),paste("Var",as.character(j),sep=""),sep="_")})}))
#'
#' # simulate p values for all the nodes in the network
#' ind.pos.pval.node<-rbinom(200,1,0.2)
#' pval.node<-(1-ind.pos.pval.node)*runif(200)+ind.pos.pval.node*rbeta(200,0.1,1)
#' names(pval.node)<-paste("Var", as.character(1:200),sep="")
#'
#' # generate the node score-NodeScore, edge scores-EdgeScore and igraph object-Network
#' network.test<-uniform.beta.node.edge.score(pval.node,pval.edge,0.05,0.05,dat1)
#'
#' network.test1<-network.test$Network
#' node.scores<-network.test$NodeScore
#' edge.scores<-network.test$EdgeScore
#'
#' # detect the optimized subnetwork in simulated network
#' subnetwork.opt<-runFastHeinz.e(network.test1, node.scores,edge.scores,weightratio.edge.node=1)
#' @export
#'
#'




runFastHeinz.e<-function(network, node.scores,edge.scores,weightratio.edge.node=1) {
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
    
    
    if(is.null(vid)){
      warning("No nodes for subnetwork")
      break
    }
    
    if(is.null(E(graph)$name)){
      warning("No edges in the network")
      
    }
    
    
    if(!is.null(V(graph)$score)){
      node.scores<-V(graph)$score
      
      
      
      
      names(node.scores)<-V(graph)$name
      node.score.sub<-node.scores[vid]
      
    }
    if(!is.null(V(graph)$weight)){
      node.weight<-V(graph)$weight
      
      names(node.weight)<-V(graph)$name
      node.weight.sub<-node.weight[vid]
      
    }
    
    
    
    
    
    if(!is.null(E(graph)$weight)){
      edge.weight<-E(graph)$weight
      names(edge.weight)<-E(graph)$name}
    
    
    
    
    if(is.null(eid)){
      warning("No edges for subnetwork")
      if(is.null(V(graph)$score) && is.null(V(graph)$weight) ){
        pos.nodes.dat<-data.frame(name=vid)}
      
      if(is.null(V(graph)$score) && !is.null(V(graph)$weight) ){
        pos.nodes.dat<-data.frame(name=vid,weight= node.weight.sub)}
      
      if(!is.null(V(graph)$score) && is.null(V(graph)$weight) ){
        pos.nodes.dat<-data.frame(name=vid,score= node.score.sub)}
      g <- graph_from_data_frame(edge.pos.nodes.dat, directed=F, vertices=pos.nodes.dat)
      return(g)
      break
    }
    
    
    if(!is.null(E(graph)$score)){
      edge.scores<-E(graph)$score
      names(edge.scores)<-E(graph)$name
      edge.score.sub<-edge.scores[eid]
      
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
    
    
    
    if(!is.null(V(graph)$score)){
      
      node.score.sub<-node.scores[vid]
      names(node.score.sub)<-vid
      
    }
    if(!is.null(V(graph)$weight)){
      
      node.weight.sub<-node.weight[vid]
      names(node.weight.sub)<-vid
    }
    
    
    
    
    
    
    
    
    names.edge.sub.nodes<-paste( from.name.edge.sub.nodes,to.name.edge.sub.nodes,sep="_")
    
    if(!is.null(E(graph)$score)){
      edge.score.sub.nodes<-edge.scores[names.edge.sub.nodes]
    }
    
    
    if(!is.null(E(graph)$weight)){
      edge.weight.sub.nodes<-edge.weight[names.edge.sub.nodes]}
    
    
    if(is.null(V(graph)$score) && is.null(V(graph)$weight) ){
      pos.nodes.dat<-data.frame(name=vid)}
    
    if(is.null(V(graph)$score) && !is.null(V(graph)$weight) ){
      pos.nodes.dat<-data.frame(name=vid,weight= node.weight.sub)}
    
    if(!is.null(V(graph)$score) && is.null(V(graph)$weight) ){
      pos.nodes.dat<-data.frame(name=vid,score= node.score.sub)}
    
    if(!is.null(V(graph)$score) && !is.null(V(graph)$weight) ){
      pos.nodes.dat<-data.frame(name=vid,score= node.score.sub,weight=node.weight.sub)}
    
    if(!is.null(E(graph)$weight) && !is.null(E(graph)$score)){
      edge.pos.nodes.dat<-data.frame(from= from.name.edge.sub.nodes,to=to.name.edge.sub.nodes,
                                     score=edge.score.sub.nodes,name= names.edge.sub.nodes,weight=edge.weight.sub.nodes)}
    
    if(is.null(E(graph)$weight) && !is.null(E(graph)$score) ){
      edge.pos.nodes.dat<-data.frame(from= from.name.edge.sub.nodes,to=to.name.edge.sub.nodes,
                                     score=edge.score.sub.nodes,name= names.edge.sub.nodes)}
    
    
    if(is.null(E(graph)$weight) && is.null(E(graph)$score) ){
      edge.pos.nodes.dat<-data.frame(from= from.name.edge.sub.nodes,to=to.name.edge.sub.nodes,
                                     name= names.edge.sub.nodes)}
    
    
    if(!is.null(E(graph)$weight) && is.null(E(graph)$score) ){
      edge.pos.nodes.dat<-data.frame(from= from.name.edge.sub.nodes,to=to.name.edge.sub.nodes,
                                     weight=edge.weight.sub.nodes,name= names.edge.sub.nodes)}
    
    
    
    g <- graph_from_data_frame(edge.pos.nodes.dat, directed=F, vertices=pos.nodes.dat)
    return(g)}
  
  
  ################################################################################################################################
  #sum the scores of the negative nodes and the adjacent positive clusters and their connectng edge scores
  #also sum the node scores and the negative edge score of two positive scoring clusters connecting with negative scoring edges
  #input:node.score-the vector of the node scores from all the nodes and the positive scoring clusters
  #input:edge.score-the edge scores connecting all the pairs among all the negative scoring nodes and positive scoring clusters
  getPathScore.e <- function(path, graph1, graph2, node.score,edge.score)
  {
    #get the edge scores connecting either the negative scoring nodes and their adjacent positive scoring clusters
    #or the negative edge scores connecting the two positive scoring clusters
    
    edge.score.mst <-c()
    #for each positive scoring (meta-) node(cluster) or negative scoring nodes, if a positive scoring cluster is adjacent,
    # sum all the exising connecting edge scores
    #Otherwise, just assign 0
    for (i in 1:length(path)) {
      if (length(V(graph1)[path[i]]$clusters[[1]]) > 0)
      { edge.score.clusters.mst<-c()
      for (j in 1:length(V(graph1)[path[i]]$clusters[[1]])){
        
        edge.score.test<-edge.score[paste(V(graph1)[path[i]]$name,
                                          
                                          V(graph2)[unlist(V(graph1)[path[i]]$clusters)]$name[j],sep="_")]
        
        edge.score.clusters.mst<-c(edge.score.clusters.mst,
                                   ifelse(is.na(edge.score.test),
                                          edge.score[paste( V(graph2)[unlist(V(graph1)[path[i]]$clusters)]$name[j],
                                                            V(graph1)[path[i]]$name,sep="_")],edge.score.test))}
      
      if(grepl("cluster",V(graph1)[path[i]]$name) && (length(intersect(V(graph1)[path[i]]$name,V(graph2)[unique(unlist(V(graph1)[path]$clusters))]$name))==0)){
        
        edge.score.mst <- c( edge.score.mst, sum(edge.score.clusters.mst,
                                                 node.score[V(graph1)[path[i]]$name]))
      }
      
      if(grepl("cluster",V(graph1)[path[i]]$name) && (length(intersect(V(graph1)[path[i]]$name,V(graph2)[unique(unlist(V(graph1)[path]$clusters))]$name))==1)){
        
        edge.score.mst <- c( edge.score.mst, sum(edge.score.clusters.mst))
      }
      if(!grepl("cluster",V(graph1)[path[i]]$name))
      {
        edge.score.mst <- c( edge.score.mst, sum(edge.score.clusters.mst,node.score[V(graph1)[path[i]]$name]))
      }
      }
      else {
        
        if(grepl("cluster",V(graph1)[path[i]]$name) && (length(intersect(V(graph1)[path[i]]$name,V(graph2)[unique(unlist(V(graph1)[path]$clusters))]$name))==0)){
          
          edge.score.mst <- c( edge.score.mst,
                               node.score[V(graph1)[path[i]]$name])
        }
        
        if(grepl("cluster",V(graph1)[path[i]]$name) && (length(intersect(V(graph1)[path[i]]$name,V(graph2)[unique(unlist(V(graph1)[path]$clusters))]$name))==1)){
          
          edge.score.mst <- c( edge.score.mst, 0)
        }
        if(!grepl("cluster",V(graph1)[path[i]]$name))
        {
          edge.score.mst <- c( edge.score.mst, node.score[V(graph1)[path[i]]$name])
        }
        
        
      }
    }
    
    
    return(sum(c(edge.score.mst, node.score[V(graph2)[unique(unlist(V(graph1)[path]$clusters))]$name])))
  }
  ################################################################################################################################
  #calculate sum of the node score and the edge score of one particular vertex's neighboring vertices in the network
  #input:network, the input graph containg the vertex name, vertex score and edge score and edge names
  #v.name, the name of the particular vertex from the network
  #output:the vector containing the combined scores of the node scores of v's neighboring vertices
  # and the edge scores of v and its neighbouring vertices in the network
  neighbor.combined.scores<-function(network,v.name){
    # library(igraph)
    if (is.null(V(network)$name)) {
      warning("Unnamed nodes in the network")
    }
    if (any(is.na(V(network)$name))) {
      warning("NA in names of nodes in the network, the node score will be removed")
      node.scores <-V(network)$score[!is.na(V(network)$name)]
      names(node.scores)<-V(network)$name[!is.na(V(network)$name)]
    }
    
    if (is.null(E(network)$name)) {
      warning("Unnamed edges in the network")
    }
    if (any(is.na(E(network)$name))) {
      warning("NA in names of edges in the network, the edge score will be removed")
      edge.scores <- E(network)$score[!is.na(E(network)$name)]
      names(edge.scores) <- E(network)$name[!is.na(E(network)$name)]
    }
    
    if (all(!is.na(V(network)$name)) && all(!is.na(E(network)$name))){
      node.scores<-V(network)$score
      names(node.scores)<-V(network)$name
      edge.scores<-E(network)$score
      names(edge.scores)<-E(network)$name}
    neighbor.pos.nodes<-matrix(neighbors(network,
                                         v = v.name),ncol=1)
    
    neighbor.edge.score<-apply(neighbor.pos.nodes,1,function(x){return(ifelse(is.na(edge.scores[paste(V(network)[x]$name,v.name,sep="_")]),
                                                                              edge.scores[paste(v.name,V(network)[x]$name,sep="_")],
                                                                              edge.scores[paste(V(network)[x]$name,v.name,sep="_")]))})
    neighbor.node.score<-node.scores[V(network)[as.vector(neighbor.pos.nodes)]$name]
    
    neighbor.score<-neighbor.edge.score+neighbor.node.score
    
    return(neighbor.score)}
  
  
  
  net.flag <- FALSE
  if (is.null(names(node.scores))) {
    warning("Unnamed node scores")
  }
  if (any(is.na(names(node.scores)))) {
    warning("NA in names of node scores, the node score will be removed")
    node.scores <- node.scores[!is.na(names(node.scores))]
  }
  
  if (is.null(names(edge.scores))) {
    warning("Unnamed edge scores")
  }
  if (any(is.na(names(edge.scores)))) {
    warning("NA in names of edge scores, the edge score will be removed")
    edge.scores <- edge.scores[!is.na(names(edge.scores))]
  }
  
  
  if (is(network, "graphNEL")) {
    network <- igraph.from.graphNEL(network)
    net.flag <- TRUE
  }
  if (is.null(E(network)$name)) {
    E(network)$name <- names(edge.scores)
  }
  
  
  network<-subnetwork.e(network,vid=names(node.scores),eid=names(edge.scores),remove.vertex=F)
  
  
  V(network)$score <- node.scores[V(network)$name]
  E(network)$name<-names(edge.scores)
  E(network)$score<-weightratio.edge.node*edge.scores
  
  #the edges with edge scores>0
  
  pos.edges<- names(edge.scores[which(edge.scores > 0)])
  #identify the names of the nodes connecting the edges  with edge scores>0
  
  pos.nodes <- names(node.scores[which(node.scores > 0)])
  
  
  
  if (length(pos.nodes) == 0 ) {
    warning("No positive nodes")
    module <- graph.empty(n = 0, directed = FALSE)
    if (net.flag) {
      nE <- ecount(module)
      module <- simplify(module, remove.multiple = TRUE, edge.attr.comb="sum")
      if (nE != ecount(module)) {
        warning("Multiple edges between two nodes had to be removed for calculation")
      }
      module <- igraph.to.graphNEL(module)
    }
    return(module)
  }
  
  if (length(pos.nodes) == 1 && length(pos.edges)==0) {
    module <- .subNetwork0(pos.nodes, network)
    if (net.flag) {
      module <- igraph.to.graphNEL(module)
    }
    return(module)
  }
  
  
  
  
  
  if (length(pos.nodes) == 1 && length(pos.edges)>=1) {
    neighbor.score<-neighbor.combined.scores(network,pos.nodes)
    if(max(neighbor.score)<=0){
      module <- .subNetwork0(pos.nodes, network)}
    if(max(neighbor.score)>0){
      neighbor.pos.nodes<-names(neighbor.score)[which.max(neighbor.score)]
      v.id.module<-c(pos.nodes,neighbor.pos.nodes)
      
      
      neighbor.score<-neighbor.score[-which.max(neighbor.score)]
      while(max(neighbor.score)>0){
        neighbor.pos.nodes<-names(neighbor.score)[which.max(neighbor.score)]
        v.id.module<-c(v.id.module,neighbor.pos.nodes)
        neighbor.score<-neighbor.score[-which.max(neighbor.score)]
        
      }
      module <-subnetwork.e(network,vid=v.id.module,eid=E(network)$name,remove.vertex = F )}
    
    if (net.flag) {
      nE <- ecount(module)
      module <- simplify(module, remove.multiple = TRUE, edge.attr.comb="sum")
      if (nE != ecount(module)) {
        warning("Multiple edges between two nodes had to be removed for calculation")
      }
      module <- igraph.to.graphNEL(module)
    }
    return(module)
  }
  
  
  if (length(pos.nodes) > 1 && length(pos.edges)==0) {
    max.pos.node <-pos.nodes[which.max(node.scores[pos.nodes])]
    neighbor.score<-neighbor.combined.scores(network,max.pos.node)
    
    if(max(neighbor.score)<=0){
      module <- .subNetwork0(max.pos.node, network)}
    
    if(max(neighbor.score)>0){
      neighbor.pos.nodes<-names(neighbor.score)[which.max(neighbor.score)]
      v.id.module<-unique(c(max.pos.node,neighbor.pos.nodes))
      
      
      neighbor.score<-neighbor.score[-which.max(neighbor.score)]
      while(max(neighbor.score)>0){
        neighbor.pos.nodes<-names(neighbor.score)[which.max(neighbor.score)]
        v.id.module<-unique(c(v.id.module,neighbor.pos.nodes))
        neighbor.score<-neighbor.score[-which.max(neighbor.score)]
        
      }
      module <-subnetwork.e(network,vid=v.id.module,eid=E(network)$name,remove.vertex = F )}
    
    if (net.flag) {
      nE <- ecount(module)
      module <- simplify(module, remove.multiple = TRUE, edge.attr.comb="sum")
      if (nE != ecount(module)) {
        warning("Multiple edges between two nodes had to be removed for calculation")
      }
      module <- igraph.to.graphNEL(module)
    }
    return(module)
  }
  #if positive nodes do not have any connecting edges in the network, no need to move forward,
  #no connecting nodes 
  
 if( length(intersect(unique(unlist(strsplit(E(network)$name,split="_"))),pos.nodes))==0){
    node.max.score<- names(node.scores)[which.max(node.scores)]
    module <-subnetwork.e(network,vid=node.max.score,eid=E(network)$name,remove.vertex = F )
    if (net.flag) {
      nE <- ecount(module)
      module <- simplify(module, remove.multiple = TRUE, edge.attr.comb="sum")
      if (nE != ecount(module)) {
        warning("Multiple edges between two nodes had to be removed for calculation")
      }
      module <- igraph.to.graphNEL(module)
    }
    return(module)
 }
  
  
  #find positive scoring nodes with the positive edges connecting them
  pos.subgraph <- subnetwork.e(network,pos.nodes,pos.edges,remove.vertex=F )
  
  #identify components of that subgraph via some algorithm (not clear to me how it works)
  #I believe the idea is that there may be components of the positive subgraph that are not
  #actually connected. These non connected clusters are the 'graph components'
  #ultimately yields list of clusters of graph nodes, with each list element a separate 'component'
  conn.comp.graph <- decompose.graph(pos.subgraph)
  
  
  #cluster.degree.within<-unlist(lapply(lapply(conn.comp.graph,function(x) igraph::degree(x)),sum))
  
  
  #cluster.degree.tot<-unlist(lapply(lapply(conn.comp.graph,function(x) igraph::degree(subnetwork.e(network,V(x)$name,E(network)$name,remove.vertex=T ))[V(x)$name]),sum))
  
  #cluster.degree<-cluster.degree.tot-cluster.degree.within
  #names(cluster.degree)<-new.names
  
  
  
  
  #get the sum of the node scores of all component elements
  node.score.comp <- unlist(lapply(lapply(conn.comp.graph, vertex_attr,
                                          "score"), sum))
  #get the sum of the edge scores of all component elements
  edge.score.comp <- unlist(lapply(lapply(conn.comp.graph, edge_attr,
                                          "score"), sum))
  #sum of total node scores and edge scores in each component
  score.comp<-edge.score.comp+node.score.comp
  
  #re-order the positive scoring clusters, and the sum of the cluster scores,
  #such that the highest scores appear first
  comp.graph <- score.comp [order(score.comp, decreasing = TRUE)]
  score.comp <- sort(score.comp, TRUE)
  
  #assign cumulative scores to each cluster
  for (i in 1:length(conn.comp.graph)) {
    conn.comp.graph[[i]]$score <- score.comp[i]
  }
  
  #now working with the network input to the function (not just the positive scoring nodes, as above)
  #get vertices in the entire network
  v.id <- seq(1, vcount(network))
  names(v.id) <- V(network)$name
  
  #get edges in network
  edgelist <- get.edgelist(network, FALSE)
  edgelist1 <- edgelist[, 1]
  edgelist2 <- edgelist[, 2]
  #now loop through and assign the same id to edge list elements
  #in the same cluster of positively scoring nodes
  for (i in 1:length(conn.comp.graph)) {
    new.id <- length(V(network)) + i
    for (j in as.character(v.id[V(conn.comp.graph[[i]])$name])) {
      edgelist1[which(edgelist1 == j)] <- new.id
      edgelist2[which(edgelist2 == j)] <- new.id
    }
  }
  
  #also give more specific names to those ids as cluster1, cluster2, etc.
  new.ids <- seq(length(V(network)) + 1, length(V(network)) +
                   length(conn.comp.graph))
  new.names <- paste("cluster", seq(1:length(conn.comp.graph)),
                     sep = "")
  names(new.ids) <- new.names
  
  #put together an object of all vertex ids, including cluster ids, and their names
  v.id <- c(v.id, new.ids)
  v.name <- names(v.id)
  names(v.name) <- v.id
  #define the combined scores for each positive scoring cluster
  cluster.score<-unlist(lapply(conn.comp.graph,graph_attr,"score"))
  names(cluster.score)<-new.names
  #put the node scores and the positive scoring cluster combined scores together
  v.score<-c(V(network)$score,cluster.score)
  names(v.score)<-c(V(network)$name,names(cluster.score))
  
  #now put together a new edgelist, with the negatively scoring nodes given unique ids,
  #and the positively scoring nodes still clustered as above
  new.edgelist <- cbind(v.name[as.character(edgelist1)], v.name[as.character(edgelist2)])
  
  
  #convert the edgelist to a graph
  interactome2 <- graph.edgelist(new.edgelist, FALSE)
  #combine the edge scores for the edges where the nodes within one cluster connected to other nodes ???
  
  
  E(interactome2)$score<-weightratio.edge.node*edge.scores
  #end part 1 identify the clusters of positive scoring nodes with positive scoring connecting edges
  #reduce the network to the clusters connecting with other nodes in the previous network
  #########################################################################################################
  
  #start part 2-define the edge weights for the network-interactome2
  ############################################################################################################################
  #now initialize weights for each edge
  E(interactome2)$weight <- rep(0, length(E(interactome2)))
  
  #remove any edges that start and end at the same node, and remove duplicate edges
  #This will leave only one member of each of the positively scoring clusters, with the cumulative score
  #for that cluster assigned to it, and the rest of the negatively scoring nodes
  #Beacuse all members of a cluster were assigned the same id and score, these all appeared as self loops
  interactome2 <- simplify(interactome2, remove.multiple = TRUE, edge.attr.comb = "sum")
  edge.name.inter.pos<-paste(get.edgelist(interactome2,T)[,1],get.edgelist(interactome2,T)[,2],sep="_")[E(interactome2)$score>0]
  E(interactome2)$name<-paste(get.edgelist(interactome2,T)[,1],get.edgelist(interactome2,T)[,2],sep="_")
  
  #define two types of degree of nodes for both meta-nodes(positive scoring nodes) and negative scoring nodes
  #The degree and newdegree of clusters is based on the subnetwork of network with edge score positive
  #to mimic the possible reality network with significant edge connectin from multiple testing of pairwise correlations
  cluster.degree.within<-unlist(lapply(lapply(conn.comp.graph,function(x) igraph::degree(x)),sum))
  
  
  cluster.degree.tot<-unlist(lapply(lapply(conn.comp.graph,function(x) igraph::degree(subnetwork.e(network,V(x)$name,pos.edges,remove.vertex=F ))[V(x)$name]),sum))
  
  cluster.degree.inter<-cluster.degree.tot-cluster.degree.within+1
  names(cluster.degree.inter)<-new.names
  
  
  cluster.newdegree.inter<-max(cluster.degree.inter)-cluster.degree.inter+1
  
  
  #find positive scoring nodes with the positive edges connecting them
  network.pos <- subnetwork.e(network,pos.nodes,pos.edges,remove.vertex=F )
  
  #make sure any least one edges in network.pos,if the nodes with positive scores not in any edges
  #meaningless to generate MST by minimizing the edge weights
  length(E(network.pos))
  
  
  
  #The degree and newdegree of clusters is based on the subnetwork of interactome2 with edge score positive
  #to mimic the possible reality network with significant edge connections from multiple testing of pairwise correlations
  negnode.name<-V(interactome2)$name[!grepl("cluster",V(interactome2)$name)]
  negnode.degree.inter<-igraph::degree(network)[negnode.name]+1
  negnode.newdegree.inter<-max(negnode.degree.inter)-negnode.degree.inter+2
  
  edgelist.name.inter<-paste(get.edgelist(interactome2, TRUE)[,1],get.edgelist(interactome2, TRUE)[,2],sep="_")
  edge.score.inter<-E(interactome2)$score
  names(edge.score.inter)<-edgelist.name.inter
  edgelist.name.inter.1<-get.edgelist(interactome2, T)[,1]
  edgelist.name.inter.2<-get.edgelist(interactome2, T)[,2]
  edgelist.name.inter.mat<-cbind(edgelist.name.inter.1,edgelist.name.inter.2)
  #the function to define the weights for the edges in interactome2-intuitive definition to set positive edge scores zero
  weight.edge.degree<-function(x){
    edge.degree.name<-paste(x[1],x[2],sep="_")
    if(grepl("cluster",x[1]) && grepl("cluster",x[2])){
      
      edge.score.degree<-ifelse(edge.score.inter[edge.degree.name]>=0,0,
                                edge.score.inter[edge.degree.name]/(mean(c(cluster.degree.inter[x[1]],cluster.degree.inter[x[2]]))))
      
      edge.degree<--edge.score.degree
      
    }
    if(grepl("cluster",x[1]) && !grepl("cluster",x[2])){
      
      edge.score.degree<-ifelse(edge.score.inter[edge.degree.name]>=0,0,
                                edge.score.inter[edge.degree.name]/(mean(c(cluster.degree.inter[x[1]],negnode.degree.inter[x[2]]))))
      
      
      edge.degree<--(node.scores[x[2]]/negnode.degree.inter[x[2]]+edge.score.degree)
      
      
    }
    if(!grepl("cluster",x[1]) && grepl("cluster",x[2])){
      
      edge.score.degree<-ifelse(edge.score.inter[edge.degree.name]>=0,0,
                                edge.score.inter[edge.degree.name]/(mean(c(cluster.degree.inter[x[2]],negnode.degree.inter[x[1]]))))
      
      
      edge.degree<--(node.scores[x[1]]/negnode.degree.inter[x[1]]+edge.score.degree)
      
    }
    
    if(!grepl("cluster",x[1]) && !grepl("cluster",x[2])){
      
      edge.score.degree<-ifelse(edge.score.inter[edge.degree.name]>=0,0,
                                edge.score.inter[edge.degree.name]/(mean(c(negnode.degree.inter[x[2]],negnode.degree.inter[x[1]]))))
      
      
      edge.degree<--(node.scores[x[1]]/negnode.degree.inter[x[1]]+node.scores[x[2]]/negnode.degree.inter[x[2]]+weightratio.edge.node*edge.score.degree)
      
    }
    
    return(edge.degree)}
  
  ##################################################################################################################################
  #the function to define the weights for the edges in interactome2-use exponential negative edge scores to set edge weights
  #normalize the exponential weights by divided the possible sum of maximum nodes score and edge score
  max.node.scores<-max(cluster.score,na.rm=T)
  max.edge.scores<-max(edge.score.inter,na.rm=T)
  sum.max.nodeedge.scores<-ifelse((2*max.node.scores+max.edge.scores)>0,(2*max.node.scores+max.edge.scores)/max.edge.scores,-(2*max.node.scores+max.edge.scores)/max.edge.scores)
  weight.edge.degree.exp<-function(x){
    edge.degree.name<-paste(x[1],x[2],sep="_")
    if(grepl("cluster",x[1]) && grepl("cluster",x[2])){
      edge.score.degree<-ifelse(edge.score.inter[edge.degree.name]>=0,edge.score.inter[edge.degree.name]/(mean(c(cluster.newdegree.inter[x[1]],cluster.newdegree.inter[x[2]]))),
                                edge.score.inter[edge.degree.name]/(mean(c(cluster.degree.inter[x[1]],cluster.degree.inter[x[2]]))))
      edge.degree<-exp((-edge.score.degree-cluster.score[x[1]]/cluster.newdegree.inter[x[1]]-cluster.score[x[2]]/cluster.newdegree.inter[x[2]])/sum.max.nodeedge.scores)
      
    }
    if(grepl("cluster",x[1]) && !grepl("cluster",x[2])){
      edge.score.degree<-ifelse(edge.score.inter[edge.degree.name]>=0,edge.score.inter[edge.degree.name]/(mean(c(cluster.newdegree.inter[x[1]],negnode.newdegree.inter[x[2]]))),edge.score.inter[edge.degree.name]/(mean(c(cluster.degree.inter[x[1]],negnode.degree.inter[x[2]]))))
      edge.degree<-exp((-cluster.score[x[1]]/cluster.newdegree.inter[x[1]]-node.scores[x[2]]/negnode.degree.inter[x[2]]-edge.score.degree)/sum.max.nodeedge.scores)
    }
    if(!grepl("cluster",x[1]) && grepl("cluster",x[2])){
      
      edge.score.degree<-ifelse(edge.score.inter[edge.degree.name]>=0,edge.score.inter[edge.degree.name]/(mean(c(cluster.newdegree.inter[x[2]],negnode.newdegree.inter[x[1]]))),
                                edge.score.inter[edge.degree.name]/(mean(c(cluster.degree.inter[x[2]],negnode.degree.inter[x[1]]))))
      
      edge.degree<-exp((-cluster.score[x[2]]/cluster.newdegree.inter[x[2]]-node.scores[x[1]]/negnode.degree.inter[x[1]]-edge.score.degree)/sum.max.nodeedge.scores)
      
    }
    
    if(!grepl("cluster",x[1]) && !grepl("cluster",x[2])){
      
      edge.score.degree<-ifelse(edge.score.inter[edge.degree.name]>=0,edge.score.inter[edge.degree.name]/(mean(c(negnode.newdegree.inter[x[2]],negnode.newdegree.inter[x[1]]))),
                                edge.score.inter[edge.degree.name]/(mean(c(negnode.degree.inter[x[2]],negnode.degree.inter[x[1]]))))
      edge.degree<-exp(-(node.scores[x[1]]/negnode.degree.inter[x[1]]+node.scores[x[2]]/negnode.degree.inter[x[2]]+edge.score.degree)/sum.max.nodeedge.scores)
      
    }
    
    return(edge.degree)}
  #intuitive edge weights to set the positive edge scores zero
  # E(interactome2)$weight<-apply( edgelist.name.inter.mat,1, weight.edge.degree)
  #exponentialized edge weights to utilize the positive edge scores
  E(interactome2)$weight<-apply( edgelist.name.inter.mat,1, weight.edge.degree.exp)
  #end of part 2
  #weight.inter<-E(interactome2)$weight
  #if(all(weight.inter>0) && all(weight.inter<Inf)){
  #weight.inter.log.new<--log(weight.inter[weight.inter<1])
  #weight.inter[weight.inter<1]<-(max(weight.inter.log.new)-weight.inter.log.new+1)/sum(weight.inter.log.new)
  #E(interactome2)$weight<-weight.inter}
  #if(all(weight.inter>0)==T && all(weight.inter<Inf)==F){
  #which.inf<-weight.inter==Inf
  #max.weight.inter<-max(weight.inter[!which.inf])
  #if(max.weight.inter>1){
  #weight.inter[which.inf]<-exp(log(max.weight.inter)+runif(1,0,(log(max.weight.inter))))
  #}
  #if(max.weight.inter==1){
  #weight.inter[which.inf]<-exp(log(max.weight.inter)+runif(1,0,0.01))
  #}
  #if(max.weight.inter<1){
  #weight.inter[which.inf]<-1
  #}}
  
  #if( all(weight.inter<Inf)==F){
  #which.inf<-weight.inter==Inf
  #max.weight.inter<-max(weight.inter[!which.inf])
  #if(max.weight.inter>1){
  #weight.inter[which.inf]<-exp(log(max.weight.inter)+runif(1,0,(log(max.weight.inter))))
  #}
  #if(max.weight.inter==1){
  #weight.inter[which.inf]<-exp(log(max.weight.inter)+runif(1,0,0.01))
  #}
  #if(max.weight.inter<1){
  #weight.inter[which.inf]<-1
  #}}
  #if(all(weight.inter>0)==F ){
  #which.zero<-weight.inter==0
  #min.weight.inter<-min(weight.inter[!which.zero])
  #if(min.weight.inter>1){
  #weight.inter[which.zero]<-exp(log(min.weight.inter)-runif(1,0,log(min.weight.inter)))
  #}
  #if(min.weight.inter==1){
  #weight.inter[which.zero]<-exp(log(min.weight.inter)-runif(1,0,0.1))
  #}
  #if(min.weight.inter<1){
  #weight.inter[which.zero]<-exp(log(min.weight.inter)+runif(1,log(min.weight.inter),0))
  #}}
  
  #E(interactome2)$weight<-weight.inter
  
  ########################################################################################################################################################
  ########################################################################################################################################################
  #start of part 3
  #now grab the node scores for all nodes in the graph-interactome2
  node.score <- v.score[V(interactome2)$name]
  
  #also make sure the elements of the vector are named equivalently to the graph-interactome2
  names(node.score) <- V(interactome2)$name
  
  #next, assign the correct scores to the positively scoring nodes
  #recall each one is represented by a single member of the cluster,
  #and the scores are currently set to zero
  #so to get the proper scores, go back to the conn.comp.graph
  #object and grab them for each cluster
  node.score.cluster <- sapply(conn.comp.graph, get.graph.attribute,
                               "score")
  names(node.score.cluster) <- new.names
  
  #now assign those scores to the proper elements of the 'node.score' vector
  node.score[grep("cluster", names(node.score))] <- node.score.cluster[names(node.score[grep("cluster",
                                                                                             names(node.score))])]
  
  #identify components of the current graph, composed of negatively scoring elements and single
  #nodes representing positively scoring clusters
  decomp.graphs <- decompose.graph(interactome2)
  
  #compute sum of the scores for the positively scoring elements of each
  #of the sub-clusters of the graph
  sum.pos <- lapply(decomp.graphs, function(x) {
    sum(node.score[names(which(node.score[V(x)$name] > 0))])
  })
  
  
  #now restrict the graph to the cluster of nodes with the highest scoring positive elements
  #note here that the which.max() command is used to return the index of the maximum value
  interactome2 <- decomp.graphs[[which.max(sum.pos)]]
  rm(decomp.graphs)
  
  #now find a minimum spanning tree for the restricted graph
  #this will likely produce a still very large graph
  mst <- minimum.spanning.tree(interactome2, weights = E(interactome2)$weight)
  
  V(mst)$name
  #end of part 3-identify the minimum spanning tree of the component with largest positive scoring meta-nodes
  #########################################################################################################################################################
  ########################################################################################################################################################
  #start of part 4 Identify the shortest path of negative nodes connected to positive (meta-)nodes
  #or possibly the two positive scoring (meta)-nodes connecting together with negative edge scores
  #identify the indices of the positively scoring clusters in the mst
  mst.cluster.id <- grep("cluster", V(mst)$name)
  
  #assign names to those indices corresponding to cluster id
  names(mst.cluster.id) <- V(mst)[mst.cluster.id]$name
  
  
  #re-order by cluster number
  mst.cluster.id <- mst.cluster.id[order(as.numeric(matrix(unlist(strsplit(names(mst.cluster.id),
                                                                           "cluster")), nrow = 2)[2, ]))]
  
  #initiate vector to find ids of all the negative scoring nodes
  # in the shortest paths connecting the positive scoring (meta-) nodes(clusters)
  all.ids <- c()
  if(length(mst.cluster.id) == 1){
    neg.node.ids.2 = c()
    
  }else
  {
    
    
    #shortest_paths(mst, from = V(mst)$name[mst.cluster.id[4]], to = V(mst)$name[mst.cluster.id[(4 + 1)]])
    
    
    #plotModule(subnetwork.e(mst,vid=V(mst)$name[mst.cluster.id],eid=E(mst)$name,remove.vertex=T))
    
    
    #this loops through the positively scoring nodes and finds shortest paths to all other
    #positively scoring nodes
    for (j in 1:(length(mst.cluster.id) - 1)) {
      path <- unlist(shortest_paths(mst, from =mst.cluster.id[j],
                                    to =mst.cluster.id[(j + 1):length(mst.cluster.id)])$vpath)
      path<-unique(path)
      all.ids <- c(all.ids, path)
    }
    
    
    
    #this somehow identifies all nodes related to shortest paths in the graph
    all.ids <- unique(all.ids)
    #now this subsets the graph to include only those ids involved in shortest paths
    sub.interactome2 <- subnetwork.e(interactome2,vid=V(mst)[all.ids]$name, eid= E(mst)$name,remove.vertex=T)
    #plotModule(sub.interactome2)
    #now identify negative scoring nodes in the graph
    neg.node.ids <- which(node.score[V(sub.interactome2)$name] < 0)
    
    ##identify the cluster id connecting directly with negative connecting edges scores
    #for each of the positive scoring cluster, and possible connecting clusters , determine if a positively scoring
    #cluster lies adjacent with negative connecting edge scores, and store the cluster ids in a list
    cluster.id<- grep("cluster",V(sub.interactome2)$name)
    for (i in cluster.id) {
      V(sub.interactome2)[i]$clusters <- list(neighbors(sub.interactome2,
                                                        v = i)[grep("cluster", V(sub.interactome2)[neighbors(sub.interactome2, v = i)]$name)])
    }
    
    V(sub.interactome2)$clusters
    ########################################################################################################################################################
    #for each of the negatively scoring nodes and possible connecting clusters with negative connecting edge scores, determine if a positively scoring
    #cluster lies adjacent, and store the cluster ids in a list
    if( length(neg.node.ids)!=0){
      for (i in neg.node.ids) {
        V(sub.interactome2)[i]$clusters <- list(neighbors(sub.interactome2, v = i)[grep("cluster", V(sub.interactome2)[neighbors(sub.interactome2, v = i)]$name)])
      }}
    #initiate vector for negative node scores
    score.neg.nodes <- c()
    if(length(neg.node.ids)!=0){
      #for each negative scoring node, if a positive scoring cluster is adjacent, sum the score
      #for that node and the positive cluster score
      #Otherwise, just assign the score of the node
      for (i in neg.node.ids) {
        if (length(V(sub.interactome2)[i]$clusters[[1]]) > 0)
        {edge.score.negnodes<-c()
        for (j in 1:length(V(sub.interactome2)[i]$clusters[[1]])){
          
          edge.score.test<-edge.score.inter[paste(V(sub.interactome2)$name[i],
                                                  
                                                  V(sub.interactome2)[i]$clusters[[1]]$name[j],sep="_")]
          
          edge.score.negnodes<-c(edge.score.negnodes,
                                 ifelse(is.na(edge.score.test),
                                        edge.score.inter[paste(V(sub.interactome2)[i]$clusters[[1]]$name[j],
                                                               V(sub.interactome2)$name[i],sep="_")],edge.score.test))}
        score.neg.nodes <- c(score.neg.nodes, sum(node.score[V(sub.interactome2)[c(i,
                                                                                   V(sub.interactome2)[i]$clusters[[1]])]$name],
                                                  edge.score.negnodes))
        }
        else {
          score.neg.nodes <- c(score.neg.nodes, node.score[V(sub.interactome2)[i]$name])
        }
      }
      names(score.neg.nodes)<-V(sub.interactome2)$name[neg.node.ids]}
    #########################################################################################################################################################
    #initiate vector for possibly connecting positive scoring (meta-) nodes(clusters) connecting with negative scoring edges
    score.clusters.conn<-c()
    if(length(cluster.id)!=0){
      #for each positive scoring (meta-) node(cluster), if a positive scoring cluster is adjacent,
      # sum the scores of two connecting positive scoring (meta-) nodes(clusters) and the connecting edge score
      #Otherwise, just assign NA
      for (i in cluster.id) {
        if (length(V(sub.interactome2)[i]$clusters[[1]]) > 0)
        { edge.score.clusters.conn<-c()
        for (j in 1:length(V(sub.interactome2)[i]$clusters[[1]])){
          
          edge.score.test<-edge.score.inter[paste(V(sub.interactome2)$name[i],
                                                  
                                                  V(sub.interactome2)[i]$clusters[[1]]$name[j],sep="_")]
          
          edge.score.clusters.conn<-c(edge.score.clusters.conn,
                                      ifelse(is.na(edge.score.test),
                                             edge.score.inter[paste(V(sub.interactome2)[i]$clusters[[1]]$name[j],
                                                                    V(sub.interactome2)$name[i],sep="_")],edge.score.test))
        }
        
        score.clusters.conn <- c( score.clusters.conn, sum(cluster.score[as.numeric(unlist(strsplit(V(sub.interactome2)$name[i],"cluster"))[2])],
                                                           cluster.score[as.numeric(unlist(strsplit(V(sub.interactome2)[i]$clusters[[1]]$name,
                                                                                                    "cluster"))[seq(2,2*length(V(sub.interactome2)[i]$clusters[[1]]),by=2)])],
                                                           edge.score.clusters.conn))
        }
        else {
          score.clusters.conn  <- c(score.clusters.conn , NA)
        }
      }
      names(score.clusters.conn )<-V(sub.interactome2)$name[cluster.id] }
    
    #########################################################################################################################################################
    score.conn<-c(na.omit(score.clusters.conn),score.neg.nodes)
    #id any negative nodes attached to positive clusters, whose sum total score is positive
    #id any postive scoring clusters attached to other positive cluster,where sum of the two cluster scores and the negative connecting edge score is positive
    neg.node.ids.2 <- names(score.conn)[score.conn > 0]
  }
  
  
  
  
  
  
  #if there are no negative nodes that connect to positive clusters such that the sum total
  #score is positive, then grab the highest scoring cluster and assign it to object 'module'
  #need to check the combined neighboring negative nodes and the connecting edge scores
  if (length(neg.node.ids.2) == 0) {
    cluster.max.nodes<-unlist(lapply(conn.comp.graph, vertex_attr,
                                     
                                     "name")[as.numeric(matrix(unlist(strsplit(names(cluster.score)[which.max(cluster.score)],
                                                                               "cluster")), nrow = 2)[2, ])])
    #the term in the square brackets is a very complicated way of getting the number of the cluster with
    #the highest score
    V(mst)$score<-v.score[V(mst)$name]
    max.cluster<-names(cluster.score)[which.max(cluster.score)]
    #under some special situation, clusters with minimum spanning tree might not the same as the cluster with highest cluster score
    if (names(mst.cluster.id)!=max.cluster){
      neighbor.score<-neighbor.combined.scores(network,cluster.max.nodes)
    }else{
      neighbor.score<-neighbor.combined.scores(mst,max.cluster)
      
    }
    
    if (length(neighbor.score) >0 ){
      
    if(max(neighbor.score)<=0){
      module <- subnetwork.e(network,vid=cluster.max.nodes,eid=names(edge.scores),remove.vertex=T )}
    if(max(neighbor.score)>0){
      neighbor.pos.nodes<-names(neighbor.score)[which.max(neighbor.score)]
      
      neighbor.pos.nodes<-neighbor.pos.nodes[!grepl("cluster",neighbor.pos.nodes)]
      
      
      v.id.module<-c(cluster.max.nodes,neighbor.pos.nodes)
      module <-subnetwork.e(network,vid=v.id.module,eid=names(edge.scores))}
    
    } else{
      module <- subnetwork.e(network,vid=cluster.max.nodes,eid=names(edge.scores),remove.vertex=F )
    }
    
    
    if (net.flag) {
      nE <- ecount(module)
      module <- simplify(module, remove.multiple = TRUE, edge.attr.comb="sum")
      if (nE != ecount(module)) {
        warning("Multiple edges between two nodes had to be removed for calculation")
      }
      module <- igraph.to.graphNEL(module)
    }
    
    return(module)
  }
  
  #Otherwise, if there are some negative nodes attached to positive scoring clusters, such that
  #the total combined score is positive, or any negative edges connecting two positive scoring clusters
  #such as the total combined node scores and edge score is positive
  #identify the largest component of the graph containing
  #those negative nodes or negative connecting edges
  subg <- largestComp(induced.subgraph(sub.interactome2, neg.node.ids.2))
  
  if(length(vcount(subg))>2){
    #find a minimum spanning tree for that graph
    mst.subg <- minimum.spanning.tree(subg, E(subg)$weight)
    ###############################################################################################################################
    
    
    #initialize maximum score
    max.score <- 0
    
    #initialize best path vector
    best.path <- c()
    
    #loop through the nodes in the mst subgraph
    for (i in 1:(length(V(mst.subg)))) {
      
      #identify shortest paths from all nodes
      path <- get.all.shortest.paths(mst.subg, from = V(mst.subg)[i])
      
      #sum the scores of the negative nodes and the adjacent positive clusters and their connectng edge scores
      #also sum the node scores and the negative edge score of two positive scoring clusters connecting with negative scoring edges
      path.score <- unlist(lapply(path$res, getPathScore.e,
                                  graph1 = mst.subg, graph2 = sub.interactome2,
                                  node.score = v.score,edge.score=edge.score.inter))
      #get the index of the highest scoring path
      best.pos <- which.max(path.score)
      
      #if the score associated with that path is higher than the current
      #high score, re-assign the high score and best path
      if (path.score[[best.pos]] > max.score) {
        best.path <- path$res[[best.pos]]
        max.score <- path.score[[best.pos]]
      }
    }
    ###############################################################################################################################
    if (length(best.path) != 1) {
      
      #for each of the elements in the highest scoring path, determine the adjacent clusters
      cluster.list <- V(mst.subg)[best.path]$clusters
      names.list <- as.character(1:length(cluster.list))
      names(cluster.list) <- names.list
      names(best.path) <- names.list
      
      #convoluted way of checking for overlap among clusters and filtering some out for some reason (I believe having to do with repeat clusters)
      for (i in names.list) {
        res <- lapply(cluster.list, intersect, cluster.list[[i]])
        if (length(intersect(unlist(cluster.list[as.character(which(as.numeric(names.list) <
                                                                    as.numeric(i)))]), unlist(cluster.list[as.character(which(as.numeric(names.list) >
                                                                                                                              as.numeric(i)))]))) > 0) {
          if (length(setdiff(res[[i]], unique(unlist(res[names(res) !=
                                                         i])))) == 0) {
            cluster.list <- cluster.list[names(cluster.list) !=
                                           i]
            names.list <- names.list[names.list != i]
          }
        }
      }
      
      
      #id vertices numbers in final best path
      best.path <- best.path[names.list]
    }
    
    
    
    
    #get names of negative scoring vertices in final best path
    module <- V(mst.subg)[best.path]$name
    
    #get names of adjacent positive scoring vertex clusters
    pos.cluster <- V(sub.interactome2)[unique(unlist(V(mst.subg)[best.path]$clusters))]$name
    
    pos.cluster.final<-unique(union(module[grepl("cluster",module)],pos.cluster))
    module<-module[!grepl("cluster",module)]
    #combine into one character vector of the final vertex names
    module <- c(module, unlist(lapply(conn.comp.graph, vertex_attr,
                                      "name")[as.numeric(matrix(unlist(strsplit(pos.cluster.final,
                                                                                "cluster")), nrow = 2)[2, ])]))
    #make final subset igraph object
    module <- subnetwork.e(network,vid=module,eid=E(network)$name,remove.vertex=F)}
  
  if(length(vcount(subg))<=2){
    
    
    
    cluster.max.nodes<-unlist(lapply(conn.comp.graph, vertex_attr,
                                     
                                     "name")[as.numeric(matrix(unlist(strsplit(names(cluster.score)[which.max(edge.score.comp+node.score.comp)],
                                                                               "cluster")), nrow = 2)[2, ])])
    #the term in the square brackets is a very complicated way of getting the number of the cluster with
    #the highest score
    V(mst)$score<-v.score[V(mst)$name]
    max.cluster<-names(cluster.score)[which.max(cluster.score)]
    neighbor.score<-neighbor.combined.scores(mst,max.cluster)
    if(max(neighbor.score)<=0){
      module <- subnetwork.e(network,vid=cluster.max.nodes,eid=names(edge.scores),remove.vertex=F )}
    if(max(neighbor.score)>0){
      neighbor.pos.nodes<-names(neighbor.score)[which.max(neighbor.score)]
      
      neighbor.pos.nodes<-neighbor.pos.nodes[!grepl("cluster",neighbor.pos.nodes)]
      
      
      v.id.module<-c(cluster.max.nodes,neighbor.pos.nodes)
      neighbor.score<-neighbor.score[-which.max(neighbor.score)]
      
      while(max(neighbor.score)>0){
        neighbor.pos.nodes<-names(neighbor.score)[which.max(neighbor.score)]
        
        neighbor.pos.nodes<-neighbor.pos.nodes[!grepl("cluster",neighbor.pos.nodes)]
        
        
        v.id.module<-c(v.id.module,neighbor.pos.nodes)
        neighbor.score<-neighbor.score[-which.max(neighbor.score)]
        
        
      }
      module <-subnetwork.e(network,vid=v.id.module,eid=names(edge.scores))}
    
  }
  
  
  #plotModule(module)
  if (net.flag) {
    nE <- ecount(module)
    module <- simplify(module, remove.multiple = TRUE,edge.attr.comb="sum")
    if (nE != ecount(module)) {
      warning("Multiple edges between two nodes had to be removed for calculation")
    }
    
    #covnert from igraph to graphNEL
    module <- igraph.to.graphNEL(module)
  }
  return(module)
  
  
  
} 
  
