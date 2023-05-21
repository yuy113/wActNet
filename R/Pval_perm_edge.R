#' @title Generate P values by multiple permulation tests on all the possible edges
#'
#'
#' @description This function calculate P values in multiple testing to check associations
#' among any two variables(nodes) in case control studies
#'  (biological network) by comparing observed Pearson's correlations
#' with those in permutated samples through nonparametric permutation tests.
#' @param dat  The observations of p variables in (matched) case control studies
#'  in matrix or dataframe format. The names of the variables must be specified.
#' @param nsim The number of permutations for observations for each variable
#'   in the dataset-dat, nsim set to a large number such as 10000
#' @param do.parallel Indicator variable(T-TRUE/F-FALSE) whether useing parallel computing
#'   with default value-F,not using parallel computing
#' @param no_cores The number of computing units(cores) for parallel computing
#' with default value-NULL when do.parallel=F,must specify the number is
#' do.parallel=T, if PC or laptop can use the value-detectCores() - 1 or smaller
#' @param MatchId A vector representing matched id for cases and controls.
#'   If the study doesn't have matching, then set to NULL, if matched studies,
#'   MatchId can't contain any missing values
#' @return A vector of p-values for Pearson's correlations
#' between any two variables, with corresponding edge names connecting with '_'
#' @examples
#' dat1<-matrix(rnorm(4000),ncol=40,nrow=100)
#' colnames(dat1)<-paste("Var", as.character(1:40),sep="")
#' pval.edge<-pval.perm.corr(dat1,nsim=1000,MatchId=NULL,do.parallel=FALSE)
#'
#' @export
#'
#'
pval.perm.corr<-function(dat,nsim,MatchId=NULL,do.parallel=FALSE,no_cores=NULL){
  if(do.parallel==T && is.null(no_cores)){
    warnings("Please specify the number of cores")
    break
  }

  #permutation test
  #sample size(the number of observations-n
  n<-dim(dat)[1]
  #the number of variales for correlatins
  p<-dim(dat)[2]

  n.submat<-ceiling(p/20)

  #some simple function used for p-values for permutation test sampling of correlations

  ##########################################################################
  #cross correlations among permutated samples for two data observations
  #input 1: dat1 observations for one set of p1 variables in dataframe or matrix format
  #input 2: dat2 observations for another set of p2 variables in dataframe or matrix format
  #MatchId: pair id group labeling for matched studies default-NA without matching
  #output: the p1*p2 cross correlations between the variables in permutated dat1
  #and those in permutated dat2
  permut_cross_corr<-function(dat1,dat2,MatchId=NULL){
    n<-nrow(dat1)
    dat<-cbind(dat1,dat2)
    #permuation random samples for further correlation
    if(is.null(MatchId)){
      var.perm<-apply(dat,2,function(x){sample(x,size=n)})}
    if(!is.null(MatchId)){
      size.b<-length(unique(MatchId))
      var.perm<-apply(dat,2,function(x){x[rep(sample(unique(MatchId),size=size.b),each=n/size.b)]})
    }
    corr.perm<-cor(var.perm[,1:ncol(dat1)],var.perm[,(1+ncol(dat1)):ncol(dat)], use="complete.obs")


    as.vector(corr.perm)

  }
  ##########################################################################

  ##########################################################################
  #correlations among the variables among permutated samples in the data-dat
  #input 1: dat observations for one set of p variables in dataframe or matrix format
  #MatchId: pair id group labeling for matched studies default-NA without matching
  #output: the p*(p-1)/2 correlations among two of all the variables in permutated dat
  permut_within_corr<-function(dat,MatchId=NULL){
    n<-nrow(dat)
    #permuation random samples for further correlation
    if(is.null(MatchId)){
      var.perm<-apply(dat,2,function(x){sample(x,size=n)})}
    if(!is.null(MatchId)){
      size.b<-length(unique(MatchId))
      var.perm<-apply(dat,2,function(x){x[rep(sample(unique(MatchId),size=size.b),each=n/size.b)]})
    }
    corr.perm<-cor(var.perm, use="complete.obs")


    corr.perm[lower.tri(corr.perm, diag = FALSE)]

  }
  ##########################################################################




  ##########################################################################
  #the function-pval.perm.corr.cross to calculate the P-values of the correlations=0
  #of the variables between the two input datasets-dat1 and dat2
  #based on permuation test sampling assuming these variables are pairwisely independent
  #equivalently the correlation=0
  #input: dat1,dat2-the two datasets containing the samples of the variables in columns
  #input:nsim, the number of sample size for permutation test sampling
  #input:MatchId-block paired id variable for paired study
  #input: do.parallel, the option whether parallel computing method=T or F
  #do.parallel=F, do not use parallel computing,T, use parallel computing
  #assuming the paired observation at each pairid group has the same size
  #MatchId can't have any nonmissing observations
  #output:the vector of p-values based on permuation resampling samples with size-p1*p2
  #where p1-the number of variables(columns) in dat1;p2-the number of variables(columns) in dat2;
  pval.perm.corr.cross<-function(dat1,dat2,nsim,MatchId=NULL){
    n<-nrow(dat1)
    if(any(is.na(colnames(dat1)))){
      warning("Names of variables in data 1 contain missing values")
    }
    if(is.null(colnames(dat1))){
      warning("Names of variables in data 1 are missing, please specify")
    }
    if(!is.null(colnames(dat1))){
    if(any(is.na(colnames(dat2)))){
      warning("Names of variables in data 2 contain missing values")
    }}
    if(is.null(colnames(dat2))){
      warning("Names of variables in data 2 are missing, please specify")
    }

    if(!is.null(MatchId)){
    if( any(is.na(MatchId)))
    {warning("NA in paired group id")
      break}
    if(!any(is.na(MatchId))){
      size.b<-length(unique(MatchId))
    }}
    #the number of variables for correlatins
    p1<-dim(dat1)[2]
    p2<-dim(dat2)[2]
    #the number of correlations from p variables
    pp<-p1*p2
    mat.corr.perm<-matrix(NA,ncol=pp,nrow=nsim)
    dat<-cbind(dat1,dat2)
    for ( i in 1:nsim){

      mat.corr.perm[i,]<-permut_cross_corr(dat1,dat2,MatchId=MatchId)

    }


    #Pearson's correlation from real data
    corr.true<-as.vector(cor(dat1,dat2, use="complete.obs"))
    #calculate the P-values from the permutation samples with the assumptions of indepedent variables in dataset-dat
    p.val.perm<-sapply(1:ncol(mat.corr.perm),function(x){(sum(abs(mat.corr.perm[,x])>=abs(corr.true[x]))+1)/(nsim+1)})
    from<-rep(colnames(dat1),ncol(dat2))
    to<-rep(colnames(dat2),each=ncol(dat1))
    names(p.val.perm)<-paste(from,to,sep="_")
    p.val.perm
  }

  ######################################################################################################################
  #calculate the P-values from the permutation samples with the assumptions of indepedent variables in dataset-dat
  #inputs:dat-the dataset containing the variables for correlations,nsim: the number of permutation samples
  #input:MatchId-block paired id variable for paired study
  #assuming the paired observation at each pairid group has the same size
  #MatchId can't have any nonmissing observations
  #input: do.parallel, the option whether parallel computing method=T or F
  #do.parallel=F, do not use parallel computing,T, use parallel computing
  #output:the vector of p-values based on permuation resampling samples with size-p(p-1)/2
  #where p-the number of variables(columns) in dat
  pval.perm.corr.submat<-function(dat,nsim,MatchId=NULL){
    n<-dim(dat)[1]
    #the number of variales for correlatins
    p<-dim(dat)[2]
    if(is.null(colnames(dat))){
      warning("colnames of variables in the data are missing, please specify")
      break
    }
    if(any(is.na(colnames(dat)))){
      warning("Names of variables in the data contain missing values")
      break
    }

    if(!is.null(MatchId)){
    if( any(is.na(MatchId)))
    {warning("NA in paired group id")
      break}
    if(!any(is.na(MatchId))){
      size.b<-length(unique(MatchId))
    }}
    #the number of correlations from p variables
    pp<-p*(p-1)/2

    mat.corr.perm<-matrix(NA,ncol=pp,nrow=nsim)

    for ( ii in 1:nsim){

      mat.corr.perm[ii,]<-permut_within_corr(dat,MatchId=MatchId)
    }

    #Pearson's correlation from real data
    corr.true<-cor(dat, use="complete.obs")[lower.tri(cor(dat, use="complete.obs"), diag = FALSE)]
    #calculate the P-values from the permutation samples with the assumptions of indepedent variables in dataset-dat
    p.val.perm<-sapply(1:ncol(mat.corr.perm),function(x){(sum(abs(mat.corr.perm[,x])>=abs(corr.true[x]))+1)/(nsim+1)})

    from<-c()
    to<-c()
    for (i in 1:(dim(dat)[2]-1)){
      from<-c(from,rep(colnames(dat)[i],length((i+1):dim(dat)[2])))
      to<-c(to,colnames(dat)[-c(1:i)])}
    names(p.val.perm)<-paste(from,to,sep="_")
    return(p.val.perm)
  }
  ############################################################################

  if(n.submat==1){
    return(pval.perm.corr.submat(dat,nsim,MatchId))}
  if(n.submat>1){

    if(!do.parallel){
      pval.perm.corr.sub<-as.vector(sapply(1:n.submat,function(i){pval.perm.corr.submat(dat[,c(((i-1)*20+1):min(p,20*i))],nsim,MatchId)}))
      pval.perm.corr.cross2<-unlist(sapply(1:(n.submat-1),function(i){sapply((i+1):n.submat,function(j){pval.perm.corr.cross(dat[,c(((i-1)*20+1):min(p,20*i))],dat[,c(((j-1)*20+1):min(p,20*j))],nsim,MatchId)})}))
    }

    if(do.parallel){
     # library(parallel)
    #  library(foreach)
    #  library(doParallel)
      # Calculate the number of cores
     # no_cores <- detectCores() - 1
      cl<-makeCluster(no_cores)
      registerDoParallel(cl)

      pval.perm.corr.sub<-foreach( i = 1:n.submat,.combine=c,multicombine=TRUE,.export=c("pval.perm.corr.submat","permut_within_corr"))%dopar%
        pval.perm.corr.submat(dat[,c(((i-1)*20+1):min(p,20*i))],nsim,MatchId)


      pval.perm.corr.cross2<-foreach(i = 1:(n.submat-1),.combine=c,multicombine=TRUE,.export="pval.perm.corr.cross")%:%foreach(j=(i+1):n.submat,.combine=c,multicombine=TRUE,.export=c("pval.perm.corr.cross","permut_cross_corr"))%dopar%
        pval.perm.corr.cross(dat[,c(((i-1)*20+1):min(p,20*i))],dat[,c(((j-1)*20+1):min(p,20*j))],nsim,MatchId)

      stopCluster(cl)

    }
  }
  c(pval.perm.corr.sub,pval.perm.corr.cross2)
}
###########################################################################################################################################################################################################################

