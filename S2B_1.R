S2B_1=function(seed_graph,index1,index2,nrep,nrep2){
  meandist=igraph::mean_distance(seed_graph, directed = TRUE)
  bt=subS2B_version1(seed_graph,index1,index2,meandist) # specify subS2B version to be used to compute S2B scores
  pbt=rep(0,gorder(seed_graph))
  nscore=rep(0,gorder(seed_graph))
  if (nrep>0){
    rbt_matrix2=matrix(nrow=length(bt$allcount),ncol=nrep)
    for (i in 1:nrep){
      rindex1=sample(igraph::gorder(seed_graph),length(index1),replace=FALSE) # select random seeds for index 1
      rindex2=sample(igraph::gorder(seed_graph),length(index2),replace=FALSE) # select random seeds for index 2
      rbt=subS2B_version1(seed_graph,rindex1,rindex2,meandist) # specify subS2B version, S2B scores after shuffle of seed identity
      nscore[rbt$allcount<bt$allcount]=nscore[rbt$allcount<bt$allcount]+1
      rbt_matrix2[,i]=rbt$allcount
    }
    nscore=nscore/nrep
  } else {
    rbt_matrix2=matrix()
  }
  
  if (nrep2>0){
    rbt_matrix=matrix(nrow=length(bt$allcount),ncol=nrep2)
    for (i in 1:nrep2){
      ee=igraph::ecount(seed_graph)
      rg=igraph::rewire(seed_graph, with=keeping_degseq(loops=FALSE, niter=ee*10)) # shuffles two randomly choosen edges at a time for a number of iterations equal to 10 times the number edges in seed_graph
      rbt=subS2B_version1(rg,index1,index2,meandist) # specify subS2B version, S2B scores after shuffle of edges
      pbt[rbt$allcount<bt$allcount]=pbt[rbt$allcount<bt$allcount]+1
      rbt_matrix[,i]=rbt$allcount
    }
    pbt=pbt/nrep2
  } else {
    rbt_matrix=matrix()
  }
  
  bigvertexlist=igraph::vertex_attr(seed_graph)
  allstat=data.frame(protein=bigvertexlist[[1]],bcount=bt$allcount,score=pbt,nscore=nscore)
  s2btable=makes2btable(allstat,seed_graph,index1,index2) # returns table with S2B scores (bt), specificity scores (pbt and nscore), node classification as candidate, seed from disease 1 or seed from disease 2, direct neighbors seeds of the nodes and crossbridges formed between seeds with corresponding specificity score
  list(s2btable=s2btable,seedmat1=bt$smat1,seedmat2=bt$smat2,maxS2B=bt$maxS2B)
}