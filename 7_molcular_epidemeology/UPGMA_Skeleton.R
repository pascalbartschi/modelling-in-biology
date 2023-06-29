#14.12.22, BIO445

# Reconstruct trees using UPGMA algorithm and the distance matrices you calculated above
# Make a function to calculate upgma tree from a given matrix with follwing inputs
# inputmatrix: matrix with distances between nodes in the tree
# internalnodestart: hold the counter for the internal nodes, it starts at 2n-1
# nodecounts: is a vector holding the counts of tips under each node with a
# corresponding entry in the distance matrix
# edge: is a matrix holding the connections (edges) between the nodes, this
# will be the output of the recursive function
# edge_length: is a vector defining the lengts of edges defined in the "edge" matrix
# it is not implemented with distance, that could be the next task

upgma<-function(inputmatrix,internalnodestart,nodecounts,edge,edge_length){
  # This function is a recursion, so we need a way to stop it,
  # if we reach smaller than a 2-by-2 matrix
  
  stopnow <- ???
  
  # which row and column do we select?
  minrow<-???
  mincol<-???
  
  #define arbitrary edge_length
  edge_length<-c(0,edge_length)
  edge_length<-c(0,edge_length)
  
  # define edge 
  edge<-rbind(c(internalnodestart,as.numeric(dimnames(inputmatrix)[[1]][minrow])),edge) 
  edge<-rbind(c(internalnodestart,as.numeric(dimnames(inputmatrix)[[2]][mincol])),edge) 
  
  # Stop the recursion if we reached a 2-by-2 matrix and want to go to 1-by-1
  if (dim(inputmatrix)[1]==2){stopnow<-???} 
  else {
    tempmatrix<-???
    # change one of the columns to hold average distance to all other entries 
    for (i in 1:nrow(inputmatrix)){
      if (i!=minrow && i!=mincol){
        inputmatrix[min(i,min(mincol,minrow)),max(i,min(mincol,minrow))]<-
          (nodecounts[mincol]*tempmatrix[min(i,mincol),max(i,mincol)]+nodecounts[minrow]*tempmatrix[min(i,minrow),max(i,minrow)])/(nodecounts[mincol]+nodecounts[minrow])}
    }
    # change the name of that column
    dimnames(inputmatrix)[[1]][min(mincol,minrow)]<-???
    dimnames(inputmatrix)[[2]][min(mincol,minrow)]<-???
    
    # change the number of nodes/tips that column holds
    nodecounts[min(mincol,minrow)]<-nodecounts[mincol]+nodecounts[minrow]
    nodecounts<-nodecounts[-max(mincol,minrow)]
    
    # remove the other row and column from the matrix
    if(!stopnow){ 
      inputmatrix<-inputmatrix[-max(mincol,minrow),]
      inputmatrix<-inputmatrix[,-max(mincol,minrow)]
    }
    
    # change the internal node counter for next recursion
    internalnodestart<- ???
    
    #The method is recursive 
    upgmaout<- ???
    edge<- ???
    edge_length<- ???
  }
  return(list(edge,edge_length)) 
}



## make a help function for passing matrices to upgma function and returning trees
reconstruct<-function(inputmatrix){
  tempmatrix<-???
  
  # for convenience, change the dimension names of the distance matrix to numbers
  dimnames(tempmatrix)[[1]]<-seq(1:nrow(inputmatrix))
  dimnames(tempmatrix)[[2]]<-seq(1:ncol(inputmatrix))
  
  # check if the matrix has NaN entries, if so do not reconstruct the tree
  # Try to figure out what does NaN entry mean? Where and when does it occur? 
  if (sum(is.nan(inputmatrix))>0){
    return("NULL") }
  
  # internalnodestart hold the counter for the internal nodes, it starts at 2n-1
  internalnodestart<-nrow(tempmatrix)*2-1
  
  # nodecounts is a vector holding the counts of tips under each node with a #corresponding entry in the distance matrix 
  nodecounts<-rep(1,nrow(tempmatrix))
  output<-upgma(tempmatrix,internalnodestart=internalnodestart,
                nodecounts=nodecounts,edge=NULL,edge_length=NULL) 
  
  # define the tree as "phylo" class
  tree<-list(edge = output[[1]]) 
  tree$tip.label<-???
  tree$edge.length<-???
  tree$Nnode<-dim(inputmatrix)[1]-1
  class(tree)<-"phylo"
  return(tree) 
}


#Plotting the tree!
dist_ha[upper.tri(dist_ha)]<-NA 
diag(dist_ha)<-NA 
tree<-reconstruct(dist_ha) 
plot(tree,use.edge.length = FALSE)

