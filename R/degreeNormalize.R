#' FUNCTION Weight edges in adjaceny matrix by the degree of interacting vertices.
#' 
#' Takes as input an adjacency matrix and return as degree normalized matrix 
#' of the same dimensions.  The adjacency matrix should be from and undirected network
#' and thus should be symmetric.
#' 
#' @param adj.mat  Adjacency matrix.  The matrix needs should be from an undirected network 
#' and should be symmetric across the diagonal.
#' 
#' @return degree normalized matrix
#' 
#' @examples
#' myMat = matrix(c(0, 1, 1, 1, 1, 0, 0, 1, 1, 0 , 0, 0, 1, 1, 0, 0), 4, 4)
#' myMat.norm = degreeNormalize(myMat)
#' myMat.norm
#' 
#' @export
degreeNormalize = function(adj.mat){
  # sum the number of interactions (degree) for columns and rows
  D_i = apply(adj.mat, 2, function(z)sum(z != 0))^0.5
  D_j = apply(adj.mat, 1, function(z)sum(z != 0))^0.5
  # divide by degree
  W = sweep(adj.mat, 2, D_i, "/")
  W = sweep(W, 1, D_i, "/")
  return(W)
}
