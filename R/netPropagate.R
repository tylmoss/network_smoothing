#' FUNCTION Propagate information across an interaction network
#' 
#' Function takes in node information in the form of a vector or matrix
#' and smooths/propagates the information across network of interacting
#' nodes represeted as an adjaceny matrix.
#' 
#' @param mutData The mutation data in vector or matrix format. Samples 
#' in rows and features/genes in columns
#' 
#' @param adj.mat The network information represented as an adjaceny matrix.
#' 
#' @param alpha Weight (0 to 1) given to interaction data.
#' 
#' @param verbose Decide whether or not print individual step results to console.
#' 
#' @param include.input Decide whether or not to include input data with results.
#' 
#' @return results List of results and input data/parameters
#' 
#' @examples \dontrun{
#' # Generate random undirected network
#' genes = paste0("g",1:150)
#' net = data.frame(from = sample(genes, 50, TRUE), to = sample(genes, 10, TRUE), stringsAsFactors = F)
#' net = net[net$from != net$to,] # remove self loops
#' # sort interacting nodes alphabetically and remove duplicate interactions
#' for(i in 1:nrow(net)){net[i,] = sort(net[i,])}
#' net = net[!duplicated(apply(net,1,function(s)paste(s,collapse = "_"))),]
#' # create igraph object
#' library(igraph)
#' net = graph.data.frame(net, directed = F)
#' # create adjacency matrix
#' net.adj = as.matrix(get.adjacency(net))
#' # weight edges by degree of connected nodes
#' W = degreeNormalize(net.adj)
#' # Generate random mutation data for single sample
#' Y = sample(c(0,1), ncol(W), replace = T, prob = c(.8,.2))
#' # smooth data over network
#' res1 = netPropagate(mutData = Y, adj.mat = W, alpha = .7)
#' myLayout = layout.fruchterman.reingold(net)
#' # plot input data
#' plot(net, vertex.label.color = "black", 
#'  vertex.color = unlist(lapply(Y, function(x)rgb(colorRamp(c("white","red"))(x), 
#'    maxColorValue = 255))), 
#'  vertex.label.cex = .7, layout = myLayout, vertex.label = NULL, main = paste("input"))
#' # plot smoothed data
#' plot(net, vertex.label.color = "black", 
#'  vertex.color = unlist(lapply(res1$result, function(x)rgb(colorRamp(c("white","red"))(x), 
#'    maxColorValue = 255))),
#'  vertex.label.cex = .7, layout = myLayout, main = paste("step",res1$steps))
#' }
#' 
#' @export
netPropagate = function(mutData, adj.mat, alpha = 0.7, verbose = F, include.input = F){
  if(alpha == 1){
    stop('Will not converge on a solution if alpha is set to "1".\nPlease pick a value: 0 < alpha < 1')  
  }
  # initialize variables
  steps = 1 # steps to similuate a 'random walk'
  F_0 = alpha*(mutData %*% adj.mat) + (1-alpha)*mutData # First step of main function
  F_t = F_0 # container matrix for loop
  convergance = 1
  while(convergance > 10^-6){
    F_tt = alpha*(F_t %*% adj.mat) + (1-alpha)*mutData
    steps = steps+1
    convergance = norm(F_tt - F_t) # calculate matrix norm of F_t - F_t-1
    
    if(verbose){
      # print results progress.
      message(paste0("range F_tt: ", paste(range(F_t), collapse = " - ")))
      message(paste0("steps: ", steps))
      message(paste0("norm: ", convergance))
    }
    # update matrix and proceed to next step until convergance
    F_t = F_tt
  }
  # list of results
  if(include.input){
    results = list(alpha = alpha, steps = steps-1, convergance = convergance, result = F_t, input = mutData)
  }else{
    results = list(alpha = alpha, steps = steps-1, convergance = convergance, result = F_t)
  }
  return(results)
}
