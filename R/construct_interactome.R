#' Construct an interaction network
#'
#' Given a list of edges, \code{construct_interactome} generates 
#' an interaction network which is used as a template network to interpret the highthrougput data.
#' 
#' @param ppi A list of edges. A \code{data.frame} composed of three columns, where each
#'  row corresponds to an edge in which the first element is a \code{head}, the second 
#'  element is a \code{tail}, and the last element represents the \code{cost} of the edge.
#'  
#' @return An interaction network as \pkg{igraph} object.
#' @import igraph
#' @export
#' 
#' @examples 
#' \dontrun{
#' library("PCSF")
#' data("STRING")
#' ppi <- construct_interactome(STRING)}
#' 
#' @author Murodzhon Akhmedov
#' 

construct_interactome <-
function(ppi){
  
  # Checking function arguments
  if (missing(ppi))
    stop("  Need to specify a list of edges to construct an interaction network.
    Provide a data.frame composed of three columns, where each row corresponds 
    to an edge in which the first element is a head node, the second element 
    is a tail node, and the last element represents the cost of the edge.")
  if (nrow(ppi)<1 || ncol(ppi) != 3 || class(ppi) != "data.frame")
    stop("  Need to provide a data.frame composed of three columns, where each row corresponds 
    to an edge in which the first element is a head node, the second element 
    is a tail node, and the last element represents the cost of the edge.")

  # Interpolate the node prizes
  node_names = unique(c(as.character(ppi[,1]),as.character(ppi[,2])))
  
  # Contruct an interaction network as igraph object
  ppi.graph = graph.data.frame(ppi[,1:2],vertices=node_names,directed=F)
  E(ppi.graph)$weight=as.numeric(ppi[,3])
  ppi.graph = simplify(ppi.graph)

  return (ppi.graph)

}
