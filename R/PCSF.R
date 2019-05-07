#' Prize-collecting Steiner Forest (PCSF)
#'
#' \code{PCSF} returns a subnetwork obtained by solving the PCSF on the given interaction network.
#'
#' @param ppi An interaction network, an \pkg{igraph} object.
#' @param terminals  A list of terminal genes with prizes to be analyzed in the PCSF context.
#' A named \code{numeric} vector, where terminal genes are named same as in the interaction network
#' and numeric values correspond to the importance of the gene within the study.
#' @param w A \code{numeric} value for tuning the number of trees in the output. A default value is 2.
#' @param b A \code{numeric} value for tuning the node prizes. A default value is 1.
#' @param mu A \code{numeric} value for a hub penalization. A default value is 0.0005.
#' @param dummies A list of nodes that are to connected to the root of the tree. If missing the root will be connected to all terminals.
#' @return The final subnetwork obtained by the PCSF.
#' It return an \pkg{igraph} object with the node prize and edge cost attributes.
#' @import igraph
#' @export
#'
#' @details
#'
#' The PCSF is a well-know problem in graph theory.
#' Given an undirected graph \emph{G = (V, E)}, where the vertices are labeled with prizes
#' \eqn{p_{v}} and the edges are labeled with costs \eqn{c_{e} > 0}, the goal is to identify
#' a subnetwork \emph{G' = (V', E')} with a forest structure. The target is to minimize
#' the total edge costs in \emph{E'}, the total node prizes left out of \emph{V'}, and the
#' number of trees in \emph{G'}. This is equivalent to  minimization of the following
#' objective function:
#'
#' \deqn{F(G')= Minimize \sum_{ e \in E'} c_{e} + \beta*\sum_{v \not\in V'} p_v + \omega*k}
#' where, \emph{k} is the number of trees in the forest, and it is regulated by parameter \eqn{\omega}.
#' The parameter \eqn{\beta} is used to tune the prizes of nodes.
#'
#' This optimization problem nicely maps onto the problem of finding differentially
#' enriched subnetworks in the cell protein-protein interaction (PPI) network.
#' The vertices of interaction network correspond to genes or proteins, and edges
#' represent the interactions among them. We can assign prizes
#' to vertices based on measurements of differential expression, copy number, or
#' mutation, and costs to edges based on confidence scores for those intra-cellular
#' interactions from experimental observation, yielding a proper input to the PCSF
#' problem. Vertices that are assigned a prize are referred to \emph{terminal} nodes,
#' whereas the vertices which are not observed in patient data are not assigned a
#' prize and are called \emph{Steiner} nodes. After scoring the interactome, the
#' PCSF is used to detect a relevant subnetwork (forest), which corresponds to a
#' portion of the interactome, where many genes are highly correlated in terms of
#' their functions and may regulate the differentially active biological process
#' of interest. The PCSF aims to identify  neighborhoods in interaction networks
#' potentially belonging to the key dysregulated pathways of a disease.

#' In order to avoid a bias towards the hub nodes of PPI networks to appear in solution
#' of PCSF, we penalize the prizes of \emph{Steiner} nodes according to their degree
#' distribution in PPI, and it is regulated by parameter \eqn{\mu}:
#'
#'  \deqn{p'_{v} = p_{v} - \mu*degree(v)}
#'
#'  The parameter \eqn{\mu} also affects the total number of \emph{Steiner} nodes in the solution.
#'  Higher the value of \eqn{\mu} smaller the number of \emph{Steiners} in the subnetwork,
#'  and vice-versa. Based on our previous analysis the recommended range of \eqn{\mu}
#'  for biological networks is between 1e-4 and 5e-2, and users can choose the values
#'  resulting subnetworks with vertex sets that have desirable \emph{Steiner/terminal}
#'  node ratio and average \emph{Steiner/terminal} in-degree ratio
#'  in the template interaction network.
#'
#' @examples
#' \dontrun{
#' library("PCSF")
#' data("STRING")
#' data("Tgfb_phospho")
#' terminals <- Tgfb_phospho
#' ppi <- construct_interactome(STRING)
#' subnet <- PCSF(ppi, terminals, w = 2, b = 1, mu = 0.0005)}
#'
#' @author Murodzhon Akhmedov
#'
#' @seealso \code{\link{PCSF_rand}}, \code{\link{plot.PCSF}}
#'
#' @references
#' Akhmedov M., LeNail A., Bertoni F., Kwee I., Fraenkel E., and Montemanni R. (2017)
#' A Fast Prize-Collecting Steiner Forest Algorithm for Functional Analyses in Biological Networks.
#' \emph{Lecture Notes in Computer Science}, to appear.


PCSF <-
function(ppi, terminals, w = 2, b = 1, mu = 0.0005, dummies){

  # Checking function arguments
  if (missing(ppi))
    stop("Need to specify an interaction network \"ppi\".")
  if (class(ppi) != "igraph")
    stop("The interaction network \"ppi\" must be an igraph object.")
  if (missing(terminals))
    stop("  Need to provide terminal nodes as a named numeric vector,
    where node names must be same as in the interaction network.")
  if(is.null(names(terminals)))
    stop("  The terminal nodes must be provided as a named numeric vector,
    where node names must be same as in the interaction network.")




  # Gather the terminal genes to be analyzed, and their scores
  terminal_names = names(terminals)
  terminal_values = as.numeric(terminals)

  # Incorporate the node prizes
  node_names = V(ppi)$name
  node_prz = vector(mode = "numeric", length = length(node_names))
  index = match(terminal_names, node_names)
  percent = signif((length(index) - sum(is.na(index)))/length(index)*100, 4)
  if (percent < 5)
    stop("  Less than 1% of your terminal nodes are matched in the interactome, check your terminals!")
  cat(paste0("  ", percent, "% of your terminal nodes are included in the interactome\n"))
  terminal_names = terminal_names[!is.na(index)]
  terminal_values = terminal_values[!is.na(index)]
  index = index[!is.na(index)]
  node_prz[index] =  terminal_values

  if(missing(dummies)||is.null(dummies)||is.na(dummies))
    dummies = terminal_names #re-assign this to allow for input

  ## Prepare input file for MST-PCSF implementation in C++

  cat("  Solving the PCSF...\n")

  # Calculate the hub penalization scores
  node_degrees = igraph::degree(ppi)
  hub_penalization = - mu*node_degrees

  # Update the node prizes
  node_prizes = b*node_prz
  index = which(node_prizes==0)
  node_prizes[index] = hub_penalization[index]


  # Construct the list of edges
  edges = ends(ppi,es = E(ppi))
  from = c(rep("DUMMY", length(dummies)), edges[,1])
  to = c(dummies, edges[,2])

  cost = c(rep(w, length(dummies)), E(ppi)$weight)

  #PCSF will faill if there are NAs in weights, this will check and fail gracefully
  if(any(is.na(E(ppi)$weight))){

  }

  ## Feed the input into the PCSF algorithm
  output = call_sr(from,to,cost,node_names,node_prizes)

  # Check the size of output subnetwork and print a warning if it is 0
  if(length(output[[1]]) != 0){

    #names(output) = c("from", "to", "cost", "terminal_names", "terminal_prizes")

    # Contruct an igraph object from the MST-PCSF output
    e = data.frame(output[[1]], output[[2]], output[[3]])
    #e = e[which(e[,1]!="DUMMY"), ]
    Ee = e[which(e[,2]!="DUMMY"), ]
    names(e) = c("from", "to", "weight")


    # Differentiate the type of nodes
    type = rep("Steiner", length(output[[4]]))
    index = match(terminal_names, output[[4]])
    index = index[!is.na(index)]
    type[index] = "Terminal"

    v = data.frame(output[[4]], output[[5]], type)
    names(v) = c("terminals", "prize", "type")
    subnet = graph.data.frame(e,vertices=v,directed=F)
    E(subnet)$weight=as.numeric(output[[3]])
    subnet = delete_vertices(subnet, "DUMMY")
    subnet = delete_vertices(subnet, names(which(degree(subnet)==0)))


    class(subnet) <- c("PCSF", "igraph")

    return (subnet)

  } else{

    stop("  Subnetwork can not be identified for a given parameter set.
    Provide a compatible b or mu value with your terminal prize list...\n\n")

  }

}
