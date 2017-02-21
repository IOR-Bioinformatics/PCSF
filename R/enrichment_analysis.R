#' Perform enrichment analysis on the subnetwork
#' 
#' \code{enrichment_analysis} performs functional enrichment analysis on the subnetwork 
#' obtained by the \code{\link{PCSF_rand}}, and returns an annotated subnetwork with top 15 
#' functional enrichments and a list of tables with a complete enrichment analysis for 
#' each cluster.
#' 
#' @param subnet A subnetwork provided by \code{\link{PCSF_rand}}, which is obtained by merging 
#' a multiple outputs of the PCSF with random noise added edge costs. An \pkg{igraph} object 
#' with edge cost and vertex prize attributes representing the total number of 
#' show ups throughout all runs.
#' 
#' @return A list composed of an interactive subnetwork and a table with enrichment 
#' analysis results. An interactive subnetwork annotated with enrichment analysis 
#' can be reached by $subnet. A full list of enrichment analysis for each cluster 
#' can be reached by $enrichment.
#' 
#' @export
#' 
#' @details 
#' An enrichment analysis of the final subnetwork obtained by multiple runs of the PCSF 
#' (with rando noise added edge costs) is performed for functional interpretation. 
#' The subnetwork is clustered using an edge betweenness clustering algorithm from 
#' the \pkg{igraph} package, and for each cluster functional enrichment is done by 
#' employing the ENRICHR API (Chen \emph{et al.}, 2013). An interactive visualization of 
#' the final subnetwork is plotted, where the node sizes and edge widths are proportional 
#' to the frequency of show ups throughout total runs. Nodes are colored according to the 
#' cluster membership, and the top 15 functional enrichment terms are displayed in tabular 
#' format during the hover-over of the node in that cluster. 
#' 
#' 
#' @examples 
#' \dontrun{
#' library("PCSF")
#' data("STRING")
#' data("Tgfb_phospho")
#' terminals <- Tgfb_phospho
#' ppi <- construct_interactome(STRING)
#' subnet <- PCSF_rand(ppi, terminals, n = 10, r = 0.1, w = 2, b = 1, mu = 0.0005)
#' res <- enrichment_analysis(subnet)}
#' 
#' \dontrun{
#' plot(res$subnet)
#' write.table(res$enrichment[[1]],file="cluster1_complete_enrichment.txt", 
#'              append = FALSE, quote = FALSE, sep ="\t", row.names=FALSE)}
#'              
#' @author Murodzhon Akhmedov
#'
#' @seealso \code{\link{PCSF_rand}}, \code{\link{plot.PCSFe}}
#' 
#' 
#' 
#' @references 
#' Chen E.Y., Christopher M.T., Yan K., Qiaonan D., Zichen W., Gabriela V.M., Neil R.C., 
#' and Avi M. (2013) "Enrichr: Interactive and Collaborative Html5 Gene List Enrichment 
#' Analysis Tool." \emph{BMC Bioinformatics} 14 (1). BioMed Central: 1.

enrichment_analysis <-function(subnet){
  
  # Checking function arguments
  if (missing(subnet))
    stop("Need to specify the subnetwork obtained from the PCSF algorithm.")
  if (class(subnet)[1] != "PCSF" || class(subnet)[2] != "igraph")
    stop("The subnetwork must be a \"PCSF\" object derived from an \"igraph\" class.")
  

  cat("  Performing enrichment analysis...\n\n")
  
  # Obtain clusters in the subnet using edge betweenness clustering algorithm from igraph package. 
  clusters = cluster_edge_betweenness(subnet)
  
  # Perform ebrichment analysis for each cluster using ENRICHR through its API.
  enrich = call_enr(clusters)
  
  
  
  enrichment = enrich[[1]]
  enrichment_complete = enrich[[2]]
  
  # Add 'group" and 'title' attributes to subnet 
  V(subnet)$group = clusters$membership
  V(subnet)$title = paste0("Cluster ",clusters$membership,": Enrichment analysis")
  for( i in 1:length(V(subnet))){
    V(subnet)$title[i] = paste0( V(subnet)$title[i], enrichment[[V(subnet)$group[i]]])
  }
  
  # Derive a "PCSFe" object from an "igraph" class. 
  class(subnet) <- c("PCSFe", "igraph")
  # Combine the subnetwork and colplete enrichment analysis tables.
  output = list(subnet, enrichment_complete)
  names(output) = c("subnet", "enrichment")
  
  return (output)
}