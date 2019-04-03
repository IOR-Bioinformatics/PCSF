#' Plot an interactive subnetwork with functional enrichment analysis
#'
#' \code{plot.PCSFe} plots an interactive figure of the subnetwork 
#' to display the functionla enrichment analysis, which is obtained by employing 
#' \code{enrichment_analysis} on the subnetwork.
#' 
#' @param x An output subnetwork provided by the \code{enrichment_analysis}. 
#' It is "PCSFe" object derived from an \pkg{igraph} class, and it has the edge 
#' cost and vertex prize attributes.
#' @param edge_width A \code{numeric} value to emphasize a maximum edge width. 
#' A default value is 5. This value must be greater than 1.
#' @param node_size A \code{numeric} value to emphasize a maximum node size. 
#' A default value is 30. This value must be greater than 10.
#' @param node_label_cex A \code{numeric} value to set a node label size. 
#' A default value is 1.
#' @param Steiner_node_legend A \code{string} to set a legend for \code{Steiner} nodes. 
#' A default legend is "Steiner".
#' @param Terminal_node_legend A \code{string} to set a legend for \code{terminal} nodes. 
#' @param extra_node_colors A \code{list} with colors of extra types of nodes added to the PCSF result, with the names of the list being the node type
#' A default legend is "Terminal".
#' @param ... Ignored.
#' @import igraph visNetwork
#' @method plot PCSFe
#' @export
#' 
#' @details 
#' 
#' An enrichment analysis of the final subnetwork obtained by multiple runs of the PCSF 
#' (with random noise added edge costs) is performed by using \code{\link{enrichment_analysis}}. 
#' The subnetwork is clustered using an edge betweenness clustering algorithm from the 
#' \pkg{igraph} package, and for each cluster functional enrichment is done by employing the 
#' ENRICHR API (Chen \emph{et al.}, 2013). An interactive visualization of the final subnetwork 
#' is plotted, where the node sizes and edge widths are proportional to the frequency of show 
#' ups in total randomised runs. Nodes are colored according to the cluster membership, and 
#' the top 15 functional enrichment terms are displayed in tabular format during the hover-over 
#' of the node in that cluster. A specific cluster can be displayed separately in the figure 
#' by selecting from the icon list at the top left side of the figure.
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
#' res <- enrichment_analysis(subnet)
#' plot(res$subnet)}
#' 
#' @author Murodzhon Akhmedov
#' 
#' @references 
#' Chen E.Y., Christopher M.T., Yan K., Qiaonan D., Zichen W., Gabriela V.M., Neil R.C., and Avi M. (2013) 
#' Enrichr: Interactive and Collaborative Html5 Gene List Enrichment Analysis Tool. \emph{BMC Bioinformatics} 14 (1). 
#' BioMed Central: 1.
#' 
#' @seealso \code{\link{enrichment_analysis}}, \code{\link{PCSF_rand}}, \code{\link{plot.PCSF}}



plot.PCSFe <-function(x, edge_width = 5, node_size = 30, node_label_cex = 1, 
                      Terminal_node_legend = "Terminal",
                      Steiner_node_legend = "Steiner",
                       extra_node_colors = list(),
                      ...){
  
  
  subnet = x
  # Checking function arguments
  if (missing(subnet))
    stop("Need to specify the subnetwork obtained from the PCSF algorithm.")
  if (class(subnet)[1] != "PCSFe" || class(subnet)[2] != "igraph")
    stop("The subnetwork must be a \"PCSFe\" object derived from an \"igraph\" class.")
  if (edge_width < 2)
    stop("The edge_width must be greater than 2.")
  if (node_size < 10)
    stop("The node_size must be greater than 10.")
  
  # Add 'label' and 'size' attributes
  V(subnet)$label.cex = node_label_cex
  prize = abs(V(subnet)$prize)
  min1 = 10
  max1 = node_size
  r1 = max1 - min1
  min2 =  min(prize)
  max2 = max(prize)
  r2 = max2 - min2
  adjusted_prize = r1*(prize - min2)/r2 + min1
  V(subnet)$size = adjusted_prize
  
  # Add edge width' attributes
  weight = E(subnet)$weight
  min1 = 2
  max1 = edge_width
  r1 = max1 - min1
  min2 =  min(weight)
  max2 = max(weight)
  r2 = max2 - min2
  adjusted_weight = r1*(weight - min2)/r2 + min1
  E(subnet)$width = adjusted_weight
  
  
  # Associate the type of nodes to shape
  shape = V(subnet)$type
  shape[which(shape=="Steiner")] = "triangle"
  shape[which(shape=="Terminal")] = "circle"
  V(subnet)$shape = shape
  
  # Visualize the subnet
  visIgraph(subnet) %>%
    visIgraphLayout(layout = "layout_with_fr") %>%
    visOptions(highlightNearest = list(enabled = T), selectedBy = "group")%>%
    visLegend(addNodes = list(
      list(label = Terminal_node_legend, shape = "dot", size = 15, label.cex = 0.3),
      list(label = Steiner_node_legend, shape = "triangle",size = 9, label.cex = 0.3)), width = 0.2,
      useGroups = FALSE)
}