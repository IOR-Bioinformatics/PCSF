#' Plot an interactive subnetwork
#'
#' \code{plot.PCSF} plots an interactive figure of the subnetwork obrained by 
#' the PCSF method.
#' 
#' @param x A subnetwork obtained by the PCSF method. It is a "PCSF" object derived 
#' from \pkg{igraph} class and it has the edge cost and vertex prize attributes.
#' @param style A \code{boolean} value to determine the visualization style of the network, 
#' where \code{0} plots the \code{static} network and \code{1} plots the \code{dynamic} 
#' network. The default valu is 0.
#' @param edge_width A \code{numeric} value to emphasize the maximum edge width. A default value is 5. 
#' This value must be greater than 1.
#' @param node_size A \code{numeric} value to emphasize the maximum node size. A default value is 40. 
#' This value must be greater than 10.
#' @param node_label_cex A \code{numeric} value to set the node label size. A default value is 30.
#' @param Steiner_node_color A \code{string} to set the color of \code{Steiner} nodes. 
#' A default value is "lightblue".
#' @param Terminal_node_color A \code{string} to set the color of \code{terminal} nodes. 
#' A default value is "lightgreen".
#' @param ... Ignored.
#' @import igraph visNetwork
#' @method plot PCSF
#' @export
#' 
#' 
#' @details 
#' This function plots an interactive subnetwork obtained by the \code{\link{PCSF}} and \code{\link{PCSF_rand}}. 
#' The node sizes and edge widths are respectively proportional to the node prizes and edge costs 
#' while plotting the subnetwork from \code{\link{PCSF}}. In contrast, the node sizes and edge widths are 
#' proportional to the total number of abondance in randomized runs while plotting the subnetwork 
#' from \code{\link{PCSF_rand}}. The node names are displayed during the hover-over. 
#' 
#' @examples 
#' \dontrun{
#' library("PCSF")
#' data("STRING")
#' data("Tgfb_phospho")
#' terminals <- Tgfb_phospho
#' ppi <- construct_interactome(STRING)
#' subnet <- PCSF(ppi, terminals, w = 2, b = 1, mu = 0.0005)
#' plot(subnet)}
#' 
#' @author Murodzhon Akhmedov
#'  
#' @seealso \code{\link{PCSF}}, \code{\link{plot.PCSFe}}

plot.PCSF <-
function(x, style = 0, edge_width=5, node_size=40, node_label_cex = 30, Steiner_node_color = "lightblue", 
         Terminal_node_color = "lightgreen", ...){
  
  subnet = x
  # Checking function arguments
  if (missing(subnet))
    stop("Need to specify the subnetwork obtained from the PCSF algorithm.")
  if (class(subnet)[1] != "PCSF" || class(subnet)[2] != "igraph")
    stop("The subnetwork must be a \"PCSF\" object derived from an \"igraph\" class.")
  if (edge_width < 1)
    stop("The edge_width must be greater than 1.")
  if (node_size < 10)
    stop("The node_size must be greater than 10.")
  
  # Calculate the adjusted node prizes
  prize = abs(V(subnet)$prize)
  min1 = 10
  max1 = node_size
  r1 = max1 - min1
  min2 =  min(prize)
  max2 = max(prize)
  r2 = max2 - min2
  adjusted_prize = r1*(prize - min2)/r2 + min1
  
  # Calculate the adjusted edge weights
  weight = E(subnet)$weight
  min1 = 1
  max1 = edge_width
  r1 = max1 - min1
  min2 =  min(weight)
  max2 = max(weight)
  r2 = max2 - min2
  adjusted_weight = r1*(weight - min2)/r2 + min1
  
  
  if(style){
    
    # List of nodes in the subnet
    nodes = data.frame(1:length(V(subnet)), V(subnet)$name)
    names(nodes) = c("id", "name")
    
    # Differentiate the type of nodes
    nodes$group = V(subnet)$type
    
    # Attach the node attributes
    nodes$size = adjusted_prize
    nodes$title = nodes$name
    nodes$label = nodes$name
    nodes$label.cex = node_label_cex
    nodes$font.size = node_label_cex
    
    # List of edges in the subnet
    edges = data.frame(ends(subnet,es = E(subnet)), adjusted_weight)
    names(edges) = c("from", "to", "width")
    edges$from = match(edges$from, nodes$name)
    edges$to = match(edges$to, nodes$name)
    
    # Visualize the subnet
    visNetwork(nodes,edges) %>%
      visNodes( shadow = list(enabled = TRUE, size = 10))  %>%
      visGroups(groupname = "Steiner", color = list(background = Steiner_node_color, border = "blue"), shape = "triangle") %>%
      visGroups(groupname = "Terminal", color = list(background = Terminal_node_color, border = "green"), shape = "dot") %>%
      visOptions(highlightNearest = list(enabled = T)) %>%
      visLegend(addNodes = list(
        list(label = "Terminal", shape = "dot", size = 15, color = list(background = Terminal_node_color, border = "green"), label.cex = 0.8),
        list(label = "Steiner", shape = "triangle",size = 10, color = list(background = Steiner_node_color, border = "blue"), label.cex = 0.8 )), width = 0.15,
        useGroups = FALSE)
    
    
  } else{
    
    # Attach the node type attribute
    V(subnet)$group = V(subnet)$type
    
    # Attach the node attributes: size, title, label, label.cex, font size
    V(subnet)$size = adjusted_prize
    V(subnet)$title = V(subnet)$name
    V(subnet)$label = V(subnet)$name
    V(subnet)$label.cex = node_label_cex/30
    V(subnet)$font.size = node_label_cex/30
    
    # Attach the edge width attribute
    E(subnet)$width = adjusted_weight
    
    
    # Visualize the subnet
    visIgraph(subnet) %>%
      visIgraphLayout(layout = "layout_with_fr") %>%
      visNodes( shadow = list(enabled = TRUE, size = 10))  %>%
      visGroups(groupname = "Steiner", color = list(background = Steiner_node_color, border = "blue"), shape = "triangle") %>%
      visGroups(groupname = "Terminal", color = list(background = Terminal_node_color, border = "green"), shape = "dot") %>%
      visOptions(highlightNearest = list(enabled = T)) %>%
      visLegend(addNodes = list(
        list(label = "Terminal", shape = "dot", size = 14, color = list(background = Terminal_node_color, border = "green"), label.cex = 0.5),
        list(label = "Steiner", shape = "triangle",size = 8, color = list(background = Steiner_node_color, border = "blue"), label.cex = 0.5)), width = 0.2,
        useGroups = FALSE)

  }

}
