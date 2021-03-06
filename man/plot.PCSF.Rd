% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.PCSF.R
\name{plot.PCSF}
\alias{plot.PCSF}
\title{Plot an interactive subnetwork}
\usage{
\method{plot}{PCSF}(x, style = 0, edge_width = 5, node_size = 40,
  node_label_cex = 30, Steiner_node_color = "lightblue",
  Terminal_node_color = "lightgreen", Terminal_node_legend = "Terminal",
  Steiner_node_legend = "Steiner", ...)
}
\arguments{
\item{x}{A subnetwork obtained by the PCSF method. It is a "PCSF" object derived 
from \pkg{igraph} class and it has the edge cost and vertex prize attributes.}

\item{style}{A \code{boolean} value to determine the visualization style of the network, 
where \code{0} plots the \code{static} network and \code{1} plots the \code{dynamic} 
network. The default valu is 0.}

\item{edge_width}{A \code{numeric} value to emphasize a maximum edge width. A default value is 5. 
This value must be greater than 1.}

\item{node_size}{A \code{numeric} value to emphasize a maximum node size. A default value is 40. 
This value must be greater than 10.}

\item{node_label_cex}{A \code{numeric} value to set a node label size. A default value is 30.}

\item{Steiner_node_color}{A \code{string} to set a color for \code{Steiner} nodes. 
A default value is "lightblue".}

\item{Terminal_node_color}{A \code{string} to set a color for \code{terminal} nodes.}

\item{Terminal_node_legend}{A \code{string} to set a legend for \code{terminal} nodes. 
A default legend is "Terminal".}

\item{Steiner_node_legend}{A \code{string} to set a legend for \code{Steiner} nodes. 
A default legend is "Steiner".}

\item{...}{Ignored.}
}
\description{
\code{plot.PCSF} plots an interactive figure of the subnetwork obrained by 
the PCSF method.
}
\details{
This function plots an interactive subnetwork obtained by the \code{\link{PCSF}} and \code{\link{PCSF_rand}}. 
The node sizes and edge widths are respectively proportional to the node prizes and edge costs 
while plotting the subnetwork from \code{\link{PCSF}}. In contrast, the node sizes and edge widths are 
proportional to the total number of abondance in randomized runs while plotting the subnetwork 
from \code{\link{PCSF_rand}}. The node names are displayed during the hover-over.
}
\examples{
\dontrun{
library("PCSF")
data("STRING")
data("Tgfb_phospho")
terminals <- Tgfb_phospho
ppi <- construct_interactome(STRING)
subnet <- PCSF(ppi, terminals, w = 2, b = 1, mu = 0.0005)
plot(subnet)}

}
\seealso{
\code{\link{PCSF}}, \code{\link{plot.PCSFe}}
}
\author{
Murodzhon Akhmedov
}
