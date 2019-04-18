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
#' @param mode A binary variable to choose the method for enrichment analysis, where 0 is for EnrichR API and 1 is for \pkg{topGO} package.
#' @param gene_universe A complete list of genes (vector of gene symbols) used as background in enrichment analysis by \pkg{topGO} package.
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
#' employing either EnrichR API (Chen \emph{et al.}, 2013) or
#' \pkg{topGO} (Alexa and Rahnenfuhrer, 2009)
#' package that is specified by the user. Important to note that EnrichR API requires
#' a working Internet connection to perform the enrichment. If the user does not
#' specify which tool to use for enrichment analysis, the package employs EnrichR
#' as a default if there is Internet connection, otherwise it uses \pkg{topGO}.
#'
#' An interactive visualization of
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
#' res <- enrichment_analysis(subnet)
#' res <- enrichment_analysis(subnet, mode=0)}
#' \dontrun{
#' library(topGO)
#' gene_universe <- V(ppi)$name
#' res <- enrichment_analysis(subnet, mode=1, gene_universe)}
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
#' and Avi M. (2013) Enrichr: Interactive and Collaborative Html5 Gene List Enrichment
#' Analysis Tool. \emph{BMC Bioinformatics} 14 (1). BioMed Central: 1.
#'
#' Alexa A. and Rahnenfuhrer J. (2009). topGO: Enrichment Analysis for Gene Ontology.
#'  R package version 2.28.0.

enrichment_analysis <-function(subnet, mode=NULL, gene_universe){

  # Checking function arguments
  if (missing(subnet))
    stop("Need to specify the subnetwork obtained from the PCSF algorithm.")
  if (class(subnet)[1] != "PCSF" || class(subnet)[2] != "igraph")
    stop("The subnetwork must be a \"PCSF\" object derived from an \"igraph\" class.")
  if (!is.null(mode)){
    if(mode==1 && missing(gene_universe))
      stop("Need to specify a list of genes (vector of gene symbols) used as background in enrichment analysis by topGO package")
  }


  cat("  Performing enrichment analysis...\n\n")

  # Obtain clusters in the subnet using edge betweenness clustering algorithm from igraph package.
  clusters = cluster_edge_betweenness(subnet)

  # Perform ebrichment analysis for each cluster using EnrichR through its API or topGO.

  havingInternet <- function() {
    if (.Platform$OS.type == "windows") {
      ipmessage <- system("ipconfig", intern = TRUE)
    } else {
      ipmessage <- system("ifconfig", intern = TRUE)
    }
    validIP <- "((25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)[.]){3}(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)"
    any(grep(validIP, ipmessage))
  }

  internet_connection <- havingInternet()

  if(!is.null(mode)){
    if(mode==0){
      if(internet_connection){
        cat("  Enrichment is being performed by EnrichR (http://amp.pharm.mssm.edu/Enrichr) API ...\n")
        enrich = call_enr(clusters, mode = 0, gene_universe)
      }
      else{
        stop("There is no working Internet connection, perform your enrichment with topGO package with mode=1 by providing background gene list ...\n")
      }
    }
    else{
      cat("  Enrichment is being performed by topGO package ...\n")
      enrich = call_enr(clusters, mode = mode, gene_universe)
    }
  }
  else
    {
    if(internet_connection){
      cat("  Enrichment is being performed by EnrichR (http://amp.pharm.mssm.edu/Enrichr) API ...\n")
      enrich = call_enr(clusters, mode = 0, gene_universe)
    }
    else{
      stop("There is no working Internet connection, perform your enrichment with topGO package with mode=1 by providing background gene list ...\n")
    }
  }

    if('Compound'%in% V(subnet)$type){##then we have drugs!
        require(dplyr)
        comps=data.frame(Drug=V(subnet)$name[which(V(subnet)$type=='Compound')],
                         Cluster=clusters$membership[which(V(subnet)$type=='Compound')])%>%
            dplyr::group_by(Cluster)%>%
            dplyr::summarise(DrugsByBetweenness=paste(Drug,collapse=';'))

    }
    else{
        comps <-NULL
        }
  enrichment = enrich[[1]]
    enrichment_complete = enrich[[2]]

    novals<-which(unlist(sapply(enrich[[2]],function(x) is.null(dim(x)))))
    if(length(novals)>0)
        enrichment_complete <- enrichment_complete[-novals]
    enrichment_tab = do.call(rbind,lapply(c(1:length(enrichment_complete)),function(x) data.frame(Cluster=x,enrichment_complete[[x]])))
    more.than.two=which(sapply(enrichment_tab$Genes,function(x) length(unlist(strsplit(x,split=';')))>2))
    if(length(more.than.two)>0)
        enrichment_tab=enrichment_tab[more.than.two,]
    if(!is.null(comps))
        enrichment_tab = enrichment_tab%>%dplyr::left_join(comps,by='Cluster')

  # Add 'group" and 'title' attributes to subnet
  V(subnet)$group = clusters$membership
  V(subnet)$title = paste0("Cluster ",clusters$membership,": Enrichment analysis")
  for( i in 1:length(V(subnet))){
    V(subnet)$title[i] = paste0( V(subnet)$title[i], enrichment[[V(subnet)$group[i]]])
  }

  # Derive a "PCSFe" object from an "igraph" class.
  class(subnet) <- c("PCSFe", "igraph")
  # Combine the subnetwork and colplete enrichment analysis tables.
  output = list(subnet, enrichment_tab)
  names(output) = c("subnet", "enrichment")

  return (output)
}
