#' Internal function \code{call_enr}
#'
#' This function is internally used to perform enrichment analysis employing ENRICHR API.
#'
#' @keywords internal
#'
#' @author Murodzhon Akhmedov
#'
#' @param clusters A subnetwork clustered using edge betweenness algorithm of \pkg{igraph} package.
#' @param mode A binary variable to choose the method for enrichment analysis, where 0 is for EnrichR API and 1 is for \pkg{topGO} package.
#' @param gene_universe A complete list of genes (vector of gene symbols) used as background in enrichment analysis by \pkg{topGO} package.

call_enr <- function(clusters, mode=0, gene_universe){

  if( mode == 0){

    # Enrichment analysis is performed with ENRICHR
    ENRICHR_ADDLIST = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    ENRICHR_EXPORT = 'http://amp.pharm.mssm.edu/Enrichr/export'

    # The list of databases to be checked in Enrichment Analysis
    database = c("GO_Biological_Process_2015","KEGG_2016", "Reactome_2016", "BioCarta_2016")

    # Enrichment results
    enrichment_result = as.list(1:length(clusters))
    enrichment_result_complete = as.list(1:length(clusters))

    # Perform Enrichment Analysis for each cluster in the forest
    for( a in 1:length(clusters)){

      # List of genes to be regusted for enrichment via ENRICHR API
      genes = clusters[[a]]
      request = list(list = paste(genes, collapse = "\n"))
      complete_request = POST(ENRICHR_ADDLIST, body = request)
      output =content(complete_request, "text", encoding = "ISO-8859-1")
      userListID  = strsplit(strsplit(output, "\n")[[1]][3], ": ")[[1]][2]

      response_collection=NULL

      # Request enrichment for each database and comnine them all
      for( b in 1:length(database)){

        # Gather an EXPORT URL and the Response
        url = paste0(ENRICHR_EXPORT, "?userListId=",userListID, "&backgroundType=", database[b])
        response = GET(url)
        response = content(response, "text",  encoding = "ISO-8859-1")
        response = strsplit(response, "\n")[[1]]
       #  if(length(unlist(response  ))==1)
       #       next
        response = lapply(response, function(x){sp = strsplit(x, "\t")[[1]]; return (sp)})

        # If the response contains some elements then combine it
        if(length(response)>1){
          x = length(response)-1
          m_resp = as.data.frame(matrix(0, nrow = x, ncol = length(response[[1]])))
          colnames(m_resp) = response[[1]]
          for(i in 1:x){
            m_resp[i,] = response[[i+1]]
          }
          response_collection = rbind(response_collection,m_resp)
        }
      }


        if(is.null(response_collection))
            next
      # Reorder the enrichment according to the "Adjusted P-value" and select the top 15 enrichments
      ordered_resp = data.frame(response_collection$`Term`, response_collection$`Adjusted P-value`, response_collection$`Combined Score`)
      ordered_resp = ordered_resp[order(ordered_resp[,2]),][1:15,]
      ordered_resp[,2] = signif(as.numeric(as.character(ordered_resp[,2])), 3)
      ordered_resp[,3] = signif(as.numeric(as.character(ordered_resp[,3])), 3)

      # Convert the enrichment table into HTML format in order to display it
      enrich = "<!DOCTYPE html> <html> <head> <style>
      table {font-family: arial, sans-serif; font-size: 10px; border-collapse: collapse;width: 100%;} td,
      th { border: 1px solid #dddddd; text-align: center; padding: 5px;}
      tr:nth-child(even) {background-color: #dddddd;}
      </style> </head> <body>
      <table> <tr>  <th>Term</th> <th>Adjusted P-value</th> <th>Combined Score</th> </tr>";
      for(i in 1:nrow(ordered_resp)){
        enrich = paste0(enrich, " <tr>")
        for(j in 1:ncol(ordered_resp)){
          enrich = paste0(enrich, "<td>",ordered_resp[i,j], "</td>")
        }
        enrich = paste0(enrich, "</tr> ")
      }
      enrich = paste0(enrich, "</table> </body> </html>")

      # Attach the Enrichment Analysis for the current cluster
      enrichment_result[[a]] = enrich
      enrichment_result_complete[[a]] = response_collection
      }

  } else{

    # Enrichment analysis is performed by topGO

    # Enrichment results
    enrichment_result = as.list(1:length(clusters))
    enrichment_result_complete = as.list(1:length(clusters))

    # Perform Enrichment Analysis for each cluster in the forest
    for( a in 1:length(clusters)){

      # List of genes to be regusted for enrichment by topGO
      genes = clusters[[a]]
      fg <- factor(as.integer(gene_universe %in% genes))
      names(fg) <- gene_universe
      tgData <- new("topGOdata", description = "simple_session", ontology = "BP",
                    allGenes=fg, nodeSize=15, annot=annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")

      resultFisher <- runTest(tgData, algorithm = "classic", statistic = "fisher")
      resultKS <- runTest(tgData, algorithm = "classic", statistic = "ks")

      res_table_top15 <- GenTable(tgData, classicFisher = resultFisher,
                         classicKS = resultKS,
                         orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 15)
      res_table_top1000 <- GenTable(tgData, classicFisher = resultFisher,
                            classicKS = resultKS,
                            orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 1000)


      # Reorder the enrichment according to the "Adjusted P-value" and select the top 15 enrichments
      #ordered_resp[,2] = signif(as.numeric(as.character(ordered_resp[,2])), 3)
      #ordered_resp[,3] = signif(as.numeric(as.character(ordered_resp[,3])), 3)

      # Convert the enrichment table into HTML format in order to display it
      enrich = "<!DOCTYPE html> <html> <head> <style>
      table {font-family: arial, sans-serif; font-size: 10px; border-collapse: collapse;width: 100%;} td,
      th { border: 1px solid #dddddd; text-align: center; padding: 5px;}
      tr:nth-child(even) {background-color: #dddddd;}
      </style> </head> <body>
      <table> <tr> <th>GO.ID</th> <th>Term</th> <th>Annotated</th> <th>Significant</th> <th>Expected</th> <th>classicFisher</th> <th>classicKS</th> </tr>";
      for(i in 1:nrow(res_table_top15)){
        enrich = paste0(enrich, " <tr>")
        for(j in 1:ncol(res_table_top15)){
          enrich = paste0(enrich, "<td>",res_table_top15[i,j], "</td>")
        }
        enrich = paste0(enrich, "</tr> ")
      }
      enrich = paste0(enrich, "</table> </body> </html>")

      # Attach the Enrichment Analysis for the current cluster
      enrichment_result[[a]] = enrich
      enrichment_result_complete[[a]] = res_table_top1000
      }

  }


  return (list(enrichment_result, enrichment_result_complete))

}
