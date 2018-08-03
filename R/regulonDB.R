#' regulonDB
#'
#' @description Management of SQLite Databases and Query of Fields tools
#' Depends of \RSQLite, \DBI, \dplyr, \dbplyr
#' @param database RSQLite object imported previously from an SQLite database file to R
#' @param mart Dataset table into the database, use the syntax "dataset"
#' @param filters A character vector with filters with the restrictions of the query restriction of the query
#' @param value A character vector with the corresponding values for filters
#' @param cond Conditional when more than one filter is used "and", "or","not
#' @examples
#' @export
#' ## import the SQLite database
#' ## path_sqlite3 <- file.choose()
#' regulon <- DBI::dbConnect(RSQLite::SQLite(), path_sqlite3)
#' listDatasets(regulon)
#' listAttributes(regulon, "TF")
#' getAttr(
#' database = regulon,
#' mart = "GENE",
#' filters = c("operon_id", "operon_name"),
#' value = c("ECK125235983","cmoM-mukFEB" ),
#' cond = "and")
#' GetRegulatedGenesByTF(regulon, "Ara")
#'  dbDisconnect(regulon)

#Install packages if needed

if (!require("RSQLite")) {
  install.packages("RSQLite", ask =FALSE)
  library(RSQLite)
}
if (!require("dplyr")) {
  install.packages("dplyr", ask =FALSE)
  library(dplyr)
}
if (!require("dbplyr")) {
  install.packages("dbplyr", ask =FALSE)
  library(dbplyr)
}
if (!require("DBI")) {
  install.packages("DBI", ask =FALSE)
  library(DBI)
}

#' @export
listDatasets <- function(database){
  dbListTables(database)
}

#' @export
listAttributes <- function(database, mart) {
  tryCatch(
    {
      regulonAttributes = DBI::dbListFields(database, mart)
      regulonAttributes
    },
    error=function(cond) {
      message("Error: Database or mart does not exist")
    }
  )
}

#' @export
GetRegulatedGenesByTF <- function(database,tf){
  genes_regulon <- tbl(database, "GENE")
  tf_table <- genes_regulon %>%
    select(name, gene_tf)
  result <- tf_table %>% filter( dplyr::sql( paste0("gene_tf LIKE '%",tf,"%' ") )) %>% select("name")
  tibble_result <- collect(select(result, name))[,1]
  as.matrix(tibble_result)

}

#' @export
getAttr <- function(database, mart, filters, value, cond ){
  query_cmd <- "SELECT * FROM "
  for (i in 1:length(filters)) {
    if (i == 1) {
      query_cmd <- paste0(query_cmd,mart," WHERE ",filters[i]," = '",value[i],"' ")
    } else {
      query_cmd <- paste0(query_cmd,cond," ",filters[i]," = '",value[i],"' ")
    }
  }
  res <- DBI::dbSendQuery(database, query_cmd)
  attributes <- DBI::dbFetch(res)
  attributes
}


