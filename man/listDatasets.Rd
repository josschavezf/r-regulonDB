% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regulonDB.R
\name{listDatasets}
\alias{listDatasets}
\title{listDatasets}
\usage{
listDatasets(database)
}
\arguments{
\item{database}{RSQLite object imported previously from an SQLite database file to R}
}
\description{
Prints all marts contained in the SQL object
}
\examples{
## import the SQLite database
## path_sqlite3 <- file.choose()
regulon <- DBI::dbConnect(RSQLite::SQLite(), path_sqlite3)
listDatasets(regulon)
listAttributes(regulon, "TF")
getAttr(
database = regulon,
mart = "GENE",
filters = c("operon_id", "operon_name"),
value = c("ECK125235983","cmoM-mukFEB" ),
cond = "and")
GetRegulatedGenesByTF(regulon, "Ara")
 dbDisconnect(regulon)
}
