#' regulonDB
#'
#' @description Management of SQLite Databases and Query of Fields tools
#' @param database RSQLite object imported previously from an SQLite database file to R
#' @param mart Dataset table into the database, use the syntax "dataset"
#' @param filters A character vector with filters with the restrictions of the query restriction of the query
#' @param value A character vector with the corresponding values for filters
#' @param cond Conditional when more than one filter is used "and", "or","not
#' @import RSQLite, DBI, dplyr, dbplyr
#' @examples ## import the SQLite database
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
	## TO-DO Check if the LIKE operator is the right syntax for extracting concrete TF info (should "=" be used instead)
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

#' @export
### TO-DO
## () improve sql request, conections, and FIND if connections needs to be closed to avoid the warnings that say: "Closing open result set, pending rows"
## () design unit test
## () dont create every posible data table format, and then print the one requested; instead, add conditionals to know which format was requested, and only execute the corresponding code blocks if needed.

## Function to retrieve data for Gene-TF pairs in regulonDB SQLite database

##########
getGeneRegulation <- function( database, genes, format, output.type  ) {
  ## Output.type has no stable definition an is not implemented yet. Near future implementation after more design thinking
  ## Create base stringfor the query request
  query_cmd <- "SELECT * FROM GENE WHERE"

  ## loop to dynamically construct the final query request
  for (i in 1:length(genes)) {
    if (i == 1) {
      query_cmd <- paste0(query_cmd," name = '",genes[i],"'")
    } else {
      query_cmd <- paste0(query_cmd," or"," name = '",genes[i],"'")
    }
  }

  ## Create the query instruction
  query_cmd <- dbSendQuery(database, query_cmd)
  ## request the query, keep only columns name and gene_tf, and store it in a temporary df
  tempdf <- dbFetch(query_cmd)[,c("name","gene_tf")]
  ## change header names
  colnames(tempdf) <- c("genes","regulators")

  ## create a third temporal column to store useful regulator ID data (this wont be showed in the user output)
  tempdf$regulator_id <- ""

  ## create another empty df to populate with one line per regulator, for concantenating transient results
  tempdf3 <- tempdf[0,]

  ## Vectorize each cell in the regulators column
  ## operate in a for loop
  for (j in 1:nrow(tempdf)) {
    ## replace tabs by spaces
    tempvec0 <- gsub("\t", " ", tempdf[j,"regulators"])
    ## split every term in from the original regulator values
    tempvec0 <- unlist(strsplit(tempvec0, split = c(" ")))

    ## Eliminate the ECKnnnnnnn values, to keep only regulator Names
    tempvec1 <- tempvec0[!grepl("ECK", x = tempvec0)]

    ## find the ECKnnnnnnn values, to keep only regulator ECK id's
    tempvec2 <- tempvec0[grepl("ECK", x = tempvec0)]

    ## create an empty df to transiently populate with one line per regulator
    tempdf2 <- tempdf[0,]

    ## For every element in the regulators vector, print one line of the original tempdf in the second tempdf
    ## n keeps count of the vector of regulators, and i keeps track of the regulated gene beign evaluated
    for (n in 1:length(tempvec1)) {
      ## populate gene names
      tempdf2[n,"genes"] <- tempdf[j,"genes"]
      ## populate one by one, the regulator factors detected for the gene
      tempdf2[n,"regulators"] <- tempvec1[n]
      ## populate one by one, the ID for regulator factors detected for the gene
      tempdf2[n,"regulator_id"] <- tempvec2[n]
    }

    tempdf3 <- rbind(tempdf3,tempdf2)

  }

  ## Recover effect data for every gene:regulator pair
  ## Said data is in NETWORK table

  ##Create empty column to store effect
  tempdf3$effect <- ""

  for (k in 1:nrow(tempdf3)) {
    ## check if regulator value is NA and assign NA to effect function
    if ( is.na(tempdf3[k,"regulators"])) {
      tempdf3[k,"effect"] <- NA
      ## else, extract the effect value for the gene-TF pair
    } else {
      ## Create full string for the query request
      ##Query should ask for the gene and its regulator
      query_cmd <- paste0("SELECT * FROM NETWORK WHERE regulator_id = '",
                          tempdf3[k,"regulator_id"],
                          "' and regulated_name = '",
                          tempdf3[k,"genes"],
                          "' and network_type = 'TF-GENE'")
      ## Create the query instruction
      query_cmd <- dbSendQuery(database, query_cmd)
      ## request the query, keep only columns effect
      tempdf3[k,"effect"] <- dbFetch(query_cmd)[,c("effect")]
    }
  }

  ## translate string values of effect to +,-,+- codes
  tempdf3[,"effect"] <- gsub("dual","+-",tempdf3[,"effect"])
  tempdf3[,"effect"] <- gsub("activator","+",tempdf3[,"effect"])
  tempdf3[,"effect"] <- gsub("repressor","+-",tempdf3[,"effect"])

  ## remove the regulator_id not originally requested by the design (but used to query the SQLite db)
  tempdf3 <- tempdf3[,c("genes","regulators","effect")]

  ###
  ## transform results to one row per gene results

  ## start an empty temp df to store transient data
  tempdf4 <- tempdf3[0,c("genes","regulators")]

  ## Define a vector with every unique gene in the data set
  unique_genes <- unique(tempdf3[,"genes"])

  ## Start a for loop to iteratively generate the data row for every gene
  for (l in 1:length(unique_genes)) {

    ## create an internal df to perform concatenations more readably
    tempdf <- tempdf3[ tempdf3$genes == unique_genes[l],]

    ## Store data per row, per gene
    tempdf4[l,"genes"] <- unique_genes[l]
    ## concatenate columns and paste collapse into a single string, comma separated
    tempdf4[l,"regulators"] <- paste0(tempdf$regulators,"(",tempdf$effect,")", collapse = ",")

  }

  ## Pass the multi row format to table format
  tempdf5 <- spread(tempdf3, regulators, effect)

  ## Up to this points, two data objects have been generated:
  ## table 3 is the multirow format requested
  ## table 4 is the onerow format requested
  ## table 5 is the table format requested
  ## depending on the requested format, one of these will be printed
  if ( format == "multirow" ) {
    tempdf3
  } else if (format == "onerow" ){
    tempdf4
  } else if (format == "table" ){
    tempdf5
  }
}
