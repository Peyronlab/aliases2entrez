#' @export
#' @import RCurl
#' @import readr
#'
#' @title Update last HGNC correspondance database
#' @description
#' This function is used to update gene symbol correspondance from HGNC database
#' @name update_symbols
#' @param url user can provide url (default is NULL)
#' @examples
#' HGNC <- update_symbols()
#' @return returns a data.frame containing gene symbols with status, previous symbols and synonyms as well as their corresponding entrezIDs
#'

update_symbols <- function(url = NULL) {
  cat("Fetching url...\n")
  if (is.null(url)) {
    url <- "https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_eg_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
  }
  if (url.exists(url)) {
    cat("Accessing data...\n")
    HGNC <- read_delim(url, "\t", escape_double = FALSE, trim_ws = TRUE)
    HGNC <- data.frame(HGNC)
    cat("Checking validity...\n")
    expected <- c("Approved.symbol", "Status", "Previous.symbols", "Synonyms", "NCBI.Gene.ID")
    if (dim(HGNC)[2] != 5) {
      warning(paste(paste(expected, collapse = ", "), "	 should be in HGNC query"))
      stop(paste("HGNC table should contain the 5 columns
"))
    } else if (!identical(colnames(HGNC), expected)) {
      warning(paste(paste(expected, collapse = ", "), "	 not found in HGNC query"))
      stop(paste("HGNC table should contain the 5 columns
"))
    }

    name.and.status <- paste(HGNC$Approved.symbol, ifelse(HGNC$Status == "Symbol Withdrawn", "~withdrawn", ""), sep = "")
    HGNC$Approved.symbol <- name.and.status
    HGNC <- HGNC[, -c(2)]
    cat("done...\n")
    return(HGNC)
  } else {
    warning("url not found. Process aborted, please refer to documentation.")
  }
}
