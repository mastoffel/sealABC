#' format microsat data from allelic to locus format
#'
#'
#'
#' @details takes allelic microsatellite data (two columns per locus, individuals in rows) and
#' transforms it into a one column per locus format, whereby the alleles are concetenated
#' with a given separator.
#'
#'
#' @param genotypes microsatellite data.frame in allelic format (two columns per locus).
#' @param sep how should the two alleles be concatenated? (defaults to "")
#'
#' @author Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'
#' @export
#'
#'
#'
allele2locus <- function(genotypes, sep = "") {

    # one column per locus
    short_geno <- matrix(nrow = nrow(genotypes), ncol = ncol(genotypes) / 2)
    short_geno <- data.frame(short_geno, stringsAsFactors = FALSE)
    length_data <- ncol(genotypes)
    genotypes <- data.frame(genotypes, stringsAsFactors = FALSE)
    col_num <- 1
    for (i in seq(from = 1, to = length_data, by = 2)) {
        short_geno[col_num] <- paste0(genotypes[[i]], sep, genotypes[[i+1]])
        col_num <- col_num + 1
    }

    short_geno

}
