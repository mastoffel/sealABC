#' format allelic microsats to genepop
#'
#' outputs a data.frame with genepop format
#'
#' @details so far, the input is the standard format
#' for the bottleneck analysis, i.e. first column "id",
#' second column "pop", third column "cluster", all
#' following columns are loci (two columns per locus).
#' "id" and "cluster" will be deleted and "pop" will
#' be kept for the formatting.
#'
#'
#' @param x genotypes
#'
#' @author  Emily Humble (emily.lhumble@gmail.com)
#'          Martin Stoffel (martin.adam.stoffel@@gmail.com)
#'
#' @export
#'
#'
make_genepop <- function(x) {
    x <- x[-c(1,3)]
    colnames(x)[1] <- "POP"
    x[1] <- "pop1,"

    # get locus names
    loci <- matrix(c(colnames(x)[-1],colnames(x[1])))
    loci <- unique(gsub("_a|_b", "", loci))

    threedigits <- function(x) {
        x[is.na(x)] <- "000"
        xchar <- as.character(x)
        stringr::str_length(xchar)
        new_x <- stringr::str_pad(xchar, 3, pad = "0")
    }

    genotypes <- apply(x[-1], 2, threedigits)

    # one column per locus
    short_geno <- matrix(nrow = nrow(genotypes), ncol = ncol(genotypes) / 2)
    short_geno <- data.frame(short_geno, stringsAsFactors = FALSE)
    length_data <- ncol(genotypes)
    genotypes <- data.frame(genotypes, stringsAsFactors = FALSE)
    col_num <- 1
    for (i in seq(from = 1, to = length_data, by = 2)) {
        short_geno[col_num] <- paste0(genotypes[[i]], genotypes[[i+1]])
        col_num <- col_num + 1
    }

    # join everything together
    genotypes <- cbind(x[1], short_geno)
    pop <- matrix(ncol = length(genotypes),
        nrow = nrow(loci) + 1)

    pop[1,1] <- "Title Line: genotypes"
    pop[2:nrow(pop),1] <- loci
    pop <- as.data.frame(pop)
    colnames(genotypes) <- colnames(pop)
    file <- rbind(pop, genotypes)

}
