#' Approximate m-ratio
#'
#' Calculate Garza-Williamson M ratio (bottleneck) statistic for microsatellite data.
#' The function is based on mRatio() from the strataG package. However, instead of requiring
#' the same repeat length for all alleles within a locus, it calculates the m-ratio by deviding
#' the allele-range by the repeat lengths that is most common. Is is therefore an approximation
#' but maximises the amount of information used for non-optimally scored microsatellite datasets.
#'
#' @param gtypes_geno microsatellite genotypes formatted as gtypes object with strataG
#'
#'
#' @export
#'
#'
#'


m_ratio <- function(gtypes_geno){

    calc_mratio <- function(freqs) {
        freqs <- freqs[, 1]

        if (length(freqs) == 1) {
            warning("only one allele")
            return(NA)
        }
        # if (!all.is.numeric(names(freqs))) {
        #     warning("allele names are non-numeric")
        #     return(NA)
        # }
        if (all(freqs == 0)) {
            warning("all frequencies are 0")
            NA
        }

        if (is.null(afs)){
            afs <- strataG::alleleFreqs(gtypes_geno)
        }

        else {
            freqs <- freqs[order(as.numeric(names(freqs)))]
            sizes <- as.numeric(names(freqs))
            size_diff <- diff(sizes)

            # instead of checking whether all alleles stick to a certain repeat size, calculate the
            # mode as a proxy
            calc_mode <- function(x) {
                ux <- unique(x)
                ux[which.max(tabulate(match(x, ux)))]
            }

            mode_af <- calc_mode(size_diff)

            smallest <- sizes[which.min(sizes[freqs > 0])]
            largest <- sizes[which.max(sizes[freqs > 0])]
            n <- (largest - smallest)/mode_af
            sum(freqs > 0)/(n + 1)
        }
    }
    freqs <- alleleFreqs(gtypes_geno, by.strata = FALSE)
    out <- sapply(freqs, calc_mratio)

}
