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


m_ratio <- function(gtypes_geno, rpt_size = 8:2){

    calc_mratio <- function(freqs) {
        freqs <- freqs[, 1]

        if (length(freqs) == 1) {
            warning("only one allele")
            return(NA)
        }

        # if (length(freqs) <= 3) {
        #     return(NA)
        # }
        # if (!all.is.numeric(names(freqs))) {
        #     warning("allele names are non-numeric")
        #     return(NA)
        # }
        if (all(freqs == 0)) {
            warning("all frequencies are 0")
            NA
        }

        if (is.null(freqs)){
            freqs <- strataG::alleleFreqs(gtypes_geno)
        }

        else {
            freqs <- freqs[order(as.numeric(names(freqs)))]
            sizes <- as.numeric(names(freqs))
            size_diff <- diff(sizes)

            rpt_found <- FALSE
            all_rpts <- sort(rpt_size, decreasing = FALSE)
            all_possible <- matrix(data = NA, nrow = length(rpt_size), ncol = length(size_diff))
            count <- 1
            for (r in all_rpts) {
                all_possible[count, ] <- size_diff%%r == 0
                count <- count + 1

                if (all(size_diff%%r == 0)) {
                    rpt_found <- TRUE
                    break
                }
            }

            # instead of checking whether all alleles stick to a certain repeat size take the most
            # common repeat size as a divident

            if (!(rpt_found)) {
                r <- all_rpts[which.max(rowSums(all_possible))]
            }

            smallest <- sizes[which.min(sizes[freqs > 0])]
            largest <- sizes[which.max(sizes[freqs > 0])]
            n <- (largest - smallest)/r
            sum(freqs > 0)/(n + 1)
        }
    }
    freqs <- strataG::alleleFreqs(gtypes_geno, by.strata = FALSE)
    out <- sapply(freqs, calc_mratio)

}
