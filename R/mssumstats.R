#' takes output from microsimr and outputs summary statistics
#'
#' @param simd_data output from microsimr
#' @param type "microsimr" for microsimr output, "microsats" for allelic microsatellite format
#'
#' @example
#'
#'
#'
#'
#'
#' @export
#'
#'


mssumstats <- function(simd_data, type = c("microsimr", "microsats")) {

    if (type == "microsimr") {
        # reshape a little bit
        simd_data <- simd_data[-c(1:2)]
        # genotypes <- splitstackshape::cSplit_f(simd_data, splitCols = names(simd_data), sep = "|")
        genotypes <- splitstackshape::cSplit(simd_data, splitCols = names(simd_data), sep = ".")
    } else {
        genotypes <- simd_data
    }


    g_types_geno <- strataG::df2gtypes(genotypes, ploidy = 2, id.col = NULL, strata.col = NULL,
                                       loc.col = 1)

    # summary statistics
    # according to DIYabc
    # number of alleles across loci
    num_alleles <- strataG::numAlleles(g_types_geno)
    num_alleles_mean <- mean(num_alleles, na.rm = TRUE)
    num_alleles_sd <- sd(num_alleles, na.rm = TRUE)

    # allele size variance (actually, sd) across loci
    # and allelic range
    allele_range_size <- sapply(seq(from = 1, to = ncol(genotypes), by = 2),
        function(x) {
            allele_sd <- sd(unlist(genotypes[x:(x+1)]), na.rm = TRUE)
            allele_range <- range(unlist(genotypes[x:(x+1)]), na.rm = TRUE)
            allele_range <- abs(allele_range[1] - allele_range[2])
            out <- c(allele_sd, allele_range)
            })

    mean_allele_size_sd <- mean(allele_range_size[1, ], na.rm = TRUE)
    sd_allele_size_sd <- sd(allele_range_size[1, ], na.rm = TRUE)

    mean_allele_range <- mean(allele_range_size[2, ], na.rm = TRUE)
    sd_allele_range <- sd(allele_range_size[2, ], na.rm = TRUE)


    # expected heterozygosity
    exp_het <- strataG::exptdHet(g_types_geno)
    exp_het_mean <- mean(exp_het, na.rm = TRUE)
    exp_het_sd<- sd(exp_het, na.rm = TRUE)

    # observed heterozygosity
    obs_het <- strataG::obsvdHet(g_types_geno)
    obs_het_mean <- mean(obs_het, na.rm = TRUE)
    obs_het_sd <- sd(obs_het, na.rm = TRUE)


    out <- data.frame(
        num_alleles_mean = num_alleles_mean, num_alleles_sd = num_alleles_sd,
        mean_allele_size_sd = mean_allele_size_sd, sd_allele_size_sd = sd_allele_size_sd,
        mean_allele_range = mean_allele_range, sd_allele_range,
        exp_het_mean = exp_het_mean,  exp_het_sd =  exp_het_sd,
        obs_het_mean = obs_het_mean,  obs_het_sd =  obs_het_sd)

    out

}
