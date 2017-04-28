#' Takes output from microsimr or an empirical microsatellite dataset and outputs summary statistics
#'
#' @param simd_data output from microsimr
#' @param type "microsimr" for microsimr output, "microsats" for allelic microsatellite format
#' @param by_pop name of population variable. If specified, all summary statistics will be calculated within populations
#' @param start_geno integer, specifying the first column with genotypes. If now specified, all column are expected to be genotypes
#' @param data_type if "empirical", the mratio will be calculated differently, see ?m_ratio
#' @param rarefaction if TRUE, calculates mean and sd number of alleles as mean of n_samp individuals and n_loc loci over n_boot bootstraps
#' @param nsamp number of samples to subssample
#' @param nloc number of loci to subsample
#' @param nboot number of bootstraps
#'
#' @examples
#' data(microsimr_data)
#' mssumstats(microsimr_data)
#'
#' data(microsat_data)
#' mssumstats(microsat_data, type = "microsats", data_type = "empirical")
#'
#' data(fur_seal)
#' mssumstats(fur_seal, type = "microsats", by_pop = "pop", start_geno = 4,
#'            data_type = "empirical")
#'
#' @export
#'
#'

mssumstats_old2 <- function(simd_data, type = c("microsimr", "microsats"), by_pop = NULL, start_geno = NULL,
    data_type = c("simulated", "empirical"), rarefaction = FALSE, nsamp = NULL, nloc = NULL, nboot = 1000) {

    # default to microsimr
    if (length(type) == 2) type <- type[1]

    # define how to calculate MRatio
    if (length(data_type) == 2) data_type  <- data_type[1]

    if (is.null(start_geno)) start_geno <- 1

    if (type == "microsimr") {
        # reshape a little bit
        simd_data <- simd_data[-c(1:2)]
        # genotypes <- splitstackshape::cSplit_f(simd_data, splitCols = names(simd_data), sep = "|")
        genotypes <- splitstackshape::cSplit(simd_data, splitCols = names(simd_data), sep = ".")
    } else if (type == "microsats") {
        genotypes <- simd_data
    }


    calculate_sumstats <- function(genotypes){
        # transform to gyptes object
        g_types_geno <- strataG::df2gtypes(as.data.frame(genotypes), ploidy = 2, id.col = NULL, strata.col = NULL, loc.col = 1)

        # summary statistics
        # number of alleles across loci
        num_alleles <- strataG::numAlleles(g_types_geno)
        num_alleles_mean <- mean(num_alleles, na.rm = TRUE)
        num_alleles_sd <- stats::sd(num_alleles, na.rm = TRUE)

        # with rarefaction

        if(rarefaction == TRUE){
            allelic_richness <- allele_rich(genotypes, nboot, nsamp, nloc)
            num_alleles_mean <- allelic_richness[, 1]
            num_alleles_sd <- allelic_richness[, 1]
        }
        # AR <- num_alleles_mean / log(nrow(genotypes))

        # mean allele size variance across loci
        calc_ranges <- function(x){
            geno <- unlist(genotypes[x:(x+1)])
            allele_size_var <- stats::sd(geno, na.rm = TRUE)
            allele_range <- range(geno, na.rm = TRUE)
            allele_range <- abs(allele_range[1] - allele_range[2])
            out <- c(allele_size_var, allele_range)
        }

        allele_stats <- vapply(seq(from = 1, to = ncol(genotypes), by = 2), calc_ranges, numeric(2))

        mean_allele_size_var <- mean(allele_stats[1, ], na.rm = TRUE)
        sd_allele_size_var <- stats::sd(allele_stats[1, ], na.rm = TRUE)

        mean_allele_range <- mean(allele_stats[2, ], na.rm = TRUE)
        sd_allele_range <- stats::sd(allele_stats[2, ], na.rm = TRUE)


        # measure similar to mRatio
        if (data_type == "simulated") {
            mratio <- strataG::mRatio(g_types_geno, by.strata = FALSE, rpt.size = 1)
            mratio_mean <- mean(mratio, na.rm = TRUE)
            mratio_sd <- stats::sd(mratio, na.rm = TRUE)
        } else if (data_type == "empirical") {
            mratio <- m_ratio(g_types_geno)
            # mratio <- mratio[mratio != 1] # not sure if makes sense
            mratio_mean <- mean(mratio, na.rm = TRUE)
            mratio_sd <- stats::sd(mratio, na.rm = TRUE)

        } else {
            stop("specify data_type = NULL or 'empirical' ")
        }

        #allelic richness
        # AR <- strataG::allelicRichness(g_types_geno)
        # AR_mean <- mean(AR, na.rm = TRUE)
        # AR_sd <- sd(AR, na.rm = TRUE)

        # expected heterozygosity
        exp_het <- strataG::exptdHet(g_types_geno)
        exp_het_mean <- mean(exp_het, na.rm = TRUE)
        exp_het_sd<- stats::sd(exp_het, na.rm = TRUE)

        # observed heterozygosity
        obs_het <- strataG::obsvdHet(g_types_geno)
        obs_het_mean <- mean(obs_het, na.rm = TRUE)
        obs_het_sd <- stats::sd(obs_het, na.rm = TRUE)

        # prop low allele freq
        afs <- strataG::alleleFreqs(g_types_geno)

        # # allele frequency variation
        # afs_sd <- mean(unlist(lapply(afs, function(x) sd(x, na.rm = TRUE))), na.rm = TRUE)

        # calculate proportion fo low frequency alleles for one locus
        prop_low_af <- function(afs){
            low_afs <- (afs[, "freq"] / sum(afs[, "freq"])) < 0.05
            prop_low <- sum(low_afs) / length(low_afs)
        }
        # and mean for all
        prop_low_afs <- unlist(lapply(afs, prop_low_af))
        prop_low_afs_mean <- mean(prop_low_afs, na.rm = TRUE)
        prop_low_afs_sd <- stats::sd(prop_low_afs, na.rm = TRUE)

        # het excess
        het_exc <- (exp_het - obs_het) < 0
        het_excess <- sum(het_exc) / length(het_exc)

        out <- data.frame(num_alleles_mean, num_alleles_sd,
            mean_allele_size_var,  sd_allele_size_var,
            mean_allele_range, sd_allele_range,
            exp_het_mean, exp_het_sd,
            obs_het_mean,  obs_het_sd,
            mratio_mean, mratio_sd,
            prop_low_afs_mean,  prop_low_afs_sd,
            het_excess)
    }

    # per population
    if (!is.null(by_pop)){
        genotypes[[by_pop]] <- as.character(genotypes[[by_pop]])
        all_sumstats <- lapply(names(table(genotypes[[by_pop]])), function(x) calculate_sumstats(genotypes[genotypes[[by_pop]] == x,
                                                                              start_geno:ncol(genotypes)]) )
        out <- as.data.frame(do.call(rbind, all_sumstats))
    # overall
    } else {
        out <- calculate_sumstats(genotypes[, start_geno:ncol(genotypes)])
    }

    out

}
