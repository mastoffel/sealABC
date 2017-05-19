#' Takes output from microsimr or an empirical microsatellite dataset and outputs summary statistics
#'
#' @param data microsatellites, whereby every locus is represented by two adjacent columns. Or microsimr output.
#' @param datatype defaults to "microsats" for allelic microsatellite format, "microsimr" for microsimr output
#' @param by_pop name of population variable. If specified, all summary statistics will be calculated within populations
#' @param start_geno integer, specifying the first column with genotypes. If now specified, all column are expected to be genotypes
#' @param mratio defaults to "strict". if "loose", the mratio will be calculated differently, see ?m_ratio
#' @param rarefaction if TRUE, calculates mean and sd number of alleles as mean of n_samp individuals and n_loc loci over n_boot bootstraps
#' @param nsamp number of samples to subssample
#' @param nloc number of loci to subsample
#' @param nboot number of bootstraps
#'
#' @examples
#' data(microsimr_data)
#' mssumstats(microsimr_data, datatype = "microsimr")
#'
#' data(microsat_data)
#' mssumstats(microsat_data, datatype = "microsats")
#'
#' data(fur_seal)
#' mssumstats(fur_seal, datatype = "microsats", by_pop = "pop", start_geno = 4)
#'
#' @export
#'
#'

mssumstats_old3 <- function(data, datatype = c("microsats", "microsimr"), by_pop = NULL, start_geno = NULL,
    mratio = c("strict", "loose"), rarefaction = FALSE, nsamp = NULL, nloc = NULL, nboot = 1000) {

    # default to microsimr
    if (length(datatype) == 2) datatype <- datatype[1]

    # define how to calculate MRatio
    if (length(mratio) == 2) mratio  <- mratio[1]

    if (is.null(start_geno)) start_geno <- 1

    if (datatype == "microsimr") {
        # reshape a little bit
        data <- data[-c(1:2)]
        # genotypes <- splitstackshape::cSplit_f(data, splitCols = names(data), sep = "|")
        genotypes <- splitstackshape::cSplit(data, splitCols = names(data), sep = ".")
    } else if (datatype == "microsats") {
        genotypes <- data
    }


    calculate_sumstats <- function(genotypes){

        # transform to gyptes object
        g_types_geno <- strataG::df2gtypes(as.data.frame(genotypes), ploidy = 2, id.col = NULL, strata.col = NULL, loc.col = 1)

        # summary statistics
        loci_summary <- as.data.frame(strataG::summarizeLoci(g_types_geno))
        summary_means <- do.call(rbind, (lapply(loci_summary, function(x) out <- data.frame(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)))))
        names(summary_means) <- c("mean", "sd")

        # number of alleles
        num_alleles_mean <- summary_means["num.alleles", "mean"]
        num_alleles_sd <- summary_means["num.alleles", "sd"]

        # expected heterozygosity
        exp_het_mean <- summary_means["exptd.heterozygosity", "mean"]
        exp_het_sd<- summary_means["exptd.heterozygosity", "sd"]

        # observed heterozygosity
        obs_het_mean <- summary_means["obsvd.heterozygosity", "mean"]
        obs_het_sd <- summary_means["obsvd.heterozygosity", "sd"]

        # allelic richness
        allel_richness_mean <- summary_means["allelic.richness", "mean"]
        allel_richness_sd <- summary_means["allelic.richness", "sd"]

        # with rarefaction
        # if(rarefaction == TRUE){
        #     allelic_richness <- allele_rich(genotypes, nboot, nsamp, nloc)
        #     num_alleles_mean <- allelic_richness[, 1]
        #     num_alleles_sd <- allelic_richness[, 1]
        # }
        # AR <- num_alleles_mean / log(nrow(genotypes))

        if (datatype == "microsats") {
            # for empirical data (with repeat sizes being not 1 as in the simulated microsimr data)
            # find the most common repeat size per locus
            rpt_size <- 2:8
            freqs <- strataG::alleleFreqs(g_types_geno)

            # what is the (most likely) repeat size?
            allele_sizes <- lapply(freqs, function(x) sort(as.numeric(row.names(x))))
            allele_size_diffs <- lapply(allele_sizes, function(x) diff(x))

            find_repeat_size <- function(size_diff, rpt_size){

                all_possible <- matrix(data = NA, nrow = length(rpt_size), ncol = length(size_diff))
                count <- 1
                for (r in rpt_size) {
                    all_possible[count, ] <- size_diff%%r == 0
                    count <- count + 1
                }
                # instead of checking whether all alleles stick to a certain repeat size take the most
                # common repeat size
                r <- rpt_size[which.max(rowSums(all_possible))]
            }
            repeat_size_per_locus <- as.numeric(lapply(allele_size_diffs, find_repeat_size, rpt_size))

        } else if (datatype == "microsimr") {
            repeat_size_per_locus <- 1 # for microsimr data
        }

        # mean allele size variance across loci
        calc_ss_per_loc <- function(x){
            geno <- unlist(genotypes[x:(x+1)])
            # allele size sd has to be divided by repeat number
            allele_size_sd <- stats::sd(geno, na.rm = TRUE)
            # kurtosis is independent of repeat number
            allel_size_kurtosis <- moments::kurtosis(geno, na.rm = TRUE) ## kurtosis is the same independent of repeat size
            # range has to be devided by repeat number too
            allele_range <- range(geno, na.rm = TRUE)
            allele_range <- abs(allele_range[1] - allele_range[2])

            out <- c(allele_size_sd, allele_range, allel_size_kurtosis)
        }

        allele_stats <- vapply(seq(from = 1, to = ncol(genotypes), by = 2), calc_ss_per_loc, numeric(3))

        ## kurtosis of allele sizes (see popabc software)
        mean_allele_size_kurtosis <- mean(allele_stats[3, ], na.rm = TRUE)
        sd_allele_size_kurtosis <- sd(allele_stats[3, ], na.rm = TRUE)
        median_allele_size_kurtosis <- stats::median(allele_stats[3, ], na.rm = TRUE)

        ## allele size variance (sd) (see popabc software), standardized by repeat length
        mean_allele_size_sd <- mean(allele_stats[1, ] / repeat_size_per_locus, na.rm = TRUE)
        sd_allele_size_sd <- stats::sd(allele_stats[1, ] / repeat_size_per_locus, na.rm = TRUE)
        median_allele_size_sd <- stats::median(allele_stats[1, ] / repeat_size_per_locus, na.rm = TRUE)

        ## allele range, standardized by repeat length
        mean_allele_range <- mean(allele_stats[2, ] / repeat_size_per_locus, na.rm = TRUE)
        sd_allele_range <- stats::sd(allele_stats[2, ] / repeat_size_per_locus, na.rm = TRUE)
        median_allele_range <- stats::median(allele_stats[2, ] / repeat_size_per_locus, na.rm = TRUE)

        # measure similar to mRatio
        if (mratio == "strict") {
            rpt_size <- 8:2
            if (datatype == "microsimr") rpt_size <- 1
            mratio <- strataG::mRatio(g_types_geno, by.strata = FALSE, rpt.size = rpt_size)
            mratio_mean <- mean(mratio, na.rm = TRUE)
            mratio_sd <- stats::sd(mratio, na.rm = TRUE)
        } else if (mratio == "loose") {
            mratio <- m_ratio(g_types_geno)
            # mratio <- mratio[mratio != 1] # not sure if makes sense
            mratio_mean <- mean(mratio, na.rm = TRUE)
            mratio_sd <- stats::sd(mratio, na.rm = TRUE)
        } else {
            stop("specify whether mratio is calculated strict or loose ")
        }


        # prop low allele freq
        afs <- strataG::alleleFreqs(g_types_geno)

        # calculate proportion fo low frequency alleles for one locus
        prop_low_af <- function(afs){
            # low_afs <- (afs[, "freq"] / sum(afs[, "freq"])) < 0.05
            low_afs <- afs[, "prop"] < 0.05
            prop_low <- sum(low_afs) / length(low_afs)
        }
        # and mean for all
        prop_low_afs <- unlist(lapply(afs, prop_low_af))
        prop_low_afs_mean <- mean(prop_low_afs, na.rm = TRUE)
        prop_low_afs_sd <- stats::sd(prop_low_afs, na.rm = TRUE)

        # # het excess
        # het_exc <- (exp_het - obs_het) < 0
        # het_excess <- sum(het_exc) / length(het_exc)

        out <- data.frame(
            num_alleles_mean, num_alleles_sd,
            exp_het_mean, exp_het_sd,
            obs_het_mean, obs_het_sd,
            allel_richness_mean, allel_richness_sd,
            mean_allele_size_kurtosis, sd_allele_size_kurtosis, median_allele_size_kurtosis,
            mean_allele_size_sd, sd_allele_size_sd, median_allele_size_sd,
            mean_allele_range, sd_allele_range, median_allele_range,
            mratio_mean, mratio_sd,
            prop_low_afs_mean,  prop_low_afs_sd)
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
