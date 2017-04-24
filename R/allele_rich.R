#' Calculates allelic richness with rarefaction
#'
#' @param genotypes microsatellite genotypes in two-column per locus format
#' @param nboot number of bootstraps for rarefaction
#' @param nsamp number of individuals to resample, usually is smaller than the smallest population to compare
#' @param nloc number of loci to resample, usually is smaller than the number of loci of the population with
#' the fewest loci
#'
#' @examples
#' data(microsat_data)
#' out <- allele_rich(microsat_data, nboot = 5, nsamp = 5, nloc = 5)
#' out
#'
#' @export
#'
#'


allele_rich <- function(genotypes, nboot = 1000, nsamp = NULL, nloc = NULL){

   genotypes <- as.data.frame(genotypes)

   if(nsamp > nrow(genotypes)) stop("nsamp has to be smaller than the number of samples")
   if(nloc > ncol(genotypes) / 2) stop("nloc has to be smaller than the number of loci")

   boot_geno <- function(niter, genotypes, nsamp, nloc){
       # subset
       inds <- sample(1:nrow(genotypes), nsamp)
       loci <- sample(seq(from = 1, to = ncol(genotypes), by = 2), nloc)
       loci <- sort(c(loci, loci + 1))
       geno <- genotypes[inds, loci]

       # transform to gyptes object
       g_types_geno <- strataG::df2gtypes(geno, ploidy = 2, id.col = NULL, strata.col = NULL, loc.col = 1)

       # summary statistics
       # number of alleles across loci
       num_alleles <- strataG::numAlleles(g_types_geno)
       num_alleles_mean <- mean(num_alleles, na.rm = TRUE)
       num_alleles_sd <- stats::sd(num_alleles, na.rm = TRUE)

       out <- data.frame(num_alleles_mean, num_alleles_sd)
   }

   all_allele_rich <- lapply(1:nboot, boot_geno, genotypes, nsamp, nloc)
   all_allele_rich_df <- as.data.frame(do.call(rbind, all_allele_rich))

   mean_num_alleles <- mean(all_allele_rich_df[[1]])
   sd_num_alleles <- mean(all_allele_rich_df[[2]])

   out <- data.frame(mean_num_alleles, sd_num_alleles)
}
