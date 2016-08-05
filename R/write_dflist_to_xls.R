#' takes list of named data.frames and writes to xls with respective names as sheet names
#'
#' @param simd_data output from microsimr
#' @param type "microsimr" for microsimr output, "microsats" for allelic microsatellite format
#'
#'
#'
#'
#'
#'
#' @export
#'
#'

write_dflist_to_xls <- function(dflist, file_name){

    envir <- environment()
    list_to_df <- function(species, dfs, envir){
        assign(species, dfs[[species]], envir)
    }

    lapply(names(all_seals_new), list_to_df, all_seals_new, envir)
    WriteXLS::WriteXLS(names(all_seals_new), ExcelFileName = "file_name")

}



