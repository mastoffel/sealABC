#' takes list of named data.frames and writes to xls with respective names as sheet names
#'
#' @param dflist  list of named data.frames
#' @param file_name file name for xls file
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

    lapply(names(dflist), list_to_df, dflist, envir)
    WriteXLS::WriteXLS(names(all_seals_new), ExcelFileName = "file_name")

}



