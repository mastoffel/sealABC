#' read multiple excel sheets
#'
#' loads multiple data frames into list from xlsx and assigns the
#' tab names from the spreadsheet to the list elements
#'
#'
#' @param data_path path to excel file with multiple data.frames
#'
#' @export


read_excel_sheets <- function(data_path) {

    # sheet numbers to load
    dataset_names <- readxl::excel_sheets(data_path)

    load_dataset <- function(dataset_names) {
        readxl::read_excel(data_path, sheet = dataset_names)
    }
    # load all datasets
    all_seals <- lapply(dataset_names, load_dataset)
    names(all_seals) <- dataset_names
    all_seals
}
