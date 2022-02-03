#' Formats numbers for text output
#'
#' @param x A numeric value.
#' @param capitalized A logical value indicating if the output should start with
#' a capital letter (\code{TRUE}) or not (\code{FALSE}).
#' @return A character value with the formatted value from \code{x}.
#' @export
format_numbers <- function(x, capitalized = FALSE) {
  stopifnot(is.numeric(x) && length(x) == 1)
  stopifnot(is.logical(capitalized) && length(capitalized) == 1)
  if(x %% 1 == 0 && x <=10) {
    x <- c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten")[[x]]
    if(capitalized) {
      x <- paste0(toupper(stringr::str_sub(x, 1, 1)), stringr::str_sub(x, 2, nchar(x)))
    }
  }
  x
}
