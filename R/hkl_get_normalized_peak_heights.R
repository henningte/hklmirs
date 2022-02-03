#' Extracts normalized peak heights for a set of peaks from mid infrared spectra

#'@param x An object of class \code{\link[ir:ir_new_ir]{ir}}.
#'@param ... Further arguments passed to
#'\code{pmird:::irp_content_klh_hodgkins_main}.
#'@param scale A logical value indicating if peak heights should be
#'z-transformed (\code{TRUE}) or not (\code{FALSE}), using
#'\code{\link{scale}}.
#'@return A data frame with a row for each spectrum in \code{x} and a column for
#'each peak extracted by \code{pmird:::irp_content_klh_hodgkins_main}.
#'@export
hkl_get_normalized_peak_heights <- function(x, ..., do_scale = FALSE) {

  stopifnot(inherits(x, "data.frame"))
  stopifnot(length(do_scale) == 1 && is.logical(do_scale))

  # extract peaks
  res <-
    irpeat:::irp_content_klh_hodgkins_main(data = x, ...)$norm.Acorr %>%
    as.data.frame()

  # define arom15arom16 peak
  res <-
    res %>%
    dplyr::mutate(arom15arom16 = arom15 + arom16)

  # scale
  if(do_scale) {
    res <-
      res %>%
      purrr::map_dfc(scale, center = TRUE, scale = TRUE)
  }

  res

}
