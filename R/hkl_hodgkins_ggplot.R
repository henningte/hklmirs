#' Plots spectra with identified peaks as ggplot2 plots

#' @details \code{hkl_hodgkins_ggplot} is a modified part of the code from
#' \url{https://github.com/shodgkins/FTIRbaselines}. It extract a defined set of
#' peaks from mid infrared spectra and creates a plot of these peaks as
#' \code{\link[ggplot2:ggplot2]{ggplot2-object}}.

#'@param x A data frame as returned by
#'\code{irpeat:::irp_content_klh_hodgkins_prepare}.
#'@return A list of \code{ggplot2-object} with a plot for each spectrum in
#'\code{x}.
#'
#'@source The code is modified from
#'\url{https://github.com/shodgkins/FTIRbaselines}.
#'
#'@export
hkl_hodgkins_ggplot <- function (x) {

    export <- NULL
    verbose <- FALSE
    make_plots <- TRUE
    data <- x

    irp_content_klh_hodgkins_minima <- function(region, default = NULL) {
      success <- TRUE
      offset.dy <- c(0, region$dy)[1:(length(region$dy))]
      just.crossed <- which(sign(region$dy * offset.dy) ==
                              -1 & region$d2y > 0)
      crossing <- c(just.crossed, just.crossed - 1)
      if (length(crossing) > 0) {
        index.min <- which(region$y == min(region$y[crossing]) &
                             region$x %in% region$x[crossing])
        index.min <- which(region$y == min(region$y[(index.min -
                                                       min(m, index.min - 1)):(index.min + min(m, nrow(region) -
                                                                                                 index.min))]) & region$x %in% region$x[(index.min -
                                                                                                                                           min(m, index.min - 1)):(index.min + min(m, nrow(region) -
                                                                                                                                                                                     index.min))])
        minimum <- region$x[index.min]
      }
      else {
        if (max(region$d2y) > 0) {
          minimum <- region$x[which(region$d2y == max(region$d2y))]
        }
        else {
          minimum <- default
          success <- FALSE
        }
      }
      if (length(minimum) > 1) {
        if (verbose) {
          message(paste(" Warning: multiple minima for",
                        default, "peak boundary:", paste(minimum, collapse = " "),
                        "\n"))
        }
      }
      return(list(minimum = minimum, success = success))
    }
    irp_content_klh_hodgkins_maxima <- function(region, default = NULL) {
      success <- TRUE
      offset.dy <- c(0, region$dy)[1:(length(region$dy))]
      just.crossed <- which(sign(region$dy * offset.dy) ==
                              -1 & region$d2y < 0)
      crossing <- c(just.crossed, just.crossed - 1)
      if (length(crossing) > 0) {
        index.max <- which(region$y == max(region$y[crossing]) &
                             region$x %in% region$x[crossing])
        index.max <- which(region$y == max(region$y[(index.max -
                                                       min(m, index.max - 1)):(index.max + min(m, nrow(region) -
                                                                                                 index.max))]) & region$x %in% region$x[(index.max -
                                                                                                                                           min(m, index.max - 1)):(index.max + min(m, nrow(region) -
                                                                                                                                                                                     index.max))])
        maximum <- region$x[index.max]
      }
      else {
        if (min(region$d2y) < 0) {
          maximum <- region$x[which(region$d2y == min(region$d2y))]
        }
        else {
          maximum <- default
          success <- FALSE
        }
      }
      if (length(maximum) > 1) {
        if (verbose) {
          message(paste(" Warning: multiple maxima for",
                        default, "peak:", paste(maximum, collapse = " "),
                        "\n"))
        }
      }
      return(list(maximum = maximum, success = success))
    }
    n <- 29
    m <- (n - 1)/2
    peaks <- c("carb", "arom15", "arom16", "trough16", "acids",
               "aliph28", "trough28", "aliph29")
    Wp <- data.frame(matrix(nrow = length(data), ncol = length(peaks),
                            dimnames = list(names(data), peaks)))
    W1 <- Wp
    W2 <- Wp
    Ap <- Wp
    A1 <- Wp
    A2 <- Wp
    Ab <- Wp
    Acorr <- Wp
    success.W1 <- Wp
    success.W2 <- Wp
    area.wholepeak <- Wp
    area.corrpeak <- Wp
    notes <- data.frame(aliph.type = rep(NA, length(data)), acids.type = rep(NA,
                                                                             length(data)), notes = rep(NA, length(data)), row.names = names(data))
    x <- as.numeric(row.names(data))
    for (i in seq_along(data)) {
      if (verbose) {
        message(paste0("Processing ", names(data)[i], " (",
                       i, " of ", length(data), ")\n"))
      }
      y <- data[, i]
      xydata <- data.frame(x = x, y = y)
      dy <- rep(NA, length(y))
      for (j in (m + 1):(length(dy) - m)) {
        dy[j] <- stats::coefficients(stats::lsfit(x[(j -
                                                       m):(j + m)], y[(j - m):(j + m)]))[2]
      }
      xydata <- cbind(xydata, dy)
      d2y <- rep(NA, length(y))
      for (j in (2 * m + 1):(length(d2y) - 2 * m)) {
        d2y[j] <- stats::coefficients(stats::lsfit(x[(j -
                                                        m):(j + m)], dy[(j - m):(j + m)]))[2]
      }
      xydata <- cbind(xydata, d2y)
      w1.result <- irp_content_klh_hodgkins_minima(xydata[x %in%
                                                            890:920, ], default = 905)
      w1 <- max(w1.result[[1]])
      a1 <- y[x == w1]
      W1[i, "carb"] <- w1
      A1[i, "carb"] <- a1
      success.W1[i, "carb"] <- w1.result[[2]]
      w2.result <- irp_content_klh_hodgkins_minima(xydata[x %in%
                                                            1150:1210, ], default = 1185)
      w2 <- min(w2.result[[1]])
      a2 <- y[x == w2]
      W2[i, "carb"] <- w2
      A2[i, "carb"] <- a2
      success.W2[i, "carb"] <- w2.result[[2]]
      w.wholepeak <- w2:w1
      a.wholepeak <- y[x %in% w.wholepeak]
      a.corrpeak <- a.wholepeak - ((a2 - a1) * (w.wholepeak -
                                                  w1)/(w2 - w1) + a1)
      area.wholepeak[i, "carb"] <- sum(a.wholepeak)
      area.corrpeak[i, "carb"] <- sum(a.corrpeak)
      Acorr[i, "carb"] <- max(a.corrpeak)
      Wp[i, "carb"] <- min(w.wholepeak[a.corrpeak == Acorr[i,
                                                           "carb"]])
      Ap[i, "carb"] <- a.wholepeak[w.wholepeak == Wp[i, "carb"]]
      Ab[i, "carb"] <- Ap[i, "carb"] - Acorr[i, "carb"]
      rm(w1.result, w2.result, w1, w2, a1, a2, w.wholepeak,
         a.wholepeak, a.corrpeak)
      w1.result <- irp_content_klh_hodgkins_minima(xydata[x %in%
                                                            1470:1500, ], default = 1485)
      w1 <- max(w1.result[[1]])
      a1 <- y[x == w1]
      W1[i, "arom15"] <- w1
      A1[i, "arom15"] <- a1
      success.W1[i, "arom15"] <- w1.result[[2]]
      w2.result <- irp_content_klh_hodgkins_minima(xydata[x %in%
                                                            1515:1540, ], default = 1535)
      w2 <- min(w2.result[[1]])
      a2 <- y[x == w2]
      W2[i, "arom15"] <- w2
      A2[i, "arom15"] <- a2
      success.W2[i, "arom15"] <- w2.result[[2]]
      w.wholepeak <- w2:w1
      a.wholepeak <- y[x %in% w.wholepeak]
      a.corrpeak <- a.wholepeak - ((a2 - a1) * (w.wholepeak -
                                                  w1)/(w2 - w1) + a1)
      area.wholepeak[i, "arom15"] <- sum(a.wholepeak)
      area.corrpeak[i, "arom15"] <- sum(a.corrpeak)
      Acorr[i, "arom15"] <- max(a.corrpeak)
      if (Acorr[i, "arom15"] == 0) {
        Wp[i, "arom15"] <- 1510
        Acorr[i, "arom15"] <- a.corrpeak[w.wholepeak == Wp[i,
                                                           "arom15"]]
        Ap[i, "arom15"] <- a.wholepeak[w.wholepeak == Wp[i,
                                                         "arom15"]]
        Ab[i, "arom15"] <- Ap[i, "arom15"] - Acorr[i, "arom15"]
      }
      else {
        Wp[i, "arom15"] <- min(w.wholepeak[a.corrpeak ==
                                             Acorr[i, "arom15"]])
        Ap[i, "arom15"] <- a.wholepeak[w.wholepeak == Wp[i,
                                                         "arom15"]]
        Ab[i, "arom15"] <- Ap[i, "arom15"] - Acorr[i, "arom15"]
      }
      rm(w1.result, w2.result, w1, w2, a1, a2, w.wholepeak,
         a.wholepeak, a.corrpeak)
      if (A2[i, "arom15"] != min(y[x %in% W2[i, "arom15"]:1590])) {
        w1.result <- list(minimum = x[y == min(y[x %in% W2[i,
                                                           "arom15"]:1590])], success = TRUE)
      }
      else {
        w1.result <- list(minimum = W2[i, "arom15"], success = success.W2[i,
                                                                          "arom15"])
      }
      w1 <- max(w1.result[[1]])
      a1 <- y[x == w1]
      W1[i, "arom16"] <- w1
      A1[i, "arom16"] <- a1
      success.W1[i, "arom16"] <- w1.result[[2]]
      W1[i, "trough16"] <- w1
      A1[i, "trough16"] <- a1
      success.W1[i, "trough16"] <- w1.result[[2]]
      W1[i, "acids"] <- w1
      A1[i, "acids"] <- a1
      success.W1[i, "acids"] <- w1.result[[2]]
      w2 <- min(x[x %in% 1760:1850 & y <= (min(y[x %in% 1760:1850]) +
                                             0.01 * max(y[x > 1800]))])
      w2.result <- list(minimum = w2, success = TRUE)
      a2 <- y[x == w2]
      W2[i, "arom16"] <- w2
      A2[i, "arom16"] <- a2
      success.W2[i, "arom16"] <- w2.result[[2]]
      W2[i, "trough16"] <- w2
      A2[i, "trough16"] <- a2
      success.W2[i, "trough16"] <- w2.result[[2]]
      W2[i, "acids"] <- w2
      A2[i, "acids"] <- a2
      success.W2[i, "acids"] <- w2.result[[2]]
      w.wholepeak <- w2:w1
      a.wholepeak <- y[x %in% w.wholepeak]
      a.corrpeak <- a.wholepeak - ((a2 - a1) * (w.wholepeak -
                                                  w1)/(w2 - w1) + a1)
      Acorr[i, "arom16"] <- max(a.corrpeak[w.wholepeak < 1660])
      if (Acorr[i, "arom16"] == 0) {
        Wp[i, "arom16"] <- 1630
        Acorr[i, "arom16"] <- a.corrpeak[w.wholepeak == Wp[i,
                                                           "arom16"]]
        Ap[i, "arom16"] <- a.wholepeak[w.wholepeak == Wp[i,
                                                         "arom16"]]
        Ab[i, "arom16"] <- Ap[i, "arom16"] - Acorr[i, "arom16"]
      }
      else {
        Wp[i, "arom16"] <- min(w.wholepeak[a.corrpeak ==
                                             Acorr[i, "arom16"]])
        Ap[i, "arom16"] <- a.wholepeak[w.wholepeak == Wp[i,
                                                         "arom16"]]
        Ab[i, "arom16"] <- Ap[i, "arom16"] - Acorr[i, "arom16"]
      }
      default.trough16 <- 1685
      default.acids <- 1725
      tr.ini <- irp_content_klh_hodgkins_minima(xydata[x %in%
                                                         max(1650, Wp[i, "arom16"] + 20):min(1740, (w2 - 20)),
      ], default = default.trough16)[[1]]
      if (all(dy[x %in% max(1650, Wp[i, "arom16"] + 20):min(1740,
                                                            (w2 - 20))] < 0)) {
        acids.result <- irp_content_klh_hodgkins_maxima(xydata[x %in%
                                                                 max(tr.ini, 1685):w2, ], default = default.acids)
        Wp[i, "acids"] <- min(acids.result[[1]])
        Acorr[i, "acids"] <- a.corrpeak[w.wholepeak == Wp[i,
                                                          "acids"]]
        Ap[i, "acids"] <- a.wholepeak[w.wholepeak == Wp[i,
                                                        "acids"]]
        Ab[i, "acids"] <- Ap[i, "acids"] - Acorr[i, "acids"]
        notes[i, "acids.type"] <- ifelse(acids.result[[2]],
                                         "shoulder", "no peak")
        rm(acids.result)
      }
      else {
        if (max(a.corrpeak[w.wholepeak > max(tr.ini, 1685)]) !=
            0) {
          Acorr[i, "acids"] <- max(a.corrpeak[w.wholepeak >
                                                max(tr.ini, 1685)])
          Wp[i, "acids"] <- min(w.wholepeak[a.corrpeak ==
                                              Acorr[i, "acids"]])
          Ap[i, "acids"] <- a.wholepeak[w.wholepeak ==
                                          Wp[i, "acids"]]
        }
        else {
          Ap[i, "acids"] <- max(a.wholepeak[w.wholepeak >
                                              max(tr.ini, 1685)])
          Wp[i, "acids"] <- min(w.wholepeak[a.wholepeak ==
                                              Ap[i, "acids"]])
          Acorr[i, "acids"] <- a.corrpeak[w.wholepeak ==
                                            Wp[i, "acids"]]
        }
        Ab[i, "acids"] <- Ap[i, "acids"] - Acorr[i, "acids"]
        notes[i, "acids.type"] <- "peak"
      }
      if (notes[i, "acids.type"] == "peak") {
        Acorr[i, "trough16"] <- min(a.corrpeak[w.wholepeak %in%
                                                 Wp[i, "arom16"]:Wp[i, "acids"]])
        Wp[i, "trough16"] <- min(w.wholepeak[a.corrpeak ==
                                               Acorr[i, "trough16"]])
      }
      else if (notes[i, "acids.type"] == "shoulder") {
        Wp[i, "trough16"] <- min(tr.ini)
        Acorr[i, "trough16"] <- a.corrpeak[w.wholepeak ==
                                             Wp[i, "trough16"]]
      }
      else {
        Wp[i, "trough16"] <- default.trough16
        Acorr[i, "trough16"] <- a.corrpeak[w.wholepeak ==
                                             Wp[i, "trough16"]]
      }
      Ap[i, "trough16"] <- a.wholepeak[w.wholepeak == Wp[i,
                                                         "trough16"]]
      Ab[i, "trough16"] <- Ap[i, "trough16"] - Acorr[i, "trough16"]
      if (Acorr[i, "trough16"] < 0 | Acorr[i, "acids"] < 0) {
        notes[i, "acids.type"] <- paste0(notes[i, "acids.type"],
                                         ", negative")
      }
      area.wholepeak[i, "arom16"] <- sum(a.wholepeak[w.wholepeak %in%
                                                       w1:Wp[i, "trough16"]])
      area.corrpeak[i, "arom16"] <- sum(a.corrpeak[w.wholepeak %in%
                                                     w1:Wp[i, "trough16"]])
      area.wholepeak[i, "acids"] <- sum(a.wholepeak[w.wholepeak %in%
                                                      (Wp[i, "trough16"] + 1):w2])
      area.corrpeak[i, "acids"] <- sum(a.corrpeak[w.wholepeak %in%
                                                    (Wp[i, "trough16"] + 1):w2])
      rm(w1.result, w2.result, w1, w2, a1, a2, w.wholepeak,
         a.wholepeak, a.corrpeak, tr.ini, default.trough16,
         default.acids)
      w1.result <- list(minimum = 2750, success = TRUE)
      w1 <- max(w1.result[[1]])
      a1 <- y[x == w1]
      W1[i, "aliph28"] <- w1
      A1[i, "aliph28"] <- a1
      success.W1[i, "aliph28"] <- w1.result[[2]]
      W1[i, "trough28"] <- w1
      A1[i, "trough28"] <- a1
      success.W1[i, "trough28"] <- w1.result[[2]]
      W1[i, "aliph29"] <- w1
      A1[i, "aliph29"] <- a1
      success.W1[i, "aliph29"] <- w1.result[[2]]
      w2.result <- irp_content_klh_hodgkins_minima(xydata[x %in%
                                                            2950:3050, ], default = 3000)
      w2 <- min(w2.result[[1]])
      a2 <- y[x == w2]
      W2[i, "aliph28"] <- w2
      A2[i, "aliph28"] <- a2
      success.W2[i, "aliph28"] <- w2.result[[2]]
      W2[i, "trough28"] <- w2
      A2[i, "trough28"] <- a2
      success.W2[i, "trough28"] <- w2.result[[2]]
      W2[i, "aliph29"] <- w2
      A2[i, "aliph29"] <- a2
      success.W2[i, "aliph29"] <- w2.result[[2]]
      w.wholepeak <- w2:w1
      a.wholepeak <- y[x %in% w.wholepeak]
      a.corrpeak <- a.wholepeak - ((a2 - a1) * (w.wholepeak -
                                                  w1)/(w2 - w1) + a1)
      if (max(w.wholepeak[a.wholepeak == max(a.wholepeak)]) %in%
          2915:w2) {
        if (!any(w.wholepeak[a.wholepeak == min(a.wholepeak[w.wholepeak %in%
                                                            2854:2913])] %in% c(2854, 2913))) {
          type <- "two separated peaks"
        }
        else if (sum(y[x %in% 2920:2850] - ((y[x == 2920] -
                                             y[x == 2850]) * ((2920:2850) - 2850)/(2920 -
                                                                                   2850) + y[x == 2850])) < 0) {
          type <- "two unseparated peaks"
        }
        else {
          type <- "one peak"
        }
      }
      else {
        type <- "one peak"
      }
      notes[i, "aliph.type"] <- type
      if (type == "one peak") {
        Acorr[i, "aliph29"] <- max(a.corrpeak)
        Wp[i, "aliph29"] <- max(w.wholepeak[a.corrpeak ==
                                              Acorr[i, "aliph29"]])
        Ap[i, "aliph29"] <- a.wholepeak[w.wholepeak == Wp[i,
                                                          "aliph29"]]
        Ab[i, "aliph29"] <- Ap[i, "aliph29"] - Acorr[i, "aliph29"]
        Acorr[i, "trough28"] <- Acorr[i, "aliph29"]
        Wp[i, "trough28"] <- Wp[i, "aliph29"]
        Ap[i, "trough28"] <- Ap[i, "aliph29"]
        Ab[i, "trough28"] <- Ab[i, "aliph29"]
        Acorr[i, "aliph28"] <- Acorr[i, "aliph29"]
        Wp[i, "aliph28"] <- Wp[i, "aliph29"]
        Ap[i, "aliph28"] <- Ap[i, "aliph29"]
        Ab[i, "aliph28"] <- Ab[i, "aliph29"]
      }
      else if (type == "two separated peaks") {
        Acorr[i, "aliph29"] <- max(a.corrpeak[w.wholepeak %in%
                                                2900:w2])
        Wp[i, "aliph29"] <- max(w.wholepeak[a.corrpeak ==
                                              Acorr[i, "aliph29"]])
        Ap[i, "aliph29"] <- a.wholepeak[w.wholepeak == Wp[i,
                                                          "aliph29"]]
        Ab[i, "aliph29"] <- Ap[i, "aliph29"] - Acorr[i, "aliph29"]
        Acorr[i, "trough28"] <- min(a.corrpeak[w.wholepeak %in%
                                                 2854:Wp[i, "aliph29"]])
        Wp[i, "trough28"] <- max(w.wholepeak[a.corrpeak ==
                                               Acorr[i, "trough28"]])
        Ap[i, "trough28"] <- a.wholepeak[w.wholepeak == Wp[i,
                                                           "trough28"]]
        Ab[i, "trough28"] <- Ap[i, "trough28"] - Acorr[i,
                                                       "trough28"]
        Acorr[i, "aliph28"] <- max(a.corrpeak[w.wholepeak %in%
                                                w1:Wp[i, "trough28"]])
        Wp[i, "aliph28"] <- min(w.wholepeak[a.corrpeak ==
                                              Acorr[i, "aliph28"]])
        Ap[i, "aliph28"] <- a.wholepeak[w.wholepeak == Wp[i,
                                                          "aliph28"]]
        Ab[i, "aliph28"] <- Ap[i, "aliph28"] - Acorr[i, "aliph28"]
      }
      else {
        Acorr[i, "aliph29"] <- max(a.corrpeak[w.wholepeak %in%
                                                2900:w2])
        Wp[i, "aliph29"] <- max(w.wholepeak[a.corrpeak ==
                                              Acorr[i, "aliph29"]])
        Ap[i, "aliph29"] <- a.wholepeak[w.wholepeak == Wp[i,
                                                          "aliph29"]]
        Ab[i, "aliph29"] <- Ap[i, "aliph29"] - Acorr[i, "aliph29"]
        Acorr[i, "trough28"] <- min(a.corrpeak[w.wholepeak %in%
                                                 2854:Wp[i, "aliph29"]])
        Wp[i, "trough28"] <- max(w.wholepeak[a.corrpeak ==
                                               Acorr[i, "trough28"]])
        Ap[i, "trough28"] <- a.wholepeak[w.wholepeak == Wp[i,
                                                           "trough28"]]
        Ab[i, "trough28"] <- Ap[i, "trough28"] - Acorr[i,
                                                       "trough28"]
        Wp[i, "aliph28"] <- 2850
        Acorr[i, "aliph28"] <- a.corrpeak[w.wholepeak ==
                                            Wp[i, "aliph28"]]
        Ap[i, "aliph28"] <- a.wholepeak[w.wholepeak == Wp[i,
                                                          "aliph28"]]
        Ab[i, "aliph28"] <- Ap[i, "aliph28"] - Acorr[i, "aliph28"]
      }
      area.wholepeak[i, "aliph28"] <- sum(a.wholepeak[w.wholepeak %in%
                                                        w1:Wp[i, "trough28"]])
      area.corrpeak[i, "aliph28"] <- sum(a.corrpeak[w.wholepeak %in%
                                                      w1:Wp[i, "trough28"]])
      area.wholepeak[i, "aliph29"] <- sum(a.wholepeak[w.wholepeak %in%
                                                        (Wp[i, "trough28"] + 1):w2])
      area.corrpeak[i, "aliph29"] <- sum(a.corrpeak[w.wholepeak %in%
                                                      (Wp[i, "trough28"] + 1):w2])
      rm(w1.result, w2.result, w1, w2, a1, a2, w.wholepeak,
         a.wholepeak, a.corrpeak, type)
      rm(xydata, y, dy, d2y)
    }
    if (verbose) {
      message(paste("Spectra with baselines and peak locations will be shown for all",
                    length(data), "samples."))
    }
    col.peaks <- rep(NA, length(peaks))
    names(col.peaks) <- peaks
    col.peaks["carb"] <- "blue"
    col.peaks["arom15"] <- "red"
    col.peaks["arom16"] <- "purple3"
    col.peaks["trough16"] <- "deepskyblue3"
    col.peaks["acids"] <- "forestgreen"
    col.peaks["aliph28"] <- "darkorange3"
    col.peaks["trough28"] <- "gold4"
    col.peaks["aliph29"] <- "orangered3"
    lty.peaks <- ifelse(peaks == "trough16" | peaks == "trough28",
                        "dashed", "solid")
    names(lty.peaks) <- peaks
    if (make_plots) {
      lapply(seq_along(data), function(i) {
        seg_data <-
          tibble::tibble(
            W1 = W1[i, ] %>% unlist(),
            W2 = W2[i, ] %>% unlist(),
            A1 = A1[i, ] %>% unlist(),
            A2 = A2[i, ] %>% unlist(),
            Wp = Wp[i, ] %>% unlist(),
            Ap = Ap[i, ] %>% unlist(),
            Ab = Ab[i, ] %>% unlist(),
            col.peaks = col.peaks,
            lty.peaks = lty.peaks
          )

        label_data <-
          tibble::tibble(
            x = as.numeric(rownames(data)),
            y = data[, i, drop = TRUE],
            label_me = y %in% seg_data$Ap,
            label = ""
          )
        label_data$label[label_data$label_me] <- rev(names(col.peaks))

        ggplot(seg_data) +
          geom_path(data = data, aes(x = as.numeric(rownames(data)),
                                     y = data[, i, drop = TRUE])) +
          geom_segment(aes(x = W1, y = A1, xend = W2, yend = A2),
                       colour = "green") +
          geom_segment(aes(x = Wp, y = Ap, xend = Wp, yend = Ab,
                           colour = names(col.peaks),
                           linetype = names(lty.peaks))) +
          scale_colour_manual(values = seg_data$col.peaks) +
          geom_segment(aes(x = W1, y = A1, xend = W1, yend = 0),
                       colour = "gray50",
                       linetype = "dotted") +
          geom_segment(aes(x = W2, y = A2, xend = W2, yend = 0),
                       colour = "gray50",
                       linetype = "dotted") +
          theme(legend.position = "none") +
          labs(x = expression("Wavenumber ["*cm^{-1}*"]"),
               y = "Absorbance") +
          geom_text_repel(data = label_data, aes(x = x, y = y, label = label), seed = seed, nudge_y = 0.00005, segment.colour = "gray50")
      })
    }
  }

# function to nicely format numbers in text
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
