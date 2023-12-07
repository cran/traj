#'@title Compute Measures for Identifying Patterns of Change in Longitudinal
#'  Data
#'
#'@description \code{Step1Measures} computes up to 26 measures for each
#'  longitudinal trajectory. See Details for the list of measures.
#'
#'@param Data a matrix or data frame in which each row contains the longitudinal
#'  data (trajectories).
#'@param Time either NULL, a vector or a matrix/data frame of the same dimension
#'  as \code{Data}. If a vector, matrix or data frame is supplied, its entries
#'  are assumed to be measured at the times of the corresponding cells in
#'  \code{Data}. When set to \code{NULL} (the default), the times are assumed
#'  equidistant.
#'@param ID logical. Set to \code{TRUE} if the first columns of \code{Data} and
#'  \code{Time} corresponds to an \code{ID} variable identifying the
#'  trajectories. Defaults to \code{FALSE}.
#'@param measures a vector containing the numerical identifiers of the measures
#'  to compute (see "Details" section below). The default, 1:23, corresponds to
#'  measures 1-23 and thus excludes the measures which require specifying a
#'  midpoint.
#'@param midpoint specifies which column of \code{Time} to use as the midpoint
#'  in measures 24-26. Can be NULL, an integer or a vector of integers of length
#'  the number of rows in \code{Time}. The default is NULL, in which case the
#'  midpoint is the time closest to the median of the Time vector specific to
#'  each trajectory.
#'@param x object of class trajMeasures.
#'@param object object of class trajMeasures.
#'@param ... further arguments passed to or from other methods.
#'
#'@return An object of class \code{trajMeasures}; a list containing the values
#'  of the measures, a table of the outliers which have been capped, as well as
#'  a curated form of the function's arguments.
#'
#'@details Each trajectory must have a minimum of 3 observations otherwise it
#'  will be omitted from the analysis.
#'
#'  The 26 measures and their numerical identifiers are listed below. Please
#'  refer to the vignette for the specific formulas used to compute them.
#'\enumerate{
#'\item  Range\cr
#'\item  Mean of the function\cr
#'\item  Functional standard deviation (SD)\cr
#'\item  Coefficient of variation (ratio of measure 3 to measure 2)\cr
#'\item  Overall change (initial value - final value)\cr
#'\item  Mean change per unit time\cr
#'\item  Overall change relative to initial value\cr
#'\item  Overall change relative to functional mean (ratio of measure 5 to measure 2)\cr
#'\item  Slope of the linear model\cr
#'\item  \eqn{R^2}: Proportion of variance explained by the linear model\cr
#'\item  Maximum value of the speed\cr
#'\item  Functional SD of the speed\cr
#'\item  Mean absolute speed\cr
#'\item  Maximum absolute speed\cr
#'\item  Maximum absolute speed relative to the functional mean (ratio of measure 14 to measure 2)\cr
#'\item  Maximum absolute speed relative to the slope (ratio of measure 14 to measure 9)\cr
#'\item  Functional SD of the speed relative to the slope\cr(ratio of measure 12 to measure 9)
#'\item  Mean acceleration\cr
#'\item  Mean absolute acceleration\cr
#'\item  Maximum of the absolute acceleration\cr
#'\item  Maximum of the absolute acceleration relative to the functional (ratio of measure 20 to measure 2)\cr
#'\item  Maximum of the absolute acceleration relative to the mean absolute speed (ratio of measure 20 to measure 13)\cr
#'\item  Mean absolute acceleration relative to the mean absolute speed (ratio of measure 19 to measure 13)\cr
#'\item  Early change relative to later change\cr
#'\item  Early change relative to overall change\cr
#'\item  Later change relative to overall change\cr
#'}
#'
#'@importFrom stats complete.cases coefficients lm median sd density
#'
#'@references Leffondre K, Abrahamowicz M, Regeasse A, Hawker GA, Badley EM,
#'  McCusker J, Belzile E. Statistical measures were proposed for identifying
#'  longitudinal patterns of change in quantitative health indicators. J Clin
#'  Epidemiol. 2004 Oct;57(10):1049-62. doi: 10.1016/j.jclinepi.2004.02.012.
#'  PMID: 15528056.
#'
#'  Nishiyama T, Improved Chebyshev inequality: new probability bounds with
#'  known supremum of PDF, arXiv:1808.10770v2 stat.ME
#'  https://doi.org/10.48550/arXiv.1808.10770
#'
#'@examples
#'\dontrun{
#'m1 = Step1Measures(trajdata, ID = TRUE, measures = 24:26, midpoint = NULL)
#'m2 = Step1Measures(trajdata, ID = TRUE, measures = 24:26, midpoint = 3)
#'
#'identical(s1$measures, s2$measures)
#'}
#'
#'@rdname Step1Measures
#'
#'@export
Step1Measures <-
  function (Data,
            Time = NULL,
            ID = FALSE,
            measures = 1:23,
            midpoint = NULL) {
    ###############################################################
    ##         Perform checks and uniformization of data         ##
    ###############################################################
    
    if (is.null(Time)) {
      if (ID == TRUE) {
        Time <- c(1:(ncol(Data) - 1))
      }
      if (ID == FALSE) {
        Time <- c(1:(ncol(Data)))
      }
    }
    
    Time.is.vector <- is.vector(Time)
    
    data <- Data
    data <- data.frame(data)
    names(data) <- NULL
    
    time <- Time
    time <- data.frame(time)
    names(time) <- NULL
    
    if (ID == TRUE) {
      IDvector <- data[, 1]
      data <- data[,-1]
      if (Time.is.vector == FALSE) {
        ID.time <- time[, 1]
        time <- time[,-1]
        if (identical(IDvector, ID.time) == FALSE) {
          stop("ID vector in Data differs from ID vector in Time.")
        }
      }
    } else{
      IDvector <- seq_len(nrow(data))
    }
    
    if (identical(dim(time), dim(data))) {
      if (identical(is.na(data), is.na(time)) == FALSE) {
        stop(
          "Data has cells with no corresponding time or Time has cells with no corresponding data."
        )
      }
      
      data2 <- data
      time2 <- time
      for (i in seq_len(nrow(data))) {
        NA.str_i <- is.na(data)[i,]
        w <- unname(which(NA.str_i == FALSE))
        if (length(w) < 3) {
          data2 <- data2[-i,]
          time2 <- time2[-i,]
          warning (
            paste(
              "Row ",
              i,
              " of Data contains less than 3 observations; it has been removed.",
              sep = ""
            )
          )
          IDvector <- IDvector[-i]
        } else if (!identical(w, seq_len(length(w)))) {
          stop(
            paste(
              "Row ",
              i,
              " of Data is not formatted correctly. Rows should be of the form X Y ... Z NA ... NA.",
              sep = ""
            )
          )
        }
      }
      data <- data2
      time <- time2
      
      for (i in seq_len(nrow(data))) {
        w2 <- unname(which(!is.na(time)[i,]))
        non.NA <- unname(unlist(time[i, w2]))
        if (!length(unique(non.NA)) == length(non.NA)) {
          stop(
            paste(
              "Line ",
              i,
              " of Time contains duplicates. The rows of Time should be strictly increasing sequences.",
              sep = ""
            )
          )
        }
        if (!identical(w2, order(non.NA))) {
          stop(
            paste(
              "The elements in row ",
              i,
              " of Time are not ordered chronologically.",
              sep = ""
            )
          )
        }
      }
      
    } else if (is.vector(Time) & (length(Time) == ncol(data))) {
      time <- unname(unlist(Time))
      if (TRUE %in% is.na(time)) {
        stop("If Time is supplied as a vector, it must not contain NA.")
      }
      if (!length(unique(time)) == length(time)) {
        stop(
          "The Time vector contains duplicates. The elements of Time should form a strictly increasing sequence."
        )
      }
      if (!identical(order(time), seq_len(length(time)))) {
        stop("The elements of Time are not ordered chronologically.")
      }
      if (length(time) < 3) {
        stop("Time must contain at least 3 elements.")
      }
      
      data2 <- data
      for (i in seq_len(nrow(data))) {
        NA.str_i <- is.na(data)[i,]
        w <- unname(which(NA.str_i == FALSE))
        if (length(w) < 3) {
          data2 <- data2[-i,]
          warning(
            paste(
              "Row ",
              i,
              " of Data contains less than 3 observations; it has been removed.",
              sep = ""
            )
          )
          if (ID == TRUE) {
            IDvector <- IDvector[-i]
          }
        }
      }
      data <- data2
      
      data2 <-
        unname(data.frame(matrix(
          NA, ncol = ncol(data), nrow = nrow(data)
        )))
      time2 <-
        unname(data.frame(matrix(
          NA, ncol = ncol(data), nrow = nrow(data)
        )))
      
      for (i in seq_len(nrow(data))) {
        w <- which(!is.na(data[i,]))
        data2[i, seq_len(length(w))] <- data[i, w]
        time2[i, seq_len(length(w))] <- Time[w]
      }
      data <- data2
      time <- time2
      
    } else{
      stop("The dimension of Data and Time are incompatible.")
    }
    
    
    ######################################################
    ##         Construct vector of "mid points"         ##
    ######################################################
    
    if (TRUE %in% (c(25:27) %in% measures)) {
      mid.position <- c()
      
      if (is.null(midpoint)) {
        median.time <- (max(time, na.rm = TRUE) + min(time, na.rm = TRUE)) / 2
        flag <- c()
        for (i in seq_len(nrow(data))) {
          v <- time[i,][!is.na(time[i,])]
          w <- which.min(abs(v - median.time))
          mid.position[i] <- w
          if (w == length(v) | w == 1) {
            flag <-
              c(flag, i)  #  the mid point can't be the first or the last point, so if this occurs, the ith row gets "flagged"
          }
        }
        if (length(flag) > 0) {
          mid.position <- mid.position[-flag]
          data <- data[-flag,]
          time <- time[-flag,]
          if (ID == TRUE) {
            Lines <- which(Data[, 1] %in% IDvector[flag])
          } else{
            Lines <- IDvector[flag]
          }
          if (length(Lines) == 1) {
            warning(
              paste(
                "When left blank, the 'midpoint' argument defaults to the observation time closests to the median time 0.5*(max(Time)-min(Time)), but this can't be either the first or last observation time. As a result, row ",
                Lines,
                " has been removed. To avoid this, consider excluding measures 24-26 from the analysis or providing custom 'midpoint' values.",
                sep = ""
              )
            )
          }
          if (length(Lines) > 1) {
            warning(
              paste(
                "When left blank, the 'midpoint' argument defaults to the observation time closests to the median time 0.5*(max(Time)-min(Time)), but this can't be either the first or last observation time. As a result, rows ",
                noquote(paste(Lines, collapse = ", ")),
                " have been removed. To avoid this, consider excluding measures 24-26 from the analysis or providing custom 'midpoint' values.",
                sep = ""
              )
            )
          }
          
          IDvector <- IDvector[-flag]
          
        }
      } else if (!is.vector(midpoint)) {
        stop(
          "'midpoint' should be either NULL or an integer/vector of integers corresponding to a column of Time."
        )
      } else{
        if (length(midpoint) == nrow(data)) {
          mid.position <- midpoint
        } else if (is.vector(Time) &
                   (length(Time) == ncol(data)) & length(midpoint) == 1) {
          mid.position <- rep(midpoint, nrow(data))
        } else{
          stop("'midpoint' does not have the correct format.")
        }
      }
      # Check to see if the midpoints are all greater than 1 but less than the number of observations
      for (i in seq_len(nrow(data))) {
        v <- time[i,][!is.na(time[i,])]
        if (mid.position[i] <= 1) {
          stop(
            paste(
              "Error in 'midpoint' for subject ",
              i,
              "; 'midpoint' must be greater than 1.",
              sep = ""
            )
          )
        } else if (mid.position[i] >= length(v)) {
          stop(
            paste(
              "Error in 'midpoint' for subject ",
              i,
              "; 'midpoint' must be less than the number of observations.",
              sep = ""
            )
          )
        }
      }
    } else{
      mid.position <- NULL
    }
    
    #################################################################################u########
    ##         Initialize main output data frame and compute the requested measures         ##
    ################################################################################u#########
    
    output <-
      data.frame(matrix(ncol = 1 + length(measures), nrow = nrow(data)))
    colnames(output) <- c("ID", paste("m", measures, sep = ""))
    output$ID <- IDvector
    
    data <- as.matrix(data)
    time <- as.matrix(time)
    
    if (1 %in% measures) {
      for (i in seq_len(nrow(data))) {
        output$m1[i] <-
          max(data[i,], na.rm = TRUE) - min(data[i,], na.rm = TRUE)
      }
    }
    
    if (sum(c(2, 3, 4, 8, 15, 21) %in% measures) > 0) {
      m2 <- c()
      for (i in seq_len(nrow(data))) {
        y <- data[i, complete.cases(data[i,])]
        x <- time[i, complete.cases(time[i,])]
        m2[i] <- FctMean(x, y)
      }
      if (2 %in% measures) {
        output$m2 <- m2
      }
    }
    
    if (sum(c(3, 4) %in% measures) > 0) {
      m3 <- c()
      for (i in seq_len(nrow(data))) {
        y <- data[i, complete.cases(data[i,])]
        x <- time[i, complete.cases(time[i,])]
        m3[i] <- sqrt(FctMean(x, (y - m2[i]) ^ 2))
      }
      if (3 %in% measures) {
        output$m3 <- m3
      }
    }
    
    if (4 %in% measures) {
      for (i in seq_len(nrow(data))) {
        output$m4[i] <- output$m3[i] / output$m2[i]
      }
    }
    
    if (sum(c(5, 6, 7, 8, 25, 26) %in% measures) > 0) {
      m5 <- c()
      for (i in seq_len(nrow(data))) {
        m5[i] <- Last(data[i,]) - First(data[i,])
      }
      if (5 %in% measures) {
        output$m5 <- m5
      }
    }
    
    if (6 %in% measures) {
      for (i in seq_len(nrow(data))) {
        output$m6[i] <- m5[i] / (Last(time[i,]) - First(time[i,]))
      }
    }
    
    if (7 %in% measures) {
      for (i in seq_len(nrow(data))) {
        output$m7[i] <- m5[i] / First(data[i,])
      }
    }
    
    if (8 %in% measures) {
      for (i in seq_len(nrow(data))) {
        output$m8[i] <- m5[i] / m2[i]
      }
    }
    
    if (sum(c(9, 16, 17) %in% measures) > 0) {
      m9 <- c()
      for (i in seq_len(nrow(data))) {
        b <- coefficients(lm(data[i,] ~ time[i,]))
        m9[i] <- b[2]
      }
      if (9 %in% measures) {
        output$m9 <- m9
      }
    }
    
    # Turn off warnings coming from summary(model) when the residuals are 0, which can happen when the trajectory is a straight line.
    defaultW <- getOption("warn")
    options(warn = -1)
    
    if (10 %in% measures) {
      for (i in seq_len(nrow(data))) {
        model <- lm(data[i,] ~ time[i,])
        output$m10[i] <- summary(model)$r.squared
      }
    }
    
    options(warn = defaultW)
    
    if (11 %in% measures) {
      for (i in seq_len(nrow(data))) {
        y <- data[i, complete.cases(data[i,])]
        x <- time[i, complete.cases(time[i,])]
        D <- Der(x, y)
        output$m11[i] <- max(D)
      }
    }
    
    if (sum(c(12, 17) %in% measures) > 0) {
      m12 <- c()
      for (i in seq_len(nrow(data))) {
        y <- data[i, complete.cases(data[i,])]
        x <- time[i, complete.cases(time[i,])]
        D <- Der(x, y)
        m12[i] <- sqrt(FctMean(x, (D - FctMean(x, D)) ^ 2))
      }
      if (12 %in% measures) {
        output$m12 <- m12
      }
    }
    
    if (sum(c(13, 22, 23) %in% measures) > 0) {
      m13 <- c()
      for (i in seq_len(nrow(data))) {
        y <- data[i, complete.cases(data[i,])]
        x <- time[i, complete.cases(time[i,])]
        D <- Der(x, y)
        m13[i] <- FctMean(x, abs(D))
      }
      if (13 %in% measures) {
        output$m13 <- m13
      }
    }
    
    if (sum(c(14, 15, 16) %in% measures) > 0) {
      m14 <- c()
      for (i in seq_len(nrow(data))) {
        y <- data[i, complete.cases(data[i,])]
        x <- time[i, complete.cases(time[i,])]
        D <- Der(x, y)
        m14[i] <- max(abs(D))
      }
      if (14 %in% measures) {
        output$m14 <- m14
      }
    }
    
    if (15 %in% measures) {
      for (i in seq_len(nrow(data))) {
        output$m15[i] <- m14[i] / m2[i]
      }
    }
    
    if (16 %in% measures) {
      for (i in seq_len(nrow(data))) {
        output$m16[i] <- m14[i] / m9[i]
      }
    }
    
    if (17 %in% measures) {
      for (i in seq_len(nrow(data))) {
        output$m17[i] <- m12[i] / m9[i]
      }
    }
    
    if (18 %in% measures) {
      for (i in seq_len(nrow(data))) {
        y <- data[i, complete.cases(data[i,])]
        x.i <- time[i, complete.cases(time[i,])]
        D <- Der(x.i, y)
        A <- Der(x.i, D)
        output$m18[i] <- FctMean(x.i, A)
      }
    }
    
    if (sum(c(19, 23) %in% measures) > 0) {
      m19 <- c()
      for (i in seq_len(nrow(data))) {
        y <- data[i, complete.cases(data[i,])]
        x.i <- time[i, complete.cases(time[i,])]
        D <- Der(x.i, y)
        A <- Der(x.i, D)
        m19[i] <- FctMean(x.i, abs(A))
      }
      if (19 %in% measures) {
        output$m19 <- m19
      }
    }
    
    if (sum(c(20, 21, 22) %in% measures) > 0) {
      m20 <- c()
      for (i in seq_len(nrow(data))) {
        y <- data[i, complete.cases(data[i,])]
        x.i <- time[i, complete.cases(time[i,])]
        D <- Der(x.i, y)
        A <- Der(x.i, D)
        m20[i] <- max(abs(A))
      }
      if (20 %in% measures) {
        output$m20 <- m20
      }
    }
    
    if (21 %in% measures) {
      for (i in seq_len(nrow(data))) {
        output$m21[i] <- m20[i] / m2[i]
      }
    }
    
    if (22 %in% measures) {
      for (i in seq_len(nrow(data))) {
        output$m22[i] <- m20[i] / m13[i]
      }
    }
    
    if (23 %in% measures) {
      for (i in seq_len(nrow(data))) {
        output$m23[i] <- m19[i] / m13[i]
      }
    }
    
    if (TRUE %in% (c(24:26) %in% measures)) {
      EC <- c() #  early change
      for (i in seq_len(nrow(data))) {
        EC[i] <-
          Last(data[i, 1:mid.position[i]]) - First(data[i, 1:mid.position[i]])
      }
      
      LC <- c() #  later change
      for (i in seq_len(nrow(data))) {
        LC[i] <-
          Last(data[i, mid.position[i]:ncol(data)]) - First(data[i, mid.position[i]:ncol(data)])
      }
    }
    
    if (24 %in% measures) {
      for (i in seq_len(nrow(data))) {
        output$m24[i] <- EC[i] / LC[i]
      }
    }
    
    if (25 %in% measures) {
      for (i in seq_len(nrow(data))) {
        output$m25[i] <- EC[i] / m5[i]
      }
    }
    
    if (26 %in% measures) {
      for (i in seq_len(nrow(data))) {
        output$m26[i] <- LC[i] / m5[i]
      }
    }
    
    output[is.na(output)] <-
      1 #  NAs correspond to quotient measures of the form 0/0
    
    
    # Remove the measures that are constant because (1) these are not useful for
    # discriminating between the trajectories (2) they cause problem with the
    # CheckCorrelation function later because their variance is 0 so division by 0
    # occurs when computing correlation.
    flag <- c()
    for (j in seq_len(ncol(output))) {
      col <- output[, j]
      if (max(col) == min(col)) {
        flag <- c(flag, j)
      }
    }
    if (length(flag) > 0) {
      output <- output[,-flag]
      warning(paste(
        "Measures ",
        paste(colnames(output)[3:6], collapse = ", "),
        " have been discarded due to being constant.",
        sep = ""
      ))
    }
    if (ncol(output) == 1) {
      stop("All the measures are constant.")
    }
    
    ######################################
    ##         Cap the outliers         ##
    ######################################
    
    outliers <-
      data.frame(matrix(nrow = nrow(output), ncol = ncol(output)))
    colnames(outliers) <- colnames(output)
    outliers[, 1] <- output$ID
    
    for (j in 2:ncol(output)) {
      y <- output[, j]
      y.TRUE <- y
      n <- length(y)
      which.inf <- which(is.infinite(y.TRUE))
      
      if (length(which.inf) > 0) {
        y.TRUE <- y.TRUE[-which.inf]
      }
      
      if (length(y.TRUE) > 2) {
        top <-
          rev(order(abs(y.TRUE - median(y.TRUE))))[1:ceiling(n * 0.01)] #  if n < 100, remove 1 element, so this is never empty
        y.TRUE <- y.TRUE[-top]
      }
      
      mu <- mean(y.TRUE)
      sigma <- sd(y.TRUE)
      
      k_Cheb <-
        sqrt(100 / 0.3) #  The classical Chebychev bound. Approximately 18.26.
      k <- seq(from = 0.1, to = 18.26, by = 0.1)
      M <- c()
      for (i in seq(length(k))) {
        max.left <-
          max(density(
            y.TRUE,
            from = (mu - 18.3 * sigma),
            to = (mu + 18.3 * sigma)
          )$y[density(y.TRUE,
                      from = (mu - 18.3 * sigma),
                      to = (mu + 18.3 * sigma))$x > mu + k[i] * sigma])
        max.right <-
          max(density(
            y.TRUE,
            from = (mu - 18.3 * sigma),
            to = (mu + 18.3 * sigma)
          )$y[density(y.TRUE,
                      from = (mu - 18.3 * sigma),
                      to = (mu + 18.3 * sigma))$x < mu - k[i] * sigma])
        M[i] <- max(max.left, max.right)
      }
      
      p <- 2 * pi * k ^ 2 * (exp(1) - 2 * pi / 3)
      q <-
        2 * (2 * pi * k) ^ 3 / 27 - (2 * pi) ^ 2 * k ^ 3 * exp(1) / 3 - 2 * pi * exp(1) /
        (sigma * M)
      
      root <-
        CubeRoot(-q / 2 + sqrt(q ^ 2 / 4 + p ^ 3 / 27)) + CubeRoot(-q / 2 - sqrt(q ^
                                                                                   2 / 4 + p ^ 3 / 27)) + 2 * pi * k / 3
      
      w <- which(root * sigma * M < 0.3 / 100)
      if (length(w) > 0) {
        k.opt <- min(k_Cheb, k[w[1]])
      } else{
        k.opt <- k_Cheb
      }
      
      cap <- which(abs(y - mu) > k.opt * sigma)
      if (length(cap) > 0) {
        outliers[cap, j] <- signif(y[cap], 3)
        
        y[cap] <- mu + sign(y[cap]) * k.opt * sigma
        output[, j] <- y
      }
    }
    
    row.rm <- which(rowSums(!is.na(outliers[,-c(1)])) == 0)
    outliers <- outliers[-row.rm,]
    col.rm <- which(colSums(!is.na(outliers)) == 0)
    outliers <- outliers[,-col.rm]
    
    ID <- IDvector
    
    trajMeasures <-
      structure(
        list(
          measures = output,
          outliers = outliers,
          mid = mid.position,
          data = cbind(ID, data),
          time = cbind(ID, time)
        ),
        class = "trajMeasures"
      )
    return(trajMeasures)
    
  }

#'@rdname Step1Measures
#'
#'@export
print.trajMeasures <- function(x, ...) {
  cat("Measures and time at midpoint (if applicable):\n")
  
  if (is.null(x$mid)) {
    print(x$measures)
  } else{
    measure.plus <- cbind(x$measures$ID, x$mid)
    measure.plus <- cbind(measure.plus, x$measures[,-1])
    colnames(measure.plus)[1:2] <- c("ID", "mid time")
    print(measure.plus, row.names = FALSE)
  }
}

#'@rdname Step1Measures
#'
#'@export
summary.trajMeasures <- function(object, ...) {
  cat("Description of the measures:\n\n")
  cat("m1: Range\n")
  cat("m2: Mean of the trajectory\n")
  cat("m3: Standard deviation (SD)\n")
  cat("m4: Coefficient of variation (m3/m2)\n")
  cat("m5: Overall change (initial value - final value)\n")
  cat("m6: Mean change per unit time\n")
  cat("m7: Overall change relative to initial value\n")
  cat("m8: Overall change relative to mean (m5/m2)\n")
  cat("m9: Slope of the linear model\n")
  cat("m10: Proportion of variance explained by the linear model (R squared)\n")
  cat("m11: Maximum value of the speed\n")
  cat("m12: Functional SD of the speed\n")
  cat("m13: Mean absolute speed\n")
  cat("m14: Maximum absolute speed\n")
  cat("m15: Maximum absolute speed relative to the mean (m14/m2)\n")
  cat("m16: Maximum absolute speed relative to the slope (m14/m9)\n")
  cat("m17: Functional SD of the speed relative to the slope (m12/m9)\n")
  cat("m18: Mean acceleration\n")
  cat("m19: Mean absolute acceleration\n")
  cat("m20: Maximum of the absolute acceleration\n")
  cat("m21: Maximum of the absolute acceleration relative to the mean (m20/m2)\n")
  cat(
    "m22: Maximum of the absolute acceleration relative to the mean absolute speed (m20/m13)\n"
  )
  cat("m23: Mean absolute acceleration relative to the mean absolute speed (m19/m13)\n")
  cat("m24: Early change relative to later change\n")
  cat("m25: Early change relative to overall change\n")
  cat("m26: Later change relative to overall change\n\n")
  
  
  
  cat("Summary of measures:\n")
  
  Q1 <- function(x) {
    return(quantile(x , probs = c(.25)))
  }
  
  Q2 <- function(x) {
    return(quantile(x , probs = c(.5)))
  }
  
  Q3 <- function(x) {
    return(quantile(x , probs = c(.75)))
  }
  
  measures.summary <-
    data.frame(matrix(nrow = 6, ncol = ncol(object$measures) - 1))
  rownames(measures.summary) <-
    c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  colnames(measures.summary) <- colnames(object$measures)[-1]
  
  measures.summary[1,] <- apply(object$measures, 2, min)[-1]
  measures.summary[2,] <- apply(object$measures, 2, Q1)[-1]
  measures.summary[3,] <- apply(object$measures, 2, Q2)[-1]
  measures.summary[4,] <- apply(object$measures, 2, mean)[-1]
  measures.summary[5,] <- apply(object$measures, 2, Q3)[-1]
  measures.summary[6,] <- apply(object$measures, 2, max)[-1]
  
  print(measures.summary)
  cat("\n")
  
  cat("Outliers before capping:\n")
  outliers <- object$outliers
  outliers.pre <- outliers
  outliers.pre[is.na(outliers.pre)] <- ""
  print(outliers.pre, row.names = FALSE)
  
  cat("Outliers after capping:\n")
  if (nrow(outliers) == 0) {
    print(outliers, row.names = FALSE)
  } else{
    outliers.post <- outliers
    for (j in seq_len(nrow(outliers))) {
      for (k in 2:ncol(outliers)) {
        if (is.na(outliers[j, k]) == FALSE) {
          outliers.post[j, k] <-
            signif(object$measures[object$measures$ID == outliers$ID[j], colnames(outliers)[k]], 3)
        }
      }
    }
    outliers.post[is.na(outliers.post)] <- ""
    print(outliers.post, row.names = FALSE)
  }
}
