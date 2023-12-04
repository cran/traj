#'@title Select a Subset of the Measures Using Factor Analysis
#'
#'@description This function applies the following dimension reduction algorithm
#'  to the measures computed by \code{\link[traj]{Step1Measures}}:
#' \enumerate{
#'   \item Use principal component analysis (PCA) on the measure to form factors summarizing the variability in the measures;
#'   \item Drop the factors whose variance is smaller than any one of the normalized measures;
#'   \item Performs a varimax rotation on the remaining factors;
#'   \item For each rotated factor, select the measure that has the highest correlation (aka factor loading) with it and that hasn't yet been selected;
#'   \item Drop the remaining measures.
#' }
#'
#'@param trajMeasures object of class \code{trajMeasures} as returned by
#'  \code{\link[traj]{Step1Measures}}.
#'@param num.select an optional positive integer indicating the number of
#'  factors to keep in the second stage of the algorithm. Defaults to NULL so
#'  that all factors with variance greater than any one of the normalized
#'  measures are selected.
#'@param discard an optional vector of positive integers corresponding to the
#'  measures to be dropped from the analysis. See
#'  \code{\link[traj]{Step1Measures}} for the list of measures. Defaults to
#'  NULL.
#'@param select an optional vector of positive integers corresponding to the
#'  measures to forcefully select. Defaults to NULL. If a vector is supplied,
#'  the five-steps selection algorithm described above is bypassed and the
#'  corresponding measures are selected instead.
#'
#'  Can be NULL or a numeric vector corresponding to the numerical identifier of
#'  measures present in \code{trajMeasures}. If a numeric vector is supplied,
#'  then four-steps selection algorithm described above is bypassed and the
#'  corresponding measures are selected instead.
#'@param x object of class trajSelection.
#'@param object object of class trajSelection.
#'@param ... further arguments passed to or from other methods.
#'
#'@return An object of class \code{trajSelection}; a list containing the values
#'  of the selected measures, the output of the principal component analysis as
#'  well as a curated form of the arguments.
#'
#'@importFrom psych principal
#'@importFrom stats cor
#'
#'@details In the presence of highly correlated measures (Pearson correlation >
#'  0.98), the function selects the highest-ranking measure on the list (see
#'  \code{\link[traj]{Step1Measures}}) and discards the others. Because the
#'  K-means algorithm is sensitive to outliers, measures which are quotients
#'  (i.e. 4, 7, 8, 15-17, 21-26) are prevented from taking extremely large or
#'  infinite values (caused by division by 0). Nishiyama's improved Chebychev
#'  bound is used to determine extreme values for each measure, corresponding to
#'  a 0.3% probability threshold. Extreme values beyond the threshold are capped
#'  to the 0.3% probability threshold. Measures corresponding to quotients which
#'  would be of the form 0/0 are set to 1. PCA is applied on the remaining
#'  measures using the \code{\link[psych]{principal}} function from the
#'  \code{psych} package.
#'
#'@references Leffondre K, Abrahamowicz M, Regeasse A, Hawker GA, Badley EM,
#'  McCusker J, Belzile E. Statistical measures were proposed for identifying
#'  longitudinal patterns of change in quantitative health indicators. J Clin
#'  Epidemiol. 2004 Oct;57(10):1049-62. doi: 10.1016/j.jclinepi.2004.02.012.
#'  PMID: 15528056.
#'
#' @examples
#' \dontrun{
#'m = Step1Measures(trajdata, ID = TRUE)
#'s = Step2Selection(m)
#'
#'s$RC$loadings
#'
#'s2 = Step2Selection(m, select = c(10, 12, 8, 4))
#'}
#'
#'
#'@seealso \code{\link[psych]{principal}} \code{\link[traj]{Step1Measures}}
#'
#'@rdname Step2Selection
#'
#'@export
Step2Selection <-
  function (trajMeasures,
            num.select = NULL,
            discard = NULL,
            select = NULL) {
    input <- list(num.select, discard, select)
    names(input) <- c("num.select", "discard", "select")
    
    if ((!is.null(select)) & (!is.numeric(select))) {
      stop("Argument 'select' must be either NULL or a numerical vector.")
    }
    
    data <- data.frame(trajMeasures$measures)
    ID <- data[, 1]
    data <- data.frame(data[,-1])
    
    if (!is.null(select)) {
      m.select <- paste("m", select, sep = "")
      if (FALSE %in% (m.select %in% colnames(data))) {
        stop("Select 'select' from the measures included in step1measure.")
      } else {
        output <- cbind(ID, data[, m.select])
        colnames(output) <- c("ID", paste("m", select, sep = ""))
      }
      
      trajSelection <-
        structure(
          list(
            selection = output,
            PC = NULL,
            RC = NULL,
            colinear.variables = NULL,
            measures = trajMeasures$measures,
            data = trajMeasures$data,
            time = trajMeasures$time,
            input = input
          ),
          class = "trajSelection"
        )
      
    } else {
      if ((!is.null(discard)) & (!is.numeric(discard))) {
        stop("Argument 'discard' must be either NULL or a numerical vector.")
      }
      if (!is.null(discard)) {
        mes.to.discard <- paste("m", discard, sep = "")
        if (FALSE %in% (mes.to.discard %in% colnames(data))) {
          stop("Can't discard a measure which was not included in step1measure.")
        }
        w <- which(colnames(data) %in% mes.to.discard)
        data <- data[,-w]
      }
      
      if (!is.null(num.select)) {
        if (!is.numeric(num.select)) {
          stop("Argument 'num.select' must be a numerical vector of length 1.")
        }
        if (!is.vector(num.select)) {
          stop("Argument 'num.select' must be a numerical vector of length 1.")
        }
        if (!(length(num.select) == 1)) {
          stop("Argument 'num.select' must be a numerical vector of length 1.")
        }
        if (!is.null(discard) & (num.select > ncol(data))) {
          stop(
            "After discarding the measures specified in 'discard', the requested number 'num.select' of measures to retain exceeds the number of measures available."
          )
        }
        if (is.null(discard) & (num.select > ncol(data))) {
          stop(
            "The requested number 'num.select' of measures to retain exceeds the number of measures included in step1measure."
          )
        }
      }
      
      corr.vars <-
        CheckCorrelation(data, verbose = FALSE, is.return = TRUE)
      
      if (!is.null(corr.vars)) {
        colinear.variables <- corr.vars
        corr.vars.pos <- which(names(data) %in% corr.vars[, 1])
        data <- data[,-corr.vars.pos]
      } else{
        colinear.variables <- NULL
      }
      
      if (num.select > ncol(data) && !is.null(num.select)) {
        stop(
          "After discarding the perfectly or almost perfectly correlated measures, there are ",
          ncol(data),
          " measures left, which is less than the number 'num.select' of requested measures to retain."
        )
      }
      
      Z <- scale(x = data,
                 center = TRUE,
                 scale = TRUE)
      
      if (is.null(num.select)) {
        eigen.values <-
          psych::principal(Z, rotate = "none", nfactors = ncol(data))$values
        num.select <- max(1, length(which(eigen.values > 1)))
      }
      
      if (ncol(data) > 1) {
        PC <- psych::principal(Z, rotate = "none", nfactors = ncol(data))
        RC <-
          psych::principal(Z, rotate = "varimax", nfactors = num.select)
        principal.factors <- RC
        
        principal.variables <- c()
        
        if (is.null(select)) {
          for (j in 1:num.select) {
            if (j == 1) {
              aux <- principal.factors$loadings
            }
            if (j > 1) {
              w <-
                which(row.names(principal.factors$loadings) %in% principal.variables)
              aux <- principal.factors$loadings[-w,]
            }
            principal.variables[j] <- names(which.max(abs(aux[, j])))
          }
          output <- cbind(ID, data[, principal.variables])
        }
      } else {
        output <- trajMeasures$measures
        PC <- NULL
        RC <- NULL
      }
      
      trajSelection <-
        structure(
          list(
            selection = output,
            PC = PC,
            RC = RC,
            colinear.variables = colinear.variables,
            measures = trajMeasures$measures,
            data = trajMeasures$data,
            time = trajMeasures$time,
            input = input
          ),
          class = "trajSelection"
        )
    }
    
    return(trajSelection)
    
  }

#'@rdname Step2Selection
#'
#'@export
print.trajSelection <- function(x, ...) {
  print(
    paste(
      x$correlations[, 1],
      " was discarded because it is perfectly or almost perfectly correlated with ",
      x$correlations[, 2],
      ".",
      sep = ""
    )
  )
  cat("\n")
  
  cat(
    paste(
      "In decreasing order of variance explained, the selected measures are ",
      paste(
        colnames(x$selection)[-1],
        collapse = ', ',
        sep = ""
      ),
      ".",
      sep = ""
    )
  )
  cat("\n")
  
  if (!is.null(x$RC)) {
    print(x$RC$loadings)
  }
}

#'@rdname Step2Selection
#'
#'@export
summary.trajSelection <- function(object, ...) {
  if (!is.null(object$input$select)) {
    cat(paste(
      "The measures ",
      paste(
        colnames(object$selection)[-1],
        collapse = ', ',
        sep = ""
      ),
      " were selected.",
      sep = ""
    ))
  } else{
    if (!is.null(object$input$num.select)) {
      if (!is.null(object$colinear.variables)) {
        dropped <- unique(object$colinear.variables[, 1])
        cat(
          paste(
            "The measures ",
            paste(dropped, collapse = ", "),
            " were discarded because they were perfectly or almost perfectly correlated with another measure. Upon forming the principal components from the remaining measures, the ",
            object$input$num.select,
            " factors that held the most variance were retained. Together, they explained ",
            round(
              100 * sum(object$PC$values[1:(length(colnames(object$selection)) - 1)]) /
                length(object$PC$values),
              1
            ) ,
            "% of the total variance. A varimax rotation was performed to maximize the correlation with the original measures without affecting the proportion of explained variance. For each rotated factor, the measure that had the highest correlation (loading) with it was selected. As a result of this procedure, the selected measures are, in decreasing order of variance explained, ",
            paste(
              colnames(object$selection)[-1],
              collapse = ', ',
              sep = ""
            ),
            ". Use print() to see more detailed informations.",
            sep = ""
          )
        )
      } else{
        cat(
          paste(
            "Upon forming the principal components from the measures, the ",
            object$input$num.select,
            " factors that held the most variance were retained. Together, they explained ",
            round(
              100 * sum(object$PC$values[1:(length(colnames(object$selection)) - 1)]) /
                length(object$PC$values),
              1
            ) ,
            "% of the total variance. A varimax rotation was performed to maximize the correlation with the original measures without affecting the proportion of explained variance. For each rotated factor, the measure that had the highest correlation (loading) with it was selected. As a result of this procedure, the selected measures are, in decreasing order of variance explained, ",
            paste(
              colnames(object$selection)[-1],
              collapse = ', ',
              sep = ""
            ),
            ". Use print() to see more detailed informations.",
            sep = ""
          )
        )
      }
    } else{
      if (!is.null(object$colinear.variables)) {
        dropped <- unique(object$colinear.variables[, 1])
        cat(
          paste(
            "The measures ",
            paste(dropped, collapse = ", "),
            " were discarded because they were perfectly or almost perfectly correlated with another measure. Upon forming the principal components from the remaining measures, ",
            length(colnames(object$selection)) - 1,
            " of them had a variance greater than any one of the normalized measures. Together, they explained ",
            round(
              100 * sum(object$PC$values[1:(length(colnames(object$selection)) - 1)]) /
                length(object$PC$values),
              1
            ) ,
            "% of the total variance. A varimax rotation was performed to maximize the correlation with the original measures without affecting the proportion of explained variance. For each rotated factor, the measure that had the highest correlation (loading) with it was selected. As a result of this procedure, the selected measures are, in decreasing order of variance explained, ",
            paste(
              colnames(object$selection)[-1],
              collapse = ', ',
              sep = ""
            ),
            ". Use print() to see more detailed informations.",
            sep = ""
          )
        )
      } else{
        cat(
          paste(
            "Upon forming the principal components from the measures, ",
            length(colnames(object$selection)) - 1,
            " of them had a variance greater than any one of the normalized measures. Together, they explained ",
            round(
              100 * sum(object$PC$values[1:(length(colnames(object$selection)) - 1)]) /
                length(object$PC$values),
              1
            ) ,
            "% of the total variance. A varimax rotation was performed to maximize the correlation with the original measures without affecting the proportion of explained variance. For each rotated factor, the measure that had the highest correlation (loading) with it was selected. As a result of this procedure, the selected measures are, in decreasing order of variance explained, ",
            paste(
              colnames(object$selection)[-1],
              collapse = ', ',
              sep = ""
            ),
            ". Use print() to see more detailed informations.",
            sep = ""
          )
        )
      }
    }
  }
}