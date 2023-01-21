 
#' @title Compute 24 Measures Describing the Features of the Trajectories
#' @description Compute 24 measures for each of the Trajectories. See details for the list of measures.
#' @param Data A \emph{n} by \emph{m} matrix or data frame containing the values of each individual trajectory. Each row corresponds to one of the \emph{n} trajectories, while the \emph{m} columns correspond to the ordered values of a given trajectory. See details.
#' @param ID Logical. Set to \code{TRUE} if the first column of \code{Data} corresponds to an \code{ID} variable. Defaults to \code{FALSE}.
#' @param verbose Logical. Set to \code{TRUE} to print information on screen. Defaults to \code{TRUE}.
#' @return trajMeasures, object containing the data used for the calculations and the 24 measures.
#' @details 
#' There must be a minimum of 4 observations for each trajectory or the trajectory will be omitted from the analysis. The
#'trajectories do not need to have the same number of observations, nor the same values of \code{Time}.
#'
#'
#'When \code{ID} is set to \code{FALSE}, a generic \code{ID} variable is created and appended as the first 
#'colunm of both the \code{Data} and \code{Time} data.frames.

#'The 24 measures are:
  
#'1. Range\cr
#'2. Mean-over-time*\cr
#'3. Standard deviation (SD)\cr
#'4. Coefficient of variation (CV)\cr
#'5. Change\cr
#'6. Mean change per unit time\cr
#'7. Change relative to the first score\cr
#'8. Change relative to the mean over time\cr
#'9. Slope of the linear model*\cr
#'10. R^2: Proportion of variance explained by the linear model\cr
#'11. Maximum of the first differences\cr
#'12. SD of the first differences\cr
#'13. SD of the first differences per time unit\cr
#'14. Mean of the absolute first differences*\cr
#'15. Maximum of the absolute first differences\cr
#'16. Ratio of the maximum absolute difference to the mean-over-time\cr
#'17. Ratio of the maximum absolute first difference to the slope\cr
#'18. Ratio of the SD of the first differences to the slope\cr
#'19. Mean of the second differences\cr
#'20. Mean of the absolute second differences\cr
#'21. Maximum of the absolute second differences\cr
#'22. Ration of the maximum absolute second difference to the mean-over-time\cr
#'23. Ratio of the maximum absolute second difference to mean absolute first difference\cr
#'24. Ratio of the mean absolute second difference to the mean absolute first difference\cr
#'
#'* If a measure is equal to zero, it will be set to the smallest, non-zero value of the same measure across the sample during further calculations.
#'If Y_1, the first observation of the trajectory of an individual, is equal to zero, it will aslo be replaced.
#'
#'For the exact equations of the measures, please go to "User guides, package vignettes and other documentation"
#'section of the "traj" package.
#'
#'
#'@author Marie-Pierre Sylvestre, Dan Vatnik
#'
#'marie-pierre.sylvestre@umontreal.ca
#'
#'
#'@examples 
#'\dontrun{
#'# Setup data
#'data = example.data$data
#'
#'# Run step1measures
#'s1 = step1measures(data, ID=TRUE)
#'
#'# Display measures
#'head(s1$measurments)
#'
#'# Plot mean trajectory of all individuals
#'plot(s1$measurments$ID, s1$measurments$m5)
#'
#'# The next step would be to run "step2factors"
#'}
#' 
#' @rdname step1measures
#' @export 
step1measures <- function (Data, ID = FALSE, verbose = TRUE) {
  data = Data
  # input.data = data
  m = ncol(data) # ICI
  n = nrow(data)  # ICI
  time = data.frame(t(matrix(rep(1:m, n), nrow = m))) #ICI
  names(time) = paste0("time.", 1:m)
  # input.time = time
  
  sample.size = dim(data)[1]
  if (ID == T) {
    data = Data
    IDvector = data[, 1]
    data = data[, -1]
    p = ncol(data)  # ICI
    time = data.frame(t(matrix(rep(1:p, n), nrow = p)))
    names(time) = paste0("time.", 1:p)
    # input.time = time
    
  }
  
  max.num.obs = dim(data)[2]
  sample.size = dim(data)[1]
  clean.data = matrix(ncol = max.num.obs, nrow = sample.size)
  clean.time = matrix(ncol = max.num.obs, nrow = sample.size)
  num.obs = rep(999, sample.size)
  less.than.4.obs = NULL
  for (i_sample in 1:sample.size) {
    real.obs.pos = which(!is.na(data[i_sample, ]))
    num.obs[i_sample] = length(real.obs.pos)
    clean.data[i_sample, ] = as.vector(c(unlist(data[i_sample,
                                                     real.obs.pos]), rep(NA, max.num.obs - num.obs[i_sample])))
    clean.time[i_sample, ] = as.vector(c(unlist(time[i_sample,
                                                     real.obs.pos]), rep(NA, max.num.obs - num.obs[i_sample])))
    if (length(real.obs.pos) < 4)
      less.than.4.obs = c(less.than.4.obs, i_sample)
    clean.data.pos = which(!is.na(clean.data[i_sample, ]))
    if (any(is.na(clean.time[i_sample, clean.data.pos])))
      stop(paste("There must be a time associated to every observation. Line: ",
                 i_sample, sep = ""))
  }
  if (!is.null(less.than.4.obs)) {
    clean.data = clean.data[-less.than.4.obs, ]
    clean.time = clean.time[-less.than.4.obs, ]
  }
  sample.size = nrow(clean.data)
  if (ID) {
    if (!is.null(less.than.4.obs))
      IDvector = IDvector[-less.than.4.obs]
  }
  else IDvector = seq(1:sample.size)
  data = clean.data
  time = clean.time
  output = data.frame(matrix(ncol = 25, nrow = sample.size))
  names = c("ID", "m1", "m2", "m3",
            "m4", "m5", "m6", "m7", "m8",
            "m9", "m10", "m11", "m12", "m13",
            "m14", "m15", "m16", "m17", "m18",
            "m19", "m20", "m21", "m22", "m23",
            "m24")
  colnames(output) = names
  output$ID = IDvector
  for (i in 1:sample.size) {
    output$m1[i] = max(data[i, ], na.rm = TRUE) - min(data[i,
    ], na.rm = TRUE)
  }
  for (i in 1:sample.size) {
    output$m2[i] = mean(data[i, ], na.rm = TRUE)
  }
  for (i in 1:sample.size) {
    output$m3[i] = sqrt(var(data[i, ], na.rm = TRUE))
  }
  for (i in 1:sample.size) {
    output$m4[i] = 100 * output$m3[i]/output$m2[i]
  }
  for (i in 1:sample.size) {
    output$m5[i] = last(data[i, ], na.rm = TRUE) - first(data[i,
    ], na.rm = TRUE)
  }
  for (i in 1:sample.size) {
    output$m6[i] = (last(data[i, ], na.rm = TRUE) - first(data[i,
    ], na.rm = TRUE))/(last(time[i, ], na.rm = TRUE) -
                         first(time[i, ], na.rm = TRUE) + 1)
  }
  for (i in 1:sample.size) {
    output$m7[i] = (last(data[i, ], na.rm = TRUE) - first(data[i,
    ], na.rm = TRUE))/first(data[i, ], na.rm = TRUE)
  }
  for (i in 1:sample.size) {
    output$m8[i] = (last(data[i, ], na.rm = TRUE) - first(data[i,
    ], na.rm = TRUE))/output$m2[i]
  }
  for (i in 1:sample.size) {
    b = coefficients(lm(data[i, ] ~ time[i, ]))
    output$m9[i] = b[2]
  }
  for (i in 1:sample.size) {
    model = lm(data[i, ] ~ time[i, ])
    r = resid(model)
    RSS = r %*% r
    Y = subset(data[i, ], is.na(data[i, ]) == FALSE)
    m = length(Y)
    SYY = Y %*% Y - (sum(Y)^2)/m
    SSREG = SYY[1] - RSS[1]
    output$m10[i] = SSREG/SYY
  }
  FD = matrix(nrow = sample.size, ncol = max.num.obs - 1)
  for (i in 1:sample.size) {
    for (j in 1:(max.num.obs - 1)) {
      FD[i, j] = data[i, (j + 1)] - data[i, j]
    }
  }
  for (i in 1:sample.size) {
    output$m11[i] = max(FD[i, ], na.rm = TRUE)
  }
  for (i in 1:sample.size) {
    output$m12[i] = sqrt(var(FD[i, ], na.rm = TRUE))
  }
  FDunit = matrix(nrow = sample.size, ncol = max.num.obs -
                    1)
  for (i in 1:sample.size) {
    for (j in 1:(max.num.obs - 1)) {
      FDunit[i, j] = FD[i, j]/(time[i, j + 1] - time[i,
                                                     j])
    }
  }
  for (i in 1:sample.size) {
    output$m13[i] = sqrt(var(FDunit[i, ], na.rm = TRUE))
  }
  for (i in 1:sample.size) {
    output$m14[i] = mean(abs(FD[i, ]), na.rm = TRUE)
  }
  for (i in 1:sample.size) {
    output$m15[i] = max(abs(FD[i, ]), na.rm = TRUE)
  }
  for (i in 1:sample.size) {
    output$m16[i] = output$m15[i]/output$m2[i]
  }
  for (i in 1:sample.size) {
    output$m17[i] = output$m15[i]/output$m9[i]
  }
  for (i in 1:sample.size) {
    output$m18[i] = output$m12[i]/output$m9[i]
  }
  SD = matrix(nrow = sample.size, ncol = max.num.obs - 2)
  for (i in 1:sample.size) {
    for (j in 1:(max.num.obs - 2)) {
      SD[i, j] = FD[i, (j + 1)] - FD[i, j]
    }
  }
  for (i in 1:sample.size) {
    output$m19[i] = mean((SD[i, ]), na.rm = TRUE)
  }
  for (i in 1:sample.size) {
    output$m20[i] = mean(abs(SD[i, ]), na.rm = TRUE)
  }
  for (i in 1:sample.size) {
    output$m21[i] = max(abs(SD[i, ]), na.rm = TRUE)
  }
  for (i in 1:sample.size) {
    output$m22[i] = output$m21[i]/output$m2[i]
  }
  for (i in 1:sample.size) {
    output$m23[i] = output$m21[i]/output$m14[i]
  }
  for (i in 1:sample.size) {
    output$m24[i] = mean(abs(SD[i, ]), na.rm = TRUE)/output$m14[i]
  }
  temp.data = output$m2[(output$m2 != 0)]
  abs.temp.data = abs(temp.data)
  if (length(temp.data) == 0)
    mean.0 = 1e-04
  else mean.0 = temp.data[which.min(abs.temp.data)]/100
  m4.na.pos = which(is.na(output$m4) | is.infinite(output$m4))
  if (length(m4.na.pos) != 0)
    output$m4[m4.na.pos] = 100 * output$m3[m4.na.pos]/mean.0
  m8.na.pos = which(is.na(output$m8) | is.infinite(output$m8))
  if (length(m8.na.pos) != 0) {
    if (length(m8.na.pos) > 1)
      output$m8[m8.na.pos] = (apply(data[m8.na.pos, ],
                                    1, last, na.rm = TRUE) - apply(data[m8.na.pos,
                                    ], 1, first, na.rm = TRUE))/mean.0
    else output$m8[m8.na.pos] = (last(data[m8.na.pos, ],
                                      na.rm = TRUE) - first(data[m8.na.pos, ], na.rm = TRUE))/mean.0
  }
  m16.na.pos = which(is.na(output$m16) | is.infinite(output$m16))
  if (length(m16.na.pos) != 0)
    output$m16[m16.na.pos] = output$m15[m16.na.pos]/mean.0
  m22.na.pos = which(is.na(output$m22) | is.infinite(output$m22))
  if (length(m22.na.pos) != 0)
    output$m22[m22.na.pos] = output$m21[m22.na.pos]/mean.0

  ########## HERE
  temp.data = data[(data[, 1] != 0), 1]# add a 1 here
  abs.temp.data = abs(temp.data)
#  if (nrow(temp.data) == 0)
#    y1.0 = 1e-04
 #else
  y1.0 = temp.data[which.min(abs.temp.data)]/100 # is zero but should not

  m7.na.pos = which(is.na(output$m7) | is.infinite(output$m7))
  if (length(m7.na.pos) != 0) {
    if (length(m7.na.pos) > 1)

      output$m7[m7.na.pos] = (apply(data[m7.na.pos, ], 1, last, na.rm = TRUE) - apply(data[m7.na.pos,], 1, first, na.rm = TRUE))/y1.0
    else output$m7[m7.na.pos] = (last(data[m7.na.pos, ], na.rm = TRUE) - first(data[m7.na.pos, ], na.rm = TRUE))/y1.0
  }
  form10 = vector(length = sample.size)
  for (i_test in 1:sample.size) {
    model = lm(data[i_test, ] ~ time[i_test, ])
    r = resid(model)
    RSS = r %*% r
    Y = subset(data[i_test, ], is.na(data[i_test, ]) == FALSE)
    m = length(Y)
    form10[i_test] = Y %*% Y - (sum(Y)^2)/m
    if (form10[i_test] == 0)
      form10[i_test] = NA
  }
  if (!is.na(min(form10)))
    syy.0 = min(form10)
  else syy.0 = 1e-04
  m10.na.pos = which(is.na(output$m10) | is.infinite(output$m10))
  if (length(m10.na.pos) != 0) {
    for (i_na in m10.na.pos) {
      model = lm(data[i_na, ] ~ time[i_na, ])
      r = resid(model)
      RSS = r %*% r
      Y = subset(data[i_na, ], is.na(data[i_na, ]) == FALSE)
      m = length(Y)
      SSREG = SYY[1] - RSS[1]
      output$m10[i_na] = SSREG/syy.0
      if(output$m10[i_na] > 1) output$m10[i_na] = 1
      if(output$m10[i_na] < 0) output$m10[i_na] = 0
    }
  }
  temp.data = output$m9[(output$m9 != 0)]
  abs.temp.data = abs(temp.data)
  if (length(temp.data) == 0)
    slope.0 = 1e-04
  else slope.0 = temp.data[which.min(abs.temp.data)]/100
  m17.na.pos = which(is.na(output$m17) | is.infinite(output$m17))
  if (length(m17.na.pos) != 0)
    output$m17[m17.na.pos] = output$m15[m17.na.pos]/slope.0
  m18.na.pos = which(is.na(output$m18) | is.infinite(output$m18))
  if (length(m18.na.pos) != 0)
    output$m18[m18.na.pos] = output$m12[m18.na.pos]/slope.0
  temp.data = output$m14[(output$m14 != 0)]
  abs.temp.data = abs(temp.data)
  if (length(temp.data) == 0)
    mean.abs.0 = 1e-04
  else mean.abs.0 = temp.data[which.min(abs.temp.data)]/100
  m23.na.pos = which(is.na(output$m23) | is.infinite(output$m23))
  if (length(m23.na.pos) != 0)
    output$m23[m23.na.pos] = output$m21[m23.na.pos]/mean.abs.0
  m24.na.pos = which(is.na(output$m24) | is.infinite(output$m24))
  if (length(m24.na.pos) != 0) {
    abs.na.data = abs(SD[m24.na.pos, ])
    if (length(m24.na.pos) > 1)
      output$m24[m24.na.pos] = apply(abs.na.data, 1, mean,
                                     na.rm = TRUE)/mean.abs.0
    else output$m24[m24.na.pos] = mean(abs.na.data, na.rm = TRUE)/mean.abs.0
  }
  check.correlation(output[, -1], verbose)
  trajMeasures = structure(list(measurments = output, data = cbind(IDvector,
                                                                   clean.data), time = cbind(IDvector, clean.time)), class = "trajMeasures")
  return(trajMeasures)
}
