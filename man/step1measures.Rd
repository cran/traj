\name{step1measures}
\alias{step1measures}
\title{Compute 24 Measures Describing the Features of the Trajectories}
\usage{
  step1measures(Data, Time, ID = FALSE, verbose = TRUE)
}
\arguments{
  \item{Data}{A \emph{n} by \emph{m} matrix or data frame containing the values of each individual trajectory. Each row corresponds to one of the \emph{n} trajectories, while the \emph{m} columns 
  correspond to the ordered values of a given trajectory. See details}

  \item{Time}{A \emph{n} by \emph{m} matrix or data frame containing the  measurement times corresponding to the values of the \code{Data} data frame. See details.}

  \item{ID}{Logical. Set to \code{TRUE} if the first column of \code{Data} corresponds to an \code{ID} variable. Defaults to \code{FALSE}}

  \item{verbose}{Logical. Set to \code{TRUE} to print information on screen. Defaults to \code{TRUE}}
}
\value{
  \item{trajMeasures}{Object containing the data used for the calculations and the 24 measures.}
}
\description{
  Computes 24 measures for each of the trajectories. See details for the list of measures.
}
\author{
  Marie-Pierre Sylvestre, Dan Vatnik\cr
  \email{dan.vatnik@gmail.com}
}
\details{
  There must be a minimum of 4 observations for each trajectory or the trajectory will be omitted from the analysis. The
  trajectories do not need to have the same number of observations, nor the same values of \code{Time}.

  The \code{Time} data.frame or matrix must have the same dimension as the \code{Data} 
  data frame or matrix and must not contain missing values or an error will be returned. The data can have missing values, \code{Time} can not.

  When \code{ID} is set to \code{FALSE}, a generic \code{ID} variable is created and appended as the first 
  colunm of both the \code{Data} and \code{Time} data.frames.
  
  The 24 measures are:
  
 1. Range\cr
 2. Mean-over-time*\cr
 3. Standard deviation (SD)\cr
 4. Coefficient of variation (CV)\cr
 5. Change\cr
 6. Mean change per unit time\cr
 7. Change relative to the first score\cr
 8. Change relative to the mean over time\cr
 9. Slope of the linear model*\cr
 10. R^2: Proportion of variance explained by the linear model\cr
 11. Maximum of the first differences\cr
 12. SD of the first differences\cr
 13. SD of the first differences per time unit\cr
 14. Mean of the absolute first differences*\cr
 15. Maximum of the absolute first differences\cr
 16. Ratio of the maximum absolute difference to the mean-over-time\cr
 17. Ratio of the maximum absolute first difference to the slope\cr
 18. Ratio of the SD of the first differences to the slope\cr
 19. Mean of the second differences\cr
 20. Mean of the absolute second differences\cr
 21. Maximum of the absolute second differences\cr
 22. Ration of the maximum absolute second difference to the mean-over-time\cr
 23. Ratio of the maximum absolute second difference to mean absolute first difference\cr
 24. Ratio of the mean absolute second difference to the mean absolute first difference\cr
 
 * If a measure is equal to zero, it will be set to the smallest, non-zero value of the same measure across the sample during further calculations.
 If Y_1, the first observation of the trajectory of an individual, is equal to zero, it will aslo be replaced.
 
 For the exact equations of the measures, please go to "User guides, package vignettes and other documentation"
 section of the "traj" package.
}
\examples{
# Setup data and time
data = example.data$data
time = example.data$time

# Run step1measures
s1 = step1measures(data,time, ID=TRUE)

# Display measures
head(s1$measurments)

# Plot mean trajectory of all individuals
plot(s1$measurments$ID, s1$measurments$m5)

# The next step would be to run "step2factors"

}
