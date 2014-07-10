\name{plotMedTraj}
\alias{plotMedTraj}

\title{
Plot Median Trajectory of \code{traj} Object
}
\description{
Plot cluster-specific median trajectory for one or all clusters provided by a \code{traj} object.
}
\usage{
plotMedTraj(x, clust.num = NULL, plot.percentile = TRUE, low.percentile = 0.10, 
         high.percentile = 0.90, \dots)
}

\arguments{
  \item{x}{
  \code{traj} object.
}

\item{clust.num}{
Integer indicating the cluster number to plot.\code{NULL} to print all clusters. Defaults to
\code{NULL}.
}
\item{plot.percentile}{
Value indicating if the function should plot percentiles. Defaults to \code{TRUE}.
}
\item{low.percentile}{
Value of the lower percentile to be plotted. Must be between 0 and 1. Defaults to 0.1.
}
\item{high.percentile}{
Value of the high percentile to be plotted. Must be between 0 and 1. Defaults to 0.9.
}
\item{\dots}{
Extra parameters used in the \code{plot} function.
}
}
\details{
  The function plots cluster specific median trajectory calculated at each time point, in addition to 10th and 90th percentiles. By setting the \code{clust.num} argument to an integer corresponding to a cluster number, one can plot the median trajectory of that cluster only. Any other plotting arguments can be added to the function.
<}

\author{
Dan Vatnik, Marie-Pierre Sylvestre\cr
\email{dan.vatnik@gmail.com}
}
\seealso{
\code{\link[graphics]{plot}}
\code{\link[stats]{median}}
\code{\link[stats]{quantile}}
}
\examples{
# Setup data and time
data = example.data$data
time = example.data$time

# Run step1measures, step2factors and step3clusters with
# a predetermined number of clusters
s1 = step1measures(data,time, ID=TRUE)
s2 = step2factors(s1)
s3.4clusters = step3clusters(s2, nclust = 4)

# Plot median trajectories
plotMedTraj(s3.4clusters)
}

