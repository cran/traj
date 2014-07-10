\name{plotMeanTraj}
\alias{plotMeanTraj}
\title{
Plot Mean Trajectory
}
\description{
Plot cluster-specific mean trajectory for one or all clusters provided by a \code{traj} object.
}
\usage{
plotMeanTraj(x, clust.num = NULL,  ...)
}

\arguments{
  \item{x}{
\code{traj} object.
}
  \item{clust.num}{
Integer indicating the cluster number to plot.\code{NULL} to print all clusters. Defaults to
\code{NULL}.
}
  \item{\dots}{
  Arguments to be passed to plot.
  }
}
\details{
  The function plots cluster specific mean trajectory calculated at each time point. By setting the \code{clust.num} argument to an integer corresponding to a cluster number, one can plot the mean trajectory of that cluster only. Any other plotting arguments can be added to the function.
}

\author{
Dan Vatnik, Marie-Pierre Sylvestre\cr
\email{dan.vatnik@gmail.com}
}


\seealso{
\code{\link[base]{mean}}
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

# Plot mean trajectories
plotMeanTraj(s3.4clusters)
}
