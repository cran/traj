\name{plotBoxplotTraj}
\alias{plotBoxplotTraj}

\title{
Plot Boxplot of \code{traj} Object
}
\description{
Produce a boxplot of the values of the trajectories from each cluster at every time points.
}
\usage{
plotBoxplotTraj(x, clust.num = NULL,  ...)
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
The function plots a boxplot of values of the trajectories in a cluster at each time point.
}

\author{
Dan Vatnik, Marie-Pierre Sylvestre\cr
\email{dan.vatnik@gmail.com}
}
\seealso{
\code{\link[graphics]{boxplot}}
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

# Plot boxplots
plotBoxplotTraj(s3.4clusters)
}

