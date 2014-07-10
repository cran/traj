\name{plotCombTraj}
\alias{plotCombTraj}
\title{
Plot Cluster-Specific Mean or Median Trajectories provided by a \code{traj} object
}
\description{
Plot cluster-specific mean or median trajectories.
}
\usage{
plotCombTraj(x, stat.type = "mean", colored = FALSE, ...)
}

\arguments{
  \item{x}{
\code{traj} object.
}
  \item{stat.type}{
Choice between "mean" or "median". The mean or the median calculated at each time point for a cluster-specific set of trajectories will be plotted.
Defaults to "mean."
}
  \item{colored}{
Boolean in dictating if the plot should use colors. If not, the trajectory lines will be distinctively patterned.
Defaults to \code{FALSE} for patterns.
}
  \item{\dots}{
Any extra parameter used by the \code{plot} function.
}
}
\details{
The function plots the mean or the median cluster-specific trajectory, calculated at each time point. A legend is generated in the
top left corner of the plot. Other plotting parameter(s) can be added to the function with the use of \dots.
}
\author{
Dan Vatnik, Marie-Pierre Sylvestre\cr
\email{dan.vatnik@gmail.com}
}

\seealso{
\code{\link[base]{mean}}
\code{\link[stats]{median}}
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

# Plot mean combination trajectories
plotCombTraj(s3.4clusters)
}

