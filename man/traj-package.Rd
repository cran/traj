\name{traj-package}
\alias{traj-package}
\alias{traj}
\docType{package}
\title{
Trajectory Analysis
}
\description{
Implements the three step procedure proposed by Leffondree et al. (2004) to identify clusters of individual longitudinal trajectories. 

The procedure involves (1) calculating 24 measures describing the features of the trajectories, (2) using factor analysis to select a subset of the 24 measures and (3) using cluster analysis to identify clusters of trajectories, and classify each individual trajectory in one of the clusters.
}
\details{
\tabular{ll}{
Package: \tab traj\cr
Version: \tab 1.1\cr
Date: \tab 2014-07-10\cr
License: \tab GPL-2\cr
LazyData: \tab true\cr
Depends: \tab R (>= 3.0.3)\cr
Imports: \tab cluster, psych, pastecs, NbClust, graphics, grDevices, stats,
utils, GPArotation\cr
Packaged: \tab 2014-10-21 17:48:21 UTC, Dan\cr
Built: \tab R 3.0.3, , 2014-10-21 17:49:32 UTC, windows\cr
}

Index:
\preformatted{
example.data            Example Data
plotBoxplotTraj         Plot Boxplot of 'traj' Object
plotCombTraj            Plot Cluster-Specific Mean or Median
                        Trajectories provided by a 'traj' object
plotMeanTraj            Plot Mean Trajectory
plotMedTraj             Plot Median Trajectory of 'traj' Object
step1measures           Compute 24 Measures Describing the Features of
                        the Trajectories
step2factors            Performs Factor Analysis to Select a Subset of
                        the 24 Measures
step3clusters           Cluster Trajectories According to the Subset of
                        Measures Selected Previously
wrapperTraj             Wrapper Function to Perform Trajectory Analysis
}

}
\author{
Marie-Pierre Sylvestre, Dan Vatnik

Maintainer: Dan Vatnik <dan.vatnik@gmail.com>

}
\references{
\enumerate{

\item Sylvestre MP, et al. (2006). Classification of patterns of delirium severity scores over time in an elderly population. International Psychogeriatrics, 18(4), 667-680. doi:10.1017/S1041610206003334.
\item Leffondree, K. et al. (2004). Statistical measures were proposed for identifying longitudinal patterns of change in quantitative health indicators. Journal of Clinical Epidemiology, 57, 1049-1062. doi : 10.1016/j.jclinepi.2004.02.012.

}
}

\keyword{ package }
\seealso{
\code{\link[NbClust]{NbClust}}
\code{\link[stats]{kmeans}}~
}
\examples{
# Setup data and time
data = example.data$data
time = example.data$time

# Run step1measures, step2factors and step3clusters
s1 = step1measures(data,time, ID=TRUE)
s2 = step2factors(s1)
s3 = step3clusters(s2)

# Print and plot "traj object"
s3
plot(s3)
}
