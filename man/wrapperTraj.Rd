\name{wrapperTraj}
\alias{wrapperTraj}

\title{
Wrapper Function to Perform Trajectory Analysis
}
\description{
Run three steps of trajectory analysis with default parameters.
}
\usage{
wrapperTraj(Data, Time, ID = FALSE)
}

\arguments{
\item{Data}{Data frame containing trajectory data.
  Each line should contain sequential observations. See details.}

  \item{Time}{Data frame containing the time of each observation in the
  \code{Data} table.See details.}

  \item{ID}{Logical. Set to \code{FALSE} if the first column of \code{Data} corresponds to an \code{ID} variable. Defaults to \code{FALSE}}

}
\details{
The function runs the full three step trajectory analysis and returns a \code{traj} object.
It will execute \code{\link[traj]{step1measures}}, \code{\link[traj]{step2factors}} and \code{\link[traj]{step3clusters}}
sequentially with their default parameters. The result of \code{step3clusters} will be returned. Details regarding 
the \code{data} and the \code{time} arguments are found in the 'Details' section of \code{\link[traj]{step1measures}}.
}
\value{
The result is a \code{traj} object. Details can be found in the 'Value' section of \code{\link[traj]{step3clusters}}.
}
\author{
Dan Vatnik, Marie-Pierre Sylvestre\cr
\email{dan.vatnik@gmail.com}
}

\examples{
# Setup data and time
data = example.data$data
time = example.data$time

# Run clustering wrapper function
wt = wrapperTraj(data, time, ID = TRUE)

# Display and plot "traj" object
wt
summary(wt)

plot(wt)
}