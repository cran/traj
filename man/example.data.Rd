\name{example.data}
\alias{example.data}
\docType{data}
\title{
Example Data
}
\description{
Example data used to display the features of the \code{traj} package.
}
\usage{example.data}
\format{

List of 2 data.frames:

  $ data:	130 subjects with 7 variables. X1 to X6 correspond to the 6 measurements obtained on the 130 subjects:\cr
     ID: num [1:130] \cr 
     X1: num [1:130] \cr
     X2: num [1:130] \cr
     X3: num [1:130] \cr
     X4: num [1:130] \cr
     X5: num [1:130] \cr
     X6: num [1:130] \cr
   
 $ time: 130 subjects with 7 variables. time1 to time6 correspond to the measurement times for the vairbales X1 to X6:\cr
     ID   : num [1:130] \cr
     time1: num [1:130] \cr
     time2: num [1:130] \cr
     time3: num [1:130] \cr
     time4: num [1:130] \cr
     time5: num [1:130] \cr
     time6: num [1:130] \cr
}

\examples{

# data and time
data = example.data$data
time = example.data$time

}
\keyword{datasets}
