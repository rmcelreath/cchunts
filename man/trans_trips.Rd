\name{trans_trips}
\alias{trans_trips}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert data into trip format}
\description{
  Takes the results of \code{prep_data} and converts them into a list of trips. 
}
\usage{
trans_trips( d )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{d}{A result from \code{prep_data}}
}
\details{
  This function takes the results of \code{prep_data} and converts them into a list of trips and their associated data. This is the format used internally by the Stan model, but the Stan model does this conversion on its own. This function is provided for debugging and transparency.
}
\value{
}
\references{}
\seealso{ \link{\code{make_joint}}, \link{\code{prep_data}} }
\examples{
dat <- make_joint()
dat_list <- prep_data( dat , dogs_miss=-1 , guns_miss=-1 )
trip_form <- trans_trips( dat_list )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

