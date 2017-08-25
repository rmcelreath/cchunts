\name{cch_predict}
\alias{cch_predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute posterior predictions}
\description{
  Computes the posterior predictive distribution for a supplied set of cases.
}
\usage{
cch_predict( post , data , func , verbose=FALSE , logh=TRUE , raw=FALSE , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{post}{List of posterior samples}
    \item{data}{Data list to structure predictions}
    \item{func}{Optional function to apply over cases of result. Otherwise returns samples from posterior predictive distribution.}
    \item{verbose}{When \code{TRUE}, shows debugging information.}
    \item{logh}{When \code{TRUE}, uses log trip duration (as in model)}
    \item{raw}{When \code{TRUE}, ignores production functions and returns raw skill}
    \item{...}{Arguments to pass to \code{func} via \code{apply}}
}
\details{
  This function uses arbitrary cases in \code{data} to feed into the statistical model and compute posterior predictions, using the posterior samples in \code{post}. The plotting functions use this function internally to compute skill and production curves.
}
\value{
}
\references{}
\seealso{  }
\examples{

# data range over 50 ages from birth to age 80 (standardized to 1)
pdat <- data.frame( age=seq(from=0,to=1,length.out=50) )

# requires posterior samples in post object
l <- cch_predict( post , pdat )
cch_predict( post , pdat , func=mean )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

