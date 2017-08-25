\name{cch_plot_curve}
\alias{cch_plot_curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot skill and production curves}
\description{
  Draws individual site skill and production curves, using posterior samples.
}
\usage{
cch_plot_curve( map_id , soc_id , ids=NULL , x=seq(from=0,to=1,length.out=50) , h=1 , col=cchpal[4] , col2=col.alpha("black",0.2) , alpha , pts=FALSE , sample_post=FALSE , legend=FALSE , jitter=FALSE  , sw=0.1 , ybuffer=0.2 , skill_only=FALSE , ra=80 , part=3 , main , show_peak=TRUE , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{map_id}{ID of site, using West-to-East map ID}
    \item{soc_id}{ID of site, using internal data ID}
    \item{ids}{IDs of individuals to show. NULL indicates all.}
    \item{x}{sequence of x-axis (standardized age) values}
    \item{h}{Labor input (trip duration). 1 means average.}
    \item{col}{Curve color}
    \item{col2}{Curve color where ages not observed}
    \item{alpha}{transparency}
    \item{pts}{Whether to show raw points (not adjusted for labor)}
    \item{sample_post}{When an integer value, plots that number of new individuals sampled from posterior rather than observed individuals}
    \item{legend}{When \code{TRUE}, plots legend}
    \item{jitter}{Adds jitter}
    \item{sw}{unused}
    \item{ybuffer}{Head room in plot}
    \item{skill_only}{When \code{TRUE}, raw skill curve rather than production}
    \item{ra}{Reference age for standardization. 80 by default.}
    \item{part}{Component of production to show. 1=success, 2=harvest, 3=expected returns.}
    \item{main}{Label for plot. If left missing, fills with site ID, name code, and sample sizes.}
    \item{show_peak}{When \code{TRUE}, draws a line to indicate average age of peak skill/production.}
    \item{...}{Additional arguments to pass to plotting functions}
}
\details{
  This functions uses posterior samples in a \code{post} to draw skill or production curves for a site.
}
\value{
}
\references{}
\seealso{  }
\examples{

par(mgp = c(1.5, 0.2, 0), mar = c(1.2, 0.5, 1.5, 0.25) + 0.1, tck = -0.02)
par(mfrow=c(2,2))

j <- 16 # Ache

# skill
cch_plot_curve( map_id=j , col=cchpal[4] , alpha=0.5 , pts=FALSE , sample_post=FALSE , lwd=1 , skill_only=TRUE , main="ACH - skill" )

# success rate
cch_plot_curve( map_id=j , col=cchpal[5] , alpha=0.5 , pts=TRUE , sample_post=FALSE , lwd=1 , skill_only=FALSE , part=1 , main="ACH - success" )
abline( h=0.5 , lty=2 , lwd=0.5 )

# harvests
cch_plot_curve( map_id=j , col=cchpal[3] , alpha=0.5 , pts=TRUE , sample_post=FALSE , lwd=1 , skill_only=FALSE , part=2 , main="ACH - harvest" )

# expected production
cch_plot_curve( map_id=j , col="#D4ACE8" , alpha=0.5 , pts=TRUE , sample_post=FALSE , lwd=1 , skill_only=FALSE , main="ACH - production" )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

