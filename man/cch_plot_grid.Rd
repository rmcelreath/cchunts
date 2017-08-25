\name{cch_plot_grid}
\alias{cch_plot_grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot multiple skill and production curves in a grid}
\description{
  Draws a grid of site skill and production curves, using posterior samples.
}
\usage{
cch_plot_grid( map_id=1:40 , nrow=7 , ncol=6 , col=cchpal[4] , col2=col.alpha("black",0.5) , fskillonly=TRUE , dosample=FALSE , draw_globalmean=TRUE , draw_legend=TRUE , alpha=0.5 , show_points=FALSE , lwd=1.5 , adj_margins=TRUE , part=3 , skip=0 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{map_id}{List of site IDs, as they appear on map, as in \code{cch_plot_curve}}
    \item{nrow}{Number of rows in grid}
    \item{ncol}{Number of columns in grid}
    \item{col}{Curve color, as in \code{cch_plot_curve}}
    \item{col2}{Color for unobserved ages, as in \code{cch_plot_curve}}
    \item{fskillonly}{When \code{TRUE}, plots only raw skill}
    \item{dosample}{When an integer, samples that number of new foragers from posterior}
    \item{draw_globalmean}{Whether to add a plot of top-level mean over all sites}
    \item{draw_legend}{Adds a legend plot}
    \item{alpha}{Transparency}
    \item{show_points}{As in \code{cch_plot_curve}}
    \item{lwd}{Line width}
    \item{adj_margins}{Shrinks margins}
    \item{part}{As in \code{cch_plot_curve}}
    \item{skip}{Number of spots to skip at start}
}
\details{
  This functions uses posterior samples in a \code{post} to draw skill or production curves for a site. It just calls \code{cch_plot_curve} for each element in \code{map_id}.
}
\value{
}
\references{}
\seealso{  }
\examples{

# for main figure
cch_plot_grid( 1:40 , nrow=7 , ncol=6 , col=cchpal[4] , fskillonly=TRUE , dosample=FALSE , draw_globalmean=TRUE , draw_legend=TRUE , alpha=0.5 , show_points=FALSE , lwd=1.5 )

# sample from skill functions to show individual variation better (SUPPLEMENTAL)
set.seed(9000)
cch_plot_grid( 1:40 , nrow=7 , ncol=6 , col=cchpal[4] , fskillonly=TRUE , dosample=20 , draw_globalmean=TRUE , draw_legend=TRUE , alpha=0.5 , show_points=FALSE , lwd=1.5 )

# success
cch_plot_grid( 1:40 , nrow=7 , ncol=6 , col=cchpal[5] , fskillonly=FALSE , dosample=FALSE , draw_globalmean=FALSE , draw_legend=FALSE , alpha=0.5 , show_points=FALSE , lwd=1.5 , part=1 , skip=2 )

# harvest
cch_plot_grid( 1:40 , nrow=7 , ncol=6 , col=cchpal[3] , fskillonly=FALSE , dosample=FALSE , draw_globalmean=FALSE , draw_legend=FALSE , alpha=0.5 , show_points=FALSE , lwd=1.5 , part=2 , skip=2 )

# production
cch_plot_grid( 1:40 , nrow=7 , ncol=6 , col="#D4ACE8" , fskillonly=FALSE , dosample=FALSE , draw_globalmean=FALSE , draw_legend=FALSE , alpha=0.5 , show_points=FALSE , lwd=1.5 , skip=2 )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

