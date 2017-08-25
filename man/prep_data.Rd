\name{prep_data}
\alias{prep_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert raw data into format model expects}
\description{
  This function takes the results of \code{make_joint} and converts them into a list of the form expected by the Stan model used in the analysis.
}
\usage{
prep_data( dat , debug=FALSE , dogs_miss=0 , guns_miss=0 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{dat}{A result from \code{make_joint}}
    \item{debug}{Flag for showing various debugging information}
    \item{dogs_miss}{What code to assign to \code{NA} values for dogs on trips. Defaults to 0, which assumes dog absent. Assigning -1 flags for marginalization in model.}
    \item{guns_miss}{What code to assign to \code{NA} values for firearms on trips. Defaults to 0, which assumes firearms absent. Assigning -1 flags for marginalization in model.} 
}
\details{
  This functions does little more than flag age, dog, and firearms values for imputation or marginalization, while constructing a matrix of participants for each trip. It returns a \code{list} appropriate for passing to the Stan model.
}
\value{
% cat(paste( "\\item{" , names(dat_list) , "}{ description }" , collapse="\n") )
A list suitable for a Stan model:
\item{ N }{ Number of harvests }
\item{ N_hunters }{ Number of individual foragers }
\item{ N_societies }{ Number of sites }
\item{ hunter_id }{ Unique identifier for each forager }
\item{ forager_female }{ Indicator of whether forager is female }
\item{ soc_id }{ Unique identifier for each site }
\item{ hours }{ Standardized duration of trip }
\item{ day_trip }{ Indicator of whether trip was under 24 hours in duration }
\item{ harvest }{ Standardized harvest amount }
\item{ age_interval }{ Indicator of when onl interval information available for age }
\item{ age }{ First age value - central estimate or lower bound }
\item{ age2 }{ Second age value - standard error or upper bound }
\item{ trip_year_offset }{ Used internally to track ages at each site }
\item{ ref_age }{ Age to use in standardizing age values internally. Defaults to 80. }
\item{ N_trips }{ Number of trips }
\item{ trip_id }{ Unique identifier for each trip }
\item{ pooled }{ Indicator of pooled harvest }
\item{ A }{ Assistant matrix for each trip }
\item{ max_foragers_per_trip }{ Largest number of foragers on any trip. Allocates memory in model. }
\item{ dogs }{ Idicator of use of dogs on trip }
\item{ firearms }{ Idicator of use of firearms on trip }
}
\references{}
\seealso{ \code{make_joint}, \code{trans_trips} }
\examples{
dat <- make_joint()
dat_list <- prep_data( dat , dogs_miss=-1 , guns_miss=-1 )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

