\name{make_joint}
\alias{make_joint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Merge individual foraging data sets}
\description{
  Merges available data sets in package cchunts.
}
\usage{
make_joint( data_sets , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{data_sets}{A character vector of named data sets to load and combine. If omitted, defaults to all data sets listed in \code{?cchunts}.}
    \item{...}{Additonal arguments, currently unused}
}
\details{
  This function provides a way to merge more than one sites data into a single data frame. The result is suitable for the analysis functions. Minimal consistency checks are performed, and all index variables are managed so that individuals and trips remain unique across sites. This is done by reconstructing index variables after merging.
}
\value{
% cat(paste( "\\item{" , names(all_data) , "}{ description }" , collapse="\n") )
A data frame in which each row is a returned harvest and containing the columns:
\item{ trip_id }{Unique index of trip, across all sites}
\item{ trip_id_orig }{Trip identifier in original data. Needed for integrity checks.}
\item{ observed }{ 0/1 indicator of whether trip was observed or is verbal report }
\item{ trip_date }{ Raw date string for trip }
\item{ julian_date_s }{ Standardized Julian date for trip (standardized within each site) }
\item{ day_id }{ Unique index for each day within each site }
\item{ trip_duration }{ Duration of trip, in hours }
\item{ day_trip }{ 0/1 indicator of trip being less than a single day (1) or exceeding more than one day (0) }
\item{ group_type }{ String indicating type of foarging group. One of: \code{Independent} (indicating independent foragers), \code{Cooperative} (cooperative foraging groups), \code{Loose} (looser but potentially cooperative groups), \code{Unknown} (Unknown group type) }
\item{ pooled }{ 0/1 indicator of a pooled harvest, not assignable to the individual }
\item{ harvest }{ Harvest amount in kilograms }
\item{ forager_id }{ Unique forager identifier, across all sites }
\item{ forager_id_orig }{ Forager identifier in original data, for consistency checks }
\item{ age_type }{ String identifying type of age measurement. One of: \code{Exact}, \code{Uncertain}, \code{Interval} }
\item{ age_dist_1 }{ Either center of age estimate (for \code{Exact} or \code{Uncertain} estimates) or bottom of \code{Interval} }
\item{ age_dist_2 }{ Either standard error of estimate (\code{Exact} or \code{Uncertain}) or top of \code{Interval} }
\item{ dogs }{ 0/1 indicator of whether dogs participated in harvest }
\item{ gun }{ 0/1 indicator of whether firearms were used }
\item{ a_n_id }{ Assistant n identifier }
\item{ a_n_sex }{ Assistant n sex }
\item{ a_n_age_dist_1 }{ Assistant n age estimate (same meaning as \code{age_dist_1}) }
\item{ a_n_age_dist_2 }{ Assistant n age estimate (same meaning as \code{age_dist_2}) }
\item{ a_n_firearms }{ Whether assistant n used firearms }
\item{ society }{ Data set label (usually researcher name or names) }
\item{ society_id }{ Index for society }
\item{ trip_id_soc }{ Index of trip within society (not unique across societies) }
\item{ day_id_soc }{ Index of day within society (not unique across societies) }
\item{ forager_id_soc }{ Index of forager within society (not unique across societies }
}
\references{}
\seealso{ \link{\code{prep_data}}, \link{\code{trans_trips}} }
\examples{
some_data <- make_joint(c("Koster","Trumble_Gurven"))

all_data <- make_joint()

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

