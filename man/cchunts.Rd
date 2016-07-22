\name{cchunts}
\alias{cchunts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Human forager hunting returns data}
\description{
  Available data sets in package cchunts.
}
\usage{
data(Alvard)
data(Beckerman)
data(Bird_Bird_Codding)
data(Coad)
data(Duda)
data(Ellen)
data(Fernandez_Llamazares)
data(Franzen)
data(Gallois)
data(Gueze)
data(Headland)
data(Healey_Nen_PNG)
data(Healey)
data(Hill_Kintigh)
data(Koster)
data(Kramer_Greaves)
data(Lupo_Schmitt)
data(Marks)
data(Napitupulu)
data(Nielsen)
data(Pacheco)
data(Pangau_Adam)
data(Ready)
data(Reyes_Garcia)
data(Sillitoe)
data(Siren)
data(Trumble_Gurven)
data(Van_Vliet_et_al_Africa_three_sites)
data(Van_Vliet_et_al_Baego)
data(Van_Vliet_et_al_Djoutou)
data(Van_Vliet_et_al_Gabon)
data(Van_Vliet_et_al_Ingolo)
data(Van_Vliet_et_al_Ngombe)
data(Van_Vliet_et_al_Ovan)
data(Van_Vliet_et_al_Phalanga)
data(Van_Vliet_et_al_South_America_sites)
data(Venkataraman_et_al)
data(Winterhalder)
data(Yu_et_al)
data(Ziker)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\format{
    \enumerate{
        \item trip_id : Identifier for trip
        \item trip_id_orig : Identifier in raw data (for reference and cross-check)
        \item observed : 0/1 indicator of whether trip was directory observed or verbally reported
        \item trip_date : Calendar date on which trip occurred
        \item julian_date_s : Standardized Julian date of trip
        \item day_id : Identifier for unique day in sample
        \item trip_duration : Duration in hours
        \item day_trip : 0/1 indicator of whether was a single day trip (1) or multiple days with camping (0)
        \item group_type : \code{Independent}, \code{Cooperative}, or \code{Pooled} foraging trips
        \item pooled : 0/1 indicator for pooled trips
        \item harvest : Harvest in kilograms
        \item forager_id : Unique indentifier for forager
        \item forager_id_orig : ID in raw data (for reference and cross-check)
        \item age_type : Either \code{Exact} age or \code{Interval} age
        \item age_dist_1 : When age_type \code{Exact}, age in years. Otherwise lower bound of interval.
        \item age_dist_2 : When age_type \code{Exact}, standard error of age. Otherwise upper bound of interval.
        \item dogs : 0/1 indicator of whether dogs were used on trip
        \item gun : 0/1 indicator of firearm use
        \item a_n_x : Variables referencing assistant number \code{n}, differing by availability in each data set
    }
}
\value{
}
\references{Hill and Kintigh. 2009. Current Anthropology 50:369-377.}
\seealso{}
\examples{
data(Alvard)
head(Alvard)
plot( harvest ~ trip_duration , Alvard )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

