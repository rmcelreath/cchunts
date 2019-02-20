\name{sim_foragers}
\alias{sim_foragers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate foraging data}
\description{
  Simulates foraging return records
}
\usage{
sim_foragers(
    N=30, N_years=10, N_days=20, N_foragers=1, soc_id=NA, 
    a=c(3,log(10),0.5), coop=c(0,0,-0.2,0.5,0,1), prob_drop_fail=0 , prob_pool=0 , 
    prob_assistant=0,
    duration_func=function(n,l=5) rpois(n,l)+1 , prob_hours_miss=0.1 ,
    Mu=c(log(0.02),log(0.04)), sigma=c(0.05,0.05), rho=0.3, 
    k,m,b,age,age_err,chance_hunt,ref_age=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{N}{Number of individuals to simulate}
    \item{N_years}{Number of years of data to simulate}
    \item{N_days}{Number of days per year to simulate}
    \item{N_foragers}{Number of foragers per trip, either a single value or a vector of values to sample from}
    \item{soc_id}{Optional ID value for society (used by \code{sim_multi})}
    \item{a}{Intercept parameter values}
    \item{coop}{Cooperative production coefficients}
    \item{prob_drop_fail}{Probability failed trips drop from data (simulates biased sampling)}
    \item{prob_pool}{Proportion of trips with pooled harvests}
    \item{prob_assistant}{Proportion trips with assistants}
    \item{duration_func}{Function that returns a distribution of trip durations. Default is 1 + Poisson with mean 5 hours.}
    \item{prob_hours_miss}{Proportion of trip with missing duration values.}
    \item{Mu}{Average log life history parameters (k and m, respectively)}
    \item{sigma}{Standard deviations of life history parameters (k and m, respectively)}
    \item{rho}{Correlation between k and m}
    \item{k,m,b}{Optional vectors of life history parameters to use}
    \item{age}{Optional vector of forager ages}
    \item{age_err}{Amount of observation error for age}
    \item{chance_hunt}{Probability a forager goes foraging on any particular sample day}
    \item{ref_age}{Reference age (used to standardize scale of life history parameters)}
}
\details{
  This function
}
\value{
A data frame in which each row is a returned harvest and containing the columns:
\item{ trip_id }{Unique index of trip, across all sites}

}
\references{}
\seealso{ \code{sim_multi} }
\examples{

dat <- sim_foragers( N_foragers=1 , prob_drop_fail=0 , prob_pool=0 , prob_assistant=0.5 , coop=c(0,0,0,0,0,1) , a=c(3,log(10),0.5) )
show_foragers(dat)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

