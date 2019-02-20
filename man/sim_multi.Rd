\name{sim_multi}
\alias{sim_multi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate foraging data from multiple sites}
\description{
  Simulates foraging return records from multiple sites
}
\usage{
sim_multi( N_soc=2 , N=30 , N_years=10 , N_days=20 , N_foragers=1 ,
    a , coop=c(0,0,-0.2,0.5,0,1) , prob_drop_fail=0 , prob_pool=0 , prob_assistant=0,
    Mu_soc=c(log(0.02),log(0.04),2) , sigma_soc=c(0.05,0.05,0.2) , Rho_soc=diag(3) ,
    sigma_i=c(0.1,0.2) , rho_km_i=0.3 , new_names=TRUE , 
    ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{N_soc}{Number of sites to simulate}
    \item{N}{Number of individuals at each site}
    \item{N_years}{Number of years of data to simulate}
    \item{N_days}{Number of days per year to simulate}
    \item{N_foragers}{Number of foragers per trip, either a single value or a vector of values to sample from}
    \item{a}{Intercept parameter values}
    \item{coop}{Cooperative production coefficients}
    \item{prob_drop_fail}{Probability failed trips drop from data (simulates biased sampling)}
    \item{prob_pool}{Proportion of trips with pooled harvests}
    \item{prob_assistant}{Proportion trips with assistants}
    \item{Mu_soc}{Average life history parameters (k and m and b respectively)}
    \item{sigma_soc}{Stddev of life history parameters}
    \item{Rho_soc}{Correlations among life history parameters}
    \item{sigma_i}{Standard deviations of life history parameters (k and m, respectively) in each site}
    \item{rho_km_i}{Correlation between k and m in each site}
    \item{new_names}{When TRUE names output so same as statistical model expects}
    \item{...}{Additonal arguments to pass to sim_foragers}
}
\details{
  This function calls sim_foragers multiple times to assemble a cross-site simulated data set. This was used to validate and design the statistical model.
}
\references{}
\seealso{ \code{sim_foragers} }
\examples{

datsim <- sim_multi( N_soc=3 , N=20 , N_years=5 , N_days=10 , N_foragers=1 , prob_drop_fail=0 , prob_pool=0 , prob_assistant=0 , prob_hours_miss=0.01 )

library(rstan)
sim_fit <- cchmodel_run( simdat=datsim , test=TRUE )



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

