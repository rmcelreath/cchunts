\name{cchmodel_run}
\alias{cchmodel_run}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compile and run statistical model of foraging returns}
\description{
  Loads data, compiles model, and samples from posterior.
}
\usage{
cchmodel_run( model_code=cchmodel_full , sites=cchunts_data_sets , simdat , dogs_miss=-1 , guns_miss=-1 , prior_scale=0.5 , n_chains=4 , the_seed=1208 , warmup=500 , iter=1000 , control=list( adapt_delta=0.99 , max_treedepth=13 ) , start , test=FALSE , ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{model_code}{A character vector specifying a Stan model. Defaults to \code{cchmodel_full}, the model used in the analysis.}
    \item{sites}{Character vector of sites to include. Defaults to full list in \code{cchunts_data_sets}}
    \item{simdat}{When present, a table of simulated data to use, produced by sim_multi.}
    \item{dogs_miss}{Value to use for \code{NA} values in dogs variable. -1 indicates marginalize in model.}
    \item{guns_miss}{Value to use for \code{NA} values in firearms variable. -1 indicates marginalize in model.}
    \item{prior_scale}{Common factor used to adjust scale of many priors in model.}
    \item{n_chains}{Number of chains to run.}
    \item{the_seed}{Random number seed to use. Default (1208) was used in described results.}
    \item{warmup}{\code{stan} warmup argument}
    \item{iter}{\code{stan} iter argument}
    \item{control}{\code{stan} control parametert list}
    \item{start}{Optional list of parameter inits. When left missing, default is zero inits for many parameters.}
    \item{test}{Flag to just test compilation and do 4 samples from a single chain.}
    \item{...}{Additonal arguments, to pass to \code{stan}}
}
\details{
  This function takes the code for a Stan model, a list of sites, and various optional arguments. It then calls \code{make_joint}, \code{prep_data}, and executes the \code{stan} model.
}
\value{
When successful, an object of class \code{stanfit}.
}
\references{}
\seealso{ \code{prep_data}, \code{make_joint} }
\examples{

# test compilation and sampling 
# for the full model on all sites
# set test=FALSE to run for real --- will take a while!
library(rstan)
test_fit <- cchmodel_run( test=TRUE )

# plot grid of skill functions
dat <- make_joint()
dat_list <- prep_data( dat , dogs_miss=-1 , guns_miss=-1 )
post <- extract(test_fit)
cch_plot_grid( 1:40 , nrow=7 , ncol=6 , col=cchpal[4] , fskillonly=TRUE , 
    dosample=FALSE , draw_globalmean=TRUE , draw_legend=TRUE , alpha=0.5 , 
    show_points=FALSE , lwd=1.5 )

# see the Stan model code
# this is long and heavily commented
cat( cchmodel_full )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

