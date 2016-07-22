# forager simulations
# N : number of foragers
# age : initial ages
# N_years : number of years to run sim
# N_days : hunting days per year
# chance_hunt : prob hunt per day per age
# k : rates of knowledge increase
# m : rates of senescence
# b : allometry parameters (importance of knowledge to production)
#   [ S = exp(-m*x)*(1-exp(-k*x))^b ]

# utility random number function
# same as rmvnorm2 in rethinking package
rmvnorm2 <- function (n, Mu = rep(0, length(sigma)), sigma = rep(1, length(Mu)), 
    Rho = diag(length(Mu)), method = "chol") 
{
    DS <- diag(sigma)
    SIGMA <- DS %*% Rho %*% DS
    rmvnorm(n = n, mean = Mu, sigma = SIGMA, method = method)
}

sim_foragers <- function(
    N=30, N_years=10, N_days=20, N_foragers=1, soc_id=NA, 
    a=c(3,log(10),0.5), coop=c(0,0,-0.2,0.5,0,1), prob_drop_fail=0 , prob_pool=0 , 
    prob_assistant=0,
    duration_func=function(n,l=5) rpois(n,l)+1 , prob_hours_miss=0.1 ,
    Mu=c(log(0.02),log(0.04)), sigma=c(0.05,0.05), rho=0.3, 
    k,m,b,age,age_err,chance_hunt,ref_age=1) {

    # format for coop argument marginal effects: 
    # ( group hunt f , group hunt h , pooled f , pooled h , assist f , assist h )

    if (missing(age)) age <- sample( 5:60 , size=N , replace=TRUE )
    if (missing(age_err)) age_err <- sample( (-5:5) , size=N , replace=TRUE )
    if (missing(chance_hunt)) chance_hunt <- inv_logit(-4+0.25*(1:100))*0.8

    if (missing(k) & missing(m)) {
        # sample individual forager deviations from means
        v <- rmvnorm2( N , Mu=Mu , sigma=sigma , Rho=matrix(c(1,rho,rho,1),nrow=2) )
        # compose total individual forager parameters
        k <- exp(v[,1])
        m <- exp(v[,2])
    }
    if (missing(b)) 
        b <- rep( 2 , N )
    else
        b <- rep_len( b , N )

    af <- a[1]
    ah <- a[2]

    dat <- data.frame(soc_id=soc_id,id=0,age=0,age_obs=0,harvest=0,year=0,day=0,trip=0,hours=0,pooled=0,A1=0,A2=0,A3=0,A4=0)

    row <- 1
    day <- 1
    trip <- 1

    # different hunt functions, to handle different types of hunts
    skillV <- Vectorize(skill,c("x","k","m","b"))
    do_hunt <- function( id , pars , coop , merge=TRUE , assist_id , hours ) {
        # assist_id should be same length as id and contain 0 where no assistant for that forager
        # or set to missing or NA when no assistants for any forager
        n <- length(id)
        if ( missing(assist_id) ) n_assist <- 0
        else {
            n_assist <- sum(assist_id > 0)
        }
        s <- skillV(age[id]/ref_age,k[id],m[id],b[id])
        hours_term <- hours/(5+hours) # michaelis-menton style
        if ( merge==TRUE ) {
            # for individual and pooled harvests
            lm_f <- exp(pars[1] + coop[1]*(n-1) + coop[3]*n_assist + hours_term)
            fail <- rbinom( 1 , size=1 , prob=prob_fail(mean(s),lm_f) )
            lm_h <- exp(pars[2] + coop[2]*(n-1) + coop[4]*n_assist + hours_term)
            h <- (1-fail)*rgamma2( 1 , (lm_h)*mean(s) , pars[3] )
        } else {
            # for group trips with individualized harvests
            has_assistant <- ifelse(assist_id>0,1,0)
            lm_f <- exp(pars[1] + coop[1]*(n-1) + coop[3]*has_assistant + hours_term)
            fail <- rbinom( n , size=1 , prob=prob_fail(s,lm_f) )
            lm_h <- exp(pars[2] + coop[2]*(n-1) + coop[4]*has_assistant + hours_term)
            h <- (1-fail)*rgamma2( n , (lm_h)*s , pars[3] )
        }
        return( h )
    }

    for ( y in 1:N_years ) {
        for ( d in 1:N_days ) {
            # initialize new list of foragers not yet foraged this day
            available <- 1:N
            for ( i in 1:N ) {
                if ( (runif(1) < chance_hunt[age[i]]) & (i %in% available) ) {
                    # do a hunt
                    # get number of hunters
                    ids <- i
                    use_a <- a
                    ####################
                    # build foraging party
                    n_party <- sample(N_foragers,size=1)
                    if ( n_party > 1 ) {
                        # add in more foragers
                        if ( length(available) >= n_party ) {
                            ids <- c( ids , sample(available[-i],size=n_party-1,prob=chance_hunt[age[available[-i]]]) )
                            available <- available[-ids]
                            use_a[1] <- a[1] + rnorm(1,0,0.1)
                            use_a[2] <- a[2] + rnorm(1,0,0.1)
                        } else {
                            # not enough foragers remaining for today
                            n_party <- 1
                        }
                    }
                    ####################
                    # build assistants to party - maybe one for each forager
                    assistant_ids <- rep(0,n_party)
                    gets_assist <- rbinom(n_party,size=1,prob=prob_assistant)
                    n_assist <- sum(gets_assist)
                    if ( n_assist > 0 ) {
                        # sample assistants
                        for ( j in 1:n_party ) {
                            if ( length(available) > 0 & gets_assist[j]>0 ) {
                                assistant_ids[j] <- sample(available,size=1,prob=1-chance_hunt[age[available]])
                                available <- available[-assistant_ids[j]]
                            }
                        }#j
                    }
                    n_assist <- sum( assistant_ids > 0 )
                    ####################
                    # generate harvests
                    pooled <- FALSE
                    coop_use <- coop[c(1,2,5,6)]
                    if ( n_party > 1 ) {
                        if ( runif(1) < prob_pool ) {
                            # pooled hunt with merged harvests
                            pooled <- TRUE
                            coop_use <- coop[c(3,4,5,6)]
                        }
                    } 
                    hours <- duration_func(1)
                    harvest <- do_hunt(ids,use_a,coop_use,merge=pooled,assist_id=assistant_ids,hours=hours)
                    ####################
                    # record data
                    record_it <- TRUE
                    # chance to drop failed (zero harvest) trips
                    if ( sum(harvest)==0 & runif(1) < prob_drop_fail ) record_it <- FALSE
                    # chance to miss hours variable
                    if ( runif(1) < prob_hours_miss ) hours <- NA
                    if ( record_it==TRUE ) {
                        assistant_ids <- ifelse( assistant_ids==0 , NA, assistant_ids)
                        if ( n_party == 1 ) {
                            dat[row,] <- c( soc_id , i , age[i] , age[i]+age_err[i] , harvest , y , day , trip , hours , 0 , assistant_ids[1] , NA , NA , NA )
                            row <- row + 1
                        } else {
                            # pooled or group hunt
                            if ( length(harvest)==1 ) harvest <- rep(harvest,n_party)
                            hi <- 1
                            for ( j in ids ) {
                                if ( pooled==TRUE ) {
                                    assistant_ids <- c(assistant_ids[!is.na(assistant_ids)],NA,NA,NA,NA,NA)
                                    dat[row,] <- c( soc_id , j , age[j] , age[j]+age_err[j] , harvest[hi] , y , day , trip , hours , as.integer(pooled) , assistant_ids[1], assistant_ids[2] , assistant_ids[3] , assistant_ids[4] )
                                } else {
                                    # unique assistant and harvest on each row
                                    aid <- assistant_ids[hi]
                                    dat[row,] <- c( soc_id , j , age[j] , age[j]+age_err[j] , harvest[hi] , y , day , trip , hours , as.integer(pooled) , aid, NA , NA , NA )
                                }
                                row <- row + 1
                                hi <- hi + 1
                            }#j
                        }
                        trip <- trip + 1
                    }#record trip
                }#if hunt
            }#i
            day <- day + 1
        }#d
        age <- age + 1
    }#y

    attr(dat,"pars") <- list(k=k,m=m,b=b,a=a,coop=coop)

    return(dat)

}#sim_foragers

# dat <- sim_foragers( N_foragers=1 , prob_drop_fail=0 , prob_pool=0 , prob_assistant=0.5 , coop=c(0,0,0,0,0,1) , a=c(3,log(10),0.5) )

show_foragers <- function( sim , cols=c(rangi2,"red","black","green") ) {
    plot( NULL , xlim=c(0,80) , ylim=c(0,0.1) , xlab="age" , ylab="skill" )
    a <- attr(sim,"pars")
    if ( is.null(a$k) ) {
        for ( s in 1:length(a) ) {
            k <- a[[s]]$k
            m <- a[[s]]$m
            b <- a[[s]]$b
            for ( i in 1:length(k) )
                curve( skill(x,k[i],m[i],b[i]) , from=0 , to=80 , add=TRUE , col=col.alpha(cols[s],0.5) )
        }#s
    } else {
        # only one society
        k <- a$k
        m <- a$m
        b <- a$b
        for ( i in 1:length(k) )
            curve( skill(x,k[i],m[i],b[i]) , from=0 , to=80 , add=TRUE , col=col.alpha(rangi2,0.5) )
    }
}
# show_foragers(dat)

# function to produce data from more than one society
# just calls sim_foragers multuple times with proper arguments
sim_multi <- function( N_soc=2 , N=30 , N_years=10 , N_days=20 , N_foragers=1 ,
    a , coop=c(0,0,-0.2,0.5,0,1) , prob_drop_fail=0 , prob_pool=0 , prob_assistant=0,
    Mu_soc=c(log(0.02),log(0.04),2) , sigma_soc=c(0.05,0.05,0.2) , Rho_soc=diag(3) ,
    sigma_i=c(0.1,0.2) , rho_km_i=0.3 ,
    ... ) {

    dat_all <- list()

    n1 <- function(x,i) {
        if ( length(x)>1 ) 
            return(x[i])
        else 
            return(x[1])
    }

    if (missing(a)) {
        a <- rmvnorm2( N_soc , Mu=c(3,log(10),0) , sigma=c(0.5,1,0.1) , Rho=diag(3) )
        a[,3] <- exp(a[,3])
        print(a)
    }

    Mu_list <- rmvnorm2( N_soc , Mu=Mu_soc , sigma=sigma_soc , Rho=Rho_soc )
    print(Mu_list)

    for ( s in 1:N_soc ) {
        dat_all[[s]] <- sim_foragers( 
            N=n1(N) , 
            soc_id=s ,
            N_years=n1(N_years) ,
            N_days=n1(N_days) ,
            N_foragers=n1(N_foragers),
            a=a[s,],
            coop=coop,
            prob_drop_fail=n1(prob_drop_fail),
            prob_pool=n1(prob_pool),
            prob_assistant=n1(prob_assistant),
            Mu=Mu_list[s,1:2],
            b=Mu_list[s,3],
            sigma=sigma_i,
            rho=rho_km_i,
            ...
        )
    }

    dat <- dat_all[[1]]
    for ( i in 2:N_soc ) dat <- rbind( dat , dat_all[[i]] )
    # fix index variables
    id <- paste(dat$soc_id,dat$id)
    dat$id <- coerce_index(id)
    tid <- paste(dat$soc_id,dat$trip)
    dat$trip_id <- coerce_index(tid)
    did <- paste(dat$soc_id,dat$day)
    dat$day_id <- coerce_index(did)

    # correct naming convention for other input code
    dat$society_id <- dat$soc_id
    dat$soc_id <- NULL
    dat$forager_id <- dat$id
    dat$id <- NULL
    dat$trip_duration <- dat$hours
    dat$hours <- NULL

    pars_list <- list()
    for ( s in 1:N_soc ) pars_list[[s]] <- attr(dat_all[[s]],"pars")
    attr(dat,"pars") <- pars_list

    return(dat)

}

# dat <- sim_multi( N_soc=3 , N_foragers=1:3 , prob_drop_fail=0 , prob_pool=0 , prob_assistant=0.1 , prob_hours_miss=0.1 )

