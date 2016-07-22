# common functions for posterior predictions

cch_predict <- function( post , data , func , verbose=FALSE , logh=TRUE , ... ) {
    d <- as.data.frame(data)
    n <- nrow(d)
    # fill in necessary but missing variables
    if ( is.null(d$soc_id) ) {
        d$soc_id <- 1
        if ( verbose==TRUE ) message("No soc_id in data. Using '1'.")
    }
    fhid <- 1 # flag to use hunter id
    if ( is.null(d$hunter_id) ) {
        d$hunter_id <- 1
        fhid <- 0
        if ( verbose==TRUE ) message("No hunter_id in data. Producing predictions for average individual.")
    }
    if ( is.null(d$hours) ) {
        d$hours <- 1
        if ( verbose==TRUE ) message("No hours in data. Using '1'.")
    }
    if ( is.null(d$age) ) {
        d$age <- 0.5
        if ( verbose==TRUE ) message("No age in data. Using '40' (standardized to '0.5').")
    }
    if ( is.null(d$n_foragers) ) d$n_foragers <- 1
    if ( is.null(d$n_assistants) ) d$n_assistants <- 0
    if ( is.null(d$dogs) ) d$dogs <- 0
    if ( is.null(d$firearms) ) d$firearms <- 0
    if ( is.null(d$pooled) ) d$pooled <- 0

    # detect non-zero intercept in skill functions
    flag_nz <- FALSE
    if ( dim(post$lifehistmeans)[2]==4 ) {
        flag_nz <- TRUE
    }

    # predict trip failures
    p_f <- sapply( 1:n, 
        function(i) {
            x <- d$age[i]
            j <- d$soc_id[i]
            h <- d$hours[i]
            hid <- d$hunter_id[i]
            n_foragers <- d$n_foragers[i]
            n_assistants <- d$n_assistants[i]
            dogs <- d$dogs[i]
            firearms <- d$firearms[i]
            fpool <- d$pooled[i]
            all_af <- post$af[,1,j] + 
                      post$af[,2+fpool,j]*(n_foragers-1) +
                      post$af[,4,j]*n_assistants +
                      post$b_dogs[,1,j]*dogs +
                      post$b_firearms[,1,j]*firearms
            labor_af <- h^post$b_hours[,1,j]
            if ( logh==FALSE ) {
                all_af <- all_af + post$b_hours[,1,j]*h
                labor_af <- 1
            }
            skill_elasticity <- 1
            if ( !is.null(post$sef) ) {
                # skill elasticity model, so pull samples out
                skill_elasticity <- exp(post$sef[,j])
            }
            if ( flag_nz==FALSE ) {
                pf <- prob_fail(
                    skill(x, #age (standardized to 80 = 1)
                        (exp(post$lifehistmeans[,1] + post$vs[,j,1] + post$vh[,hid,1]*fhid)), #k
                        (exp(post$lifehistmeans[,2] + post$vs[,j,2] + post$vh[,hid,2]*fhid)), #m
                        (exp(post$lifehistmeans[,3] + post$vs[,j,3])))^skill_elasticity, #b
                    exp(all_af) * labor_af #intercept
                )
            } else {
                # non-zero intercept function
                pf <- prob_fail(
                    skillnz(x, #age (standardized to 80 = 1)
                        (exp(post$lifehistmeans[,1] + post$vs[,j,1] + post$vh[,hid,1]*fhid)), #k
                        (exp(post$lifehistmeans[,2] + post$vs[,j,2] + post$vh[,hid,2]*fhid)), #m
                        (exp(post$lifehistmeans[,3] + post$vs[,j,3])), #b
                        (exp(post$lifehistmeans[,4] + post$vs[,j,4])) )^skill_elasticity, 
                    exp(all_af) * labor_af #intercept
                )
            }
            return(pf)
        })
    
    # predict trip harvests
    p_h <- sapply( 1:n,
        function(i) {
            x <- d$age[i]
            j <- d$soc_id[i]
            h <- d$hours[i]
            hid <- d$hunter_id[i]
            n_foragers <- d$n_foragers[i]
            n_assistants <- d$n_assistants[i]
            dogs <- d$dogs[i]
            firearms <- d$firearms[i]
            fpool <- d$pooled[i]
            all_ah <- post$ah[,1,j] + 
                      post$ah[,2+fpool,j]*(n_foragers-1) +
                      post$ah[,4,j]*n_assistants +
                      post$b_dogs[,2,j]*dogs +
                      post$b_firearms[,2,j]*firearms
            labor_ah <- h^post$b_hours[,2,j]
            if ( logh==FALSE ) {
                all_ah <- all_ah + post$b_hours[,2,j]*h
                labor_ah <- 1
            }
            skill_elasticity <- 1
            if ( !is.null(post$seh) ) {
                # skill elasticity model, so pull samples out
                skill_elasticity <- exp(post$seh[,j])
            }
            if ( flag_nz==FALSE ) {
                ph <- exp(all_ah)*labor_ah* #intercept
                    skill(x,
                        (exp(post$lifehistmeans[,1] + post$vs[,j,1] + post$vh[,hid,1]*fhid)), #k
                        (exp(post$lifehistmeans[,2] + post$vs[,j,2] + post$vh[,hid,2]*fhid)), #m
                        (exp(post$lifehistmeans[,3] + post$vs[,j,3])) #b
                    )^skill_elasticity
            } else {
                ph <- exp(all_ah)*labor_ah* #intercept
                    skillnz(x,
                        (exp(post$lifehistmeans[,1] + post$vs[,j,1] + post$vh[,hid,1]*fhid)), #k
                        (exp(post$lifehistmeans[,2] + post$vs[,j,2] + post$vh[,hid,2]*fhid)), #m
                        (exp(post$lifehistmeans[,3] + post$vs[,j,3])) , #b
                        (exp(post$lifehistmeans[,4] + post$vs[,j,4])) #a
                    )^skill_elasticity
            }
            return(ph)
        })
    # result
    result <- list( failure=p_f , harvest=p_h )
    if ( !missing(func) ) {
        r2 <- result
        for ( i in 1:2 ) r2[[i]] <- apply(r2[[i]],2,func,...)
        result <- r2
    }
    return(result)
}

# pdat <- data.frame( age=seq(from=0,to=1,length.out=50) )
# l <- cch_predict( post , pdat )
# cch_predict( post , pdat , func=mean )
