if_else <- function(x,y,z) ifelse(x,y,z)
int_step <- function(x) ifelse(x>0,TRUE,FALSE)

purgeNA <- function(x,by=x) x[ !is.na(by) ]

checkNA <- function(x) {
    if ( any(is.na(x)) ) {
        i <- which(is.na(x))
        print(i)
    }
}

na_to_x <- function(x,y=0) ifelse(is.na(x),y,x)
skill <- function(x,k,m,b,a=1) a*( exp(-m*x)*(1-exp(-k*x))^b )

# skill function with non-zero intercept 'a'
skillnz <- function(x,k,m,b,a=0) {
    s1 <- exp(-m*x)*(1-exp(-k*x))^b
    s2 <- s1 * inv_logit(100*(x-a))
    return(s2)
}
# curve( skillnz(x,0.05,0.01,2,10) , from=0, to=80)

skill_std <- function(x,k,m,b) exp(-m*x)*(1-exp(-k*x))^b * (1+b*k/m)^(m/k) * (b*k/(b*k+m))^(-b)
prob_fail <- function(s,a) 2*(1 - inv_logit( s * a ))

hid_in_soc <- function(soc_id) {
    message(unique(dat$society[dat_list$soc_id==soc_id]))
    sort(unique(dat_list$hunter_id[dat$society_id==soc_id]))
}

cchunts_data_sets <- c(
'Alvard',
'Beckerman',
'Bird_Bird_Codding',
'Coad',
'Duda',
'Ellen',
'Fernandez_Llamazares',
'Franzen',
'Gallois',
'Gueze',
'Headland',
'Healey_Nen_PNG',
'Healey',
'Hill_Kintigh',
'Koster',
'Kramer_Greaves',
'Lupo_Schmitt',
'Marks',
'Napitupulu',
'Nielsen',
'Pacheco',
'Pangau_Adam',
'Ready',
'Reyes_Garcia',
'Sillitoe',
'Siren',
'Trumble_Gurven',
#'Van_Vliet_et_al_Africa_three_sites',
'Van_Vliet_et_al_Baego',
'Van_Vliet_et_al_Djoutou',
'Van_Vliet_et_al_Gabon',
'Van_Vliet_et_al_Ingolo',
'Van_Vliet_et_al_Ngombe',
'Van_Vliet_et_al_Ovan',
'Van_Vliet_et_al_Phalanga',
'Van_Vliet_et_al_South_America_sites',
'Venkataraman_et_al',
'Winterhalder',
'Yu_et_al',
'Ziker',
'Ross' # added most recently
)

do_counts <- function( data_sets , ... ) {

    if ( missing(data_sets) ) {
        data_sets <- cchunts_data_sets
    }
    
    # initial merge
    for ( i in 1:length(data_sets) ) {
        do.call( data , list(data_sets[i]) )
        # get copy of it from global environment
        d <- get( data_sets[i] , envir=.GlobalEnv )
        print(data_sets[i])
        print('foragers:')
        try(check_index(d$forager_id))
        print('trips:')
        try(check_index(d$trip_id))
    }

}

make_joint <- function( data_sets , ... ) {

    if ( missing(data_sets) ) {
        data_sets <- cchunts_data_sets
    }
    
    # pass through individual data tables and gather column names
    all_d <- list()
    varlist <- c()
    for ( i in 1:length(data_sets) ) {
        # load
        do.call( data , list(data_sets[i]) )
        # get copy of it from global environment
        d <- get( data_sets[i] , envir=.GlobalEnv )
        d$society <- as.character(rep( data_sets[i] , nrow(d) ))
        all_d[[i]] <- d
        # union of column names
        varlist <- union( varlist , names(d) )
    }#i
    # now pass over again and row bind, filling in missing columns with NA
    for ( i in 1:length(data_sets) ) {
        for ( j in 1:length(varlist) ) {
            avar <- varlist[j]
            if ( is.null(all_d[[i]][[avar]]) ) {
                all_d[[i]][[avar]] <- rep(NA,nrow(all_d[[i]]))
            }
        }#j
        if ( i==1 ) {
            joint_data <- all_d[[i]]
        } else {
            joint_data <- rbind(joint_data,all_d[[i]])
        }
    }#i

    # now housekeeping
    # need to delete some rows with missing values

    # find NA harvests
    idx <- which( is.na(joint_data$harvest) )
    if ( length(idx)>0 ) {
        message(concat("Found ",length(idx)," NA harvests:"))
        print( data.frame( idx , harvest=joint_data$harvest[idx] , society=joint_data$society[idx] , forager=joint_data$forager_id[idx] ) )
        message("Removing these from the joint data.")
        joint_data <- joint_data[ -idx , ]
        idx2 <- which( is.na(joint_data$harvest) )
        if ( length(idx2)>0 ) stop("Still NA harvests somehow! Panic!")
    }

    # find NA forager IDs
    idx <- which( is.na(joint_data$forager_id) )
    if ( length(idx)>0 ) {
        message(concat("Found ",length(idx), " NA forager_id:"))
        print(idx)
        joint_data <- joint_data[ -idx , ]
        message("Removing these rows from joint data.")
    }

    # make society index
    joint_data$society_id <- coerce_index(joint_data$society)
    
    # make new unique trip_id and hunter_id
    tid <- paste( joint_data$trip_id , joint_data$society_id )
    tid <- coerce_index(tid)
    joint_data$trip_id_soc <- joint_data$trip_id
    joint_data$trip_id <- tid

    did <- paste( joint_data$day_id , joint_data$society_id )
    did <- coerce_index(did)
    joint_data$day_id_soc <- joint_data$day_id
    joint_data$day_id <- did
    
    # have to merge all the assistant ids as well
    idx <- grep("a_._id",names(joint_data))
    fid <- list()
    fid[[1]] <- joint_data$forager_id + joint_data$society_id * 1000
    joint_data$forager_id_soc <- joint_data$forager_id
    if ( length(idx)>0 ) {
        # for each assistant ID variable, merge with society ID
        for ( i in 1:length(idx) ) {
            fid[[i+1]] <- joint_data[[idx[i]]] + joint_data$society_id * 1000
        }
        hid <- do.call(coerce_index,fid)
        joint_data$forager_id <- hid[[1]]
        for ( i in 1:length(idx) ) {
            joint_data[[idx[i]]] <- hid[[i+1]]
        }
    } else {
        hid <- coerce_index(fid[[1]])
        joint_data$forager_id <- hid
    }
    
    message("Trips:")
    check_index(joint_data$trip_id)
    message("Unique foragers (incl. assistants):")
    check_index(c(joint_data$forager_id,unlist(joint_data[,idx])))
    message("Societies:")
    check_index(joint_data$society_id)
    message(paste(unique(joint_data$society),collapse="\n"))
    
    return( invisible( joint_data ) )
}

# preps joined data sets so ready to feed into Stan model
prep_data <- function( dat , debug=FALSE , dogs_miss=0 , guns_miss=0 ) {

    N_soc <- length(unique(dat$society_id)) 

    # find NA pooled
    idx <- which( is.na(dat$pooled) )
    if ( length(idx)>0 ) {
        message(concat("Found ",length(idx), " NA pooled values:"))
        print(data.frame(idx,society=dat$society[idx],forager=dat$forager_id[idx]))
        message("Marking these as pooled==0.")
        dat$pooled[idx] <- rep(0,length(idx))
    }

    # find NA dogs
    if ( !is.null(dat$dogs) ) {
        if ( dogs_miss==0 ) {
            # just replace all NA with 0, as if dogs absent
            idx <- which( is.na(dat$dogs) )
            if ( length(idx)>0 ) {
                message(concat("Found ",length(idx), " NA dogs values:"))
                #print(data.frame(idx,society=dat$society[idx],forager=dat$forager_id[idx]))
                message("Marking these as dogs==0.")
                dat$dogs[idx] <- rep(0,length(idx))
            }
        }
        if ( dogs_miss!=0 ) {
            # mark as missing with -1, but only if site has some non-NA values
            for ( i in 1:N_soc ) {
                # fetch dogs in society i
                dogsi <- dat$dogs[dat$society_id==i]
                idx <- which( is.na(dat$dogs) & dat$society_id==i )
                if ( any(is.na(dogsi)) ) {
                    if ( all(is.na(dogsi)) ) {
                        # all NA, so just set all to zero
                        # if all missing, nothing to gain by marginalizing
                        dat$dogs[idx] <- rep(0,length(idx))
                    } else {
                        # replace NAs with dogs_miss
                        dat$dogs[idx] <- rep(dogs_miss,length(idx))
                    }
                }
            }#i
        }
    }

    # find NA guns
    if ( !is.null(dat$gun) ) {
        if ( guns_miss==0 ) {
            # just replace all NA with 0, as if guns absent
            idx <- which( is.na(dat$gun) )
            if ( length(idx)>0 ) {
                message(concat("Found ",length(idx), " NA firearms values:"))
                #print(data.frame(idx,society=dat$society[idx],forager=dat$forager_id[idx]))
                message("Marking these as 0.")
                dat$gun[idx] <- rep(0,length(idx))
            }
        }
        if ( guns_miss!=0 ) {
            # mark as missing with -1, but only if site has some non-NA values
            for ( i in 1:N_soc ) {
                # fetch guns in society i
                gunsi <- dat$gun[dat$society_id==i]
                idx <- which( is.na(dat$gun) & dat$society_id==i )
                if ( any(is.na(gunsi)) ) {
                    if ( all(is.na(gunsi)) ) {
                        # all NA, so just set all to zero
                        # if all missing, nothing to gain by marginalizing
                        dat$gun[idx] <- rep(0,length(idx))
                    } else {
                        # replace NAs with dogs_miss
                        dat$gun[idx] <- rep(guns_miss,length(idx))
                    }
                }
            }#i
        }
    }

    # calculate year offset within each society
    # need this to correctly impute ages that are of type "Uniform"
    # calculation is just:
    # (1) find earliest date in society
    # (2) for each harvest, enter harvest year - earliest year
    # then in model, we add this offset to uniform imputation for each forager
    # this procedure obeys uniform interval on age while also aging forager by a year each year
    dat$trip_year_offset <- rep(NA,nrow(dat))
    for ( j in 1:N_soc ) {
        soc_name <- unique(dat$society[dat$society_id==j])
        min_year <- 0
        #trip_year <- as.integer(substr(dat$trip_date[dat$society_id==j],1,4))
        trip_year <- as.integer( dat$trip_year[dat$society_id==j] )
        min_year <- min( trip_year , na.rm=TRUE )
        if ( min_year==Inf ) {
            message(concat("BEWARE: No trip years recorded for ",soc_name,". Using 1970 for min_year."))
            min_year <- 1970
        }
        if ( debug==TRUE ) print(paste(soc_name," - ",min_year))
        dat$trip_year_offset[dat$society_id==j] <- trip_year - min_year
    }
    dat$trip_year_offset <- ifelse(is.na(dat$trip_year_offset),0,dat$trip_year_offset)

    # prep age data input for Stan
    age_obs <- ifelse( dat$age_type %in% c("Exact","Uncertain") , dat$age_dist_1 , (dat$age_dist_1+dat$age_dist_2)/2 )
    age2 <- dat$age_dist_2
    sdunif <- function(a,b) sqrt((1/12)*(a-b)^2)
    age2 <- ifelse( dat$age_type %in% c("Exact","Uncertain") , age2 , sdunif(dat$age_dist_1,dat$age_dist_2) )
    age2 <- ifelse( dat$age_type == "Uncertain" , age2/2 , age2 ) # Uncertain coded as 95% interval, not standard error
    # formula above inserts stddev of uniform interval

    # count max number of foragers on a single trip (needed internally for matrix declaration)
    max_n <- 0
    tid_list <- sort(unique(dat$trip_id))
    for ( j in tid_list ) {
        n <- length(which(dat$trip_id==j))
        if ( n > max_n ) max_n <- n
    }

    # build assistant matrix
    idx <- grep( "a_._id" , names(dat) )
    if ( length(idx)==0 ) idx <- grep( "A" , names(dat) ) # simulated data uses A# format
    Assist_matrix <- matrix( 0 , nrow=nrow(dat) , ncol=9 )
    for ( i in 1:length(idx) ) Assist_matrix[,i] <- na_to_x(dat[[idx[i]]])

    # count number of societies and unique foragers (including assistants)
    Total_foragers <- max( c( dat$forager_id , unlist(dat[,idx]) ) , na.rm=TRUE )

    # Lupo_Schmitt rows with missing durations but should not be
    if ( FALSE ) {

    idx <- which(dat$society=="Lupo_Schmitt" & dat$trip_id_soc==1 & is.na(dat$trip_duration) )
    if ( length(idx) > 0 )
        dat$trip_duration[idx] <- 4.3 # duration from previous
    idx <- which(dat$society=="Lupo_Schmitt" & dat$trip_id_soc==22 & is.na(dat$trip_duration) )
    if ( length(idx) > 0 )
        dat$trip_duration[idx] <- 5.3
    idx <- which(dat$society=="Lupo_Schmitt" & dat$trip_id_soc==23 & is.na(dat$trip_duration) )
    if ( length(idx) > 0 )
        dat$trip_duration[idx] <- 4

    }

    # standardize harvests and trip durations
    hstd <- dat$harvest
    duration_std <- dat$trip_duration
    for ( j in 1:N_soc ) {
        idx <- dat$society_id==j
        hmean <- mean( dat$harvest[idx] )
        hstddev <- sd( dat$harvest[idx] )
        hstd[idx] <- ( hstd[idx] / hmean )#/hstddev
        dmean <- mean( dat$trip_duration[idx] , na.rm=TRUE )
        dstddev <- sd( dat$trip_duration[idx] , na.rm=TRUE )
        duration_std[idx] <- ( duration_std[idx] / dmean )#/dstddev
        if ( debug==TRUE ) {
            print(concat(j," H /",hmean,"/",hstddev))
            print(concat(j," TD /",dmean,"/",dstddev))
        }#debug
    }
    dat$hstd <- hstd
    dat$duration_std <- duration_std

    # single-day indicator variable
    # for now, set missing to zero
    dat$day_trip <- ifelse( is.na(dat$day_trip) , 1 , dat$day_trip )

    # check for empty dogs and gun columns
    if ( is.null(dat$dogs) ) dat$dogs <- 0
    if ( is.null(dat$gun) ) dat$gun <- 0

    # build data to pass to Stan
    dat_list <- list(
        N = nrow(dat),
        N_hunters = Total_foragers,
        N_societies = N_soc,
        hunter_id = dat$forager_id,
        forager_female = ifelse( dat$sex=="F" , 1 , 0 ),
        soc_id = dat$society_id,
        hours = na_to_x( dat$duration_std , -1 ), # -1 is mark for imputation in Stan code
        day_trip = dat$day_trip,
        harvest = hstd, #dat$harvest,
        age_interval = ifelse(dat$age_type=="Uniform",1,0)*0,
        age = age_obs,
        age2 = age2,
        trip_year_offset = dat$trip_year_offset,
        ref_age = 80,
        N_trips = max(dat$trip_id),
        trip_id = dat$trip_id,
        pooled = dat$pooled,
        A = Assist_matrix,
        max_foragers_per_trip = max_n,
        dogs = na_to_x(dat$dogs,0),
        firearms = na_to_x(dat$gun,0)
    )

    # check for zero trip durations
    for ( i in 1:length(dat_list$hours) ) {
        y <- dat_list$hours[i]
        if ( y > -1 ) {
            if ( log(y)==-Inf ) {
                message( concat("duration ZERO found at index ",i," - patching to exp(-5)") )
                dat_list$hours[i] <- exp(-5)
            }
        }
    }

    # find all NA age values and specify for imputation
    idx <- which( is.na(dat_list$age) )
    # cbind( idx , dat$society[idx] , dat$forager_id[idx] )
    if ( length(idx)>0 ) {
        message(concat("Found ",length(idx)," NA age values. Marking for imputation as: age ~ normal(42.5,10.8)"))
        #print( data.frame( idx , society=dat$society[idx] , forager=dat$forager_id[idx] ) )
        for ( i in idx ) {
            # specify as gaussian
            dat_list$age[i] <- 85/2
            #dat_list$age2[i] <- 0 + 0*sdunif(5,80)
            dat_list$age2[i] <- sdunif(5,80)/2
            dat_list$age_interval[i] <- 0
        }
    }

    invisible(dat_list)

}
# dat_list <- prep_data(dat)

# function to take prepped data and transform into trip records
# the Stan code does this internally
trans_trips <- function( d ) {

    N <- d$N
    N_hunters <- d$N_hunters
    N_trips <- d$N_trips

    n_foragers <- rep(NA,N_trips)
    n_assistants <- matrix(NA,nrow=N_trips,ncol=d$max_foragers_per_trip)
    n_dogs <- matrix(NA,nrow=N_trips,ncol=d$max_foragers_per_trip)
    n_firearms <- matrix(NA,nrow=N_trips,ncol=d$max_foragers_per_trip)
    trip_pooled <- rep(NA,N_trips)
    trip_soc_id <- rep(NA,N_trips)
    forager_ids <- matrix(NA,nrow=N_trips,ncol=d$max_foragers_per_trip)
    forager_age <- matrix(NA,nrow=N_trips,ncol=d$max_foragers_per_trip)
    trip_harvests <- matrix(NA,nrow=N_trips,ncol=d$max_foragers_per_trip)
    trip_hours <- rep(NA,N_trips)
    forager_soc_id <- rep(NA,N_hunters)
    trip_counted <- rep(NA,N_trips)
    hours_miss_idx <- rep(NA,N_trips)
    age_impute_idx <- rep(NA,N_hunters)
    age_impute_table <- matrix(0,nrow=N_hunters,ncol=3)
    
    # process age data
    # need to identify imputed ages
    # Exact ages: each unique forager gets an estimated age error
    # Interval ages: each forager gets an imputed age
    for ( i in 1:N_hunters ) age_impute_idx[i] <- (-1)
    N_ages_impute <- 0
    for ( i in 1:N ) {
        if ( age_impute_idx[d$hunter_id[i]] < 0 ) {
            # still need to process this forager
            if ( d$age_interval[i]==1 ) {
                # interval age, so always need imputation
                N_ages_impute <- N_ages_impute + 1
                age_impute_idx[d$hunter_id[i]] <- N_ages_impute;
                age_impute_table[d$hunter_id[i],1] <- 2;  # uniform
                age_impute_table[d$hunter_id[i],2] <- d$age[i]; # lower bound
                age_impute_table[d$hunter_id[i],3] <- d$age2[i]; # upper bound
            } else {
                # exact age, so maybe need imputation
                if ( d$age2[i]>0 ) {
                    N_ages_impute <- N_ages_impute + 1;
                    age_impute_idx[d$hunter_id[i]] <- N_ages_impute;
                    age_impute_table[d$hunter_id[i],1] <- 1;  # normal
                    age_impute_table[d$hunter_id[i],2] <- d$age[i]; # mu
                    age_impute_table[d$hunter_id[i],3] <- d$age2[i]; # sigma
                } else {
                    age_impute_idx[d$hunter_id[i]] <- 0;
                }
            }
        }
    }#i
    # loop over rows and store soc_id of each unique hunter_id
    for ( i in 1:N ) {
        forager_soc_id[d$hunter_id[i]] <- d$soc_id[i];
        # and check any assistants as well
        for ( a in 1:9 )
            if ( d$A[i,a]>0 ) forager_soc_id[d$A[i,a]] <- d$soc_id[i];
    }#i
    # some per-trip information
    N_hours_missing <- 0;
    for ( i in 1:N_trips ) {
        n_foragers[i] <- 0;
        hours_miss_idx[i] <- 0;
        for ( j in 1:d$max_foragers_per_trip ) {
            n_assistants[i,j] <- 0;
            n_dogs[i,j] <- 0;
            n_firearms[i,j] <- 0;
        }
        trip_counted[i] <- 0;
    }
    # loop over rows and process trip information
    for ( i in 1:N ) {
        j <- d$trip_id[i];
        trip_soc_id[j] <- d$soc_id[i];
        n_foragers[j] <- n_foragers[j] + 1;
        k <- n_foragers[j];
        forager_ids[j,k] <- d$hunter_id[i];
        forager_age[j,k] <- d$age[i];
        # if uniform age imputation, store forager age as midpoint + trip year offset
        if ( age_impute_table[d$hunter_id[i],1]==2 )
            forager_age[j,k] <- (d$age[i] + d$trip_year_offset[i]);
        trip_harvests[j,k] <- d$harvest[i];
        trip_pooled[j] <- d$pooled[i];
        trip_hours[j] <- d$hours[i];    # should already be standardized to mean 1
        if ( trip_hours[j] < 0 & hours_miss_idx[j]==0 ) {
            # missing trip duration
            N_hours_missing <- N_hours_missing + 1;
            # hours_miss_idx holds position of corresponding imputation parameter
            hours_miss_idx[j] <- N_hours_missing;
        }
        # count assistants, dogs, firearms
        if ( d$pooled[i]==0 ) {
            # unique assistants for each row i, which is forager k in this trip
            for ( a in 1:9 ) if ( d$A[i,a]>0 ) n_assistants[j,k] <- n_assistants[j,k] + 1;
            # dogs on row i?
            n_dogs[j,k] <- d$dogs[i];
            # firearms on row i?
            n_firearms[j,k] <- d$firearms[i];
        } else {
            # repeated assistants on all rows for this trip j
            # so be sure not to repeat count them
            # use trip_counted flag
            if ( trip_counted[j]==0 ) {
                for ( a in 1:9 )
                    if ( d$A[i,a]>0 ) n_assistants[j,1] <- n_assistants[j,1] + 1;
                trip_counted[j] <- 1;
            }
            n_dogs[j,1] <- d$dogs[i];
            n_firearms[j,1] <- d$firearms[i];
        }
    }#i over rows

    r <- list(
            n_foragers = n_foragers,
            n_assistants = n_assistants,
            n_dogs = n_dogs,
            n_firearms = n_firearms,
            trip_pooled = trip_pooled,
            trip_soc_id = trip_soc_id,
            forager_ids = forager_ids,
            forager_age = forager_age,
            trip_harvests = trip_harvests,
            trip_hours = trip_hours,
            forager_soc_id = forager_soc_id,
            trip_counted = trip_counted,
            hours_miss_idx = hours_miss_idx,
            age_impute_idx = age_impute_idx,
            age_impute_table = age_impute_table,
            N_hours_missing = N_hours_missing,
            N_ages_impute = N_ages_impute
        )

    return(r)

}

# td <- trans_trips( dat_list )
