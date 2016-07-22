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
skill <- function(x,k,m,b) exp(-m*x)*(1-exp(-k*x))^b

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
'Van_Vliet_et_al_Africa_three_sites',
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
'Ziker'
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

# function to take raw data and transform into trip records
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
