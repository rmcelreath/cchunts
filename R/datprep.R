# converts csv 
datprepper <- function( file , objectname , save=TRUE , path="data_rda/" , redact=TRUE , debug=FALSE ) {

    # file e.g.: "data/Alvard.csv"
    # objectname e.g.: "Alvard" --- is what data() will load it as

    require(chron)
    require(rethinking)
    yn2bin <- function(x) ifelse(is.na(x),NA,ifelse(x %in% c("Yes","yes","Y","y"),1,0))

    d <- read.csv( file , na.strings = c("","NA") , stringsAsFactors=FALSE )

    # target columns
    if ( redact==FALSE )
    resultlist <- c(
        "trip_id", "trip_id_orig", "observed", "trip_date", "julian_date_s", "day_id", 
        "trip_duration", "day_trip", "group_type", "pooled", "harvest", 
        "forager_id", "forager_id_orig", "age_type" , "age_dist_1", "age_dist_2", "sex",
        "dogs", "gun", "trip_year"
    )
    else
    resultlist <- c(
        "trip_id", "trip_id_orig", "observed", "trip_date", "julian_date_s", "day_id", 
        "trip_duration", "day_trip", "group_type", "pooled", "harvest", 
        "forager_id", "age_type" , "age_dist_1", "age_dist_2", "sex",
        "dogs", "gun", "trip_year"
    )

    # make unique trip ID index
    d$trip_id <- coerce_index(d$Trip.ID)
    d$trip_id_orig <- d$Trip.ID

    # make binary variable to denote if trip was observed
    d$observed <- yn2bin(d$Observed)

    # Julian date, standardized
    if ( !is.null(d$Hunt.date) ) {
        d$trip_date <- as.Date( d$Hunt.date , format="%m/%d/%y" ) # short year format mm/dd/yy
    }
    if ( !is.null(d$Hunt.Date) )
        d$trip_date <- as.Date( d$Hunt.Date , format="%m/%d/%Y" ) # long year format mm/dd/yyyy
    d$julian_date_s <- scale( julian( d$trip_date ) )

    # trip year
    if ( !is.null(d$Hunt.year) ) {
        d$trip_year <- as.integer( d$Hunt.year )
    } else {
        if ( !is.null(d$year) ) {
            d$trip_year <- as.integer( d$year )
        } else {
            warning(concat("No Hunt.year for ",file))
        }
    }

    # unique day indicator, for clustering by specific date
    d$day_id <- coerce_index( d$trip_date )

    # Retain information on group type and harvest indicator
    d$group_type <- d$Group.type
    d$pooled <- yn2bin(d$Pooled)

    # Rename harvest outcomes, the response variable
    d$harvest <- d$Harvest

    # check for assistant columns and merge them into forager IDs when present
    num_assist <- 0
    fstop <- FALSE
    arg_list <- list(d$Hunter)
    while ( !fstop ) {
        colname <- concat( "Assistant.",num_assist+1 )
        if ( is.null(d[[colname]]) ) {
            fstop <- TRUE
        } else {
            num_assist <- num_assist + 1
            arg_list[[length(arg_list)+1]] <- d[[colname]]
        }
    }
    if ( num_assist > 0 ) {
        hidx <- do.call( coerce_index , arg_list )
        d$forager_id <- hidx[[1]]
        for ( j in 1:num_assist ) {
            varname <- concat( "a_",j,"_id" )
            d[[varname]] <- hidx[[j+1]]
            resultlist <- c(resultlist,varname)
        }
        total_foragers <- max(unlist(hidx),na.rm=TRUE)
    } else {
        d$forager_id <- coerce_index( d$Hunter )
        total_foragers <- max(d$forager_id)
    }

    # only keep original hunter name when redact argument FALSE
    if ( redact==FALSE )d$forager_id_orig <- d$Hunter

    #### forager sex
    if ( !is.null(d$Hunter.Sex) )
        d$sex <- d$Hunter.Sex
    else
        d$sex <- "M"

    ########################
    #### Age
    ########################
    d$age_type <- as.character(d$Age.type)
    d$age_dist_1 <- d$Age_dist_1
    d$age_dist_2 <- d$Age_dist_2

    # Denotes presence of gun on trip
    d$gun <- yn2bin(d$Firearms)

    # Hunting time in hours
    d$trip_duration <- d$Hours
    if ( class(d$trip_duration)=="character" ) {
        # Coad coding error? 5:30 e.g.
        # get user input
        ids <- grep(":",d$trip_duration)
        if ( length(ids)>0 ) {
            for ( i in ids )
                d$trip_duration[i] <- readline(concat("Duration '",d$trip_duration[i],"' found. Enter correct value:"))
            d$trip_duration <- as.numeric(d$trip_duration)
        }
    }

    # Denotes whether hunters spend an overnight away
    # from primary residence or camp.
    if ( is.numeric(d$Single.day) )
        d$day_trip <- (d$Single.day)
    else
        d$day_trip <- yn2bin(d$Single.day)

    # dogs variable to denote presence of dogs
    d$dogs <- yn2bin(d$Dogs)

    # process assistant covariates
    covlist <- c("sex","Age_dist_1","Age_dist_2","age_dist_1","age_dist_2","firearms")
    ynlist <- c(0,0,0,0,0,1) # whether to convert to binary indicator
    if ( num_assist > 0 ) {
        for ( j in 1:num_assist ) {
            inprefix <- concat("Assistant.",j,".")
            outprefix <- concat("a_",j,"_")
            for ( k in 1:length(covlist) ) {
                inname <- concat(inprefix,covlist[k])
                if ( !is.null(d[[inname]]) ) {
                    result <- d[[inname]]
                    outname <- concat( outprefix , tolower(covlist[k]) )
                    if ( ynlist[k]==1 ) result <- yn2bin(result)
                    d[[outname]] <- result
                    # add to resultlist
                    resultlist <- c(resultlist,outname)
                }# not is null
            }#k
        }#j
    }

    if ( debug==TRUE ) {
        print(names(d))
        print(resultlist)
    }

    # pare down
    dsave <- d[ , resultlist ]

    # run some checks
    # all pooled trips have more than one harvest (row)?
    x <- d[ d$pooled==1 , ]
    for ( i in unique(x$trip_id) ) {
        n_foragers <- nrow(x[x$trip_id==i,])
        if ( n_foragers<2 ) message(concat("Trip ",i," is pooled but has only ",n_foragers," harvest."))
    }

    # save as binary object for loading in package
    assign( objectname , dsave )

    if ( save==TRUE ) {
        destpath <- path
        destpath <- concat(destpath,objectname,".rda")
        save( list=objectname , file=destpath )
        message(concat("Result saved in binary at ",destpath))
    }

    message(concat("**",file,"**"))
    message(concat("Trips:",max(d$trip_id)))
    message(concat("Foragers:",total_foragers))

    return(invisible(dsave))

}

if ( FALSE ) {

# do <- read.csv("data/Hill_Kintigh.csv", na.strings = c("","NA"), stringsAsFactors=FALSE )

doSave <- FALSE
doDebug <- TRUE

d <- datprepper( "data/Alvard.csv" , "Alvard" , save=doSave , debug=doDebug )
d <- datprepper( "data/Beckerman.csv" , "Beckerman" , save=doSave , debug=doDebug )
d <- datprepper( "data/Bird_Bird_Codding.csv" , "Bird_Bird_Codding" , save=doSave , debug=doDebug )
d <- datprepper( "data/Coad.csv" , "Coad" , save=doSave , debug=doDebug )
d <- datprepper( "data/Ellen.csv" , "Ellen" , save=doSave , debug=doDebug )
d <- datprepper( "data/Franzen.csv" , "Franzen" , save=doSave , debug=doDebug )
d <- datprepper( "data/Healey_Nen_PNG.csv" , "Healey_Nen_PNG" , save=doSave , debug=doDebug )
d <- datprepper( "data/Healey.csv" , "Healey" , save=doSave , debug=doDebug )
d <- datprepper( "data/Koster.csv" , "Koster" , save=doSave , debug=doDebug )
d <- datprepper( "data/Kramer_Greaves.csv" , "Kramer_Greaves" , save=doSave , debug=doDebug )
d <- datprepper( "data/Lupo_Schmitt.csv" , "Lupo_Schmitt" , save=doSave , debug=doDebug )
d <- datprepper( "data/Marks.csv" , "Marks" , save=doSave , debug=doDebug )
d <- datprepper( "data/Nielsen.csv" , "Nielsen" , save=doSave , debug=doDebug )
d <- datprepper( "data/Pacheco.csv" , "Pacheco" , save=doSave , debug=doDebug )
d <- datprepper( "data/Pangau_Adam.csv" , "Pangau_Adam" , save=doSave , debug=doDebug )
d <- datprepper( "data/Ready.csv" , "Ready" , save=doSave , debug=doDebug )
d <- datprepper( "data/Sillitoe.csv" , "Sillitoe" , save=doSave , debug=doDebug )
d <- datprepper( "data/Siren.csv" , "Siren" , save=doSave , debug=doDebug )
d <- datprepper( "data/Trumble_Gurven.csv" , "Trumble_Gurven" , save=doSave , debug=doDebug )
d <- datprepper( "data/Van_Vliet_et_al_Africa_three_sites.csv" , "Van_Vliet_et_al_Africa_three_sites" , save=doSave , debug=doDebug )
d <- datprepper( "data/Van_Vliet_et_al_Gabon.csv" , "Van_Vliet_et_al_Gabon" , save=doSave , debug=doDebug )
d <- datprepper( "data/Van_Vliet_et_al_South_America_sites.csv" , "Van_Vliet_et_al_South_America_sites" , save=doSave , debug=doDebug )
d <- datprepper( "data/Yu_et_al.csv" , "Yu_et_al" , save=doSave , debug=doDebug )
d <- datprepper( "data/Ziker.csv" , "Ziker" , save=doSave , debug=doDebug )
d <- datprepper( "data/Ross.csv" , "Ross" , save=doSave , debug=doDebug )

d <- datprepper( "data/Hill_Kintigh.csv" , "Hill_Kintigh" , save=doSave , debug=doDebug )

########
# Automated script to process all csv files in data/ directory

doSave <- TRUE
doDebug <- FALSE
files <- dir("data",pattern=".csv")
for ( f in files ) {
    filein <- concat("data/",f)
    nameout <- substr(f,1,nchar(f)-4)
    datprepper( filein , nameout , save=doSave , debug=doDebug )
}

}
