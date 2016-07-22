# common utility functions for displaying model fits

shade_post <- function( x , y , col="black" , alpha=0.05 , probs=seq(0.1,0.99,by=0.025) ) {
    for ( p in probs ) {
        ci <- apply( x , 2 , PCI , prob=p )
        shade( ci , y , col=col.alpha(col,alpha) )
    }
}

plot.failure <- function( fit , hid=0 , n.ages=30 , xlab="age" , ylab="probability failure" , new=TRUE , ... ) {

    age.range <- range( fit@data$age )
    age.seq <- seq( from=age.range[1] , to=age.range[2] , length.out=n.ages )
    pred.dat <- list(
        harvest = rep( 0 , n.ages ),
        age = age.seq,
        hid = rep( hid , n.ages )
    )
    
    # check for flag to replace varying effects with zeros
    replace <- list()
    if ( hid==0 ) {
        pred.dat$hid = rep( 1 , n.ages )
        post <- extract.samples( fit )
        n <- length(post$a2)
        replace[['a2hid']] <- matrix( 0 , nrow=n , ncol=ncol(post$a2hid) )
    }
    
    # check for flag to marginalize over individual intercepts
    if ( hid==-1 ) {
        pred.dat$hid = rep( 1 , n.ages )
        post <- extract.samples( fit )
        n <- length(post$a2)
        replace[['a2hid']] <- matrix( rnorm( n , 0 , post$sigma[,2] ) , nrow=n , ncol=ncol(post$a2hid) )
    }
    
    pred <- link( fit , data=pred.dat , replace=replace )
    p <- apply( pred$p , 2 , mean )
    z <- ifelse( fit@data$harvest==0 , 1 , 0 )
    
    if ( new==TRUE ) blank()
    plot( fit@data$age , jitter(z,0.075) , col=col.alpha("black",0.35) , xlab=xlab , ylab=ylab , yaxp=c(0,1,4) )
    shade_post( pred$p , pred.dat$age , col="slateblue" )
    lines( age.seq , p , col="orange" , lwd=2 )
    
}

plot.harvest <- function( fit , hid=0 , n.ages=30 , xlab="age" , ylab="harvest (kg)" , new=TRUE , ... ) {

    age.range <- range( fit@data$age )
    age.seq <- seq( from=age.range[1] , to=age.range[2] , length.out=n.ages )
    pred.dat <- list(
        harvest = rep( 0 , n.ages ),
        age = age.seq,
        hid = rep( hid , n.ages )
    )
    
    # check for flag to replace varying effects with zeros
    replace <- list()
    if ( hid==0 ) {
        pred.dat$hid = rep( 1 , n.ages )
        post <- extract.samples( fit )
        n <- length(post$a1)
        replace[['a1hid']] <- matrix( 0 , nrow=n , ncol=ncol(post$a1hid) )
    }
    
    # check for flag to marginalize over individual intercepts
    if ( hid==-1 ) {
        pred.dat$hid = rep( 1 , n.ages )
        post <- extract.samples( fit )
        n <- length(post$a1)
        replace[['a1hid']] <- matrix( rnorm( n , 0 , post$sigma[,1] ) , nrow=n , ncol=ncol(post$a1hid) )
    }
    
    pred <- link( fit , data=pred.dat , replace=replace )
    mu <- apply( pred$mu , 2 , mean )
    z <- fit@data$harvest
    z2 <- z[z>0]
    
    if ( new==TRUE ) blank()
    plot( fit@data$age[z>0] , z2 , col=col.alpha("black",0.35) , xlab=xlab , ylab=ylab )
    shade_post( pred$mu , pred.dat$age , col="slateblue" )
    lines( age.seq , mu , col="orange" , lwd=2 )
    
}

cch_dens_intercept <- function( post , soc_id , new_plot=TRUE , lwd=2 , fcol="orange" , hcol="slateblue" , pars , ... ) {

    if ( missing(pars) ) {
        pars <- c("af")
    }

    if ( new_plot==TRUE ) blank(h=1.6,w=5.5)
    par(mfrow=c(2,10))

    j <- soc_id
    soc_name <- dat$society[dat$society_id==j][1]

    # failures
    dens( post$af[,1,j] , col=fcol , lwd=lwd , ... )
    mtext( concat("successes"," (",j,") ",soc_name) )
    dens( post$af[,2,j] , show.zero=TRUE , col=fcol , lwd=lwd , ... )
    mtext( concat("group"," (",j,")") )
    dens( post$af[,3,j] , show.zero=TRUE , col=fcol , lwd=lwd , ... )
    mtext( concat("pooled"," (",j,")") )
    dens( post$af[,4,j] , show.zero=TRUE , col=fcol , lwd=lwd , ... )
    mtext( concat("assist"," (",j,")") )
    dens( post$sef[,j] , show.zero=TRUE , col=fcol , lwd=lwd , ... )
    mtext( concat("f_elast"," (",j,")") )

    dens( post$b_hours[,1,j] , show.zero=TRUE , col=fcol , lwd=lwd , ... )
    mtext( concat("duration"," (",j,")") )
    dens( post$b_dogs[,1,j] , show.zero=TRUE , col=fcol , lwd=lwd , ... )
    mtext( concat("b_dogs"," (",j,")") )
    dens( post$se_dogs[,1,j] , show.zero=TRUE , col=fcol , lwd=lwd , ... )
    mtext( concat("se_dogs"," (",j,")") )
    dens( post$b_firearms[,1,j] , show.zero=TRUE , col=fcol , lwd=lwd , ... )
    mtext( concat("b_firearms"," (",j,")") )
    dens( post$se_firearms[,1,j] , show.zero=TRUE , col=fcol , lwd=lwd , ... )
    mtext( concat("se_firearms"," (",j,")") )

    # harvests
    dens( post$ah[,1,j] , col=hcol , lwd=lwd , ... )
    mtext( concat("harvests"," (",j,")") )
    dens( post$ah[,2,j] , show.zero=TRUE , col=hcol , lwd=lwd , ... )
    mtext( concat("group"," (",j,")") )
    dens( post$ah[,3,j] , show.zero=TRUE , col=hcol , lwd=lwd , ... )
    mtext( concat("pooled"," (",j,")") )
    dens( post$ah[,4,j] , show.zero=TRUE , col=hcol , lwd=lwd , ... )
    mtext( concat("assist"," (",j,")") )
    dens( post$seh[,j] , show.zero=TRUE , col=hcol , lwd=lwd , ... )
    mtext( concat("h_elast"," (",j,")") )

    dens( post$b_hours[,2,j] , show.zero=TRUE , col=hcol , lwd=lwd , ... )
    mtext( concat("duration"," (",j,")") )
    dens( post$b_dogs[,2,j] , show.zero=TRUE , col=hcol , lwd=lwd , ... )
    mtext( concat("b_dogs"," (",j,")") )
    dens( post$se_dogs[,2,j] , show.zero=TRUE , col=hcol , lwd=lwd , ... )
    mtext( concat("se_dogs"," (",j,")") )
    dens( post$b_firearms[,2,j] , show.zero=TRUE , col=hcol , lwd=lwd , ... )
    mtext( concat("b_firearms"," (",j,")") )
    dens( post$se_firearms[,2,j] , show.zero=TRUE , col=hcol , lwd=lwd , ... )
    mtext( concat("se_firearms"," (",j,")") )

}

# plot individual foraging functions across age
cch_plot_forager <- function( post , id , soc_id , data , n_points=100 , new_plot=TRUE , ylim , ... ) {
    x_seq <- seq(from=0,to=1,length.out=n_points)
    if ( missing(data) ) {
        data <- data.frame(age=x_seq,soc_id=soc_id,hunter_id=id)
    } else {
        data <- as.data.frame(data)
        data$soc_id <- soc_id
        data$hunter_id <- id
        if ( is.null(data$age) ) {
            data$age <- x_seq
        }
    }
    p <- cch_predict( post , data=data , ... )

    if ( new_plot==TRUE ) blank(w=2)
    par(mfrow=c(1,2))

    # failures
    pf_median <- apply(p[[1]],2,median)
    pf50 <- apply(p[[1]],2,PI,prob=0.5)
    pf80 <- apply(p[[1]],2,PI,prob=0.8)

    plot( NULL , xlim=c(0,1) , ylim=c(0,1) , xlab="age (standardized)" , ylab="probability failure" )
    lines( data$age , pf_median , lwd=2 )
    shade( pf50 , data$age )
    shade( pf80 , data$age )

    mtext( concat("Forager ",id,"(",soc_id,")") )

    # harvests
    ph_median <- apply(p[[2]],2,median)
    ph50 <- apply(p[[2]],2,PI,prob=0.5)
    ph80 <- apply(p[[2]],2,PI,prob=0.8)

    if ( missing(ylim) ) ylim <- range(ph80)
    plot( NULL , xlim=c(0,1) , ylim=ylim , xlab="age (standardized)" , ylab="harvest (standardized)" )
    lines( data$age , ph_median , lwd=2 )
    shade( ph50 , data$age )
    shade( ph80 , data$age )

}

# cch_plot_forager( post , 2 , 1 )
# cch_plot_forager( post , hid_in_soc(9)[3] , 9 , new_plot=FALSE )