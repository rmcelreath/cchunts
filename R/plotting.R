# common utility functions for displaying model fits

shade_post <- function( x , y , col="black" , alpha=0.05 , probs=seq(0.1,0.99,by=0.025) ) {
    for ( p in probs ) {
        ci <- apply( x , 2 , PCI , prob=p )
        shade( ci , y , col=col.alpha(col,alpha) )
    }
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

# function to convert internal soc_id to appropriate map_id as in paper
cch_soc2map <- function(x) {
    data(cch_mapkey)
    return( cch_mapkey$Map.label[cch_mapkey$Site.number==x] )
}

cch_map2soc <- function(x) {
    data(cch_mapkey)
    return( cch_mapkey$Site.number[cch_mapkey$Map.label==x] )
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

cchpal <- c("#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20")

# function shows skill, success, harvest, expected returns for individual sites
cch_plot_curve <- function( map_id , soc_id , ids=NULL , x=seq(from=0,to=1,length.out=50) , h=1 , col=cchpal[4] , col2=col.alpha("black",0.2) , alpha , pts=FALSE , sample_post=FALSE , legend=FALSE , jitter=FALSE  , sw=0.1 , ybuffer=0.2 , skill_only=FALSE , ra=80 , part=3 , main , show_peak=TRUE , ... ) {

    if ( !missing(alpha) ) col=col.alpha(col,alpha)

    data(cch_mapkey)
    key <- cch_mapkey

    if ( !missing(map_id) ) j <- key$Site.number[key$Map.label==map_id]
    if ( !missing(soc_id) ) {
        j <- soc_id
        map_id <- key$Map.label[key$Site.number==j]
    }

    x_seq <- x

    i_soc <- which( dat_list$soc_id==j )
    age_min <- min( dat_list$age[i_soc] , na.rm=TRUE ) / ra
    age_max <- max( dat_list$age[i_soc] , na.rm=TRUE ) / ra
    x_seq_obs <- seq(from=age_min,to=age_max,length.out=40)

    # individual foragers in society j
    if ( is.null(ids) ) ids <- unique(dat_list$hunter_id[i_soc])
    if ( legend==TRUE ) ids <- sample(ids,1)
    n <- length(ids)
    if ( sample_post!=FALSE ) n <- sample_post
    p_med <- matrix(NA,nrow=n,ncol=50)
    for ( i in 1:n ) {
        h_use <- h
        if ( jitter==TRUE ) h_use <- h + runif(1,-0.5,2)
        if ( TRUE ) {
            if ( sample_post==FALSE )
                p <- cch_predict( post , data=data.frame(age=x_seq,soc_id=j,hours=h_use,hunter_id=ids[i]) , raw=skill_only )
            else
                p <- cch_predict( post , data=data.frame(age=x_seq,soc_id=j,hours=h_use) , , raw=skill_only )
            if ( skill_only==FALSE ) {
                if ( part==3 )
                    Ep <- (1-p$failure) * p$harvest
                else {
                    if ( part==1 )
                        Ep <- (1-p$failure)
                    if ( part==2 )
                        Ep <- p$harvest
                }
            } else
                Ep <- p$skill
            p_med[i,] <- apply(Ep,2,mean)
        } else {
            # skill curves raw
            p_med[i,] <- sapply( x_seq , function(x) exp(-m_est[ids[i]]*x)*(1-exp(-k_est[ids[i]]*x))^b_est[ids[i]] )
        }
        if (sample_post!=FALSE) {
            # instead of average for each forager, show samples from average forager
            p_med[i,] <- Ep[ sample(1:1000,1) ,  ]
        }
    }
    # now plot
    ymax <- max(p_med)
    ymax <- ymax*(1+ybuffer)
    if ( legend==TRUE ) ymax <- ymax*1.33
    if ( part==1 ) ymax <- 1
    plot( NULL , xlim=c(0,1) , xlab="" , ylim=c(0,ymax) , ylab="" , xaxt="n" , yaxt="n" , bty="n" )
    xat <- c(0,0.5,1)
    xlabs <- xat*ra
    axis( 1 , at=xat , labels=xlabs , cex=0.8 , padj=-0.1 )
    inobs <- which( x_seq >= age_min & x_seq <= age_max )
    lowobs <- which( x_seq < age_min )
    upobs <- which( x_seq > age_max )
    outcol <- col2
    for ( i in 1:n ) {
        lines( x_seq[inobs] , p_med[i,inobs] , col=col , ... )
        lines( x_seq[lowobs] , p_med[i,lowobs] , col=outcol , ... )
        lines( x_seq[upobs] , p_med[i,upobs] , col=outcol , ... )
    }#i

    # compute age at peak for average individual
    if ( show_peak==TRUE ) {
        peak <- sapply( 1:nrow(p_med) , function(i) x_seq[which.max(p_med[i,])] )
        mu <- mean(peak)
        abline(v=mu,lwd=0.5,lty=2)
        text(mu,0,round(ra*mu),cex=0.9,pos=2,offset=0.1,col=gray(0.4))
    }

    if ( pts==TRUE ) {
        for ( i in unique(round(dat_list$age[i_soc])) ) {
            if ( !is.na(i) ) {
                if ( part==3 )
                    obs <- dat_list$harvest[i_soc][ round(dat_list$age[i_soc])==i & dat_list$pooled[i_soc]==0 & dat_list$A[i_soc,1]==0 ]
                if ( part==2 )
                    obs <- dat_list$harvest[i_soc][ round(dat_list$age[i_soc])==i & dat_list$pooled[i_soc]==0 & dat_list$A[i_soc,1]==0 & dat_list$harvest[i_soc]>0 ]
                if ( part==1 )
                    obs <- (dat_list$harvest[i_soc]>0)[ round(dat_list$age[i_soc])==i & dat_list$pooled[i_soc]==0 & dat_list$A[i_soc,1]==0 ]
                p <- mean( obs , na.rm=TRUE )
                points( i/ra , p , cex=0.5 )
            }
        }#i
    }

    num_hunters <- length(unique(dat_list$hunter_id[i_soc]))
    num_h <- length( dat_list$harvest[i_soc] )

    if ( missing(main) )
        sn <- key$CODE[key$Site.number==j]
    else
        sn <- main
    if ( legend==TRUE ) {
        mtext( concat(0," ","--KEY--") , adj=0 )
        mtext( concat("ids (obs)") , adj=1 , cex=0.8 )
    } else {
        mtext( concat(map_id," ",sn) , adj=0 )
        mtext( concat(num_hunters," (",num_h,")") , adj=1 , cex=0.8 )
    }

    return(invisible(list(age_peak=mu,ymax=ymax)))
}

cch_plot_avgskill <- function( ymax=0.9 , x_seq=seq(from=0,to=1,length.out=60) , shading=TRUE , postdraws=FALSE , shade_levels=c(0.5,0.6,0.7,0.8,0.9,0.99) , shade_factor=0.8 ) {

    k <- exp( post$lifehistmeans[,1] )
    m <- exp( post$lifehistmeans[,2] )
    b <- exp( post$lifehistmeans[,3] )
    skill_at_x <- sapply( x_seq , function(x) skill(x,k,m,b) )

    plot( x_seq , apply(skill_at_x,2,mean) , xlab="" , xaxt="n" , ylab="" , lwd=1 , col="black" , yaxt="n" , ylim=c(0,ymax) , bty="n" , type="l" )
    xat <- c(0,18/ra,31/ra,0.5,55/ra,1)
    xlabs <- xat*ra
    axis( 1 , at=xat , labels=xlabs , cex=0.8 , padj=-0.1 )

    if ( postdraws==TRUE ) {
        # uncertainty as draws from posterior
        for ( i in 1:50 ) {
            j <- sample( 1:2000 , 1 )
            curve( skill(x,k[j],m[j],b[j]) , add=TRUE , lwd=1 , col=col.alpha(pal[4],0.3) )
        }
        curve( skill(x,mean(k),mean(m),mean(b)) , add=TRUE , lwd=2 )
    }

    if ( shading==TRUE ) {
        # uncertainty as shaded region
        for ( sl in c(0.5,0.6,0.7,0.8,0.9,0.99) ) {
            shcol <- col.alpha(pal[3], (1-sl)*shade_factor )
            ci <- apply(skill_at_x,2,PI,prob=sl)
            shade( ci , x_seq , col=shcol )
        }
        lines( x_seq , apply(skill_at_x,2,mean) , lwd=2 )
    }
    

    # lines for illustration

    do_line_to <- function(x) {
        y <- x/ra
        sk <- skill(y,mean(k),mean(m),mean(b))
        lines( c(y,y) , c(0,sk) , lty=2 , lwd=0.5 )
    }

    do_line_to(18)
    do_line_to(31)
    do_line_to(55)
    sk18 <- mean(skill(18/ra,(k),(m),(b)))
    sk55 <- mean(skill(55/ra,(k),(m),(b)))
    lines( c(18/ra,55/ra) , c(sk18,sk55) , lty=2 )

}

cch_plot_grid <- function( map_id=1:40 , nrow=7 , ncol=6 , col=cchpal[4] , col2=col.alpha("black",0.5) , fskillonly=TRUE , dosample=FALSE , draw_globalmean=TRUE , draw_legend=TRUE , alpha=0.5 , show_points=FALSE , lwd=1.5 , adj_margins=TRUE , part=3 , skip=0 ) {

    if ( adj_margins==TRUE )
        par(mgp = c(1.5, 0.2, 0), mar = c(1.2, 0.5, 1.5, 0.25) + 0.1, tck = -0.02)

    par(mfrow=c(nrow,ncol))

    fcol <- col

    if ( skip > 0 ) for ( z in 1:skip ) plot.new()

    if ( draw_globalmean==TRUE ) {
        cch_plot_avgskill()
        mtext( concat("Global mean") , adj=0 )
    }

    # legend
    if ( draw_legend==TRUE ) {
        xp <- cch_plot_curve( map_id=8 , col=fcol , alpha=1 , pts=FALSE , lwd=lwd , legend=TRUE , col2=col2 , skill_only=fskillonly , part=part )
        text( xp$age_peak , xp$ymax/4 , "peak" , cex=1.2 , col=gray(0.2) , pos=4 , offset=0.1 )
        text( xp$age_peak+0.2 , xp$ymax/1.33 , "obs range" , col=fcol , cex=1.2 , pos=4 , offset=0.1 )
    }

    for ( i in map_id ) cch_plot_curve( map_id=i , col=fcol , alpha=alpha , pts=show_points , lwd=lwd , jitter=FALSE , sample_post=dosample , skill_only=fskillonly , part=part )
}


# cch_plot_forager( post , 2 , 1 )
# cch_plot_forager( post , hid_in_soc(9)[3] , 9 , new_plot=FALSE )