cchmodel_full <- "
data{
    int N;
    int N_hunters;
    int N_societies;
    int age_interval[N];    // 0/1 indicator for interval age data
    real age[N];            // observed age or (when interval) lower bound
    real age2[N];           // stderr of age or (when interval) upper bound
    real trip_year_offset[N];
    real hours[N];
    int day_trip[N];        // 0/1 indicator of whether trip exceeded a single day
    real harvest[N];
    int hunter_id[N];
    int forager_female[N];
    int soc_id[N];
    real ref_age;
    int N_trips;
    int trip_id[N];
    int pooled[N];
    int A[N,9];             // assistant id matrix. up to 9 assistants per harvest.
    int dogs[N];
    int firearms[N];
    int max_foragers_per_trip;
    real prior_scale;       // input scale for regularizing priors
}
transformed data{
    // loop over rows and compute trip information
    real RHO_PRIOR_REGGER;
    int n_foragers[N_trips];
    int n_assistants[N_trips,max_foragers_per_trip];
    int n_dogs[N_trips,max_foragers_per_trip];
    int n_firearms[N_trips,max_foragers_per_trip];
    int trip_pooled[N_trips];
    int trip_soc_id[N_trips];
    int forager_ids[N_trips,max_foragers_per_trip];
    vector[max_foragers_per_trip] forager_age[N_trips];
    vector[max_foragers_per_trip] trip_forager_female[N_trips];
    vector[max_foragers_per_trip] trip_harvests[N_trips];
    vector[2] zeros2;
    vector[3] zeros3;
    vector[N_trips] trip_hours;
    vector[N_trips] trip_xday;          // trip exceeded a single day (camping)
    int forager_soc_id[N_hunters];      // useful for precomputing individual skill parameters
    int trip_counted[N_trips];
    // duration imputation variables
    int N_hours_missing;
    int hours_miss_idx[N_trips];
    // age error variables
    int N_ages_impute;
    int age_impute_idx[N_hunters];
    real age_impute_table[N_hunters,3];
    // firearm imputation variables
    int N_firearm_missing;
    int firearm_miss_idx[N];    // possibly for every observed harvest

    RHO_PRIOR_REGGER = 40;

    // process firearm missingness
    N_firearm_missing = 0;
    for ( i in 1:N ) {
        if ( firearms[i] < 0 ) {
            N_firearm_missing = N_firearm_missing + 1;
            firearm_miss_idx[i] = N_firearm_missing;
        } else {
            firearm_miss_idx[i] = 0;
        }
    }
    // process age data
    // need to identify imputed ages
    // Exact ages: each unique forager gets an estimated age error
    // Interval ages: each forager gets an imputed age
    for ( i in 1:N_hunters ) age_impute_idx[i] = (-1); // just flag to process this forager
    N_ages_impute = 0;
    for ( i in 1:N ) {
        if ( age_impute_idx[hunter_id[i]]<0 ) {
            // still need to process this forager
            if ( age_interval[i]==1 ) {
                // interval age, so always need imputation
                N_ages_impute = N_ages_impute + 1;
                age_impute_idx[hunter_id[i]] = N_ages_impute;
                age_impute_table[hunter_id[i],1] = 2;  // uniform
                age_impute_table[hunter_id[i],2] = age[i]; // lower bound
                age_impute_table[hunter_id[i],3] = age2[i]; // upper bound
            } else {
                // exact age, so maybe need imputation
                if ( age2[i]>0 ) {
                    N_ages_impute = N_ages_impute + 1;
                    age_impute_idx[hunter_id[i]] = N_ages_impute;
                    age_impute_table[hunter_id[i],1] = 1;  // normal
                    age_impute_table[hunter_id[i],2] = age[i]; // mu
                    age_impute_table[hunter_id[i],3] = age2[i]; // sigma
                } else {
                    age_impute_idx[hunter_id[i]] = 0;
                }
            }
        }
    }
    // loop over rows and store soc_id of each unique hunter_id
    for ( i in 1:N ) {
        forager_soc_id[hunter_id[i]] = soc_id[i];
        // and check any assistants as well
        for ( a in 1:9 )
            if ( A[i,a]>0 ) forager_soc_id[A[i,a]] = soc_id[i];
    }
    // some per-trip information
    N_hours_missing = 0;
    for ( i in 1:2 ) zeros2[i] = 0;
    for ( i in 1:3 ) zeros3[i] = 0;
    for ( i in 1:N_trips ) {
        n_foragers[i] = 0;
        hours_miss_idx[i] = 0;
        for ( j in 1:max_foragers_per_trip ) {
            n_assistants[i,j] = 0;
            n_dogs[i,j] = 0;
            n_firearms[i,j] = 0;
        }
        trip_counted[i] = 0;
        trip_xday[i] = 0;
    }
    // loop over rows and process trip information
    for ( i in 1:N ) {
        int j;
        int k;
        real ai;
        j = trip_id[i];
        trip_soc_id[j] = soc_id[i];
        n_foragers[j] = n_foragers[j] + 1;
        k = n_foragers[j];
        forager_ids[j,k] = hunter_id[i];
        forager_age[j,k] = age[i];
        trip_forager_female[j,k] = forager_female[i];
        // if uniform age imputation, store forager age as midpoint + trip year offset
        if ( age_impute_table[hunter_id[i],1]==2 )
            forager_age[j,k] = (age[i] + trip_year_offset[i]);
        trip_harvests[j,k] = harvest[i];
        trip_pooled[j] = pooled[i];
        trip_hours[j] = hours[i];    // should already be standardized to mean 1
        trip_xday[j] = ( 1 - day_trip[i] );
        if ( trip_hours[j] < 0 && hours_miss_idx[j]==0 ) {
            // missing trip duration
            N_hours_missing = N_hours_missing + 1;
            // hours_miss_idx holds position of corresponding imputation parameter
            hours_miss_idx[j] = N_hours_missing;
        }
        // count assistants, dogs, firearms
        if ( pooled[i]==0 ) {
            // unique assistants for each row i, which is forager k in this trip
            for ( a in 1:9 ) if ( A[i,a]>0 ) n_assistants[j,k] = n_assistants[j,k] + 1;
            // dogs on row i?
            n_dogs[j,k] = dogs[i];
            // firearms on row i?
            n_firearms[j,k] = firearms[i];
        } else {
            // repeated assistants on all rows for this trip j
            // so be sure not to repeat count them
            // use trip_counted flag
            if ( trip_counted[j]==0 ) {
                for ( a in 1:9 )
                    if ( A[i,a]>0 ) n_assistants[j,1] = n_assistants[j,1] + 1;
                trip_counted[j] = 1;
            }
            n_dogs[j,1] = dogs[i];
            n_firearms[j,1] = firearms[i];
        }
    }//i over rows
    if ( N_hours_missing==0 ) N_hours_missing = 1;
}
parameters{
    // forager life history parameters
    vector[6] lifehistmeans;                // averages across societies (k,m,b,rho_km,sk,sm)
    // parameters in vector are on latent scales:
    // (1) log k - average log(k) across sites
    // (2) log m - average log(m) across sites
    // (3) log b - average log(b) across sites
    // next three parameters are second level variation:
    //   variation among individuals WITHIN sites
    //   these parameters are the global averages of WITHIN site distributions
    //   note that log(b) is absent because b doesn't vary across individuals within sites
    // (4) logit(rho_km) - average correlation btw k,m among individuals across sites
    // (5) log sk - average log stddev of k among individuals across sites
    // (6) log sm - average log stddev of m among individuals across sites
    vector<lower=0>[6] sigma_societies;       // stddev of varying effects
    // (1) stddev of log k across sites
    // (2) stddev of log m across sites
    // (3) stddev of log b across sites
    // (4) stddev of latent-scale rho_km
    //     so variation across sites in correlation btw km within sites 
    // (5) stddev across sites of variation in k within sites
    // (6) stddev across sites of variation in m within sites
    matrix[6,N_societies] zs;                 // z-scores for each society
    // indices here are same as in lifehistmeans parent vector
    // but each site gets its own vector drawn from multi_normal defined by lifehistmeans and sigma_societies
    // correlation matrix fixed to identity matrix

    // intercept components of outcomes
    real afbar;                 // average af[1,] (for pooling)
    real ahbar;                 // average ah[1,]
    real<lower=0> sigma_af;     // stddev af[1,] (for pooling)
    real<lower=0> sigma_ah;     // stddev ah[1,] (for pooling)
    real af[4,N_societies];     // intercepts, failures
    real ah[4,N_societies];     // intercepts, harvests
                                // [1] main failure intercept
                                // [2] marginal group trips
                                // [3] marginal pooled trips
                                // [4] marginal assistant effects
    real sef[N_societies]; // skill elasticity failures - log scale
    real seh[N_societies]; // skull elasticity harvests - log scale

    // harvest scale parameters (gamma distribution scale)
    real<lower=0> hscale[N_societies];

    // hunter effects
    matrix[2,N_hunters] zh;         // z-scores for each hunter
    // [1] k, [2] m
    // variance-covariance matrix for these effects comes from vs vector
    // see construction in translated parameters block

    // various regression parameters
    // row 1: failures
    // row 2: harvests
    real b_hours[2,N_societies];
    real b_dogs[2,N_societies];
    real b_firearms[2,N_societies];
    real se_dogs[2,N_societies];
    real se_firearms[2,N_societies];
    real b_xday[N_societies,2];     // trip exceeded a day (not day_trip==1)

    // imputation
    // trip durations
    real trip_duration_imputed[N_hours_missing];
    real trip_duration_mu[N_societies];         // mean and stddev of prior
    real<lower=0> trip_duration_sigma[N_societies];
    // ages
    real age_err[N_ages_impute];

    // hack to take Ache out of correlation rho_km pooling
    real ache_fix_rho;

    vector<lower=0,upper=1>[N_societies] dogs_mu; // prob trip has dogs for each society
    vector<lower=0,upper=1>[N_societies] guns_mu; // prob trip has firearms for each society
}
transformed parameters{
    // code to transform non-centered prior to centered effects
    matrix[N_hunters,2] vh;
    matrix[N_societies,6] vs;
    matrix[N_societies,2] sigmas_hunters;
    vector[N_societies] rho_km;
    //vs = (diag_pre_multiply(sigma_societies,L_Rho_societies) * zs)';
    for ( s in 1:N_societies )
        for ( i in 1:6 )
            vs[s,i] = sigma_societies[i] * zs[i,s];
    // loop over individuals and sample individual hunter k,m parameters
    // need to build different covariance matrix for each society,
    // as implied by society specific parameters (as in vs matrix)
    // cholesky decomp implies:
    // x1 = z1 * sigma1
    // x2 = (rho*z1 + sqrt(1-rho^2)*z2)*sigma2
    for ( s in 1:N_societies ) {
        sigmas_hunters[s,1] = exp( lifehistmeans[4] + vs[s,4] ); // sigma on k
        sigmas_hunters[s,2] = exp( lifehistmeans[5] + vs[s,5] ); // sigma on m
        rho_km[s] = 2*inv_logit( lifehistmeans[6] + vs[s,6] ) - 1; // rho for k,m ; 0 on latent scale is 0 
        if ( s==14 ) rho_km[s] = 2*inv_logit( lifehistmeans[6] + ache_fix_rho ) - 1;
    }
    for ( i in 1:N_hunters ) {
        // currently very inefficient, because repeat calcs
        // should chunk by society
        int s;
        s = forager_soc_id[i];
        vh[i,1] = zh[1,i] * sigmas_hunters[s,1];
        vh[i,2] = (rho_km[s]*zh[1,i] + sqrt(1-rho_km[s]^2)*zh[2,i])*sigmas_hunters[s,2];
    }
}
model{
    // temp variables
    real k[N_hunters];
    real m[N_hunters];
    real b[N_societies];
    vector[N_trips] lm_f;
    vector[N_trips] lm_h;
    real p;
    real mu;
    matrix[2,2] Sigma;
    vector[N_trips] trip_duration_merge;
    //real prior_scale;

    //prior_scale = 0.5;
    
    // priors
    // society-level life history means --- centered on global means
    // equivalent to:
    //vs ~ multi_normal( lifehistmeans , quad_form_diag(Rho_societies,sigma_societies) );
    to_vector(zs) ~ normal(0,1);

    lifehistmeans[1:2] ~ normal( 1, prior_scale ); // log k,m
    lifehistmeans[3] ~ normal( 1, prior_scale ); // log b
    lifehistmeans[6] ~ normal( 0, prior_scale ); // shifted logit rho_km
    // do prior for stddev k,m between [4,5] as normal on transformed scale
    // this allows us to define same prior for sigma_societies[1:2]
    exp(lifehistmeans[4]) ~ normal( 0 , prior_scale );
    exp(lifehistmeans[5]) ~ normal( 0 , prior_scale );
    // need Jacobian adjustments for these priors
    // log|d/dy exp y| = log|exp y| = y
    // see also section 33.2 of Stan reference manual
    target += lifehistmeans[4];
    target += lifehistmeans[5];

    sigma_societies ~ normal( 0 , prior_scale );

    dogs_mu ~ beta(2,4); // weighted towards low values to stop mode switching in site 8
    guns_mu ~ beta(2,4);

    ache_fix_rho ~ normal(0, prior_scale );

    afbar ~ normal(0, prior_scale );
    ahbar ~ normal(0, prior_scale );
    sigma_af ~ normal(0, prior_scale );
    sigma_ah ~ normal(0, prior_scale );

    for ( s in 1:N_societies ) {
        af[1,s] ~ normal(afbar,sigma_af);
        ah[1,s] ~ normal(ahbar,sigma_ah);
        for ( i in 2:4 ) {
            af[i,s] ~ normal(0,prior_scale);
            ah[i,s] ~ normal(0,prior_scale);
        }
        sef[s] ~ normal(0,prior_scale);
        seh[s] ~ normal(0,prior_scale);
        for ( i in 1:2 ) {
            b_hours[i,s] ~ normal(0,prior_scale);
            b_dogs[i,s] ~ normal(0,prior_scale);
            b_firearms[i,s] ~ normal(0,prior_scale);
            se_dogs[i,s] ~ normal(0,prior_scale);
            se_firearms[i,s] ~ normal(0,prior_scale);

            b_xday[s,i] ~ normal(0,prior_scale);
        }
    }//s
    hscale ~ normal( 1 , prior_scale );

    // varying effects
    // foragers --- these are zero-centered
    // see translation to vh in transformed parameters block
    to_vector(zh) ~ normal(0,1);

    // age imputation
    for ( i in 1:N_hunters ) {
        if ( age_impute_idx[i] > 0 ) {
            if ( age_impute_table[i,1]==1 )
                age_err[age_impute_idx[i]] ~ 
                        normal( 0 , age_impute_table[i,3] );
            //else
                //age_err[age_impute_idx[i]] ~ 
                //        uniform( 0 , age_impute_table[i,3]-age_impute_table[i,2] );
        }
    }

    // trip durations
    for ( j in 1:N_societies ) trip_duration_mu[j] ~ normal(0,1);
    trip_duration_sigma ~ exponential(1);
    for ( i in 1:N_trips ) {
        if ( trip_hours[i]<0 ) {
            // missing
            trip_duration_merge[i] = trip_duration_imputed[hours_miss_idx[i]];
        } else {
            // observed
            trip_duration_merge[i] = log(trip_hours[i]);
        }
        // prior (when missing) or likelihood (when observed)
        trip_duration_merge[i] ~ normal( trip_duration_mu[trip_soc_id[i]] , 
                                         trip_duration_sigma[trip_soc_id[i]] );
    }//i
    
    // prep hunter effects so can re-use
    for ( j in 1:N_hunters ) {
        if ( forager_soc_id[j]==14 ) {
            k[j] = exp( lifehistmeans[1] + vs[forager_soc_id[j],1] + vh[j,1] );
            m[j] = exp( lifehistmeans[2] + vs[forager_soc_id[j],2] + vh[j,2] );
        } else {
            k[j] = exp( lifehistmeans[1] + vs[forager_soc_id[j],1] + vh[j,1] );
            m[j] = exp( lifehistmeans[2] + vs[forager_soc_id[j],2] + vh[j,2] );
        }
    }
    // prep b for each society, so only have to compute once
    for ( s in 1:N_societies ) {
        b[s] = exp( lifehistmeans[3] + vs[s,3] );     // ensure positive with log link
    }
    
    // likelihoods
    lm_f = rep_vector(0,N_trips);
    lm_h = lm_f;
    // loop over trips and compute likelihoods
    for ( i in 1:N_trips ) {
        real skillj;
        real sefx;
        real sehx;
        real ai;
        int hid;
        real avg_skill;
        vector[2] LLterms;
        vector[4] LL4terms;
        int xdogs;
        int xguns;
        int n_foragers_index;
        int coopidx;
        int xdogsvec[4];
        int xgunsvec[4];
        xdogsvec[1] = 1;
        xdogsvec[2] = 1;
        xdogsvec[3] = 0;
        xdogsvec[4] = 0;
        xgunsvec[1] = 1;
        xgunsvec[2] = 0;
        xgunsvec[3] = 1;
        xgunsvec[4] = 0;
        avg_skill = 0;
        if ( trip_pooled[i]==1 ) {
            // pooled harvest
            // compute average skill in foraging group
            for ( j in 1:n_foragers[i] ) {
                // ai = forager_age[i,j];

                hid = forager_ids[i,j];
                if ( age_impute_idx[hid]==0 ) {
                    // simple case, just fetch observed age
                    ai = forager_age[i,j]; // from trip variables
                } else {
                    // need some kind of imputation
                    ai = forager_age[i,j] + age_err[age_impute_idx[hid]];
                }
                ai = ai/ref_age;

                skillj = exp(-m[hid]*ai)*pow(1-exp(-k[hid]*ai),b[trip_soc_id[i]]);
                avg_skill = avg_skill + skillj;
            }//j
            avg_skill = avg_skill/n_foragers[i] + 0.001;
            n_foragers_index = 1; // loop over just one forager
            coopidx = 3;
        } else {
            // independent harvests
            n_foragers_index = n_foragers[i];
            coopidx = 2;
        }

        for ( j in 1:n_foragers_index ) {

            if ( trip_pooled[i]==1 ) {
                skillj = avg_skill;
            } else {
                hid = forager_ids[i,j];
                if ( age_impute_idx[hid]==0 ) {
                    // simple case, just fetch observed age
                    ai = forager_age[i,j]; // from trip variables
                } else {
                    // need some kind of imputation
                    ai = forager_age[i,j] + age_err[age_impute_idx[hid]];
                }
                ai = ai/ref_age;
                skillj = exp(-m[hid]*ai)*pow(1-exp(-k[hid]*ai),b[trip_soc_id[i]]) + 0.001;
            }

            // prep common parts of linear models
            // this just removes the dogs and firearms terms
            // those get added below, depending upon missingness
            lm_f[i] = exp( af[1,trip_soc_id[i]] + 
                                af[coopidx,trip_soc_id[i]]*(n_foragers[i]-1) + 
                                af[4,trip_soc_id[i]]*n_assistants[i,1] + 
                                //b_dogs[1,trip_soc_id[i]]*n_dogs[i,1] +
                                //b_firearms[1,trip_soc_id[i]]*n_firearms[i,1] +
                                b_xday[trip_soc_id[i],1]*trip_xday[i]
                            ) *
                            exp(trip_duration_merge[i])^b_hours[1,trip_soc_id[i]];
            lm_h[i] = exp( ah[1,trip_soc_id[i]] + 
                                ah[coopidx,trip_soc_id[i]]*(n_foragers[i]-1) + 
                                ah[4,trip_soc_id[i]]*n_assistants[i,1] +
                                //b_dogs[2,trip_soc_id[i]]*n_dogs[i,1] +
                                //b_firearms[2,trip_soc_id[i]]*n_firearms[i,1] +
                                b_xday[trip_soc_id[i],2]*trip_xday[i]
                            )*
                            exp(trip_duration_merge[i])^b_hours[2,trip_soc_id[i]];
            sefx = exp( sef[trip_soc_id[i]]  
                                //se_dogs[1,trip_soc_id[i]]*n_dogs[i,1] +
                                //se_firearms[1,trip_soc_id[i]]*n_firearms[i,1]
                            );
            sehx = exp( seh[trip_soc_id[i]]  
                                //se_dogs[2,trip_soc_id[i]]*n_dogs[i,1] +
                                //se_firearms[2,trip_soc_id[i]]*n_firearms[i,1]
                            );

            // compute linear components of production functions
            // and do target update
            if ( n_dogs[i,1] != -1 && n_firearms[i,1] != -1 ) {
                // dogs and guns observed
                // estimate base rates of dogs and guns
                n_dogs[i,1] ~ bernoulli(dogs_mu[trip_soc_id[i]]);
                n_firearms[i,1] ~ bernoulli(guns_mu[trip_soc_id[i]]);
                // build production functions
                lm_f[i] = lm_f[i] * exp( b_dogs[1,trip_soc_id[i]]*n_dogs[i,1] +
                                         b_firearms[1,trip_soc_id[i]]*n_firearms[i,1] );
                lm_h[i] = lm_h[i] * exp( b_dogs[2,trip_soc_id[i]]*n_dogs[i,1] +
                                         b_firearms[2,trip_soc_id[i]]*n_firearms[i,1] );
                sefx = sefx * exp( se_dogs[1,trip_soc_id[i]]*n_dogs[i,1] +
                                   se_firearms[1,trip_soc_id[i]]*n_firearms[i,1] );
                sehx = sehx * exp( se_dogs[2,trip_soc_id[i]]*n_dogs[i,1] +
                                se_firearms[2,trip_soc_id[i]]*n_firearms[i,1] );
                p = 2*(1 - inv_logit( skillj^sefx * lm_f[i] ));
                mu = lm_h[i] * skillj^sehx;
                if ( trip_harvests[i,1]==0 )
                    //increment_log_prob(bernoulli_log(1,p));
                    1 ~ bernoulli(p);
                else {
                    //increment_log_prob(bernoulli_log(0,p) + gamma_log(trip_harvests[i,1],mu/hscale[trip_soc_id[i]],1/hscale[trip_soc_id[i]]));
                    0 ~ bernoulli(p);
                    trip_harvests[i,1] ~ gamma(mu/hscale[trip_soc_id[i]],1/hscale[trip_soc_id[i]]);
                }
            } // no dogs/guns missing
            if ( n_dogs[i,1] == -1 && n_firearms[i,1] != -1 ) {
                // dogs missing
                n_firearms[i,1] ~ bernoulli(guns_mu[trip_soc_id[i]]);
                // average over missingness
                for ( nterm in 1:2 ) {
                    xdogs = nterm-1;
                    lm_f[i] = lm_f[i] * exp( b_dogs[1,trip_soc_id[i]]*xdogs +
                                             b_firearms[1,trip_soc_id[i]]*n_firearms[i,1] );
                    lm_h[i] = lm_h[i] * exp( b_dogs[2,trip_soc_id[i]]*xdogs +
                                             b_firearms[2,trip_soc_id[i]]*n_firearms[i,1] );
                    sefx = sefx * exp( se_dogs[1,trip_soc_id[i]]*xdogs +
                                       se_firearms[1,trip_soc_id[i]]*n_firearms[i,1] );
                    sehx = sehx * exp( se_dogs[2,trip_soc_id[i]]*xdogs +
                                       se_firearms[2,trip_soc_id[i]]*n_firearms[i,1] );
                    p = 2*(1 - inv_logit( skillj^sefx * lm_f[i] ));
                    mu = lm_h[i] * skillj^sehx;
                    LLterms[nterm] = 0;
                    if ( trip_harvests[i,1]==0 ) {
                        LLterms[nterm] = LLterms[nterm] + log(p);
                    } else {
                        LLterms[nterm] = LLterms[nterm] + log1m(p);
                        LLterms[nterm] = LLterms[nterm] + gamma_lpdf( trip_harvests[i,1] | mu/hscale[trip_soc_id[i]] , 1/hscale[trip_soc_id[i]] );
                    }
                    
                }//k
                // do the mixture
                target += log_mix( dogs_mu[trip_soc_id[i]] , LLterms[2] , LLterms[1] );
            }// !dogs √firearms
            if ( n_dogs[i,1] != -1 && n_firearms[i,1] == -1 ) {
                // dogs observed but firearms missing
                n_dogs[i,1] ~ bernoulli(dogs_mu[trip_soc_id[i]]);
                // average over missingness
                for ( nterm in 1:2 ) {
                    xguns = nterm-1;
                    lm_f[i] = lm_f[i] * exp( b_dogs[1,trip_soc_id[i]]*n_dogs[i,1] +
                                             b_firearms[1,trip_soc_id[i]]*xguns );
                    lm_h[i] = lm_h[i] * exp( b_dogs[2,trip_soc_id[i]]*n_dogs[i,1] +
                                             b_firearms[2,trip_soc_id[i]]*xguns );
                    sefx = sefx * exp( se_dogs[1,trip_soc_id[i]]*n_dogs[i,1] +
                                       se_firearms[1,trip_soc_id[i]]*xguns );
                    sehx = sehx * exp( se_dogs[2,trip_soc_id[i]]*n_dogs[i,1] +
                                       se_firearms[2,trip_soc_id[i]]*xguns );
                    p = 2*(1 - inv_logit( skillj^sefx * lm_f[i] ));
                    mu = lm_h[i] * skillj^sehx;
                    LLterms[nterm] = 0;
                    if ( trip_harvests[i,1]==0 ) {
                        LLterms[nterm] = LLterms[nterm] + log(p);
                    } else {
                        LLterms[nterm] = LLterms[nterm] + log1m(p);
                        LLterms[nterm] = LLterms[nterm] + gamma_lpdf( trip_harvests[i,1] | mu/hscale[trip_soc_id[i]] , 1/hscale[trip_soc_id[i]] );
                    }
                    
                }//k
                // do the mixture
                target += log_mix( guns_mu[trip_soc_id[i]] , LLterms[2] , LLterms[1] );
            }// √dogs !firearms
            if ( n_dogs[i,1] == -1 && n_firearms[i,1] == -1 ) {
                // miss both dogs and guns
                for ( nterm in 1:4 ) {
                    xdogs = xdogsvec[nterm];
                    xguns = xgunsvec[nterm];
                    lm_f[i] = lm_f[i] * exp( b_dogs[1,trip_soc_id[i]]*xdogs +
                                             b_firearms[1,trip_soc_id[i]]*xguns );
                    lm_h[i] = lm_h[i] * exp( b_dogs[2,trip_soc_id[i]]*xdogs +
                                             b_firearms[2,trip_soc_id[i]]*xguns );
                    sefx = sefx * exp( se_dogs[1,trip_soc_id[i]]*xdogs +
                                       se_firearms[1,trip_soc_id[i]]*xguns );
                    sehx = sehx * exp( se_dogs[2,trip_soc_id[i]]*xdogs +
                                       se_firearms[2,trip_soc_id[i]]*xguns );
                    p = 2*(1 - inv_logit( skillj^sefx * lm_f[i] ));
                    mu = lm_h[i] * skillj^sehx;
                    LL4terms[nterm] = 0;
                    if ( trip_harvests[i,1]==0 ) {
                        LL4terms[nterm] = LL4terms[nterm] + log(p);
                    } else {
                        LL4terms[nterm] = LL4terms[nterm] + log1m(p);
                        LL4terms[nterm] = LL4terms[nterm] + gamma_lpdf( trip_harvests[i,1] | mu/hscale[trip_soc_id[i]] , 1/hscale[trip_soc_id[i]] );
                    }
                    // add leading factor for probability of combination of missingness
                    if ( xdogs==1 )
                        LL4terms[nterm] = LL4terms[nterm] + log(dogs_mu[trip_soc_id[i]]);
                    else
                        LL4terms[nterm] = LL4terms[nterm] + log1m(dogs_mu[trip_soc_id[i]]);
                    if ( xguns==1 )
                        LL4terms[nterm] = LL4terms[nterm] + log(guns_mu[trip_soc_id[i]]);
                    else
                        LL4terms[nterm] = LL4terms[nterm] + log1m(guns_mu[trip_soc_id[i]]);
                }//k
                // do the mixture
                target += log_sum_exp( LL4terms );
            }// miss both dogs and guns

        } //j
    }//i over trips
}
generated quantities{
    real k[N_hunters];
    real m[N_hunters];
    real b[N_hunters];
    //matrix[6,6] Rho_societies;
    // build individual forager effects
    // useful for posterior predictions
    for ( j in 1:N_hunters ) {
        k[j] = exp( lifehistmeans[1] + vs[forager_soc_id[j],1] + vh[j,1] );
        m[j] = exp( lifehistmeans[2] + vs[forager_soc_id[j],2] + vh[j,2] );
        b[j] = exp( lifehistmeans[3] + vs[forager_soc_id[j],3] );
    }
    // build ordinary correlation matrices from cholesky factors
    // useful during interpretation
    //Rho_societies = L_Rho_societies * L_Rho_societies';
}
"

concat <- function (...) 
{
    paste(..., collapse = "", sep = "")
}

coerce_index <- function (...) 
{
    L <- list(...)
    if (is.list(L[[1]]) && length(L) == 1) 
        L <- L[[1]]
    if (length(L) == 1) {
        x <- as.integer(L[[1]])
        if (any(is.na(x))) 
            x <- as.integer(as.factor(as.character(L[[1]])))
        return(x)
    }
    else {
        vnames <- match.call()
        vnames <- as.character(vnames)[2:(length(L) + 1)]
        M <- L
        for (i in 1:length(L)) M[[i]] <- as.character(L[[i]])
        Mall <- M[[1]]
        for (i in 2:length(L)) Mall <- c(Mall, M[[i]])
        Mall <- unique(Mall)
        new_levels <- levels(as.factor(Mall))
        for (i in 1:length(L)) {
            M[[i]] <- factor(M[[i]], levels = new_levels)
            M[[i]] <- as.integer(M[[i]])
        }
        names(M) <- paste(vnames, "_idx", sep = "")
        return(M)
    }
}

check_index <- function (x) 
{
    y <- sort(unique(x))
    n <- length(y)
    message(concat("Length: ", n))
    message(concat("Range: ", min(y), " / ", max(y)))
    if (max(y) != n) 
        message("Maximum index different than number of unique values")
    diffs <- sapply(2:n, function(i) y[i] - y[i - 1])
    if (any(diffs) != 1) 
        message("At least one gap in consecutive values")
}

cchmodel_run <- function( model_code=cchmodel_full , sites=cchunts_data_sets , dogs_miss=-1 , guns_miss=-1 , prior_scale=0.5 , n_chains=4 , the_seed=1208 , warmup=500 , iter=1000 , control=list( adapt_delta=0.99 , max_treedepth=13 ) , start , test=FALSE , ... ) {

    require(parallel)

    dat <- make_joint(sites)

    N_soc <- length(unique(dat$society_id))
    message(concat(N_soc," sites joined."))
    # which societies assigned which IDs
    for ( i in 1:N_soc ) {
        an_id <- dat$society[dat$society_id==i][1]
        print( concat(i," - ",an_id) )
    }

    dat_list <- prep_data( dat , dogs_miss=dogs_miss , guns_miss=guns_miss )
    trip_form <- trans_trips( dat_list )

    n_age_miss <- trip_form$N_ages_impute
    n_hour_miss <- trip_form$N_hours_missing
    N_soc <- dat_list$N_societies

    if ( missing(start) ) {
    start <- list(
        zh = matrix(0,nrow=2,ncol=dat_list$N_hunters),
        afbar = 0,
        ahbar = 0,
        sigma_af = 0.5,
        sigma_ah = 0.5,
        af = matrix( c( 0 , 0,0,0 ), nrow=4 , ncol=N_soc ),
        ah = matrix( c( 0 , 0,0,0 ), nrow=4 , ncol=N_soc ),
        hscale=rep( 0.5 , N_soc ),
        zs = matrix(0,nrow=6,ncol=N_soc),
        lifehistmeans = c(2,-1,-0.5,-1,-0.5,0),
        sigma_societies = rep(0.5,6),
        b_hours = matrix( 0 , nrow=2 , ncol=N_soc ),
        b_dogs = matrix( 0 , nrow=2 , ncol=N_soc ),
        b_firearms = matrix( 0 , nrow=2 , ncol=N_soc ),
        b_xday = matrix( 0 , nrow=N_soc , ncol=2 ),
        age_err = rep( 0 , n_age_miss ),
        trip_duration_imputed = rep( 0 , n_hour_miss ),
        sef = rep(0,N_soc),
        seh = rep(0,N_soc),
        se_dogs = matrix( 0 , nrow=2 , ncol=N_soc ),
        se_firearms = matrix( 0 , nrow=2 , ncol=N_soc ),
        ache_fix_rho = 0,
        dogs_mu = rep(0.5,N_soc),
        guns_mu = rep(0.5,N_soc)
    )
    }

    # define scale for many regularizing priors in model
    # 0.5 gives efficient mixing
    dat_list$prior_scale <- prior_scale

    init <- function() return(start)
    max_cores <- detectCores()

    set.seed(the_seed)

    if ( test==TRUE ) {
        n_chains <- 1
        warmup <- 2
        iter <- 4
    }

    mfit <- stan( 
        model_code=model_code , 
        data=dat_list , 
        warmup=warmup , iter=iter, 
        chains=n_chains , cores=min(n_chains,max_cores) , 
        control=control , 
        init=init , ... )

    return(mfit)

}
