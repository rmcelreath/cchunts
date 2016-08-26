ccmodel_full <- "
// multisociety varying effects model
// varying intercepts on forager life history parameters
// imputation of missing trip durations
// adds cobb-douglas style elasticity to skill components
// TO DO:
// add sex of forager
// potential to add: observed (indicator of anthropologist presence)
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
}
transformed data{
    // loop over rows and compute trip information
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
    print("Number of trip durations to impute: ",N_hours_missing);
    print("Number of ages to impute: ",N_ages_impute);
}
parameters{
    // forager life history parameters
    vector[3] lifehistmeans;                // averages across societies
    // ** centered
    //vector[3] vs[N_societies];              // log k, log m, b for each society
    //vector<lower=0>[3] sigma_societies;     // stddev across societies
    //corr_matrix[3] Rho_societies;           // correlation matrix
    // ** non-centered
    matrix[3,N_societies] zs;         // z-scores for each hunter
    vector<lower=0>[3] sigma_societies;       // stddev of varying effects
    cholesky_factor_corr[3] L_Rho_societies;  // cholesky factor for corr matrix

    // intercept components of outcomes
    real af[4,N_societies];     // intercepts, failures
    real ah[4,N_societies];     // intercepts, harvests
                                // [1] main failure intercept
                                // [2] marginal group trips
                                // [3] marginal pooled trips
                                // [4] marginal assistant effects
    real<lower=0> hscale[N_societies];
    real sef[N_societies]; // skill elasticity failures - log scale
    real seh[N_societies]; // skull elasticity harvests - log scale

    // hunter effects
    // ** centered parameterization
    //vector[2] vh[N_hunters];        // hunter specific parameters: k,m
    //vector<lower=0>[2] sigma;       // stddev of varying effects
    //corr_matrix[2] Rho;             // correlation among varying effects
    // ** non-centered parameterization
    matrix[2,N_hunters] zh;         // z-scores for each hunter
    vector<lower=0>[2] sigma;       // stddev of varying effects
    cholesky_factor_corr[2] L_Rho;  // cholesky factor for corr matrix

    // various regression parameters
    // row 1: failures
    // row 2: harvests
    real b_hours[2,N_societies];
    real b_dogs[2,N_societies];
    real b_firearms[2,N_societies];
    real se_dogs[2,N_societies];
    real se_firearms[2,N_societies];

    // need pooling for xday parameters, on account of some societies having no variation
    vector[2] b_xday[N_societies];     // trip exceeded a day (not day_trip==1)
    vector<lower=0>[2] sigma_xday;
    corr_matrix[2] Rho_xday;

    // imputation
    // trip durations
    real trip_duration_imputed[N_hours_missing];
    real trip_duration_mu[N_societies];         // mean and stddev of prior
    real<lower=0> trip_duration_sigma[N_societies];
    // ages
    real age_err[N_ages_impute];
}
transformed parameters{
    // code to transform non-centered prior to centered effects
    matrix[N_hunters,2] vh;
    matrix[N_societies,3] vs;
    vh = (diag_pre_multiply(sigma,L_Rho) * zh)';
    vs = (diag_pre_multiply(sigma_societies,L_Rho_societies) * zs)';
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
    
    // priors
    // society-level life history means --- centered on global means
    //vs ~ multi_normal( lifehistmeans , quad_form_diag(Rho_societies,sigma_societies) );
    //vs ~ multi_normal( zeros3 , quad_form_diag(Rho_societies,sigma_societies) );
    to_vector(zs) ~ normal(0,1);
    lifehistmeans ~ normal(1,1);
    sigma_societies ~ exponential(10);
    //Rho_societies ~ lkj_corr(4);
    L_Rho_societies ~ lkj_corr_cholesky(4);

    for ( s in 1:N_societies ) {
        for ( i in 1:4 ) {
            af[i,s] ~ normal(0,1);
            ah[i,s] ~ normal(0,1);
        }
        sef[s] ~ normal(0,1);
        seh[s] ~ normal(0,1);
        for ( i in 1:2 ) {
            b_hours[i,s] ~ normal(0,1);
            b_dogs[i,s] ~ normal(0,1);
            b_firearms[i,s] ~ normal(0,1);
            se_dogs[i,s] ~ normal(0,1);
            se_firearms[i,s] ~ normal(0,1);
        }
    }//s
    hscale ~ exponential(1);

    // xday parameters
    b_xday ~ multi_normal( zeros2 , quad_form_diag(Rho_xday,sigma_xday) );
    sigma_xday ~ exponential(1);
    Rho_xday ~ lkj_corr(2);

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

    // varying effects
    // foragers --- these are zero-centered
    //Sigma = quad_form_diag(Rho,sigma);
    //vh ~ multi_normal( zeros2 , Sigma );
    to_vector(zh) ~ normal(0,1);
    sigma ~ exponential(2);
    //Rho ~ lkj_corr(4);
    L_Rho ~ lkj_corr_cholesky(4);

    // trip durations
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
    for ( j in 1:N_societies ) trip_duration_mu[j] ~ normal(0,1);
    trip_duration_sigma ~ exponential(1);
    
    // prep hunter effects so can re-use
    for ( j in 1:N_hunters ) {
        k[j] = exp( lifehistmeans[1] + vs[forager_soc_id[j],1] + vh[j,1] );
        m[j] = exp( lifehistmeans[2] + vs[forager_soc_id[j],2] + vh[j,2] );
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
        if ( trip_pooled[i]==1 ) {
            real avg_skill;
            // pooled harvest
            avg_skill = 0;
            lm_f[i] = exp( af[1,trip_soc_id[i]] + 
                            af[3,trip_soc_id[i]]*(n_foragers[i]-1) + 
                            af[4,trip_soc_id[i]]*n_assistants[i,1] + 
                            b_dogs[1,trip_soc_id[i]]*n_dogs[i,1] +
                            b_firearms[1,trip_soc_id[i]]*n_firearms[i,1] +
                            b_xday[trip_soc_id[i],1]*trip_xday[i]
                        ) *
                        exp(trip_duration_merge[i])^b_hours[1,trip_soc_id[i]];
            lm_h[i] = exp( ah[1,trip_soc_id[i]] + 
                            ah[3,trip_soc_id[i]]*(n_foragers[i]-1) + 
                            ah[4,trip_soc_id[i]]*n_assistants[i,1] +
                            b_dogs[2,trip_soc_id[i]]*n_dogs[i,1] +
                            b_firearms[2,trip_soc_id[i]]*n_firearms[i,1] +
                            b_xday[trip_soc_id[i],2]*trip_xday[i]
                        )*
                        exp(trip_duration_merge[i])^b_hours[2,trip_soc_id[i]];
            sefx = exp( sef[trip_soc_id[i]] + 
                            se_dogs[1,trip_soc_id[i]]*n_dogs[i,1] +
                            se_firearms[1,trip_soc_id[i]]*n_firearms[i,1]
                        );
            sehx = exp( seh[trip_soc_id[i]] + 
                            se_dogs[2,trip_soc_id[i]]*n_dogs[i,1] +
                            se_firearms[2,trip_soc_id[i]]*n_firearms[i,1]
                        );
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
            }#j
            avg_skill = avg_skill/n_foragers[i] + 0.001;
            p = 2*(1 - inv_logit( avg_skill^sefx * lm_f[i] ));
            mu = lm_h[i] * avg_skill^sehx;
            if ( trip_harvests[i,1]==0 )
                #increment_log_prob(bernoulli_log(1,p));
                1 ~ bernoulli(p);
            else {
                #increment_log_prob(bernoulli_log(0,p) + gamma_log(trip_harvests[i,1],mu/hscale[trip_soc_id[i]],1/hscale[trip_soc_id[i]]));
                0 ~ bernoulli(p);
                trip_harvests[i,1] ~ gamma(mu/hscale[trip_soc_id[i]],1/hscale[trip_soc_id[i]]);
            }
        } else {
            // independent harvests to predict
            for ( j in 1:n_foragers[i] ) {
                lm_f[i] = exp( af[1,trip_soc_id[i]] + 
                                af[2,trip_soc_id[i]]*(n_foragers[i]-1) + 
                                af[4,trip_soc_id[i]]*n_assistants[i,j] +
                                b_dogs[1,trip_soc_id[i]]*n_dogs[i,j] +
                                b_firearms[1,trip_soc_id[i]]*n_firearms[i,j] +
                                b_xday[trip_soc_id[i],1]*trip_xday[i]
                            )*
                            exp(trip_duration_merge[i])^b_hours[1,trip_soc_id[i]];
                lm_h[i] = exp( ah[1,trip_soc_id[i]] + 
                                ah[2,trip_soc_id[i]]*(n_foragers[i]-1) + 
                                ah[4,trip_soc_id[i]]*n_assistants[i,j] +
                                b_dogs[2,trip_soc_id[i]]*n_dogs[i,j] +
                                b_firearms[2,trip_soc_id[i]]*n_firearms[i,j] +
                                b_xday[trip_soc_id[i],1]*trip_xday[i]
                            )*
                            exp(trip_duration_merge[i])^b_hours[2,trip_soc_id[i]];
                sefx = exp( sef[trip_soc_id[i]] + 
                            se_dogs[1,trip_soc_id[i]]*n_dogs[i,j] +
                            se_firearms[1,trip_soc_id[i]]*n_firearms[i,j]
                        );
                sehx = exp( seh[trip_soc_id[i]] + 
                            se_dogs[2,trip_soc_id[i]]*n_dogs[i,j] +
                            se_firearms[2,trip_soc_id[i]]*n_firearms[i,j]
                        );
                //ai = forager_age[i,j];

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
                p = 2*(1 - inv_logit( skillj^sefx * lm_f[i] ));
                mu = lm_h[i] * skillj^sehx;
                if ( trip_harvests[i,j] == 0 )
                    #increment_log_prob(bernoulli_log(1,p));
                    1 ~ bernoulli(p);
                else {
                    #increment_log_prob(bernoulli_log(0,p) + gamma_log(trip_harvests[i,j],mu/hscale[trip_soc_id[i]],1/hscale[trip_soc_id[i]]));
                    0 ~ bernoulli(p);
                    trip_harvests[i,j] ~ gamma(mu/hscale[trip_soc_id[i]],1/hscale[trip_soc_id[i]]);
                }
            }#j
        }
    }#i over trips
}
generated quantities{
    real k[N_hunters];
    real m[N_hunters];
    real b[N_hunters];
    matrix[2,2] Rho;
    matrix[3,3] Rho_societies;
    // build individual forager effects
    // useful for posterior predictions
    for ( j in 1:N_hunters ) {
        k[j] = exp( lifehistmeans[1] + vs[forager_soc_id[j],1] + vh[j,1] );
        m[j] = exp( lifehistmeans[2] + vs[forager_soc_id[j],2] + vh[j,2] );
        b[j] = exp( lifehistmeans[3] + vs[forager_soc_id[j],3] );
    }
    // build ordinary correlation matrices from cholesky factors
    // useful during interpretation
    Rho = L_Rho * L_Rho';
    Rho_societies = L_Rho_societies * L_Rho_societies';
}
"