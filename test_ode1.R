#R
require(PMXStan)
simple.gof = function()
{
    fs = summary(fit)
    ix = grep("y_pred", dimnames(fs$summary)[[1]])
    yp = fs$summary[ix,c("50%")]
    plot(yp, dat$conc); abline(0,1,col="red")
}

#-- generate ODE patch
ode = "
   #C2 = centr/V;
   d/dt(depot) =-ka*depot;
   d/dt(centr) = ka*depot - ke*centr;
"
instant.stan.extension(ode)

#-- stan code for generic_ode_interface
stan_src = "
data{
    int<lower=0> NOBS;                 // number of observations for each patient
    int<lower=0> NEVTS;                // number of observations for each patient
    vector[NOBS] conc;                 // observations of concentration
    vector[2] inits;                   // 
    vector<lower=0>[NEVTS] obs_time;   // observation time for all patients
    vector<lower=0>[NEVTS] evid;       // 
    vector<lower=0>[NEVTS] amt;        // 
}
parameters{
    vector<lower=-3.0, upper=3.0>[3] theta; //fixed effect
    real<lower=0> sigma2;
}
transformed parameters{
    real<lower=0> ka;
    real<lower=0> ke;
    real<lower=0> V;
    real<lower=0> sigma;
    vector[NOBS] y_pred;
    ka = exp( theta[1] );
    ke = exp( theta[2] );
    V  = exp( theta[3] );
    sigma = sqrt(sigma2);
    
    {
        vector[NOBS] g;
        vector[3] params;
        params[1] = ka;
        params[2] = ke;
        params[3] = V;
        g = generic_ode_interface(
            params,
            inits,
            obs_time,
            evid,
            amt,  
            1E-4,
            1E-4,
            NOBS,
            2);
        for(j in 1:NOBS)
            y_pred[j] = g[j]/V;
    } //end of local variable
}//end of transformed parameters block
model{
    for(k in 1:3){
        theta[k] ~ normal(0.,1000.);
    }
    sigma2 ~ inv_gamma(.01,.01);
    conc ~ normal(y_pred, sigma);
}
"
cat(stan_src, file='tmp.stan')

#-- stan data
dat = list(
NOBS = 11,
NEVTS = 12,
conc = c(-0.3796411, 4.9439890, 9.3984885, 12.2915293, 14.6313447, 
         16.2968084, 17.1462266, 17.6130543, 17.3000089, 18.1211557, 16.6071588),
inits = c(0, 0),
obs_time = c(0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
evid   = c(101, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
amt    = c(4,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
)
str(as.list(dat))

sm = stan_model("tmp.stan")
fit = sampling(sm, data=dat, chains=1, iter=400)
simple.gof()
