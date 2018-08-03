make_fake_data <- function(nSamples=20,
                           nTotal = 20,
                           nu = 0.95,
                           m0 = c(0, rep(1,15)), #initial condition
                           th = c(6.8,132.8),
                           sig = 1,
                           phi = 0.289) {
print('generating fake data')
#############################################################
source('get_nc_transition_matrix.R')
B1 = get_nc_transition_matrix(0) %>% as.vector
B2 = get_nc_transition_matrix(1) %>% as.vector
source('extract_times_and_scaling.R')
times <- extract_times_and_scaling(nSamples,0,optional_plot=FALSE,test_on_mutant_data=TRUE)

mc <- stan_model('model_comparison_at_stst2.stan')
expose_stan_functions(mc)
#currently will not be using the absorbing oocyte fn
B = construct_matrix(nu) %>% as.vector
#############################################################
samples <- stan(file = 'mrna_transport6.stan',
                data = list (
                  T  = nTotal,
                  y0 = m0,
                  t0 = (times$t0),
                  ts = (times$ts1),
                  theta = array(th, dim = 2),
                  sigma = sig,
                  phi = phi,
                  B = B
                ),
                algorithm="Fixed_param",
                seed = 42,
                chains = 1,
                iter =100, 
                refresh = -1
)  
s <- rstan::extract(samples,permuted=FALSE)
plot(s[1,1,seq(from=1, to=(nTotal*16-1), by=nTotal)])
plot(s[1,1,seq(from=nSamples, to=(nTotal*16), by=nTotal)])
boxplot(s[,1,seq(from=nSamples, to=(nTotal*16), by=nTotal)])
full_data = matrix(s[1,1,1:(16*(nTotal))],nrow=(nTotal),byrow=FALSE) #this is our fake data

return(full_data)
}

