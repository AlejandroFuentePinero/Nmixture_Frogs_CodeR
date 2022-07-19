# CodaR talk on N-mixture model
# Name: Alejandro de la Fuente
# Date:


# Libraries ---------------------------------------------------------------

library(jagsUI)

# Simulate some data ------------------------------------------------------


###
# STEP 1 : SIMULATION OF THE TRUE STATE OF NATURE
###

set.seed(123) # so we all get the same dataset

Nsites <- 30

Nsurveys <- 3

C <- matrix(NA, nrow = Nsites, ncol = Nsurveys) # this is the matrix use for the observed counts

(Chytrid <- sort(runif(Nsites, -1, 1))) # create a predictor

a <- 2 # log-scale intercept

chytrid_severity <- -3 # log-scale slope for "Chytrid"

lambda <- exp(a + chytrid_severity * Chytrid) # Expected abundance this is just a linear model (y ~ a + b*x)

# Lets do the first plot to see how things are looking

plot(Chytrid, lambda, type = "l", lwd = 3, xlab = "Chytrid severity", ylab = "Expected abundance", xlim = c(-1,1))

# Local abundance simulation

N <- rpois(Nsites, lambda) # rpois for count data

points(Chytrid, N) # Lets add the points we just estimated to the plot. Those points represent the "realised" abundance

table(N) # in 10 sites, this species is absent

# Lets visualise the true system state

plot(Chytrid, lambda, xlab = "Chytrid severity", ylab = "True abundance (N)", frame = F, cex = 1.5, xlim = c(-1,1))
lines(seq(-1,1,,100), exp(a + chytrid_severity * seq(-1,1,,100)), lwd = 3, col = "red")


###
# STEP 2 : SIMULATION OF THE OBSERVATION PROCESS
###

Env_wet <- array(runif(Nsites * Nsurveys, -1, 1), dim= c(Nsites, Nsurveys)) # environmental wetness is expected to affect frog detection

a2 <- -2 #logit-scale intercept

b_wet <- 3 # logit-scale slope for night temperature

p <- plogis(a2 + b_wet * Env_wet) # detection probability

plot(p~Env_wet, ylim=c(0,1), xlab = "Wet", ylab = "Detection probability", frame = F)

# Observation process simulation

for(survey in 1:Nsurveys){
  
  C[,survey] <- rbinom(Nsites, N, p[,survey])
  
} # loop survey

plot(Env_wet, C/max(C), xlab = "Wet", ylab = "Observed counts (scaled)", frame = F, cex = 1.5, xlim = c(-1,1))
lines(seq(-1,1,,100), plogis(a2 + b_wet * seq(-1,1,,100)), lwd = 3, col = "red")


# Data exploration --------------------------------------------------------


cbind(lambda = round(lambda,2), N = N, C1 = C[,1], C2 = C[,2], C3 = C[,3])


# Write the model ---------------------------------------------------------

cat(file = "frog_model.txt", "
    model{
    # 1. Priors
    
    a_abun ~ dunif(-10,10)
    a_det ~ dunif(-10,10)
    chytrid_severity ~ dunif(-10,10)
    b_wet ~ dunif(-10,10)
    
    # 2. Likelihood
    
      # 2.1. Ecological model for true abundance
      
      for(s in 1:Nsites){
      
        N[s] ~ dpois(lambda[s])
      
        log(lambda[s]) <- a_abun + chytrid_severity * Chytrid[s]
      
      # 2.2. Observation model for replicate counts
      
        for(j in 1:Nsurveys){
          
         C[s,j] ~ dbin(p[s,j], N[s])
         
         logit(p[s,j]) <- a_det + b_wet * Env_wet[s,j]
      
       } # loop j
      } # loop s
     } # close model
    ")


# Run model ---------------------------------------------------------------

# Compile data

str(bdata <- list(C = C, 
                  Nsites = Nsites, 
                  Nsurveys = Nsurveys, 
                  Chytrid = Chytrid, 
                  Env_wet = Env_wet))

# Initial values

Nst <- apply(C, 1, max) + 1

inits <- function() list(N = Nst,
                         a_abun = rnorm(1),
                         a_det = rnorm(1),
                         chytrid_severity = rnorm(1),
                         b_wet = rnorm(1))

# Parameters to monitor

params <- c("chytrid_severity", "b_wet",
            "N")

# MCMC set up

nc <- 3 # number of chains

ni <- 22000 # length of the chains

nb <- 2000 # burn-in

nt <- 10 # thinning

# Run model

out <- jags(data = bdata,
            inits = inits,
            parameters.to.save = params,
            model.file = "frog_model.txt",
            n.chains = nc,
            n.thin = nt,
            n.iter = ni,
            n.burnin = nb,
            parallel = TRUE) 

print(out,2)


# Let's compare the estimations to the true abundance


plot(c(1:30), N, xlab = "Site", ylab = "Abundance", frame = F, cex = 1.5, pch=16, xlim = c(1,30), col = "black", ylim = c(0,230))
points(c(1:30), out$mean$N, cex = 1.5, col = "red", pch=16)
points(c(1:30),apply(C, 1, max), cex = 1.5, col = "blue", pch=16)
arrows(x0=c(1:30), y0=out$q2.5$N, x1=c(1:30), y1=out$q97.5$N, code=3, col="red", length = 0.05, angle = 90)
legend(15, 150, legend=c("True", "Estimated", "Observed"),
       col=c("black", "red", "blue"), pch=16, cex=0.8)

