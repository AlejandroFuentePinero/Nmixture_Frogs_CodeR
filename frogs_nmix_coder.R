# CodaR talk on N-mixture model
# Name: Alejandro de la Fuente
# Date:


# Libraries ---------------------------------------------------------------

library(jagsUI)
library(tidyverse)


# Simulate some data ------------------------------------------------------


###
# STEP 1 : SIMULATION OF THE TRUE STATE OF NATURE
###

set.seed(123) # so we all get the same dataset

Nsites <- 30

Nsurveys <- 3

C <- matrix(NA, nrow = Nsites, ncol = Nsurveys) # this is the matrix use for the observed counts

(Rainfall <- sort(runif(Nsites, -1, 1))) # create a predictor

a <- 2 # log-scale intercept

chytrid_severity <- -3 # log-scale slope for "Rainfall"

lambda <- exp(a + chytrid_severity * Rainfall) # Expected abundance this is just a linear model (y ~ a + b*x)

# Lets do the first plot to see how things are looking

plot(Rainfall, lambda, type = "l", lwd = 3)

# Local abundance simulation

N <- rpois(Nsites, lambda) # rpois for count data

points(Rainfall, N) # Lets add the points we just estimated to the plot. Those points represent the "realised" abundance

table(N) # in 10 sites, this species is absent

# Lets visualise the true system state

plot(Rainfall, lambda, xlab = "Rainfall", ylab = "True abundance (N)", frame = F, cex = 1.5, xlim = c(-1,1))
lines(seq(-1,1,,100), exp(a + chytrid_severity * seq(-1,1,,100)), lwd = 3, col = "red")


###
# STEP 2 : SIMULATION OF THE OBSERVATION PROCESS
###

Night_temp <- array(runif(Nsites * Nsurveys, -1, 1), dim= c(Nsites, Nsurveys)) # night temperature is expected to affect frog detection

a2 <- -2 #logit-scale intercept

b_temp <- 3 # logit-scale slope for night temperature

p <- plogis(a2 + b_temp * Night_temp) # detection probability

plot(p~Night_temp, ylim=c(0,1))

# Observation process simulation

for(survey in 1:Nsurveys){
  
  C[,survey] <- rbinom(Nsites, N, p[,survey])
  
} # loop survey

plot(Night_temp, C/max(C), xlab = "Night temperature (°C)", ylab = "Observed counts (scaled)", frame = F, cex = 1.5)
lines(seq(-1,1,,100), plogis(a2 + b_temp * seq(-1,1,,100)), lwd = 3, col = "red")


# Data exploration --------------------------------------------------------


cbind(lambda = round(lambda,2), N = N, C1 = C[,1], C2 = C[,2], C3 = C[,3])


# Write the model ---------------------------------------------------------

cat(file = "frog_model.txt", "
    model{
    # 1. Priors
    
    a_abun ~ dunif(-10,10)
    a_det ~ dunif(-10,10)
    chytrid_severity ~ dunif(-10,10)
    b_temp ~ dunif(-10,10)
    
    # 2. Likelihood
    
      # 2.1. Ecological model for true abundance
      
      for(s in 1:Nsites){
      
        N[s] ~ dpois(lambda[s])
      
        log(lambda[s]) <- a_abun + chytrid_severity * Rainfall[s]
      
      # 2.2. Observation model for replicate counts
      
        for(j in 1:Nsurveys){
          
         C[s,j] ~ dbin(p[s,j], N[s])
         
         logit(p[s,j]) <- a_det + b_temp * Night_temp[s,j]
      
       } # loop j
      } # loop s
     } # close model
    ")


# Run model ---------------------------------------------------------------

# Compile data

str(bdata <- list(C = C, 
                  Nsites = Nsites, 
                  Nsurveys = Nsurveys, 
                  Rainfall = Rainfall, 
                  Night_temp = Night_temp))

# Initial values

Nst <- apply(C, 1, max) + 1

inits <- function() list(N = Nst,
                         a_abun = rnorm(1),
                         a_det = rnorm(1),
                         chytrid_severity = rnorm(1),
                         b_temp = rnorm(1))

# Parameters to monitor

params <- c("chytrid_severity", "b_temp",
            "N")

# MCMC set up

nc <- 3 # number of chains

ni <- 22000 # lenght of the chains

nb <- 2000 # burnin

nt <- 10 # thining

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

comp <- tibble(true  = N,
               estimated = out$mean$N,
               lo = out$q2.5$N,
               hi = out$q97.5$N)

comp %>% ggplot(aes(x = (1:30), y = true), col = "black")+
  geom_pointrange(aes(y = estimated, ymin = lo, ymax = hi), col = "red", size = 1)+
  geom_point(size = 4)+
  labs(x = "\nSite", y = "Abundance\n")+
  theme_classic()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))

