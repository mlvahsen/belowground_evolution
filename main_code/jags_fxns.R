# Vahsen et al. Belowground Trait Evolution
# JAGS models and functions to source for analysis of common garden trait data

## Write JAGS models ####
# Create Gaussian random intercept model
sink("main_code/Gaussian_RI_Model.txt")
cat("
model{
    sigma.int ~ dunif(0, 100)
    tau.int <- 1/(sigma.int^2)
    
    for (j in 1:nalpha) {
    alpha[j] ~ dnorm(mu.alpha, tau.int) # Random by-subject deflections to the intercept
    }
    
    mu.alpha ~ dnorm(0, 0.001) # The mean intercept
    
    for (k in 1:nbeta){
    beta[k] ~ dnorm(0, 0.001) # Common slopes
    }
    
    sigma.res ~ dunif(0, 100) # Residual standard deviation
    tau.res <- 1/(sigma.res^2)
    
    # Likelihood
    for (i in 1:length(y)) {
    mu[i] <- alpha[genotype[i]] + beta%*%X[i,] 
    y[i] ~ dnorm(mu[i], tau.res) # The actual (random) responses
    

    }
    
    }",fill = TRUE)
sink()

# Create Poisson random intercept model
sink("main_code/Poisson_RI_Model.txt")
cat("
model{
  sigma.int ~ dunif(0.001, 10)
  tau.int <- 1/(sigma.int^2)
  
  for (j in 1:nalpha) {
    alpha[j] ~ dnorm(mu.alpha, tau.int) # Random by-subject deflections to the intercept
  }
  
  mu.alpha ~ dnorm(0, 0.0001) # The mean intercept
  
  for (k in 1:nbeta){
    beta[k] ~ dnorm(0, 0.0001) # Common slopes
  }
  
  # Likelihood
  for (i in 1:length(y)) {
    mu[i] <- alpha[genotype[i]] + beta%*%X[i,] # Expectation
    log(lambda[i]) <- mu[i]
    y[i] ~ dpois(lambda[i]) # The actual (random) responses
  }
}",fill = TRUE)
sink()

## Create functions to extract JAGS samples ####

# Create seed
seed <- 1234
# Set burn-in
n.burn <- 1000

# Define get_samples functions (runs jags and extracts jags.samples and
# dic.samples) for three types of models (Gaussian w/ genotype random intercept,
# Poisson w/ genotype random intercept, Gaussian without genotype random intercept)
get_samples_fixed <- function(model, data, trait, seed){
  # Get model matrices from linear models
  M_all <- as.matrix(model.matrix(model)[ , -1])
  
  # Scale (and center) continuous covariates
  M_all_sc <- M_all
  M_all_sc[,1] <- scale(M_all[,1], scale = T)
  M_all_sc[,2] <- scale(M_all[,2], scale = T)
  
  # Set up data structure
  nbeta <- ncol(M_all_sc)
  y <- pull(data[,trait])

  data_jags = list(
    y = y,
    X = M_all_sc,
    nbeta = nbeta
  )
  
  # Set initial values for chains (2)
  inits = list(
    list(
      alpha = 1,
      beta = rep(0, nbeta),
      sigma.res = 20,
      .RNG.name = "base::Mersenne-Twister",
      .RNG.seed = seed
    ),
    list(
      alpha = 2,
      beta = rep(5, nbeta),
      sigma.res = 2,
      .RNG.name = "base::Mersenne-Twister",
      .RNG.seed = seed
    )
  )
  
  # Fit JAGS model with gaussian likelihood and genotype as random intercept
  jm_mono = jags.model("main_code/Gaussian_Model", data=data_jags, 
                       n.adapt = n.adapt, inits=inits, n.chains=length(inits))
  # Add burn-in
  update(jm_mono, n.burn)
  # Get DIC samples
  dic_mono = dic.samples(jm_mono, n.iter = n.iter)
  
  # Return both JAGS samples and DIC samples
  return(list(dic_mono = dic_mono, jm_mono = jm_mono))}

get_samples <- function(model, data, trait, seed){
  # Get model matrices from linear models
  M_all <- as.matrix(model.matrix(model)[ , -1])
  
  # Scale (and center) continuous covariates
  M_all_sc <- M_all
  M_all_sc[,1] <- scale(M_all[,1], scale = T)
  M_all_sc[,2] <- scale(M_all[,2], scale = T)
  
  # Set up data structure
  nbeta <- ncol(M_all_sc)
  genotype <- as.numeric(droplevels(as.factor(data$genotype)))
  y <- pull(data[,trait])
  nalpha <- length(unique(genotype))
  
  data_jags = list(
    y = y,
    X = M_all_sc,
    genotype = genotype,
    nalpha = nalpha,
    nbeta = nbeta
  )
  
  # Set initial values for chains (2)
  inits = list(
    list(
      alpha = rep(0, nalpha),
      beta = rep(0, nbeta),
      sigma.res = 20,
      sigma.int = 10,
      mu.alpha = 8,
      .RNG.name = "base::Mersenne-Twister",
      .RNG.seed = seed
    ),
    list(
      alpha = rep(5, nalpha),
      beta = rep(5, nbeta),
      sigma.res = 2,
      sigma.int = 12,
      mu.alpha = 8,
      .RNG.name = "base::Mersenne-Twister",
      .RNG.seed = seed
    )
  )
  
  # Fit JAGS model with gaussian likelihood and genotype as random intercept
  jm_mono = jags.model("main_code/Gaussian_RI_Model.txt", data=data_jags, 
                       n.adapt = n.adapt, inits=inits, n.chains=length(inits))
  # Remove burn-in
  update(jm_mono, n.burn)
  # Get DIC samples
  dic_mono = dic.samples(jm_mono, n.iter = n.iter)
  
  # Return both JAGS samples and DIC samples
  return(list(dic_mono = dic_mono, jm_mono = jm_mono))}

get_samples_pois <- function(model, data, seed){# Get model matrices from linear models
  
  M_all <- as.matrix(model.matrix(model)[ , -1])
  
  # Scale continuous covariates
  M_all_sc <- M_all
  M_all_sc[,1] <- scale(M_all[,1], scale = T)
  M_all_sc[,2] <- scale(M_all[,2], scale = T)
  
  # Set up data structure
  nbeta <- ncol(M_all_sc)
  genotype <- as.numeric(droplevels(as.factor(data$genotype)))
  y <- data$density
  nalpha <- length(unique(genotype))
  
  data_jags = list(
    y = y,
    X = M_all_sc,
    genotype = genotype,
    nalpha = nalpha,
    nbeta = nbeta
  )
  
  # Set initial values for chains (2)
  inits = list(
    list(
      alpha = rep(1, nalpha),
      beta = rep(1, nbeta),
      sigma.int = 2,
      mu.alpha = 8,
      .RNG.name = "base::Mersenne-Twister",
      .RNG.seed = seed
    ),
    list(
      alpha = rep(5, nalpha),
      beta = rep(5, nbeta),
      sigma.int = 1,
      mu.alpha = 8,
      .RNG.name = "base::Mersenne-Twister",
      .RNG.seed = seed
    )
  )
  
  jm_mono = jags.model("main_code/Poisson_RI_Model.txt", data=data_jags, 
                       n.adapt = n.adapt, inits=inits, n.chains=length(inits))
  # Remove burn-in
  update(jm_mono, n.burn)
  dic_mono = dic.samples(jm_mono, n.iter = n.iter)
  
  return(list(dic_mono = dic_mono, jm_mono = jm_mono))}

## Create function to run all models ####

# This will run the model that will be either Gaussian or Poisson based on the
# trait specified and will show convergence diagnostic plots if diag_plot is set
# to TRUE
run_models <- function(trait_data, trait, model_template, diag_plot = F, seed){
  # Set up data structure
  if(trait == "density"){
    jags_out <- get_samples_pois(model_template, trait_data, seed)
  }else{
    jags_out <- get_samples(model_template, trait_data, trait, seed)
  }
  model_coda <- coda.samples(jags_out$jm_mono, variable.names = c("beta", "alpha","sigma.res","sigma.int", "mu.alpha"), n.iter=n.iter, thin = thin)  
  
  ##
  # Model diagnostics
  ##
  if(diag_plot == T){
    # Autocorrelation
    ggs_autocorrelation(samples_alphas)
    ggs_autocorrelation(samples_betas)
    ggs_autocorrelation(samples_sigmas)
    
    # Density plots
    ggs_density(samples_alphas) + facet_wrap(~Parameter, nrow = 3) +
      theme_classic()
    ggs_density(samples_betas) + facet_wrap(~Parameter, nrow = 1) +
      theme_classic()
    ggs_density(samples_sigmas) + facet_wrap(~Parameter, nrow = 1) +
      theme_classic()
    
    # Trace plots
    ggs_traceplot(samples_alphas) + facet_wrap(~Parameter, nrow = 6) +
      theme_classic()
    ggs_traceplot(samples_betas) + facet_wrap(~Parameter, nrow = 2) +
      theme_classic()
    ggs_traceplot(samples_sigmas) + facet_wrap(~Parameter, nrow = 2) +
      theme_classic() + scale_color_manual(values = c("blue", "green"))
    
    # Gelman-Rubin
    gelman.diag(samples_coda) 
    
  }
  
  return(model_coda)
  
}

run_models_fixed <- function(trait_data, trait, model_template, diag_plot = F, seed){
  # Set up data structure
  jags_out <- get_samples_fixed(model_template, trait_data, trait, seed)
  model_coda <- coda.samples(jags_out$jm_mono, variable.names = c("beta", "alpha","sigma.res","sigma.int", "mu.alpha"), n.iter=n.iter, thin = thin)  
  
  return(model_coda)
  
}
