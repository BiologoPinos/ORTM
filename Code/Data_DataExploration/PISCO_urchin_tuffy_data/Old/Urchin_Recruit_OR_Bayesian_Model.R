# TITLE: Urchin_Recruit_OR_Bayesian_Model
# Author: Andres Pinos-Sanchez

# Objective:
# 1. estimate the normalized year to year variance (SD) of incoming recruits
#     using a log normal distribution of recruitment counts

# I'd be building a hierarchical structure where my incoming settlers 
# data informs a log-normal distributed annual means, and those annual means 
# inform a normally distributed year to year variance (SD)

#------------------------------- TOP LEVEL STUFF ------------------------------

# Clear all
rm(list = ls()) #environment

# Load libraries
library(nimble)
library(tidyverse)
library(mcmcplots)
library(MCMCvis)
library(coda)

set.seed(123)

# Set working directory (Why not)
setwd("C:/Users/pinosa/OneDrive - Oregon State University/Ed_OSU_Thesis_GradSchool/OR_Trophic_Model/Data/PISCO_Recruitment_Data")

# Nimble
source("attach.nimble.R")

# PISCO data
pisco <- read_csv("PISCO_all_years.csv") # source("PISCO_data_prep.R")

# - - - - - - - - - - - - - - - End of section - - - - - - - - - - - - - - - -



#---------------------------- Let's do some modeling --------------------------

# Filter for taxa 66,67,70 and drop unwanted sites
urchins_raw <- pisco %>%
  filter(count_classcode %in% c("66","67","70"), # data 2001 - 2011
         !site_code %in% c("CMEN00","CMES00","IFBRXX","KHLX00",
                           "PSGX00","SHCX00","TRHX00","IPTAXX")) # 14 sites

# Group (transects) by year, site_code, zone, replicate, 
# then sum within transects, and compute annual mean and sd
urchins <- urchins_raw %>%
  group_by(year, site_code, zone, replicate) %>%
  summarise(count = sum(count))

# Log-transform counts (add 1 to avoid log(0))
urchins <- urchins %>%
  mutate(log_count = log(count + 1))

# Nimble model

model_code <- nimbleCode({
  
  # Hyperpriors
  mu         ~ dnorm(0, 1)           # Overall mean
  sigma_year ~ T(dt(0, 1, 1), 0, )   # Half-Cauchy for year SD
  sigma_obs  ~ T(dt(0, 1, 1), 0, )   # Half-Cauchy for residual SD
  
  # Random year effects
  for(j in 1:n_years){
    year_effect[j] ~ dnorm(mu, sd = sigma_year)
  }#j
  
  # Observation model
  for(i in 1:n_obs){
    
    log_count[i] ~ dnorm(mean = year_effect[year_index[i]], sd = sigma_obs)
    
  }#i
  
})#model_code

# ---- 3. PREP FOR MCMC ----

# Constants & data
parameters <- c("mu", "sigma_year", "sigma_obs", "year_effect")

nimble_constants <- list(n_obs = nrow(urchins),
                         n_years = length(sort(unique(urchins$year))),
                         year_index = match(urchins$year, sort(unique(urchins$year))))

nimble_data <- list(log_count = urchins$log_count)

# MCMC settings
ni <- 50000
nb <- 10000
nt <- 40
nc <- 3

# ---- 4. RUN MCMC ----

mcmc_output <- nimbleMCMC(code      = model_code,
                          data      = nimble_data,
                          constants = nimble_constants,
                          monitors  = parameters,
                          niter     = ni,
                          nburnin   = nb,
                          nchains   = nc,
                          thin      = nt,
                          summary   = TRUE,
                          samplesAsCodaMCMC = TRUE)

# ---- 5. PROCESS OUTPUT ----
attach.nimble(mcmc_output$samples)
summary(mcmc_output)
print(mcmc_output$summary)
mcmcplot(mcmc_output$samples)

# Normalize the year effects (back-transform from log-scale)
year_mean_recruit <- exp(year_effect)  # back to count scale
year_mean_norm <- sweep(year_mean_recruit, 1, rowMeans(year_mean_recruit), "/")
norm_sd <- apply(year_mean_norm, 1, sd)

# Summary of normalized SD
hist(norm_sd, main = "Posterior of Normalized Year-to-Year SD",
     xlab = "Normalized SD")
mean(norm_sd)
quantile(norm_sd, probs = c(0.025, 0.5, 0.975))


#Plot posterior density of CV
norm_sd_mcmc <- as.mcmc(norm_sd)
plot(density(norm_sd_mcmc),
     main = "Posterior density of recruitment norm_sd",
     xlab = "Norm_sd")

densplot(mcmc_output$samples[, "sigma_obs"],
         main = "Posterior density of normalized recruitment variance\n(exp(σ²) - 1)",
         xlab = "Normalized variance")

# - - - - - - - - - - - - - - - End of section - - - - - - - - - - - - - - - -











# #3. log-transform (with +1 pseudocount)
# urchins <- urchins %>%
#   mutate(log_count = log(count_mean + 1)) %>%
#   mutate(log_sd = log(count_sd + 1))
# 
# 
# #calc yr2yr var ----
# urchins %>% 
#   ungroup %>% 
#   summarise(across(log_count:log_sd, list(mean = mean, sd = sd)))
# 
# #4. Nimble code (Let's make Josh proud)
# 
# code <- nimbleCode({
#   
#   #Priors
#   mu           ~ dnorm(mean = 0, sd = 10) # Intercept on log scale
#   beta         ~ dnorm(mean = 0, sd = 1) # Any linear trend
#   sigma        ~ T(dt(mu = 0, sigma = 1, df = 1), 0,) #Half Cauchy on log-scale SD 
#   
#   mean.recruitment ~
#     
#   annual.sigma ~
#   
#   for(n in 1:N.year){
#     year.randomeff[n] ~ dnorm(mean = mean.recruitment, sd = annual.sigma)
#     
#     survey.year.eff[n] ~ #specify a prior
#   }
# 
#   #Observation/process model
#   for(i in 1:N.survey) {
#     
#     # exp.urchins[i] <- mu + beta * year[i]
#     # 
#     # log_y[i] ~ dnorm(mean=exp.urchins[i], sd=sigma) ###ASK JOSH
#     # 
#     
#     log(exp.count[i]) <- year.randomeff[year[i]] + log(hectares.surveyed[i]) + survey.random.eff[i]
#     
#     survey.obs[i] ~ dpois(exp.count[i])
#     
#   }#i
#   
# })#code
# 
# #5. Data prep + nimble constants
# parameters  <- c("mu","beta","sigma")
# 
# nimble.constants <- list(N.year = nrow(urchins), 
#                          year = urchins$year)
# 
# nimble.data      <- list(log_y = urchins$log_n)
# 
# ni <- 40000
# nb <- 20000
# nc <- 3
# nt <- 40
# 
# #6. Run MCMC
# 
# mcmc.out <- nimbleMCMC(code       = code,
#                        data       = nimble.data,
#                        constants  = nimble.constants,
#                        monitors   = parameters,
#                        niter      = ni,
#                        nburnin    = nb,
#                        thin       = nt,
#                        nchains    = nc,
#                        summary    = TRUE,
#                        samplesAsCodaMCMC = TRUE)
# 
# #Attach sample object
# attach.nimble(mcmc.out$samples)
# 
# #7. Results and plotting
# summary(mcmc.out)
# print(mcmc.out$summary)
# mcmcplot(mcmc.out$samples)
# 
# #Trace log sigma
# plot(mcmc.out$samples[, "sigma"], 
#      main = "Trace of σ (yearly SD, log‐scale)")
# 
# #Density of log sigma
# dens <- density(as.numeric(sigma))
# plot(dens, main = "Posterior density of sigma (σ)")
# 
# # - - - - - - - - - - - - - - - End of section - - - - - - - - - - - - - - - -
# 
# 
# 
# #-------- Normalized year to year variance (SD) of incoming recruits --------
# 
# #8. Extract all sigma (σ) posterior samples (across chains)
# sigma_samples <- do.call(c, lapply(mcmc.out$samples, function(chain) {
#   as.numeric(chain[,"sigma"])
#   }))
# 
# #Transform to CV (coefficient of variation) on the original (count) scale
# cv_samples <- sqrt(exp(sigma_samples^2) - 1) # CV = sqrt(exp(σ^2) - 1)
# 
# #Posterior summaries for CV
# mean_cv <- mean(cv_samples)
# ci_cv   <- quantile(cv_samples, probs = c(0.025, 0.975))
# 
# print(c(mean_cv, ci_cv))
# 
# #Plot posterior density of CV
# cv_mcmc <- as.mcmc(cv_samples)
# plot(density(cv_mcmc),
#      main = "Posterior density of recruitment CV",
#      xlab = "Coefficient of Variation")

# ANOTEHR APPORACH
# 
# # 2. group by year & sum within transects
# # here we treat each (year, site, zone, replicate) as one transect
# transect_totals <- urchins_raw %>%
#   group_by(year, site_code, zone, replicate) %>%
#   summarise(count = sum(count), .groups="drop")
# 
# # 3. Compute annual mean recruitment
# annual_stats <- transect_totals %>%
#   group_by(year) %>%
#   summarise(mean_count = mean(count),
#             sd_count   = sd(count),
#             .groups    = "drop") %>%
#   # log‐transform with +1 pseudocount
#   mutate(log_y  = log(mean_count + 1), # center year for numerical stability
#          year_c = year - mean(year))
# 
# # 4. Nimble code (Let's make Josh proud)
# nim_code <- nimbleCode({
#   # priors for trend intercept and slope
#   alpha ~ dnorm(0, sd = 10)
#   beta  ~ dnorm(0, sd =  1)
#   
#   # prior on residual precision -> convert to SD
#   tau   ~ dgamma(0.001, 0.001)
#   sigma <- 1 / sqrt(tau)
#   
#   # likelihood: log‐mean around trend
#   for(i in 1:N.year) {
#     
#     # process model
#     mu[i]    <- alpha + beta * year_c[i]
#     
#     # observation model
#     log_y[i] ~ dnorm(mu[i], tau = tau)
#     
#   }#i
#   
#   # derived: normalized variance on original scale
#   var_norm <- exp(sigma * sigma) - 1 #Var(X)/E(X)^2 = exp(σ^2) - 1 for a lognormal
#   
# })#nim_code
# 
# # 5. Data prep + nimble constants
# parameters <- c("alpha", "beta", "sigma", "var_norm")
# 
# nimble.data <- list(log_y  = annual_stats$log_y)
# 
# nimble.constants <- list(N.year = nrow(annual_stats),
#                          year_c = annual_stats$year_c)
# 
# n_iter    <- 20000
# n_burnin  <- 2000
# n_chains  <- 3
# n_thin    <- 10
# 
# # 6. Run MCMC
# mcmc.output <- nimbleMCMC(code       = nim_code,
#                           data       = nimble.data,
#                           constants  = nimble.constants,
#                           monitors   = parameters,
#                           niter      = n_iter,
#                           nburnin    = n_burnin,
#                           thin       = n_thin,
#                           nchains    = n_chains,
#                           summary    = TRUE,
#                           samplesAsCodaMCMC = TRUE)
# 
# # attach sample objects (alpha, beta, sigma, var_norm, etc.)
# attach.nimble(mcmc.output$samples)
# 
# # 7. Summarize & Plot
# summary(mcmc.output)
# mcmcplot(mcmc.output$samples)
# print(mcmc.output$summary)
# 
# # Plot posterior density of var_norm
# densplot(mcmc.output$samples[, "var_norm"],
#          main = "Posterior density of normalized recruitment variance\n(exp(σ²) - 1)",
#          xlab = "Normalized variance")
# 
# # RESULT: normalize variance = 0.111295




##### ANOTHER APPROACH:
# # 1. Filter for taxa 66,67,70 and drop unwanted sites
# urchins_raw <- pisco %>%
#   filter(count_classcode %in% c("66","67","70"), # data 2001 - 2011
#          !site_code %in% c("CMEN00","CMES00","IFBRXX","KHLX00",
#                            "PSGX00","SHCX00","TRHX00","IPTAXX")) # 14 sites
# 
# # Compute annual mean counts
# annual_counts <- urchins_raw %>%
#   group_by(year) %>%
#   summarise(count = mean(count, na.rm = TRUE)) %>%
#   arrange(year)
#     # # Plot it
#     # ggplot(annual_counts, aes(x = as.character(year), y = count)) +
#     #   geom_point() +
#     #   labs(title = "Annual Sea Urchin Recruitment",
#     #        y = "Recruitment Count", x = "Year")
# 
# # Add small offset ε to avoid log(0)
# epsilon <- 1e-6
# annual_counts <- annual_counts %>%
#   mutate(count_off = count + epsilon)
#     # # Plot it
#     # ggplot(annual_counts, aes(x = as.character(year), y = count)) +
#     #   geom_point() +
#     #   labs(title = "Annual Sea Urchin Recruitment",
#     #        y = "Recruitment Count_off", x = "Year")
# 
# 
# # 2. Define nimble model
# urchin_code <- nimbleCode({
#   
#   # Hyperpriors
#   mu_global    ~ dnorm(0, sd = 10)             # prior on log-scale global mean
#   sigma_year   ~ T(dt(0, 1, df=1), 0,)         # half-Cauchy
#   sigma_obs    ~ T(dt(0, 1, df=1), 0,)         # half-Cauchy
#   
#   # Precision parameters
#   tau_year <- pow(sigma_year, -2)
#   tau_obs  <- pow(sigma_obs, -2)
#   
#   # Process model - Year‐level means
#   for(j in 1:n_years) {
#     log_mu[j] ~ dnorm(mu_global, tau = tau_year)
#   }#j
#   
#   # Observation model
#   for(i in 1:n_counts) {
#     log(count_off[i]) ~ dnorm(log_mu[year_id[i]], tau = tau_obs)
#   }#i
# 
# # Deterministic transformations to natural scale
# mean_recruit <- exp(mu_global + 0.5 * sigma_year^2)
# sd_recruit   <- sqrt((exp(sigma_year^2) - 1) *
#                        exp(2 * mu_global + sigma_year^2))
# 
# })#urchin_code
# 
# # 3. Data, constants, and inits for nimble
# 
# nimble_data <- list(count_off = annual_counts$count_off)
# nimble_constants <- list(n_counts = nrow(annual_counts),
#                          n_years  = nrow(annual_counts),
#                          year_id  = 1:nrow(annual_counts))
# inits <- function() {
#   list(mu_global  = log(mean(annual_counts$count_off)),
#        sigma_year = sd(log(annual_counts$count_off)),
#        sigma_obs  = sd(log(annual_counts$count_off)),
#        log_mu     = log(annual_counts$count_off))
# }#inits
# parameters <- c("mu_global", "sigma_year", "sigma_obs", "mean_recruit", "sd_recruit")
# 
# # 4. MCMC settings & run
# n_iter    <- 20000
# n_burnin  <- 2000
# n_chains  <- 3
# n_thin    <- 10
# 
# mcmc_out <- nimbleMCMC(code       = urchin_code,
#                        data       = nimble_data,
#                        constants  = nimble_constants,
#                        inits      = inits,
#                        monitors   = parameters,
#                        niter      = n_iter,
#                        nburnin    = n_burnin,
#                        nchains    = n_chains,
#                        thin       = n_thin,
#                        summary    = TRUE,
#                        samplesAsCodaMCMC = TRUE)
# 
# # 5. Quick diagnostics
# print(mcmc_out$summary)        # posterior summaries
# mcmcplot(mcmc_out$samples)     # trace & density plots
# 
# # 9. (Optional) Save posterior summary for downstream use
# post_summary <- as_tibble(mcmc_out$summary, rownames = "Parameter")
# write_csv(post_summary, "recruitment_posterior_summary_v2.csv")
# 
# # 6. Post‐processing to natural scale
# # Extract posterior samples
# samps <- as.matrix(mcmc_out$samples)
# 
# # Overall mean on natural scale
# post_mu      <- exp(samps[,"mu_global"])
# mean_mu      <- mean(post_mu)
# ci_mu        <- quantile(post_mu, c(0.025, 0.975))
# 
# # Normalized variance and coefficient of variation
# # For log-normal: Var_norm = exp(sigma^2) - 1; CV = sqrt(Var_norm)
# post_var_norm <- exp(samps[,"sigma_year"]^2) - 1
# post_cv_norm  <- sqrt(post_var_norm)
# 
# res <- tibble(
#   parameter   = c("Mean recruitment", "Normalized var", "Normalized CV"),
#   median      = c(median(post_mu), median(post_var_norm), median(post_cv_norm)),
#   CI_lower    = c(ci_mu[1], quantile(post_var_norm,0.025), quantile(post_cv_norm,0.025)),
#   CI_upper    = c(ci_mu[2], quantile(post_var_norm,0.975), quantile(post_cv_norm,0.975)))
# 
# print(res)
# 
# # # Optionally save results
# # write_csv(res, "urchin_recruitment_posterior_summary.csv")



##### ONE APPROACH:
# # 1. Prepare the PISCO data
# summary(pisco)
# 
# # Filter for sea urchin species (classcodes 66, 67, 70), and Oregon sites
# urchins_raw <- pisco %>%
#   filter(count_classcode %in% c("66", "67", "70"), # data 2001 - 2011
#          !site_code %in% c("CMEN00", "CMES00", "IFBRXX", "KHLX00", 
#                            "PSGX00", "SHCX00", "TRHX00", "IPTAXX")) # 14 sites
# summary(urchins_raw)
# summary(urchins_raw$count)
# 
# # Aggregate counts by year (main model dataset)
# annual_counts <- urchins_raw %>%
#   group_by(year) %>%
#   summarise(count = mean(count, na.rm = TRUE)) %>% #is this SUM or MEAN?
#   arrange(year)
# 
# # Plot counts by year
# ggplot(annual_counts, aes(x = as.character(year), y = count)) +
#   geom_point() +
#   labs(title = "Annual Sea Urchin Recruitment",
#        y = "Recruitment Count", x = "Year")
# 
# 
# # 3. Define the hierarchical log-normal model using robust priors
# urchin_code <- nimbleCode({
#   
#   # Priors for hyperparameters
#   mu_global   ~ dt(mu = 0, sigma = 1, df = 1)
#   sigma_year  ~ T(dt(mu=0, sigma=1, df=1), 0,) #Half Cauchy
#   sigma_obs   ~ T(dt(mu=0, sigma=1, df=1), 0,) #Half Cauchy
#   
#   tau_year   <- pow(sigma_year, -2)
#   tau_obs    <- pow(sigma_obs, -2)
#   
#   # Year-specific log-means (Proces model)
#   for (j in 1:n_years) {
#     log_mu[j] ~ dnorm(mean = mu_global, sd = sigma_year)
#     }#j
#   
#   # Observation model
#   for (i in 1:n_counts) {
#     log(count[i]) ~ dnorm(mean = log_mu[year_id[i]], sd = sigma_obs)
#     }#i
#   
# })#urchin_code
# 
# # 4. Prepare parameters, data, and constants for nimble
# parameters <- c("mu_global", "sigma_year", "sigma_obs", "log_mu")
# 
# nimble.data <- list(count = annual_counts$count)
# 
# nimble.constants <- list(n_counts = length(annual_counts$count),
#                          n_years  = length(annual_counts$year),
#                          year_id  = match(annual_counts$year, sort(unique(annual_counts$year))))
# 
# # 5. Initial values
# gen_inits <- function() {
#   list(mu_global   = log(mean(annual_counts$count)),
#        sigma_year  = sd(log(annual_counts$count)),
#        sigma_obs   = sd(log(annual_counts$count)),
#        log_mu      = log(annual_counts$count))
# }
# 
# # 7. MCMC settings
# ni <- 40000
# nb <- 2000
# nt <- 40
# nc <- 3
# 
# # 8. Run the model
# mcmc.output <- nimbleMCMC(code       = urchin_code,
#                           data       = nimble.data,
#                           constants  = nimble.constants,
#                           inits      = gen_inits(),
#                           monitors   = parameters,
#                           niter      = ni,
#                           nburnin    = nb,
#                           nchains    = nc,
#                           thin       = nt,
#                           summary    = TRUE,
#                           samplesAsCodaMCMC = TRUE)
# 
# # 9. Summarize and visualize
# library(coda)
# summary(mcmc.output)
# mcmcplot(mcmc.output$samples)
# 
# # Posterior densities
# plot(density(as.matrix(mcmc.output$samples)[, "mu_global"]), main = "Posterior of mu_global", xlab = "mu_global")
# plot(density(as.matrix(mcmc.output$samples)[, "sigma_year"]), main = "Posterior of sigma_year", xlab = "sigma_year")
# plot(density(as.matrix(mcmc.output$samples)[, "sigma_obs"]), main = "Posterior of sigma_obs", xlab = "sigma_obs")
# 
# # 10. Extract annual means
# log_mu_post <- as.matrix(mcmc.output$samples)[, grep("log_mu", colnames(as.matrix(mcmc.output$samples)))]
# mu_years <- exp(apply(log_mu_post, 2, mean))
# 
# annual_recruitment <- data.frame(
#   year = sort(unique(years)),
#   mean_recruitment = mu_years
# )

