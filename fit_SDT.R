
## this script fits the models discussed in Greene & Rhodes
## "A Tutorial on Cognitive Modeling for Cognitive Aging Research"
## make sure the working directory is set to the folder containing
## the 'models' and 'data' folders

LOAD = T 
# if you have fitted model objects saved, set this to T to load them and not re-run models
# if you don't have fitted models or want to re-run, set to F
# you can download fitted models from: https://drive.google.com/drive/folders/14gmtoYXKHMtZL7yjIzdmjKEhjIrmsGaq

# load the packages we use
library(rstan)
library(bridgesampling)
library(bayesplot)
library(HDInterval)

# rstan settings
rstan_options(auto_write = TRUE) # write the models so we don't have to recompile
options(mc.cores = 4) # can sample chains in parallel 
# use parallel::detectCores() to get number of cores available

### read the data ----
rdat = read.csv("data/example-data.csv")

head(rdat)

# plot frequency of each rating for signal and noise trials
with(rdat, barplot(table(signal, rating), beside = T))

# separate plots for each group
par(mfrow=c(1,2))
with(subset(rdat, group=="Y"), barplot(table(signal, rating), beside = T, main = "Y"))
with(subset(rdat, group=="O"), barplot(table(signal, rating), beside = T, main = "O"))
par(mfrow=c(1,1))

# code so older = 1 in the design matrix, X
rdat$group = as.factor(rdat$group)
contrasts(rdat$group) = c(1,0)

# put the data into list form for rstan
data_list = list(
  N = nrow(rdat), # number of observations (trials)
  y = rdat$rating, # response
  J = length(unique(rdat$id)), # number of participants
  id = rdat$id, # participant ids
  X = model.matrix(~ 1 + group, data = rdat), # predictors
  sig_trial = rdat$signal, # was this a signal trial? (1 = yes, 0 = no) 
  K = 6 # number of rating categories
)

### ### ### ### ### ### ### ### ### ### 
### fit the signal detection models ###
### ### ### ### ### ### ### ### ### ###

# specify things here so we can change them later
nchains = 4
nwarm = 500
niter = 2000

### model 1: ----
# age differences allowed for all parameters

if (!LOAD | !file.exists("models/SDT_m1_fit.rds")){
  SDT_m1_fit <- stan(
    file = "models/SDT_m1.stan",
    data = data_list,
    chains = nchains,
    warmup = nwarm,
    iter = niter,
    # the lines below exclude the stan parameters listed from 
    # being saved. Otherwise the model object is v large
    pars=c("d", "s", "a", "b", "c", "theta"),
    include=F
  )
  # save the model object for later
  saveRDS(SDT_m1_fit, file = "models/SDT_m1_fit.rds")
} else {
  # load
  SDT_m1_fit = readRDS("models/SDT_m1_fit.rds")
}

# these values should be close to 1
# if not the model has not converged
# https://mc-stan.org/docs/2_24/reference-manual/convergence.html
rhat(SDT_m1_fit, pars=c("B_d", "B_a", "B_b", "B_s",
                        "tau_d", "tau_a", "tau_b", "tau_s", "lp__"))

# plot the posterior samples
# medians and credible intervals
plot(SDT_m1_fit, pars=c("B_d", "B_a", "B_b", "B_s",
                    "tau_d", "tau_a", "tau_b", "tau_s"))

# steps of the chains
traceplot(SDT_m1_fit, pars=c("B_d", "B_a", "B_b", "B_s",
                         "tau_d", "tau_a", "tau_b", "tau_s"), 
          inc_warmup=T) # this will show the warm up period

# extract particular parameters and plot histograms
hist(extract(SDT_m1_fit, pars="B_s[1]")[[1]])
hist(exp(extract(SDT_m1_fit, pars="B_s[1]")[[1]]))

# extract age parameter for d
age_d = extract(SDT_m1_fit, pars="B_d[2]")[[1]]

hist(age_d, breaks = 30, main=bquote(Beta[1]^"(d)"), xlab="")

# convert back to d scale to compare groups
d_young = exp(extract(SDT_m1_fit, pars="B_d[1]")[[1]])
d_old = exp(apply(extract(SDT_m1_fit, pars="B_d")[[1]], 1, sum))
hist(d_young, xlim=c(.5,3), ylim=c(0, 2000), col="lightblue", main="Younger=blue; Older=pink")
hist(d_old, col="pink", add=T)

hist(d_old-d_young, main="Age difference in units of d")

# posterior predictive checks
yrep <- extract(SDT_m1_fit, pars = "y_rep")[[1]] # posterior simulations from the generated quantities block
y <- rdat$rating

# plot only 100 simulations from the posterior as takes a long time to do all
ppc_bars_grouped(y, yrep[sample(1:nrow(yrep), size = 100),], group = rdat$group)
ppc_bars_grouped(y, yrep[sample(1:nrow(yrep), size = 100),], 
                 group = paste0("group = ", rdat$group, 
                                "; trial = ", rdat$signal))

# by condition, A/B
# remember condition isn't in this model (see model 3)
# so predictions don't differ between conditions, which leads to mis-fit
ppc_bars_grouped(y, yrep[sample(1:nrow(yrep), size = 100),], 
                 group = paste0("condition = ", rdat$cond,
                                "; group = ", rdat$group, 
                                "; trial = ", rdat$signal))

### model 2: ----
# fix d so there is no age difference

if (!LOAD | !file.exists("models/SDT_m2_fit.rds")){
  SDT_m2_fit <- stan(
    file = "models/SDT_m2.stan",
    data = data_list,
    chains = nchains,
    warmup = nwarm,
    iter = niter,
    pars=c("d", "s", "a", "b", "c", "theta"),
    include=F
  )
  
  # save
  saveRDS(SDT_m2_fit, file = "models/SDT_m2_fit.rds")
} else{
  SDT_m2_fit = readRDS("models/SDT_m2_fit.rds")
}

plot(SDT_m2_fit, pars=c("B_d", "B_a", "B_b", "B_s",
                    "tau_d", "tau_a", "tau_b", "tau_s"))

# posterior predictions from model 2
yrep2 <- extract(SDT_m2_fit, pars = "y_rep")[[1]]

ppc_bars_grouped(y, yrep2[sample(1:nrow(yrep2), size = 100),], 
                 group = paste0("group = ", rdat$group, 
                                "; trial = ", rdat$signal))

ppc_bars_grouped(y, yrep2[sample(1:nrow(yrep), size = 100),], 
                 group = paste0("condition = ", rdat$cond,
                                "; group = ", rdat$group, 
                                "; trial = ", rdat$signal))

# compare models 1 and 2
# use the bridgesampling package to estimate log likelihood
if (!LOAD | !file.exists("models/SDT_m1_ll.rds")){
  SDT_m1_ll = bridge_sampler(SDT_m1_fit)
  SDT_m2_ll = bridge_sampler(SDT_m2_fit)
  
  saveRDS(SDT_m1_ll, file = "models/SDT_m1_ll.rds")
  saveRDS(SDT_m2_ll, file = "models/SDT_m2_ll.rds")
} else{
  SDT_m1_ll = readRDS("models/SDT_m1_ll.rds")
  SDT_m2_ll = readRDS("models/SDT_m2_ll.rds")
}

bayes_factor(SDT_m1_ll, SDT_m2_ll)

### ### ### ### ### ### ### ### ### ### 
### other modifications of model 1  ###
### ### ### ### ### ### ### ### ### ###

# go back to the original design matrix for m1
data_list$X = model.matrix(~ 1 + group, data = rdat)

### model 1.2 ----
# model the effect of items on d

data_list$item = rdat$item
data_list$M = length(unique(rdat$item))

if (!LOAD | !file.exists("models/SDT_m1.2_fit.rds")){
  SDT_m1.2_fit <- stan(
    file = "models/SDT_m1.2.stan",
    data = data_list,
    chains = nchains,
    warmup = nwarm,
    iter = niter,
    pars=c("d", "s", "a", "b", "c", "theta"),
    include=F
  )
  
  # save 
  saveRDS(SDT_m1.2_fit, file = "models/SDT_m1.2_fit.rds")
} else {
  SDT_m1.2_fit = readRDS("models/SDT_m1.2_fit.rds")
}

item_d = as.array(SDT_m1.2_fit, pars="alpha_d")

# plot the item effects (ordered by median)
mcmc_intervals(item_d[,,order(apply(item_d, 3, median), decreasing = T)]) # medians and error bars
mcmc_areas(item_d[,,order(apply(item_d, 3, median), decreasing = T)]) # density plots

# extract the coefficient for age difference in d (log scale)
age_d1.2 = extract(SDT_m1.2_fit, pars="B_d[2]")[[1]]

# compare to original model
plot(density(age_d), main=bquote(Beta[1]^"(d)"), xlab="", lwd=2, xlim=c(-1,.2), ylim=c(0,4))
lines(density(age_d1.2), lwd=2, col="red")

h1 = hdi(age_d)
h2 = hdi(age_d1.2)
segments(x0 = h1["lower"], x1 = h1["upper"], y0 = .25, y1 = .25, lwd = 3)
segments(x0 = h2["lower"], x1 = h2["upper"], y0 = .5, y1 = .5, lwd = 3, col="red")
points(x = c(median(h1), median(h2)), y = c(.25, .5), pch=16, cex=2, col=c("black", "red"))

legend("topright", legend = c("no item effect", "item effect"), text.col = c("black", "red"), bty='n')

### model 1.3 ----
# allow the groups to differ in between-participant variability for d

# get the group of each participant
grps = rdat$group[!duplicated(rdat$id)]
# convert to 0s (Y) and 1s (O)
data_list$group = as.integer(grps == "O")

if (!LOAD | !file.exists("models/SDT_m1.3_fit.rds")){
  SDT_m1.3_fit <- stan(
    file = "models/SDT_m1.3.stan",
    data = data_list,
    chains = nchains,
    warmup = nwarm,
    iter = niter,
    pars=c("d", "s", "a", "b", "c", "theta"),
    include=F
  )
  
  # save
  saveRDS(SDT_m1.3_fit, file = "models/SDT_m1.3_fit.rds")
} else{
  SDT_m1.3_fit = readRDS("models/SDT_m1.3_fit.rds")
}

plot(SDT_m1.3_fit, pars=c("B_d", "B_a", "B_b", "B_s",
                          "tau_d", "tau_a", "tau_b", "tau_s"))

plot(SDT_m1.3_fit, pars=c("tau_d"))

# extract estimates of between-participant variability (SD)
tau_d = extract(SDT_m1.3_fit, pars=c("tau_d"))[[1]]
# column 1 = younger, 2 = older

# plot the estimates and their difference
# par(mfrow=c(1,2))
plot(density(tau_d[,1]), lwd=2, main=bquote(tau^"(d)"), xlab="", xlim=c(0, .9))
lines(density(tau_d[,2]), lwd=2, col="blue")
legend("topright", legend=c("younger", "older"), text.col=c("black", "blue"), bty="n")

plot(density(apply(tau_d, 1, diff)), lwd=2, 
     main=bquote(tau[old]^"(d)" - tau[young]^"(d)"), xlab="")

### model 1.4 ----
# correlate another measure with individual differences (i.e. random effects) in d

# for this we have to read in the other measure
scores = read.csv("data/cor-scores.csv") # one score for each person

data_list$score = scores$score

if (!LOAD | !file.exists("models/SDT_m1.4_fit.rds")){
  SDT_m1.4_fit <- stan(
    file = "models/SDT_m1.4.stan",
    data = data_list,
    chains = nchains,
    warmup = nwarm,
    iter = niter,
    pars=c("d", "s", "a", "b", "c", "theta"),
    include=F
  )
  
  # save
  saveRDS(SDT_m1.4_fit, file = "models/SDT_m1.4_fit.rds")
} else{
  SDT_m1.4_fit = readRDS("models/SDT_m1.4_fit.rds")
}

plot(SDT_m1.4_fit, pars="Sigma")

# extract the correlation samples
rho = extract(SDT_m1.4_fit, pars="Sigma[1,2]")[[1]]

# posterior mean and median
mean(rho); median(rho)
# and 95% credible intervals
quantile(rho, probs = c(0.025, .975))

cor.test(x = apply(extract(SDT_m1_fit, pars="b_d")[[1]], 2, mean), 
         y = scores$score)

### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### model 3: ----
# model the effect of condition (+ individual level effect)
# and interaction between age group and condition for d

# modify the group level design matrix (X) and add one for individual level (Z)
data_list$X = model.matrix(~ group + cond + group:cond, data = rdat)
data_list$Z = model.matrix(~ cond, data = rdat)

SDT_m3_fit <- stan(
  file = "models/SDT_m3.stan",
  data = data_list,
  chains = nchains,
  warmup = nwarm,
  iter = niter,
  pars=c("d", "s", "a", "b", "c", "theta"),
  include=F
)

# save
saveRDS(SDT_m3_fit, file = "models/SDT_m3_fit.rds")
#SDT_m3_fit = readRDS("models/SDT_m3_fit.rds")

plot(SDT_m3_fit, pars=c("B_d", "B_a", "B_b", "B_s",
                        "tau_d", "tau_a", "tau_b", "tau_s"))

traceplot(SDT_m3_fit, pars=c("B_d", "B_a", "B_b", "B_s",
                        "tau_d", "tau_a", "tau_b", "tau_s"))

# extract the group level effects for d
B_d = extract(SDT_m3_fit, pars="B_d")[[1]]

# convert back to d scale for both groups and conditions
# condition A was coded 0 and B coded 1
# group Y was coded 0 and group O coded 1
# therefore, we can multiply the matrix of B_d samples
# by particular vectors that code for group and condition

d_youngA = exp( B_d %*% c(1,0,0,0) )
d_youngB = exp( B_d %*% c(1,0,1,0) )
d_oldA = exp( B_d %*% c(1,1,0,0) )
d_oldB = exp( B_d %*% c(1,1,1,1) )

# par(mfrow=c(1,2))
plot(density(d_youngA), lwd=2, xlim=c(.5,3.5), ylim=c(0,5), main="d", xlab="")
lines(density(d_youngB), lwd=2, lty=3)
lines(density(d_oldA), lwd=2, col="blue")
lines(density(d_oldB), lwd=2, col="blue", lty=3)
legend("topright", legend = c("Younger", "Older", "A", "B"), 
       lty = c(NA, NA, 1, 3), lwd=c(NA, NA, 2, 2),
       text.col = c("black", "blue", "black", "black"), bty="n")

plot(density(d_youngA - d_oldA), lwd=2, xlim=c(-1.5,2.5), ylim=c(0,3), main="d - group differences\nby condition", xlab="")
lines(density(d_youngB - d_oldB), lwd=2, lty=3)
lines(density((d_youngA - d_oldA) - (d_youngB - d_oldB)), lwd=2, col='red')
legend("topleft", legend = c("A", "B", "Difference (A - B)"), 
       lty = c(1, 3, 1), lwd=c(2, 2, 2), col=c("black", "black", "red"), bty="n")

