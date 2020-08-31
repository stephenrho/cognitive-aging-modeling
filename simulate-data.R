
## this script simulates the data set used in Greene & Rhodes
## "A Tutorial on Cognitive Modeling for Cognitive Aging Research"
## the simulated data conforms to the parsimonious SDT model of 
## Selker et al. (2019): https://link.springer.com/article/10.3758/s13428-019-01231-3

# settings
J = 48 # total sample (each group = J/2)
L = 80 # total trials (each condition = L/2)
K = 6 # number of rating categories

stopifnot(L%%4 == 0)
stopifnot(J%%2 == 0)

gen_pars = function(J, mu = c(1,1,1,1), tau=c(1,1,1,1), rho = c(1,1,1,1,1,1)){
  # This function can be used to generate individual level parameter values 
  # from a multivariate normal.  Effects at the group level only can be specified
  # by setting tau to 0 and associated rhos to 0. Only allows for 4 coefs (e.g., group by condition)
  
  cor_mat = diag(length(tau))
  
  cor_mat[1,2] = cor_mat[2,1] = rho[1]
  cor_mat[1,3] = cor_mat[3,1] = rho[2]
  cor_mat[1,4] = cor_mat[4,1] = rho[3]
  cor_mat[2,3] = cor_mat[3,2] = rho[4]
  cor_mat[2,4] = cor_mat[4,2] = rho[5]
  cor_mat[3,4] = cor_mat[4,3] = rho[6]
  
  cov_mat = diag(tau) %*% cor_mat %*% diag(tau)
  
  p = MASS::mvrnorm(n = J, mu = mu, Sigma = cov_mat)
  
  return(p)
}

ratingModel <- function(d, s, a, b, K, ratingProbs = T){
  # the parsimonius model from Selker et al
  unb = -log((1 - 1:K/K)/(1:K/K))
  
  c_k = a*unb + b
  
  f = (diff(pnorm(c(-Inf, c_k),0,1)))
  h = (diff(pnorm(c(-Inf, c_k),d,s)))
  if (!ratingProbs){
    # for plotting ROCs
    f = cumsum(rev(f))
    h = cumsum(rev(h))
    f = f[1:(K-1)]
    h = h[1:(K-1)]
  }
  
  return(cbind(f, h))
}

# make data frame
dat = expand.grid(trial = 1:(L/2), cond = c(-1,1), id = 1:J)

# code trial type, item, and group
dat$signal = ifelse(dat$trial > L/4, 1, 0)
dat$item = rep(sample(x = 1:(L/4)), J*2)
dat$group = ifelse(dat$id > J/2, -1, 1)

## parameter settings

# d
d_i = rnorm(1:L, mean = 0, sd = .2) # item effect
d_j = gen_pars(J, 
               mu = c(.5, -.3, .2,.2), # c(intercept, group, cond, group:cond)
               tau = c(.3, 0, .3, 0), 
               rho = c(0, .4, 0, 0, 0, 0)) # (1,2) (1,3) (1,4) (2,3) (2,4) (3,4)

# s
s_i = rnorm(1:L, mean = 0, sd = .1) # item effect
s_j = gen_pars(J, 
               mu = c(.18, 0, .05, 0), # av ~= 1.2
               tau = c(.1, 0, .1, 0), 
               rho = c(0, 0, 0, 0, 0, 0))

# a (scale)
a_i = rnorm(1:L, mean = 0, sd = .1) # item effect
a_j = gen_pars(J, 
               mu = c(0, 0, 0, 0), # av ~= 1
               tau = c(.3, 0, 0, 0), 
               rho = c(0, 0, 0, 0, 0, 0))

# b (shift)
b_i = rnorm(1:L, mean = 0, sd = .3) # item effect
b_j = gen_pars(J, 
               mu = c(1, .3, 0, 0), 
               tau = c(.1, 0, 0, 0), 
               rho = c(0, 0, 0, 0, 0, 0))

# to multiple individual coefs by 
X = model.matrix(~ group + cond + cond:group, data = dat)

dat$d = exp( rowSums(d_j[dat$id,] * X) + d_i[dat$item] )
dat$s = exp( rowSums(s_j[dat$id,] * X) + s_i[dat$item] )
dat$a = exp( rowSums(a_j[dat$id,] * X) + a_i[dat$item] )
dat$b = ( rowSums(b_j[dat$id,] * X) + b_i[dat$item] )

with(dat, hist(d))
with(dat, hist(s))
with(dat, hist(a))
with(dat, hist(b))

# simulate data...
dat$rating = NA
for (i in 1:nrow(dat)){
  p = with(dat[i,], ratingModel(d = d, s = s, a = a, b = b, K = K))[, dat$signal[i]+1 ]
  dat$rating[i] = sample(1:K, size = 1, replace = F, prob = p)
}

barplot(table(dat$rating))
barplot(table(dat$signal, dat$rating), beside = T)

par(mfrow=c(1,2))
with(subset(dat, group==1), barplot(table(signal, rating), beside = T) )
with(subset(dat, group==-1), barplot(table(signal, rating), beside = T) )
par(mfrow=c(1,1))

# save data
# the full simulation
write.csv(x = dat, file = "data/sim.csv", row.names = F)

# clean data for example
dat$cond[dat$cond==-1] = "A"
dat$cond[dat$cond==1] = "B"

dat$group[dat$group==-1] = "Y"
dat$group[dat$group==1] = "O"

dat=dat[, c("id", "trial", "group", "cond", "item", "signal", "rating")]

write.csv(x = dat, file = "data/example-data.csv", row.names = F)
