# simulate covariate data
n <- 100
sdx <- 5
pop_taux <- 1 / (sdx * sdx)
sdobs <- 10
taux <- 1 / (sdobs * sdobs)
truex <- rnorm(n, 0, sdx)
errorx <- rnorm(n, 0, sdobs)
obsx <- truex + errorx

# simulate response data
alpha <- -60
beta <- -5

# process error
p_sdy <- 10
p_tauy <- 1 / (p_sdy * p_sdy)
p_errory <- rnorm(n, 0, p_sdy)

# measurement error
obs_sdy <- 20
obs_tauy <- 1 / (obs_sdy * obs_sdy)
obs_errory <- rnorm(n,0,obs_sdy)

truey <- alpha + beta*truex + p_errory
obsy <- truey + obs_errory

par(mfrow=c(1,1))
plot(obsx,obsy,col="red")
points(truex,truey,col="blue")

# bundle data
jags_d <- list(x = obsx, y = obsy, n = length(obsx))

#### Ordinary linear regression. ####

# write model
cat("
    model{
    ## Priors
    alpha ~ dnorm(0, .001)
    beta ~ dnorm(0, .001)
    sdy ~ dunif(0, 100)
    tauy <- 1 / (sdy * sdy)
    
    ## Likelihood
    for (i in 1:n){
    mu[i] <- alpha + beta * x[i]
    y[i] ~ dnorm(mu[i], tauy)
    }
    }
    ",
    fill=TRUE, file="yerror.txt")

require(rjags)
# initiate model
mod1 <- jags.model("yerror.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)

# simulate posterior
out1 <- coda.samples(mod1, n.iter=1000, thin=1,
                    variable.names=c("alpha", "beta", "sdy"))

# store parameter estimates
require(ggmcmc)
ggd <- ggs(out1)
a <- ggd$value[which(ggd$Parameter == "alpha")]
b <- ggd$value[which(ggd$Parameter == "beta")]
d <- data.frame(a, b)

#### Model with observation error in x.

# specify model
cat("
    model {
    ## Priors
    alpha ~ dnorm(0, .001)
    beta ~ dnorm(0, .001)
    sdy ~ dnorm(0,0.1)T(0,)
    tauy <- 1 / (sdy * sdy)
    sdx ~ dnorm(pop_sdx,0.1)T(0,)
    taux <- 1 / (sdx * sdx)
    
    ## Likelihood
    for (i in 1:n){
    xtrue[i] ~ dnorm(0,taux)
    x[i] ~ dnorm(xtrue[i],obs_taux)
    y[i] ~ dnorm(mu[i], tauy)
    mu[i] <- alpha + beta * xtrue[i]
    }
    }
    ", fill=T, file="xyerror.txt")

# bundle data
jags_d <- list(x = obsx, y = obsy, obs_taux = taux, pop_sdx = sdx, 
               n = length(obsx))

# initiate model
mod2 <- jags.model("xyerror.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)

# simulate posterior
out2 <- coda.samples(mod2, n.iter=30000, thin=30,
                    variable.names=c("alpha", "beta", "sdy","taux","xtrue[1]"))
# store parameter estimates
ggd <- ggs(out2)
a2 <- ggd$value[which(ggd$Parameter == "alpha")]
b2 <- ggd$value[which(ggd$Parameter == "beta")]
d2 <- data.frame(a2, b2)

#### Model with measurement error in both X and y. ####

# specify model
cat("
    model {
    ## Priors
    alpha ~ dnorm(0, .001)
    beta ~ dnorm(0, .001)
    sdy ~ dnorm(pop_sdy,0.1)T(0,)
    tauy <- 1 / (sdy * sdy)
    sdx ~ dnorm(pop_sdx,0.1)T(0,)
    taux <- 1 / (sdx * sdx)
    
    ## Likelihood
    for (i in 1:n){
    xtrue[i] ~ dnorm(0,taux)
    x[i] ~ dnorm(xtrue[i], obs_taux)
    y[i] ~ dnorm(ytrue[i], obs_tauy)
    ytrue[i] ~ dnorm(mu[i],tauy)
    mu[i] <- alpha + beta * xtrue[i]
    }
    }
    ", fill=T, file="xyerror_ymeas.txt")

# bundle data
jags_d <- list(x = obsx, y = obsy, obs_taux = taux, pop_sdy=p_sdy,
               pop_sdx=sdx,obs_tauy=obs_tauy,n = length(obsx))

# initiate model
mod3 <- jags.model("xyerror_ymeas.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)

# simulate posterior
out3 <- coda.samples(mod3, n.iter=30000, thin=30,
                    variable.names=c("alpha", "beta","sdx","sdy"))
# store parameter estimates
ggd <- ggs(out3)
a3 <- ggd$value[which(ggd$Parameter == "alpha")]
b3 <- ggd$value[which(ggd$Parameter == "beta")]
d3 <- data.frame(a3, b3)


ggplot(d, aes(x=obsx, obsy)) +
  geom_point(size=2,color="grey60") +
  geom_abline(aes(intercept=a, slope=b), data=d, color="red", alpha=0.01) +
  geom_abline(aes(intercept=a2, slope=b2), data=d2, color="blue", alpha=0.01) +
  geom_abline(aes(intercept=a3, slope=b3), data=d3, color="green", alpha=0.01) +
  geom_abline(aes(intercept=alpha, slope=beta),
              color="black", size=1.5, linetype="dashed") +
  theme_bw() +
  geom_point(aes(x=truex,y=truey),size=2,color="black")+
  xlab("X") + ylab("Y") +
  guides(color=guide_legend(title="Legend"))+
  ggtitle("Model results with and without modeling error in X and Y")