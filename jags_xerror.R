# simulate covariate data
n <- 100
sdx <- 10
pop_taux <- 1 / (sdx * sdx)
sdobs <- 8
taux <- 1 / (sdobs * sdobs)
truex <- rnorm(n, 0, sdx)
errorx <- rnorm(n, 0, sdobs)
obsx <- truex + errorx

# simulate response data
alpha <- 0
beta <- 10

sdy <- 10
tauy <- 1 / (sdy * sdy)
errory <- rnorm(n, 0, sdy)
obsy <- alpha + beta*truex + errory
parms <- data.frame(alpha, beta)

par(mfrow=c(1,1))
plot(obsx,obsy,col="red")
points(truex,obsy,col="blue")

# bundle data
jags_d <- list(x = obsx, y = obsy, n = length(obsx))

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
out <- coda.samples(mod1, n.iter=1000, thin=1,
                    variable.names=c("alpha", "beta", "sdy"))

# store parameter estimates
require(ggmcmc)
ggd <- ggs(out)
a <- ggd$value[which(ggd$Parameter == "alpha")]
b <- ggd$value[which(ggd$Parameter == "beta")]
d <- data.frame(a, b)

# specify model
cat("
    model {
    ## Priors
    alpha ~ dnorm(0, .001)
    beta ~ dnorm(0, .001)
    sdy ~ dunif(0,100)
    tauy <- 1 / (sdy * sdy)
    taux ~ dnorm(pop_taux,0.01)
    
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
jags_d <- list(x = obsx, y = obsy, obs_taux = taux, pop_taux = pop_taux, 
               n = length(obsx))

# initiate model
mod2 <- jags.model("xyerror.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)

# simulate posterior
out <- coda.samples(mod2, n.iter=30000, thin=30,
                    variable.names=c("alpha", "beta", "sdy","taux","xtrue[1]"))
# store parameter estimates
ggd <- ggs(out)
a2 <- ggd$value[which(ggd$Parameter == "alpha")]
b2 <- ggd$value[which(ggd$Parameter == "beta")]
d2 <- data.frame(a2, b2)


ggplot(d, aes(x=obsx, obsy)) +
  geom_abline(aes(intercept=a, slope=b), data=d, color="red", alpha=0.05) +
  geom_abline(aes(intercept=a2, slope=b2), data=d2, color="blue", alpha=0.05) +
  geom_abline(aes(intercept=alpha, slope=beta),
              data=parms, color="green", size=1.5, linetype="dashed") +
  theme_bw() +
  geom_point(shape=1, size=3,color="red") +
  geom_point(aes(x=truex,y=obsy),shape=1,size=3,color="blue")+
  xlab("X values") + ylab("Observed Y values") +
  ggtitle("Model results with and without modeling error in X")