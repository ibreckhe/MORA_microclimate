##Functions to support community regression against environmental data.
##Author: Ian Breckheimer.

####Function to fit a Bayesian mixed-effects model in JAGS.####
fit.jags.mixed.tree <- function(x,y,species,plot,nsamples=10000,thin=50){
  require(rjags)
  
  ##Removes negative and NA values from x, y and group vectors.
  data <- data.frame(x,y,species,plot)
  data_complete <- data[complete.cases(data),]
  
  ##Prepares input.
  group <- as.numeric(as.factor(data_complete$species))
  ngroups <- length(unique(group))
  plot <- as.numeric(data_complete$plot)
  nplots <- length(unique(plot))
  
  x <- data_complete$x
  y <- as.numeric(data_complete$y > 0)
  n <- length(x)
  
  ##Fits model.
  cat(
    "
    model{
    # priors
    height.mu ~ dnorm(0,0.01)
    height.sigma ~ dunif(0,5)
    height.tau <-  pow(height.sigma, -2)
    width.mu ~ dnorm(-1.5,100)
    width.sigma ~ dunif(0.001,1)
    width.tau <- pow(width.sigma,-2)
    opt.mu ~ dnorm(0,0.001)
    opt.sigma ~ dunif(0.01,20)
    opt.tau <- pow(opt.sigma,-2)
    plot.sigma ~ dunif(0.01,20)
    plot.tau <- pow(opt.sigma,-2)
    
    for (j in 1:ngroups){
    opt.g[j] ~ dnorm(opt.mu, opt.tau)T(-3.5,3.5)
    width.g[j] ~ dnorm(width.mu, width.tau)T(,-0.2)
    height.g[j] ~ dnorm(height.mu,height.tau)
    }
    
    for (k in 1:nplots){
    height.plot[k] ~ dnorm(0,plot.tau)
    }
    
    # likelihood
    for (i in 1:n){
    y[i] ~ dbern(p[i])
    p[i] <- 1 / (1 + exp(-z[i]))
    z[i] <- width.g[group[i]] * (x[i] - opt.g[group[i]])^2 + height.g[group[i]] + height.plot[plot[i]]
    }
    
    
    }  
    ", file="jagsmodel_tree_vertexform.txt"
  )
  #   inits <-expression(list(height=rnorm(1,0,1),
  #                           width=rnorm(1,-10,1),
  #                           opt.mu=rnorm(1,0,1),
  #                           opt.sigma=runif(1,0.01,20),
  #                           b.opt=rnorm(1,0,1)))
  #   init.list <- list(eval(inits),eval(inits),eval(inits))
  jd <- list(x=x, y=y, n=n,group=group,ngroups=ngroups,plot=plot,nplots=nplots)
  mod <- jags.model("jagsmodel_tree_vertexform.txt", data= jd, n.chains=3, n.adapt=1000)
  update(mod,n.iter=nsamples)
  out <- coda.samples(mod, c("height.mu","opt.mu","opt.sigma","width.mu","width.sigma","opt.g","width.g","height.g"),
                      n.iter=nsamples,thin=thin)
  #diags <- gelman.diag(out)
  return(list(mod=mod,out=out))
  }

####Function to fit a Bayesian mixed-effects model in JAGS.####
fit.jags.mixed.2var <- function(x1,x2,y,species,plot,nsamples=10000,thin=50){
  require(rjags)
  
  ##Removes negative and NA values from x, y and group vectors.
  data <- data.frame(x1,x2,y,species,plot)
  data_complete <- data[complete.cases(data),]
  
  ##Prepares input.
  group <- as.numeric(as.factor(data_complete$species))
  ngroups <- length(unique(group))
  plot <- as.numeric(data_complete$plot)
  nplots <- length(unique(plot))
  
  x1 <- data_complete$x1
  x2 <- data_complete$x2
  y <- as.numeric(data_complete$y > 0)
  n <- length(x1)
  
  ##Fits model.
  cat(
    "
    model{
    # priors
    height.mu ~ dnorm(0,0.01)
    height.sigma ~ dunif(0,5)
    height.tau <-  pow(height.sigma, -2)
    width.mu ~ dnorm(-1.5,0.1)T(,-0.1)
    width.sigma ~ dunif(0.001,10)
    width.tau <- pow(width.sigma,-2)
    opt.mu ~ dnorm(0,0.001)
    opt.sigma ~ dunif(0.01,100)
    opt.tau <- pow(opt.sigma,-2)

    width.mu.x2 ~ dnorm(-0.2,0.1)T(,-0.01)
    width.sigma.x2 ~ dunif(0.1,10)
    width.tau.x2 <- pow(width.sigma.x2,-2)
    opt.mu.x2 ~ dnorm(0,0.001)
    opt.sigma.x2 ~ dunif(0.001,4)
    opt.tau.x2 <- pow(opt.sigma.x2,-2)

    int.mu ~ dnorm(0,0.001)
    int.sigma ~ dunif(0.1,10)
    int.tau <- pow(int.sigma,-2)

    plot.sigma ~ dunif(0.01,20)
    plot.tau <- pow(opt.sigma,-2)
    
    for (j in 1:ngroups){
      opt.g[j] ~ dnorm(opt.mu, opt.tau)T(-4.5,4.5)
      width.g[j] ~ dnorm(width.mu, width.tau)T(,-0.1)
      height.g[j] ~ dnorm(height.mu,height.tau)
      opt.g.x2[j] ~ dnorm(opt.mu.x2, opt.tau.x2)
      width.g.x2[j] ~ dnorm(width.mu.x2, width.tau.x2)T(,-0.01)
      int.g[j] ~ dnorm(int.mu,int.tau)
    }
    
    for (k in 1:nplots){
      height.plot[k] ~ dnorm(0,plot.tau)
    }
    
    # likelihood
    for (i in 1:n){
      y[i] ~ dbern(z[i])
      logit(z[i]) <- width.g[group[i]] * (x1[i] - opt.g[group[i]])^2 + 
                     width.g.x2[group[i]] * (x2[i] - opt.g.x2[group[i]])^2  + 
                     int.g[group[i]] * x1[i] * x2[i] +
                     height.g[group[i]] + height.plot[plot[i]]
    }
    
    }  
    ", file="jagsmodel_2var_vertexform.txt"
  )
  #   inits <-expression(list(height=rnorm(1,0,1),
  #                           width=rnorm(1,-10,1),
  #                           opt.mu=rnorm(1,0,1),
  #                           opt.sigma=runif(1,0.01,20),
  #                           b.opt=rnorm(1,0,1)))
  #   init.list <- list(eval(inits),eval(inits),eval(inits))
  jd <- list(x1=x1, x2=x2,y=y, n=n,group=group,ngroups=ngroups,plot=plot,nplots=nplots)
  mod <- jags.model("jagsmodel_2var_vertexform.txt", data= jd, n.chains=3, n.adapt=1000)
  update(mod,n.iter=nsamples)
  out <- coda.samples(mod, c("height.mu","opt.mu","opt.sigma","width.mu","width.sigma","opt.mu.x2","opt.sigma.x2","width.mu.x2",
                             "width.sigma.x2","int.mu","int.sigma","opt.g","opt.g.x2","width.g","width.g.x2","int.g","height.g"),
                      n.iter=nsamples,thin=thin)
  #diags <- gelman.diag(out)
  return(list(mod=mod,out=out))
  }

####Gets credible intervals for the parameters of interest.####
get.param.creds <- function(out,params=c("width.mu","opt.mu","height","opt.b1"),intervals=c(0.025,0.5,0.975)){
  gg_out <- ggmcmc::ggs(out)
  samples <- matrix(NA,ncol=length(params),nrow=(length(out[[1]])*length(out)))
  for (i in 1:length(params)){
    samples[,i] <- subset(gg_out,Parameter==params[i])$value
  }
  cred.fun <- function(x){quantile(x,probs=intervals)}
  creds <- apply(samples,FUN=cred.fun,MARGIN=2)
  colnames(creds) <- params
  return(creds)
}

##Logit and antilogit function.####
antilogit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}

##Function to parse JAGS output into a matrix of functions representing response curves
##for each species.
make_fun_matrix <- function(jags.out,par.names=c("width.g","opt.g","height.g"),
                            n.samples=1000){
  gg.out <- ggmcmc::ggs(jags.out)
  ngroups <- sum(grepl(paste(par.names[1],"[",sep=""),unique(gg.out$Parameter),fixed=TRUE))
  samples <- sample(1:dim(jags.out[[1]])[1],size=n.samples,replace=FALSE)
  chains <- sample(1:3,size=n.samples,replace=TRUE)
  gg.sample <- gg.out[gg.out$Iteration %in% samples & gg.out$Chain == chains,]
  gg.tbl <- dplyr::as.tbl(gg.sample)
  fun.list <- replicate(ngroups*n.samples,function(x){x})
  fun.matrix <- matrix(fun.list,ncol=ngroups)
  print("Making function matrix for species ")
  for(j in 1:ngroups){
    print(paste(j,"of ",ngroups))
    pname1 <- paste(par.names[1],"[",j,"]",sep="")
    pname2 <- paste(par.names[2],"[",j,"]",sep="")
    pname3 <- paste(par.names[3],"[",j,"]",sep="")
    p1.gr <- gg.sample[gg.sample$Parameter == pname1,4]
    p2.gr <- gg.sample[gg.sample$Parameter == pname2,4]
    p3.gr <- gg.sample[gg.sample$Parameter == pname3,4]
    np <- length(p3.gr)
    fun.body <- paste(rep("antilogit(",np),p1.gr,rep(" * (x1 - ",np),p2.gr,rep(")^2 + ",np),
                      p3.gr,rep(")",np),sep="")
    for (i in 1:n.samples){
      body(fun.matrix[i,j][[1]]) <- parse(text=fun.body[i])
    }
  }
  return(fun.matrix)
}

make_fun_matrix_2var <- function(jags.out,par.names=c("width.g","opt.g","height.g","width.g.x2","opt.g.x2","int.g"),
                            n.samples=1000){
  gg.out <- ggmcmc::ggs(jags.out)
  ngroups <- sum(grepl(paste(par.names[1],"[",sep=""),unique(gg.out$Parameter),fixed=TRUE))
  samples <- sample(1:dim(jags.out[[1]])[1],size=n.samples,replace=FALSE)
  chains <- sample(1:3,size=n.samples,replace=TRUE)
  gg.sample <- gg.out[gg.out$Iteration %in% samples & gg.out$Chain == chains,]
  fun.list <- replicate(ngroups*n.samples,function(x1,x2){x1})
  fun.matrix <- matrix(fun.list,ncol=ngroups)
  print("Making function matrix for species ")
  for(j in 1:ngroups){
    print(paste(j,"of ",ngroups))
    pname1 <- paste(par.names[1],"[",j,"]",sep="")
    pname2 <- paste(par.names[2],"[",j,"]",sep="")
    pname3 <- paste(par.names[3],"[",j,"]",sep="")
    pname4 <- paste(par.names[4],"[",j,"]",sep="")
    pname5 <- paste(par.names[5],"[",j,"]",sep="")
    pname6 <- paste(par.names[6],"[",j,"]",sep="")
    p1.gr <- gg.sample[gg.sample$Parameter == pname1,4]
    p2.gr <- gg.sample[gg.sample$Parameter == pname2,4]
    p3.gr <- gg.sample[gg.sample$Parameter == pname3,4]
    p4.gr <- gg.sample[gg.sample$Parameter == pname4,4]
    p5.gr <- gg.sample[gg.sample$Parameter == pname5,4]
    p6.gr <- gg.sample[gg.sample$Parameter == pname6,4]
    np <- length(p6.gr)
    fun.body <- paste(rep("antilogit(",np),p1.gr,rep(" * (x1 - ",np),p2.gr,rep(")^2 + ",np),
                      p4.gr,rep(" * (x2 - ",np),p5.gr,rep(")^2 + ",np),
                      p6.gr,rep(" * x1 * x2 + ",np),p3.gr,rep(")",np),sep="")
    for (i in 1:n.samples){
      body(fun.matrix[i,j][[1]]) <- parse(text=fun.body[i])
    }
  }
  return(fun.matrix)
}

make_med_fun_matrix_2var <- function(jags.out,par.names=c("width.g","opt.g","height.g",
                                                          "width.g.x2","opt.g.x2","int.g")){
  gg.out <- ggmcmc::ggs(jags.out)
  ngroups <- sum(grepl(paste(par.names[1],"[",sep=""),unique(gg.out$Parameter),fixed=TRUE))
  fun.list <- replicate(ngroups,function(x1,x2){x1})
  fun.matrix <- matrix(fun.list,ncol=ngroups)
  print("Making function matrix for species ")
  for(j in 1:ngroups){
    print(paste(j,"of ",ngroups))
    pname1 <- paste(par.names[1],"[",j,"]",sep="")
    pname2 <- paste(par.names[2],"[",j,"]",sep="")
    pname3 <- paste(par.names[3],"[",j,"]",sep="")
    pname4 <- paste(par.names[4],"[",j,"]",sep="")
    pname5 <- paste(par.names[5],"[",j,"]",sep="")
    pname6 <- paste(par.names[6],"[",j,"]",sep="")
    p1.gr <- quantile(gg.out[gg.out$Parameter == pname1,4],probs=c(0.5))
    p2.gr <- quantile(gg.out[gg.out$Parameter == pname2,4],probs=c(0.5))
    p3.gr <- quantile(gg.out[gg.out$Parameter == pname3,4],probs=c(0.5))
    p4.gr <- quantile(gg.out[gg.out$Parameter == pname4,4],probs=c(0.5))
    p5.gr <- quantile(gg.out[gg.out$Parameter == pname5,4],probs=c(0.5))
    p6.gr <- quantile(gg.out[gg.out$Parameter == pname6,4],probs=c(0.5))
    np <- length(p6.gr)
    fun.body <- paste(rep("antilogit(",np),p1.gr,rep(" * (x1 - ",np),p2.gr,rep(")^2 + ",np),
                      p4.gr,rep(" * (x2 - ",np),p5.gr,rep(")^2 + ",np),
                      p6.gr,rep(" * x1 * x2 + ",np),p3.gr,rep(")",np),sep="")
    body(fun.matrix[1,j][[1]]) <- parse(text=fun.body)
  }
  return(fun.matrix)
}



##Builds a function to multiply all of the functions in a given column based on a vector.
make_combofun <- function(fun.matrix,iteration=1){
  ncol <- dim(fun.matrix)[2]
  combo_fun <- function(x, nx=length(x),pres_vec=rep(FALSE,ncol)) {x}
  body <- c("(")
  for(i in 1:(dim(fun.matrix)[2]-1)){
    body <- paste(body, "ifelse(rep(pres_vec[",i,"],nx),fun.matrix[",iteration,",",i,
                  "][[1]](x),1-fun.matrix[",iteration,",",i,"][[1]](x))*",sep="")
  }
  body <- paste(body,"ifelse(rep(pres_vec[",i+1,"],nx),fun.matrix[",iteration,",",i+1,
                "][[1]](x),1-fun.matrix[",iteration,",",i+1,"][[1]](x)))",sep="")
  body(combo_fun) <- parse(text=body)
  return(combo_fun)
}

##Builds a function to multiply all of the functions in a given column based on a vector.
make_combofun_2var <- function(fun.matrix,iteration=1){
  ncol <- dim(fun.matrix)[2]
  combo_fun <- function(x1,x2, nx=length(x1),pres_vec=rep(FALSE,ncol)) {x1}
  body <- c("(")
  for(i in 1:(dim(fun.matrix)[2]-1)){
    body <- paste(body, "ifelse(rep(pres_vec[",i,"],nx),fun.matrix[",iteration,",",i,
                  "][[1]](x1,x2),1-fun.matrix[",iteration,",",i,"][[1]](x1,x2))*",sep="")
  }
  body <- paste(body,"ifelse(rep(pres_vec[",i+1,"],nx),fun.matrix[",iteration,",",i+1,
                "][[1]](x1,x2),1-fun.matrix[",iteration,",",i+1,"][[1]](x1,x2)))",sep="")
  body(combo_fun) <- parse(text=body)
  return(combo_fun)
}

##Numerically integrates the functions so that they are probability densities.
combo_dens <- function(x,combo_fun,pres_vec) { 
  y <- combo_fun(x,nx=length(x),pres_vec=pres_vec)
  yi <- integrate(combo_fun, -Inf, +Inf, pres_vec=pres_vec)
  return(y/yi[[1]])
}

##Samples from the joint probability functions using rejection sampling.
sample_combo <- function(newx1,newx2,combo_fun,pres_vec){
  pred.prob <- combo_fun(newx1,newx2,pres_vec,nx=length(newx))
  pred.sample <- sample(newx,size=1,prob=pred.prob)
  return(pred.sample)
}

##Samples from the joint probability functions using rejection sampling.
sample_combo_2var <- function(newx1,newx2,combo_fun,pres_vec){
  
  ## Binds newx1 and newx2.
  newx <- cbind(newx1,newx2)
  
  ## Generates predicted probability.
  pred.prob <- combo_fun(newx1,newx2,pres_vec,nx=length(newx1))
  
  ## Samples rows with probability p.
  newx_accept <- newx[sample(1:dim(newx)[1],size=1,prob=pred.prob,replace=FALSE),]

  return(newx_accept)
}

##Plot parameter estimates
plot_cit_jags_out <- function(data,jags.out,scaled_elev,scaled_cit,cat_factor,color_pal){
  ggd <- ggmcmc::ggs(jags.out)
  a1 <- ggd$value[which(ggd$Parameter == "alpha[1]")]
  a2 <- ggd$value[which(ggd$Parameter == "alpha[2]")]
  a3 <- ggd$value[which(ggd$Parameter == "alpha[3]")]
  a4 <- ggd$value[which(ggd$Parameter == "alpha[4]")]
  b1 <- ggd$value[which(ggd$Parameter == "beta.elev[1]")]
  b2 <- ggd$value[which(ggd$Parameter == "beta.elev[2]")]
  b3 <- ggd$value[which(ggd$Parameter == "beta.elev[3]")]
  b4 <- ggd$value[which(ggd$Parameter == "beta.elev[4]")]
  d <- data.frame(a1, a2, a3, a4, b1, b2, b3, b4)
  
  cat_numeric <- as.numeric(cat_factor)
  
  ####Predicts based on the fit model.
  d.x.pred1 <- seq(min(elev[cat_numeric==1]),max(elev[cat_numeric==1]), length.out=dim(d)[1])
  d.x.pred2 <- seq(min(elev[cat_numeric==2]),max(elev[cat_numeric==2]), length.out=dim(d)[1])
  d.x.pred3 <- seq(min(elev[cat_numeric==3]),max(elev[cat_numeric==3]), length.out=dim(d)[1])
  d.x.pred4 <- seq(min(elev[cat_numeric==4]),max(elev[cat_numeric==4]), length.out=dim(d)[1])
  
  d.predfun.1 <- function(x){quantile(d$a1 + d$b1*x,probs=c(0.025,0.25,0.5,0.75,0.975))}
  d.predfun.2 <- function(x){quantile(d$a2 + d$b2*x,probs=c(0.025,0.25,0.5,0.75,0.975))}
  d.predfun.3 <- function(x){quantile(d$a3 + d$b3*x,probs=c(0.025,0.25,0.5,0.75,0.975))}
  d.predfun.4 <- function(x){quantile(d$a4 + d$b4*x,probs=c(0.025,0.25,0.5,0.75,0.975))}
  
  
  d.pred.1 <- sapply(d.x.pred1,FUN=d.predfun.1)
  d.pred.2 <- sapply(d.x.pred2,FUN=d.predfun.2)
  d.pred.3 <- sapply(d.x.pred3,FUN=d.predfun.3)
  d.pred.4 <- sapply(d.x.pred4,FUN=d.predfun.4)
  
  ####Unscales predictions.
  elev_scale <- attr(scaled_elev,"scaled:scale")
  elev_center <- attr(scaled_elev,"scaled:center")
  cit_scale <- attr(scaled_cit,"scaled:scale")
  cit_center <- attr(scaled_cit,"scaled:center")
  
  d.x.pred1.un <- d.x.pred1 * elev_scale + elev_center
  d.x.pred2.un <- d.x.pred2 * elev_scale + elev_center
  d.x.pred3.un <- d.x.pred3 * elev_scale + elev_center
  d.x.pred4.un <- d.x.pred4 * elev_scale + elev_center
  
  d.pred1.un <- d.pred.1 * cit_scale + cit_center
  d.pred2.un <- d.pred.2 * cit_scale + cit_center
  d.pred3.un <- d.pred.3 * cit_scale + cit_center
  d.pred4.un <- d.pred.4 * cit_scale + cit_center
  
  cols <- color_pal
  cols_bg <- adjustcolor(cols,alpha.f=0.2)
  cols_ln <- adjustcolor(cols,alpha.f=0.6)
  
  plot(CIT_mean~MORA_elev_3m,data=data,col=cols[dry_snow_cat],xlab="Elevation (m)",
       ylab="Community Inferred Temp. (C)",type="n",xlim=c(500,1700),ylim=c(3.5,9))
  points(data$MORA_elev_3m,data$CIT_mean,pch=20,col=cols_bg[cat_numeric],cex=0.5)
  arrows(x0=data$MORA_elev_3m,x1=data$MORA_elev_3m,
         y0=data$CIT_lwr_quart.25.,y1=data$CIT_upr_quart.75.,
         code=3,angle=90,length=0,col=cols_bg[cat_numeric],lwd=0.5)
  polygon(x=c(d.x.pred1.un,rev(d.x.pred1.un)),y=c(d.pred1.un[1,],rev(d.pred1.un[5,])),border=NA,
          col=cols_bg[1])
  polygon(x=c(d.x.pred2.un,rev(d.x.pred2.un)),y=c(d.pred2.un[1,],rev(d.pred2.un[5,])),border=NA,
          col=cols_bg[2])
  polygon(x=c(d.x.pred3.un,rev(d.x.pred3.un)),y=c(d.pred3.un[1,],rev(d.pred3.un[5,])),border=NA,
          col=cols_bg[3])
  polygon(x=c(d.x.pred4.un,rev(d.x.pred4.un)),y=c(d.pred4.un[1,],rev(d.pred4.un[5,])),border=NA,
          col=cols_bg[4])
  
  polygon(x=c(d.x.pred1.un,rev(d.x.pred1.un)),y=c(d.pred1.un[2,],rev(d.pred1.un[4,])),border=NA,
          col=cols_ln[1])
  polygon(x=c(d.x.pred2.un,rev(d.x.pred2.un)),y=c(d.pred2.un[2,],rev(d.pred2.un[4,])),border=NA,
          col=cols_ln[2])
  polygon(x=c(d.x.pred3.un,rev(d.x.pred3.un)),y=c(d.pred3.un[2,],rev(d.pred3.un[4,])),border=NA,
          col=cols_ln[3])
  polygon(x=c(d.x.pred4.un,rev(d.x.pred4.un)),y=c(d.pred4.un[2,],rev(d.pred4.un[4,])),border=NA,
          col=cols_ln[4])
  legend("bottomleft",title="Site Class",bty="n",legend=levels(cat_factor),
         pch=20,col=cols)
}

####Function to find the vertices of a rotated ellipse####
ellipse_majoraxis <- function(ellipse){
  center <- ellipse$loc
  cov <- ellipse$cov
  s1 <- sqrt(cov[1,1])
  s2 <- sqrt(cov[2,2])
  rho <- cov[1,2]/(s1*s2)
  phi <- if(s1==s2){
           (pi/4)*sin(rho)
        }else if(s1>s2){
           0.5*atan((2*rho*s1*s2)/(s1^2 - s2^2))
        }else if(s1<s2){
           0.5*atan((2*rho*s1*s2)/(s1^2 - s2^2)) - pi/2
        }else{print("S1 or S2 missing")}
  slope <- tan(ifelse(phi < 0, pi + phi, phi))
  int <- slope * (0 - center[1]) + center[2]
  return(c(intercept=int,slope=slope))
}
