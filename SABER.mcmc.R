#===========================================================================
# SABER.mcmc.R uses MCMC sampling of parameter of interest to retrieve the water components 
#along with bathymetry and bottom reflectance (for shallow water) given the initial vals and Rrs.

#The estimates are generated as:

#1. Create the prior distribution of parameter of interest, i.e. P(theta).
#2. Create the conditional log-likelihood (P(Data|Theta) of the model noise.
#3. Estimate Posterior; P(theta|Data) ~= log(P(theta)) + log(P(Data|Theta)
#4. Optimize the posterior function using MCMC in Metropolis-Hastings framework
#5. Extract estimates from the full P(theta|Data) using central tendancies.

#Author: Mr. Soham Mukherjee - PhD Student, Aquatel, UQAR
#===========================================================================
#Calcluate mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#-----------------------------------------------------------------------------------------------
## Find Prior distribution
#Fit distributions to prior parameters
#-----------------------------------------------------------------------------------------------
bgc.data <- read.csv("S:/Work/UQAR/Datasets/WISE/in-Situ/L3/biogeochemistry_wiseman.csv",
                     header = T)
acdom.data <- read.csv("S:/Work/UQAR/Datasets/WISE/in-Situ/L3/ag_long_wiseman.csv",
                       header = T)
anap.data <- read.csv("S:/Work/UQAR/Datasets/WISE/in-Situ/L3/ad_long_wiseman.csv",
                      header = T)
chl.sample <- bgc.data$Chl
#hist(chl.sample, probability = T)

acdom.440.sample <- acdom.data$ag[acdom.data$wavelength == "440"]
#hist(acdom.440.sample, probability = T)

anap.440.sample <- anap.data$ad[anap.data$wavelength == "440"]
#hist(anap.440.sample, probability = T)

fit.chl.norm <- fitdistrplus::fitdist(chl.sample, "weibull") #try to fit dist
#plot(fit.chl.norm) #Infer results of the fit

fit.acdom440.norm <- fitdistrplus::fitdist(acdom.440.sample, "weibull") #try to fit dist
#plot(fit.acdom440.norm) #Infer results of the fit

fit.anap440.norm <- fitdistrplus::fitdist(anap.440.sample, "weibull") #try to fit dist
#plot(fit.anap440.norm) #Infer results of the fit


#-----------------------------------------------------------------------------------------------------
#Function for prior
#-----------------------------------------------------------------------------------------------------
prior = function(param){
  chl = param[1]
  acdom.440 = param[2]
  anap.440 = param[3]
  chl.prior = dweibull(x=chl,shape = fit.chl.norm$estimate[1], 
                       scale =fit.chl.norm$estimate[2] , log = T)
  
  acdom440.prior = dweibull(x=acdom.440,shape = fit.acdom440.norm$estimate[1], 
                            scale =fit.acdom440.norm$estimate[2] , log = T)
  
  anap440.prior = dweibull(x=anap.440,shape = fit.anap440.norm$estimate[1], 
                           scale =fit.anap440.norm$estimate[2] , log = T)
  
  return(chl.prior+acdom440.prior+anap440.prior)
}

#-----------------------------------------------------------------------------------------------------
#Function for Likelihood
#-----------------------------------------------------------------------------------------------------
likelihood = function(pars, dataobs, batch){
  
  # Values predicted by the forward model for single RUN
  if (batch == FALSE) {
    Gpred = Saber_forward(chl = pars[1], acdom440 = pars[2],
                          anap440 =pars[3],
                          #bbp.550 = Fit.input$bbp.550,
                          verbose=F ,realdata = dataobs)
  }
  # Gpred = Saber_forward(chl = pars[1], acdom440 = pars[2], 
  #                       anap440 =pars[3], 
                          #bbp.550 = Fit.input$bbp.550,
  #                       verbose=F ,realdata = dataobs)
  
  #Values predicted by the forward model for batch RUN
  if (batch == TRUE) {
    Gpred = Saber_forward(chl = pars[1], acdom440 = pars[2], 
                          anap440 =pars[3], bbp.550 = HL.deep.iop$bbp550[j],
                          verbose=F ,realdata = dataobs)
  }
  
  # log-likelihood 
  smull = sum(dnorm(x = dataobs, mean = Gpred[[1]]$Rrs, sd = 0.0001, log = TRUE)) #Single RUN
  #smull = sum(dnorm(x = dataobs, mean = Gpred[[1]]$Rrs, sd = 0.13129354/100, log = TRUE)) #IOCCG RUN
  return(smull)
}
#-----------------------------------------------------------------------------------------------------
#Function for Posterior
#-----------------------------------------------------------------------------------------------------
posterior = function(param, rrs.actual, lklhood.max, maxpar, batch.run){
  if (lklhood.max == TRUE) {
    return (likelihood(pars = maxpar, dataobs = rrs.actual,batch = batch.run ) + prior(param))
  } else {
    return (likelihood(pars = param, dataobs = rrs.actual, batch = batch.run ) + prior(param))
  }
  
}
#-----------------------------------------------------------------------------------------------------
#Function for MCMC using Metropolis-Hastings
#-----------------------------------------------------------------------------------------------------
# mse = startvalue/(41-3)
# proposalfunction = function(param){
#   return(abs(rnorm(3,mean = param, 
#                sd= c(sd(chl.sample),sd(acdom.440.sample),sd(anap.440.sample))
#                #sd= c(mse[1],mse[2],mse[3])
#                #sd= c(0.5,0.25,0.01)
#                ))
#   
#   )
# }

proposalfunction = function(param){
  priopart <- rweibull(param,
           shape = as.numeric(c(fit.chl.norm$estimate[1], 
                                fit.acdom440.norm$estimate[1],
                                fit.anap440.norm$estimate[1])),
           scale = as.numeric(c(fit.chl.norm$estimate[2], 
                                fit.acdom440.norm$estimate[2],
                                fit.anap440.norm$estimate[2]))
  )
  lklpart <- rnorm(3,mean = param, 
                   sd= c(sd(chl.sample),sd(acdom.440.sample),sd(anap.440.sample))
                   #sd= c(mse[1],mse[2],mse[3])
                   #sd= c(0.5,0.25,0.01)
  )
  #target.distribution <- (priopart+lklpart)/2
  target.distribution <- priopart
  
  return(target.distribution)
}
#Create the Progress Bar
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = 5000, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")

run_metropolis_MCMC = function(startvalue, iterations, rrs.actual, lklhood.max, maxpar, batchrun){
  chain = array(dim = c(iterations+1,3))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    probab = exp(posterior(proposal,rrs.actual = rrs.actual, lklhood.max, maxpar,batch.run = batchrun) - posterior(chain[i,],rrs.actual = rrs.actual, lklhood.max, maxpar, batch.run = batchrun))
    if (runif(1) < probab){
      chain[i+1,] = proposal
      #print(paste0(i," no. iteration: P(h+1)/P(h) > 1"))
    }else{
      chain[i+1,] = chain[i,]
      #print(paste0(i," no. iteration:P(h+1)/P(h) < 1"))
    }
    setTxtProgressBar(pb, i)
    #cat(paste0("\033[0;42m",i," iterations completed, ", (iterations - i), " remaining","\033[0m","\n"))
  }
  return(chain)
}
#-------------------------------------------------------------------------------------------------------
#Create function for applying the MCMC
#--------------------------------------------------------------------------------------------------------
require(coda)
mcmc.run  = function(rrs.actual, initials, iterations, burn_In, lklhood.max=FALSE, maxpar, batchrun)
{
  trueSd = 0.000001
  sampleSize = length(wavelength)
  
  # create independent x-values
  x = wavelength
  # create dependent values according to ax + b + N(0,sd)
  y = rrs.actual + rnorm(n=sampleSize,mean=0,sd=trueSd)
  #plot(x,y, main="Observed Data", type="l", col="red")
  
  startvalue = initials
  #startvalue = c(5,1,0.05)
  start.time <- Sys.time()
  chain = run_metropolis_MCMC(startvalue, iterations, rrs.actual = rrs.actual, lklhood.max, maxpar, batchrun = batchrun)
  end.time <- Sys.time()
  burnIn = burn_In
  acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
  return(list(data.frame(chain[-(1:burnIn),]),
              "acceptance"=acceptance, "runtime"=end.time - start.time))
}

#-------------------------------------------------------------------------------------------------------
#Apply MCMC to sample the posterior
#--------------------------------------------------------------------------------------------------------
startvalue = as.numeric(Fit.optimized.ssobj[1,1:3])
chain.op <- mcmc.run(rrs.actual = rrs.forward.am, iterations = 5000, 
                   initials = startvalue, burn_In = 2500,
                  maxpar = startvalue, lklhood.max = TRUE, batchrun = FALSE  )
chain1 <- chain.op[[1]]


startvalue = c(3,0.8,0.05)
chain.op <- mcmc.run(rrs.actual = rrs.forward.am, iterations = 5000, 
                     initials = startvalue, burn_In = 2500,
                     maxpar = startvalue, lklhood.max = TRUE, batchrun = FALSE  )
chain2 <- chain.op[[1]]

#-------------------------------------------------------------------------------------------------------
#Create plot to visualize result of MCMC
#--------------------------------------------------------------------------------------------------------
#png(filename = "Posterior.SABER.png", units = "in", width = 6.5, height = 4.5, res=300)
dev.new()
par(mfrow = c(2,3))
hist(chain[,1],nclass=30, main="Posterior of chl[mg^1m-3]", 
     xlab="Bayes = green; MLE =blue; actual=red", probability = T 
     #xlim=c(min(chl.sample), max(chl.sample)) 
)
lines(density(chain[,1]), col="blue", lwd=2)
abline(v = mean(chain[,1]), col="seagreen", lwd=2)
abline(v = startvalue[1], col="blue", lwd=2 )
abline(v = Fit.input$chl, col="red", lwd=2 )

hist(chain[,2],nclass=30, main="Posterior of acdom440[m^-1]",
     xlab="Bayes = green; MLE =blue; actual=red", probability = T 
     #xlim=c(min(acdom.440.sample), max(acdom.440.sample))
)
lines(density(chain[,2]), col="blue", lwd=2)
abline(v = mean(chain[,2]),col="seagreen", lwd=2)
abline(v = startvalue[2], col="blue", lwd=2 )
abline(v = Fit.input$acdom.440, col="red", lwd=2 )

hist(chain[,3],nclass=30, main="Posterior of anap440[m^-1]", 
     xlab="Bayes = green; MLE =blue; actual=red", probability = T
     #xlim=c(min(anap.440.sample), max(anap.440.sample))
)
lines(density(chain[,3]), col="blue", lwd=2)
abline(v = mean(chain[,3]) ,col="seagreen", lwd=2)
abline(v = startvalue[3], col="blue", lwd=2 )
abline(v = Fit.input$anap.440, col="red", lwd=2 )

plot(chain1[,1], type = "l", xlab="Bayes = green; actual = Red" , 
     main = "Chain of chl[mg^1m-3]")
lines(chain2[,2], col="red")
#abline(h = startvalue[1], col="red", lwd=2 )
abline(h = Fit.input$chl, col="red", lwd=2 )
abline(h =  mean(chain[,1]), col="seagreen", lwd=2 )

plot(chain[,2], type = "l", xlab="Bayes = green; actual = Red" ,
     main = "Chain of acdom440[m^-1]")
#abline(h = startvalue[2], col="red", lwd=2)
abline(h = Fit.input$acdom.440, col="red", lwd=2 )
abline(h =  mean(chain[,2]), col="seagreen" , lwd=2)

plot(chain[,3], type = "l", xlab="Bayes = green; actual = Red" , 
     main = "Chain of anap440[m^-1]")
#abline(h = startvalue[3], col="red", lwd=2 )
abline(h = Fit.input$anap.440, col="red", lwd=2 )
abline(h =  mean(chain[,3]), col="seagreen", lwd=2 )
#dev.off()

#Save results from MCMC chains 
Fit.optimized.mcmc <- data.frame("chl"=mean(chain[,1]), 
                                 "acdom.440"=mean(chain[,2]),
                                 "anap.440"=mean(chain[,3]))

# Example: plot the likelihood profile of a single parameter among the fit vector
slopevalues = function(x){return(likelihood(pars = c(x,as.numeric(Fit.optimized.ssobj[-1])), dataobs = rrs.forward.am))}
slopelikelihoods = lapply(seq(1, 10, by=.05), slopevalues )
plot (seq(1, 10, by=.05), slopelikelihoods , type="l", xlab = "values of chl", ylab = "Log likelihood")
abline(v=seq(1, 10, by=.05)[which.max(slopelikelihoods)])

#===============================================================================================




log_like_graph <- function(chl, acdom440){
  loglik = likelihood(pars = c(chl,acdom440,as.numeric(Fit.optimized.ssobj[-(1:2)])), dataobs = rrs.forward.am)
  return(loglik)
}
log_like_graph <- Vectorize(log_like_graph)

# Set grid of beta and sigma2 values 
chl <- seq(0.1,10, by =0.5)
acdom440 <- seq(0.1,10, by =0.5)
log_vals <- outer(chl, acdom440, log_like_graph)

persp(chl, acdom440, log_vals, theta = 7, phi =8 , r= 500)

x <- chl; y <- acdom440; z <- log_vals

# top/bot are coordinates(x,y,z) of the axis; pmat is output of persp()
# pos means the label is c(left(-1) or right(1), down(-1) or up(1)) of the axis

persp_lab_f <- function(top, bot, pmat, space, pos = c(-1, -1))
{
  coords_3d <- list(top = top, bot = bot)
  coords_2d <- lapply(coords_3d, function(v, pmat) unlist(trans3d(v[1], v[2], v[3], pmat)), pmat = pmat)
  coords_2d$mid <- (coords_2d$top + coords_2d$bot)/2  # mid is calculated from 2d-coordinates
  # coords_2d$mid <- unlist(trans3d(((top + bot)/2)[1], ((top + bot)/2)[2], ((top + bot)/2)[3], pmat)) if use mid in 3d-coordinates
  tb_diff <- coords_2d$top - coords_2d$bot
  angle <- 180/pi * atan2(tb_diff[2], tb_diff[1])
  names(angle) <- "angle"
  center <-  coords_2d$mid + sqrt(space^2 / sum(tb_diff^2)) * rev(abs(tb_diff)) * pos
  out <- list(angle = angle, center = as.data.frame(t(center)))
  return(out)
}

x_top <- c(max(x), max(y), min(z))
x_bot <- c(min(x), max(y), min(z))
y_top <- c(max(x), max(y), min(z))
y_bot <- c(max(x), min(y), min(z))
z_top <- c(max(x), min(y), max(z))
z_bot <- c(max(x), min(y), min(z))

pmat = persp(chl, acdom440, log_vals, theta = 140, phi = 40)

xlab_param <- persp_lab_f(x_top, x_bot, pmat, 0.1, pos = c(1, -1))
ylab_param <- persp_lab_f(y_top, y_bot, pmat, 0.1)
zlab_param <- persp_lab_f(z_top, z_bot, pmat, 0.1)

par(mar=c(1,1,1,1))


nbcol = 20
color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
zcol  = cut(z, nbcol)


persp.mat <- persp(x, y, z,
                   theta = 320, phi = 30,
                   cex.lab=1.3,
                   nticks = 5,
                   #smooth=T, 
                   col = color[zcol], shade = 0.4,
                   border="grey80", axes = T,
                   box=TRUE)

text(cex=1.3,xlab_param$center, srt = xlab_param$angle + 180, "[chl]")
text(cex=1.3,ylab_param$center, srt = ylab_param$angle, "acdom.440")
text(cex = 1.3,zlab_param$center, srt = zlab_param$angle + 180, labels = paste("logL(","Rrs","|","[chl]",",","acdom.440",")"))

# The coords at which we want ticks
x.ticks <- seq(-10,10,5)
# Transform them in 3D
x.3d <- trans3d(x.ticks, 0, 0, persp.mat)
x.3d.1 <- trans3d(x.ticks, 0, -2, persp.mat)
# The coordinates for the text
x.3d.labels <- trans3d(x.ticks, -60, -3, persp.mat)
# Draw the axis ticks
segments(x.3d$x, x.3d$y, x.3d.1$x, x.3d.1$y)
# Write the labels
text(x.3d.labels$x, x.3d.labels$y, x.ticks, cex=0.8)



pars <- as.numeric(Fit.optimized.ssobj)
startvalue = pars
startvalue = c(5,1,0.05)

chain = run_metropolis_MCMC(startvalue, 2000)

burnIn = 1000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))

