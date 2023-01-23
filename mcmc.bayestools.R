require(mcmc)
require(BayesianTools)
require(plot3D)

#Create Prior density function
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
  lklhood.prior = dunif(x = x, min = 0.0001, max = 0.01, log = T)
  
  return(chl.prior+acdom440.prior+anap440.prior+lklhood.prior)
}

#Create Prior sampling function
sampler = function(n=1){
  
  chl.prior = rweibull(n, shape = fit.chl.norm$estimate[1], 
                       scale =fit.chl.norm$estimate[2])
  
  acdom440.prior = rweibull(n,shape = fit.acdom440.norm$estimate[1], 
                            scale =fit.acdom440.norm$estimate[2])
  
  anap440.prior = rweibull(n,shape = fit.anap440.norm$estimate[1], 
                           scale =fit.anap440.norm$estimate[2])
  lklhood.prior = runif(n, min = 0.0001, max = 0.01)
  
  return(cbind(chl.prior,acdom440.prior,anap440.prior,lklhood.prior))
}

#Set actual observation
#obs.data <- rrs.HL[j,]
obs.data <- obsdata

#Create prior density and sampling class
prior.actual <- BayesianTools::createPrior(density = prior, sampler = sampler,
                                           #lower = c(0,0,0,0.0001), upper = c(30,5,0.5, 0.01),
                                           best = NULL)
                                           #best = as.numeric(Fit.optimized.ssobj))

#Create likelihood function
ll <-function(param){
  Gpred = Lee_forward(chl = param[1], acdom440 = param[2],
                        anap440 =param[3],
                        #bbp.550 = HL.deep.iop$bbp550[j],
                        bbp.550 = Fit.input$bbp.550,
                        verbose=F ,realdata = obs.data)
  #smull = sum(dnorm(obs.data, mean = Gpred[[1]]$Rrs, sd = 0.0006327431, log = TRUE), na.rm = T)
  smull = sum(dnorm(obs.data, mean = Gpred[[1]]$Rrs, sd = param[4], log = TRUE))
  #smull = sum(dnorm(obs.data, mean = Gpred[[1]]$Rrs, sd = 0.13129354/100, log = TRUE))
  
  return(smull)
}

#Create Bayesian setup for MCMC
#With prior
bayessetup <- createBayesianSetup(prior = prior.actual,
                                  likelihood = ll,#lower = c(0,0,0), upper = c(30,5,0.5),
                                      names = c("chl","acdom440","anap440", "pop.sd"), parallel = T)
#Only likelihood
bayessetup <- createBayesianSetup(prior = NULL,
                                  likelihood = ll,
                                  lower = lower.bound, upper = upper.bound,
                                  names = c("chl","acdom440","anap440", "pop.sd"), parallel = F)

checkBayesianSetup(bayessetup) #Test if the setup is initiated for theta pars

settings = list(iterations = 10000, message = TRUE, nrChains = 1, burnin=2000) #Set MCMC config
samplerlist <-c("Metropolis", "AM", "DR", "DRAM", "DE", "DEzs", "DREAM", "DREAMzs", "SMC")

#Run MCMC
out <- runMCMC(bayesianSetup = bayessetup, settings = settings, sampler = samplerlist[6] )
summary(out)

#Create Diag plots
plot(out, start = 1000)
#correlationPlot(out, start = 1000)

marginalPlot(out, start = 1000)
gelmanDiagnostics(out, plot = T)
MAP(out)
DIC(out)


# getSample(out, coda = F)

# summary(out)
chain.one = as.data.frame(out[["codaChain"]][[1]])[-1,]
chain.two = as.data.frame(out[["codaChain"]][[2]])[-1,]
chain.three = as.data.frame(out[["codaChain"]][[3]])[-1,]

chain.mean = cbind(as.data.frame(out[["codaChain"]][[1]]),
              as.data.frame(out[["codaChain"]][[2]]),
              as.data.frame(out[["codaChain"]][[3]]))

chain = as.data.frame(t(apply(chain.mean,1, function(x) tapply(x,colnames(chain.mean),mean))))
chain = chain[-1,]

#------------------------------------------------------------------------------------------------
#Create MCMC plots (histograms)
#------------------------------------------------------------------------------------------------
png(filename = paste0("Posterior.SABER.chl=",Fit.input$chl,"_acdom=",Fit.input$acdom.440,"_anap=",Fit.input$anap.440,"_bbp=",Fit.input$bbp.550,".png"), units = "in", width = 6.5, height = 4.5, res=300)
#dev.new()
par(mfrow = c(2,3))
#Posterior of [chl]
hist(chain$chl,nclass=30, main=expression(paste("Posterior of ", italic("[chl]")," [","mg"^1,"m"^-3,"]")) , 
     xlab="MAP=green; MLE=blue; observed=black", probability = T, col = "grey" 
     #xlim=c(min(chl.sample), max(chl.sample)) 
)
lines(density(chain$chl), col="purple", lwd=2)
abline(v = MAP(out)[[1]][1], col="seagreen", lwd=2)
abline(v = Fit.optimized.ssobj$chl, col="blue", lwd=2 )
abline(v = Fit.input$chl, col="black", lwd=2 )

#Posterior of acdom.440
hist(chain$acdom440,nclass=30, main=expression(paste("Posterior of ", italic(a)["CDOM"](440)," [","m"^-1,"]")) , 
     xlab="MAP=green; MLE=blue; observed=black", probability = T, col = "grey" 
     #xlim=c(min(chl.sample), max(chl.sample)) 
)
lines(density(chain$acdom440), col="purple", lwd=2)
abline(v = MAP(out)[[1]][2], col="seagreen", lwd=2)
abline(v = Fit.optimized.ssobj$acdom.440, col="blue", lwd=2 )
abline(v = Fit.input$acdom.440, col="black", lwd=2 )

#Posterior of anap.440
hist(chain$anap440,nclass=30, main=expression(paste("Posterior of ", italic(a)["NAP"](440)," [","m"^-1,"]")) , 
     xlab="MAP=green; MLE=blue; observed=black", probability = T, col = "grey" 
     #xlim=c(min(chl.sample), max(chl.sample)) 
)
lines(density(chain$anap440), col="purple", lwd=2)
abline(v = MAP(out)[[1]][3], col="seagreen", lwd=2)
abline(v = Fit.optimized.ssobj$anap.440, col="blue", lwd=2 )
abline(v = Fit.input$anap.440, col="black", lwd=2 )

#------------------------------------------------------------------------------------------------
#Create MCMC plots (chains)
#------------------------------------------------------------------------------------------------
#Chain of [chl]
plot(chain.one$chl, type = "l", xlab="iterations", col="orange", ylab= expression(paste(italic("[chl]"))),
     main = expression(paste("Chain of ", italic("[chl]")," [","mg"^1,"m"^-3,"]")))
lines(chain.two$chl, col="purple")
lines(chain.three$chl, col="grey")
abline(h = MAP(out)[[1]][1], col="seagreen", lwd=2)
abline(h = Fit.optimized.ssobj$chl, col="blue", lwd=2 )
abline(h = Fit.input$chl, col="black", lwd=2 )

#Chain of acdom.440
plot(chain.one$acdom440, type = "l", xlab="iterations", col="orange", ylab= expression(paste(italic(a)["CDOM"](440))),
     main = expression(paste("Chain of ", italic(a)["CDOM"](440)," [","m"^-1,"]")))
lines(chain.two$acdom440, col="purple")
lines(chain.three$acdom440, col="grey")
abline(h = MAP(out)[[1]][2], col="seagreen", lwd=2)
abline(h = Fit.optimized.ssobj$acdom.440, col="blue", lwd=2 )
abline(h = Fit.input$acdom.440, col="black", lwd=2 )

#Chain of anap.440
plot(chain.one$anap440, type = "l", xlab="iterations", col="orange", ylab= expression(paste(italic(a)["NAP"](440))),
     main = expression(paste("Chain of ", italic(a)["NAP"](440)," [","m"^-1,"]")))
lines(chain.two$anap440, col="purple")
lines(chain.three$anap440, col="grey")
abline(h = MAP(out)[[1]][3], col="seagreen", lwd=2)
abline(h = Fit.optimized.ssobj$anap.440, col="blue", lwd=2 )
abline(h = Fit.input$anap.440, col="black", lwd=2 )
dev.off()

#------------------------------------------------------------------------------------------------
#Create MCMC plots (Random Walk in 3-D parameter space)
#------------------------------------------------------------------------------------------------

#@@@@@@@@@@@@@@@@@@-Create 3D parameter space move-#@@@@@@@@@@@@@@@@@@
dir.create("gfx/2D")  # sub-directory to store individual plot files

iterationlist <- c(seq(1,100,5),seq(101,dim(chain.mean)[1],by = 50))
for(i in iterationlist) {
  #for(i in 1:dim(chain.mean)[1]) {
  png(file = paste0("gfx/2D/randomWalk",i,".png"),  units = "in", width = 4.5, height = 4.5, res=300)
  
  plot(0.5, 0.5, xlim=c(min(chain.mean$chl), max(chain.mean$chl)), 
       ylim=c(min(chain.mean$acdom440), max(chain.mean$acdom440)), type='n', asp=1, las=1,
       ylab=expression(paste(italic(a)["CDOM"](440)," [","m"^-1,"]")), 
       xlab=expression(paste("[",italic("chl"),"]"," [","mg"^1,"m"^-3,"]")),
       main=paste("Iteration #", i-1))
  legend(x = "topright",          # Position
         legend = c("chain1", "chain2", "chain3"),  # Legend texts
         lty = c(1, 1, 1),           # Line types
         col = c("orange", "purple", "blue"),           # Line colors
         lwd = 2)
  
  points(chain.one$chl[1], chain.one$acdom440[1], pch=23, col='orange')
  lines(chain.one$chl[1:i], chain.one$acdom440[1:i], col='grey')
  points(chain.one$chl[2:i], chain.one$acdom440[2:i], pch=16, cex=0.5, col="orange")
  
  points(chain.two$chl[1], chain.two$acdom440[1], pch=23, col='purple')
  lines(chain.two$chl[1:i], chain.two$acdom440[1:i], col='grey')
  points(chain.two$chl[2:i], chain.two$acdom440[2:i], pch=16, cex=0.5, col="purple")
  
  points(chain.three$chl[1], chain.three$acdom440[1], pch=23, col='blue')
  lines(chain.three$chl[1:i], chain.three$acdom440[1:i], col='grey')
  points(chain.three$chl[2:i], chain.three$acdom440[2:i], pch=16, cex=0.5, col="blue")
  
  if (i == max(iterationlist)) {
    points(MAP(out)[[1]][1],MAP(out)[[1]][2], pch=17, col="green", cex=1.2)
  }

  dev.off() # Do this to write plots for animation
  print(paste0(i,"th iteration saved"))
}

#@@@@@@@@@@@@@@@@@@-Create 3D parameter space move-#@@@@@@@@@@@@@@@@@@
dir.create("gfx/3D")  # sub-directory to store individual plot files

iterationlist <- c(seq(1,100,5),seq(101,dim(chain.mean)[1],by = 50))
for(i in iterationlist) {
  #for(i in 1:dim(chain.mean)[1]) {
  png(file = paste0("gfx/3D/randomWalk",i,".png"),  units = "in", width = 8, height = 6, res=300)
  
  par(mar=c(1.2,1.2,1.2,1.2))
  #Start value
  plot3D::scatter3D(x = chain$chl[1],
                    y =  chain$acdom440[1], 
                    #z = log(chain$LP[1]),
                    z =  chain$anap440[1],
                    phi = 0, bty = "g", #theta = -45,
                    theta = 145,
                    xlim=c(min(chain$chl), max(chain$chl)),
                    #xlim=c(4,6), ylim=c(0.8,1), zlim=c(0.02,0.10),
                    ylim=c(min(chain$acdom440), max(chain$acdom440)),
                    #zlim=c(min(log(chain$LP)), max(log(chain$LP))),
                    zlim=c(min(chain$anap440), max(chain$anap440)),
                    #xlab="[chl]", 
                    #ylab="acdom(440)",
                    #zlab="anap(440)",
                    ticktype="detailed", 
                    #clab = c("log(LL)"),
                    main=paste("Iteration #", i-1))
  
  c = chain$LP[1:(i+1)] #Load posterior likelihood
  
  if (i == max(iterationlist)) {
    
    #@@@@ Add the surfaces of actual and predicted@@@@@@@@@ 
    rect3D(x0 = min(chain$chl), y0 = MAP(out)[[1]][2], z0 = min(chain$anap440), x1 = MAP(out)[[1]][1],z1 = MAP(out)[[1]][3], 
           #ylim = c(0, 1), 
           bty = "g", facets = TRUE, 
           border = "blue", col ="red", alpha=0.5,
           lwd = 2, phi = 0, theta = 145,add = T)
    
    rect3D(x0 = min(chain$chl), y0 = Fit.input$acdom.440, z0 = min(chain$anap440), x1 =  Fit.input$chl,z1 = Fit.input$anap.440, 
           #ylim = c(0, 1), 
           bty = "g", facets = TRUE, 
           border = "blue", col ="darkgreen", alpha=0.5,
           lwd = 2, phi = 0, theta = 145, add = T)
    
    #@@@@ Add the legend @@@@@@@@@ 
    predict <- signif(as.numeric(MAP(out)[[1]]), digits = 2); actual <- signif(as.numeric(Fit.input[1:3]),digits = 2)
    predicted <- paste("Predicted:{", predict[1],",",predict[2],",",predict[3],"}")
    actual <- paste0("Actual:{", actual[1],",",actual[2],",",actual[3],"}")
    
    legend(x = "topleft",          # Position
           legend = c(predicted, actual),  # Legend texts
           lty = c(NA, NA),           # Line types
           pch=19,
           col = c("red", "darkgreen"),           # Line colors
           box.col="white",
           box.lwd=0,
           lwd = 2,
           #inset = c(0,1),
           cex=0.9)
  }

  #Plot chains
  plot3D::scatter3D(x = chain$chl[1:(i+1)],y =  chain$acdom440[1:(i+1)],
                    #z = log(chain$LP[1:(i+1)]),
                    z=chain$anap440[1:(i+1)],
                    add = T, phi = 0, bty = "g", type = "b", theta = 145,
                    ticktype = "detailed", pch = 20,
                    #col = cols,
                    clab = c("Posterior LL"),
                    colkey = T, 
                    colvar =c,
                    cex = c(0.5, 1, 1.5))
  
  
  #Add axis labels 
  text(cex=1.3,x = -0.40, y = -0.45, srt = 0 , "acdom(440)")
  text(cex=1.3,x = 0.32, y = -0.45, srt = 0,  "[chl]")
  text(cex = 1.3,10*zlab_param$center, srt = zlab_param$angle + 180, labels = "anap(440)")
  
  
  dev.off() # Do this to write plots for animation
  print(paste0(i,"th iteration saved"))
}

#Create plot with Rotation
plot3Drgl::plotrgl()

#Save MAP outputs of MCMC
Fit.optimized.mcmc <- data.frame("chl"=MAP(out)[[1]][1], 
                                 "acdom.440"=MAP(out)[[1]][2],
                                 "anap.440"=MAP(out)[[1]][3])







