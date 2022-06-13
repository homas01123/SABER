#===========================================================================
# main.R performs the inversion to retrieve the water components along with bathymetry and 
# bottom reflectance (for shallow water) given the wavelengths,

#Author: Mr. Soham Mukherjee - PhD Student, Aquatel, UQAR
#===========================================================================
rm(list=ls(all=TRUE))
setwd("S:/Work/UQAR/Analysis/WISE/R_inverse_wasi")
source("./saber_forward.R")
source("./saber_inverse.R")
source("./snell_law.R")
source("./solve.objective.inverse.R")
source("./lee_forward.R")
#----------------------------------------------------------------------------------------------------
require(dplyr)
require(readxl)
require(stats4)
require(MASS)
require(dglm)
library(fitdistrplus)
require(Riops)
require(Cops)
require(ggplot2)
require(rho)
require(marqLevAlg)
require(BayesianTools)
require(coda)
#--------------------------------------------------------------------------
## Inputs
#Water type specs
type_case_water = 2
type_Rrs_below = "deep"
type_Rrs_water = "below_surface"

# Desired Wavelength for the simulation
wavelength <- seq(400,800,10)

# Observed/in situ Rrs for following parameters to be input in RT
Rrs.cops <- read.csv("S:/Work/UQAR/Datasets/WISE/L2_SRIMGSAT/20190818_StationOUT-F18/IOP/hl_simulation/Cops_table_T-F18.csv", 
                    header = T)
Rrs.albert <- read.csv("S:/Work/UQAR/Datasets/WISE/L2_SRIMGSAT/20190818_StationOUT-F18/IOP/hl_simulation/Albert_table_T-F18.csv",
                       header = T)
Rrs.HL <- read.csv("S:/Work/UQAR/Datasets/WISE/L2_SRIMGSAT/20190818_StationOUT-F18/IOP/hl_simulation/HL_el_table_T-F18.csv",
                       header = T)

Rrs_obs <- as.numeric(Rrs.albert[1,-1])
Rrs_obs_wl <- c(320, 330, 340, 380, 412, 443, 465, 490, 510, 
                532, 555, 589, 625, 665, 683, 694, 710, 780, 875)

Rrs_obs.interp = Hmisc::approxExtrap(Rrs_obs_wl, Rrs_obs, xout = wavelength,
                                     method = "linear")$y

Rrs_obs.cops <- as.numeric(Rrs.cops[1,-1])
Rrs_obs.cops.interp = Hmisc::approxExtrap(Rrs_obs_wl, Rrs_obs.cops, xout = wavelength,
                                     method = "linear")$y

# Viewing geometry in degrees
view = 0
sun  = 60

#bottom depth
zB=0.5

# areal fraction of bottom surface
fA0=0; # constant 
fA1=0; # sand
fA2=0.5; # sediment
fA3=0; # Chara contraria
fA4=0.5; # Potamogeton perfoliatus
fA5=0; # Potamogeton pectinatus
fA.set= c(fA0,fA1,fA2,fA3,fA4,fA5)

## Atmospheric conditions

# Irradiance intensities [1/sr]
g_dd=0.05; g_dsr=0; g_dsa=0;

# Intensities of light sources 
f_dd= 1; f_ds= 1;

# Angstrom exponent
alpha = 1.317;

# Atmospheric pressure 
P = 1013.25; # [mbar]

# Relative Humidity
RH = 0.60;

# Scale height for ozone
Hoz=0.300; # [cm]

# Scale height of the precipitable water in the atmosphere
WV= 2.500; # [cm]

#-----------------------------------------------------------------------------------------------------
#Perform forward modeling
#-----------------------------------------------------------------------------------------------------
##Single FORWARD RUN

batch=TRUE #Set TRUE for IOCCG dataset duplication; Set FALSE for user wanted random inputs 
insitu.present=TRUE #Set TRUE if actual in situ observation or simulated observation exist; else set FALSE
plot=FALSE #Set TRUE when the ggplot2 outputs needed to be saved onto disk
j=376 # No. of IOCCG data point among 500

if (insitu.present == TRUE) {#Set Rrs observation data
  
  insitu.data <-rrs.HL[j,] #Set the observed Rrs manually here
}

if (batch == TRUE) { #For IOCCG
  Fit.input <- data.frame("chl"=  HL.deep.iop$chl[j],
                          "acdom.440"=  HL.deep.iop$acdom440[j],
                          "anap.440"= HL.deep.iop$anap440[j],
                          "bbp.550" = HL.deep.iop$bbp550[j])
} else { #Manual entries
  Fit.input <- data.frame("chl"=  4.77, 
                          "acdom.440"=  0.97,
                          "anap.440"= 0.017,
                          "bbp.550"=0.005)
}

forward.op.am <- Saber_forward(chl = Fit.input$chl, acdom440 = Fit.input$acdom.440, 
                            anap440 =Fit.input$anap.440 , bbp.550 = Fit.input$bbp.550, 
                            realdata = Rrs_obs.interp, verbose = T,z = zB,
                            rb.fraction = fA.set)

rrs.forward.am <- forward.op.am[[1]]$Rrs #Extract AM03 modelled Rrs

forward.op.lee <- Lee_forward(chl = Fit.input$chl, acdom440 = Fit.input$acdom.440, 
                              anap440 =Fit.input$anap.440 , 
                              bbp.550 = Fit.input$bbp.550, 
                              realdata = Rrs_obs.interp, verbose = T,z = zB,
                              rb.fraction = fA.set)

rrs.forward.lee <- forward.op.lee[[1]]$Rrs #Extract Lee98 modelled Rrs

#Plot the forward simulated Rrs (AM03 & Lee98 with optional in situ)
if (insitu.present == TRUE) {
  forward.rrs <- data.frame("wave"=wavelength, "rrs.obs"=insitu.data,
                              "rrs.am03"=rrs.forward.am,
                              "rrs.lee98"=rrs.forward.lee)
  
  xmin = min(forward.rrs$wave); xmax= max(forward.rrs$wave); xstp=100
  ymin= 0; ymax=max(forward.rrs$rrs.am03)+0.20*max(forward.rrs$rrs.am03);ystp= signif(ymax/5, digits = 1)
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  g1 <- ggplot()  + geom_line(data = forward.rrs,aes(y=rrs.obs,color="xx1",x = wave),
                              linetype="dashed",size=1.3,show.legend = TRUE) +
    geom_line(data = forward.rrs,aes(y=rrs.am03,color="xx2",x = wave),
              size=1.3,show.legend = TRUE) +
    geom_line(data = forward.rrs,aes(y=rrs.lee98,x = wave,color="xx3"), 
              size=1.3,show.legend = TRUE)+
    
    scale_colour_manual(labels = c(expression(paste(italic("R")["rs,model,actual"])),
                                   expression(paste(italic("R")["rs,model,AM03"])),
                                   expression(paste(italic("R")["rs,model,Lee98"]))), 
                        values = c("red","green","purple")) +
    
    scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), limits = c(xmin, xmax), 
                       breaks = seq(xmin, xmax, xstp))  +
    scale_y_continuous(name =expression(paste(italic("R"),{}[rs],"(",lambda,",", 0^"-",")[", sr^-1,"]")) , limits = c(ymin, ymax),
                       breaks = seq(ymin, ymax, ystp))+ 
    coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
                ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
    
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
          axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          axis.ticks.length = unit(.25, "cm"),
          legend.position=c(0.55, 0.9),
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.text = element_text(colour = "black", size = 20, face = "plain"),
          legend.background = element_rect(fill = NA, size = 0.5, 
                                           linetype = "solid", colour = 0),
          legend.key = element_blank(),
          legend.justification = c("left", "top"),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey", 
                                          size = 0.5, linetype = "dotted"), 
          panel.grid.minor = element_blank(),
          #legend.spacing.y = unit(2.0, 'cm'),
          plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
          legend.text.align = 0,
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  g1 
  
  if (plot == "TRUE") {
    ggsave(paste0("./Rrs_forward_chl=",Fit.input$chl,"_acdom=",Fit.input$acdom.440,"_anap=",Fit.input$anap.440,"_bbp=",Fit.input$bbp.550,".png"), plot = g1,
           scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
  }
} else {
  forward.rrs <- data.frame("wave"=wavelength, #"rrs.obs"=insitu.data,
                            "rrs.am03"=rrs.forward.am,
                            "rrs.lee98"=rrs.forward.lee)
  
  xmin = min(forward.rrs$wave); xmax= max(forward.rrs$wave); xstp=100
  ymin= 0; ymax=max(forward.rrs$rrs.am03)+0.20*max(forward.rrs$rrs.am03);ystp= signif(ymax/5, digits = 1)
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  g1 <- ggplot()  + 
    geom_line(data = forward.rrs,aes(y=rrs.am03,color="xx2",x = wave),
              size=1.3,show.legend = TRUE) +
    geom_line(data = forward.rrs,aes(y=rrs.lee98,x = wave,color="xx3"), 
              size=1.3,show.legend = TRUE)+
    
    scale_colour_manual(labels = c(expression(paste(italic("R")["rs,model,AM03"])),
                                   expression(paste(italic("R")["rs,model,Lee98"]))), 
                        values = c("green","purple")) +
    
    scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), limits = c(xmin, xmax), 
                       breaks = seq(xmin, xmax, xstp))  +
    scale_y_continuous(name =expression(paste(italic("R"),{}[rs],"(",lambda,",", 0^"-",")[", sr^-1,"]")) , limits = c(ymin, ymax),
                       breaks = seq(ymin, ymax, ystp))+ 
    coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
                ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
    
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
          axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          axis.ticks.length = unit(.25, "cm"),
          legend.position=c(0.55, 0.9),
          legend.direction = "vertical",
          legend.title = element_blank(),
          legend.text = element_text(colour = "black", size = 20, face = "plain"),
          legend.background = element_rect(fill = NA, size = 0.5, 
                                           linetype = "solid", colour = 0),
          legend.key = element_blank(),
          legend.justification = c("left", "top"),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey", 
                                          size = 0.5, linetype = "dotted"), 
          panel.grid.minor = element_blank(),
          #legend.spacing.y = unit(2.0, 'cm'),
          plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
          legend.text.align = 0,
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  g1 
  
  if (plot == "TRUE") {
    ggsave(paste0("./Rrs_forward_chl=",Fit.input$chl,"_acdom=",Fit.input$acdom.440,"_anap=",Fit.input$anap.440,"_bbp=",Fit.input$bbp.550,".png"), plot = g1,
           scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
  }
}
#---------------------------------------------------------------------------------------------------
#Batch FORWARD RUN
Fit <- data.frame("C_ph"=seq(1,10,0.5),
                  "a_cdom.440"=seq(0.5,5,0.25),
                  "a.nap.440"=seq(0.01,0.1,0.005))

Fit.input.LUT <- expand.grid(Fit) #Create Fit params LUT

rrs.forward.LUT <- matrix(nrow = length(Fit.input.LUT$C_ph), ncol = length(wavelength),0)

for (i in 1:length(Fit.input.LUT$C_ph)) { #Create Rrs LUT
  temp1 <- as.numeric(Fit.input.LUT[i,])
  temp2 <- Saber_forward(chl = temp1[1], acdom440 = temp1[2], 
                         anap440 =temp1[3],bbp.550 = Fit.input$bbp.550, realdata = rrs.forward.am )
  rrs.forward.LUT[i,] <- temp2[[1]]$Rrs
  cat(paste0("\033[0;42m",i," iterations over, ", (nrow(rrs.forward.LUT) - i), " remaining","\033[0m","\n"))
  
}

#----------------------------------------------------------------------------------------------------
#Perform inverse modeling using optimization of inverse cost function
#-----------------------------------------------------------------------------------------------------
#Set Rrs observation data
obsdata <-insitu.data#manual set obsdata to begin inversion

#Pre-FIT of initial values
pre.Fit <- data.frame("C_ph"=seq(1,10,0.5),
                  "a_cdom.440"=seq(0.5,5,0.25),
                  "a.nap.440"=seq(0.01,0.1,0.005))

pre.Fit.input.LUT <- expand.grid(pre.Fit) #Create pre-Fit params LUT

preFIT.rrs.forward.LUT <- matrix(nrow = length(pre.Fit.input.LUT$C_ph),
                                 ncol = length(wavelength),0)
#Create the Progress Bar
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = length(pre.Fit.input.LUT$C_ph), # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")

reslist = vector()
for (i in 1:length(pre.Fit.input.LUT$C_ph)) { #Create Rrs LUT
  temp1 <- as.numeric(pre.Fit.input.LUT[i,])
  temp2 <- Saber_forward(chl = temp1[1], acdom440 = temp1[2], 
                         anap440 =temp1[3], bbp.550 = Fit.input$bbp.550,
                         realdata = obsdata,verbose=F )
  
  preFIT.rrs.forward.LUT[i,] <- temp2[[1]]$Rrs
  reslist[i] = temp2[[2]]
  #cat(paste0("\033[0;43m",i," iterations over, ", (nrow(preFIT.rrs.forward.LUT) - i), " remaining","\033[0m","\n"))
  setTxtProgressBar(pb, i)
  if (i == length(pre.Fit.input.LUT$C_ph)) {
    cat(paste0("\033[0;42m","###############PRE-FIT FINISHED################","\033[0m","\n"))
  }
}

prefit.best <- pre.Fit.input.LUT[which.min(reslist),] #retrieve best initial values using C.R.I.S.T.A.L.[mimimizing SSR]

rrs.prefit <- Saber_forward(chl = prefit.best$C_ph, acdom440 = prefit.best$a_cdom.440 , 
              anap440 =prefit.best$a.nap.440, bbp.550 = Fit.input$bbp.550,realdata = obsdata, verbose = T)[[1]]$Rrs

#plot(wavelength, Rrs_obs.interp, type="l", col="red")
#lines(wavelength, prefit.Rrs, col="green")

#-----------------------------------------------------------------------------------------
##Do the optimization 
#Initial values from pre-fit 
pop.sd = "unknown"
if (pop.sd == "known") {
  par0 = c(chl = prefit.best$C_ph, acdom440 = prefit.best$a_cdom.440, #@@population sigma KNOWN
           anap440 = prefit.best$a.nap.440)#, population.sd = 0.0006327431)
} else {
  par0 = c(chl = prefit.best$C_ph, acdom440 = prefit.best$a_cdom.440, #@@population sigma UNKNOWN
           anap440 = prefit.best$a.nap.440, population.sd = 0.001)
}

#Initial values from USER
# par0 = c(chl = 2, acdom440 = 0.8, 
#          anap440 = 0.05, population.sd = 0.1) 

increament.scale <- 1
lower.bound <- par0 - 0.8*par0
upper.bound <- par0 + 5*par0

##Single RUN
obj = c("log-LL", "SSR"); obj.run <- obj[1]
methods.opt <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent","levenberg-marqardt")

inverse_output <- solve.objective.inverse(initial = par0, obsdata = obsdata,
                                          sa.model = "lee99", obj.fn =obj.run , 
                                          method.opt = methods.opt[4],
                                          lower.b = lower.bound,
                                          upper.b = upper.bound, 
                                          batch = TRUE,pop.sd = FALSE)
if (obj.run == "log-LL") {
  Fit.optimized.ssobj <- data.frame("chl"=inverse_output[[1]]$estimates[1], 
                                    "acdom.440"=inverse_output[[1]]$estimates[2],
                                    "anap.440"=inverse_output[[1]]$estimates[3])
} else{
  if (obj.run == "SSR") {
    
    Fit.optimized.ssobj <- data.frame("chl"=inverse_output[[1]]$chl,
                                      "acdom.440"=inverse_output[[1]]$acdom.440,
                                "anap.440"=inverse_output[[1]]$anap.440)
  }
}
#------------------------------------------------------------------------------------
#Implement MCMC optimization

#Create prior density and sampling class
prior.actual <- BayesianTools::createPrior(density = prior, sampler = sampler,
                                           lower = c(0,0,0), upper = c(30,5,0.5), 
                                           best = as.numeric(Fit.optimized.ssobj))
#Create Bayesian setup for MCMC
bayessetup <- createBayesianSetup(prior = prior.actual,
                                  likelihood = ll,#lower = c(0,0,0), upper = c(30,5,0.5),
                                  names = c("chl","acdom440","anap440"), parallel = T)

checkBayesianSetup(bayessetup) #Test if the setup is inititated for theta pars

settings = list(iterations = 10000, message = TRUE, nrChains = 1, burnin=2000) #Set MCMC config
samplerlist <-c("Metropolis", "AM", "DR", "DRAM", "DE", "DEzs", "DREAM", "DREAMzs", "SMC")

#Run MCMC
out <- runMCMC(bayesianSetup = bayessetup, settings = settings, sampler = samplerlist[6] )
summary(out)

#MCMC diagnostics
plot(out, start = 1000) #chain and parameter density
correlationPlot(out, start = 1000) #correlation plot among parameters
marginalPlot(out, start = 1000) #Variation in marginal prob density of prior and posterior

MAP.mcmc <- MAP(out) #Store MAP
DIC.mcmc <- DIC(out) #Store DIC

Fit.optimized.mcmc <- data.frame("chl"=MAP(out)[[1]][1], #Save MCMC MAP outputs
                                 "acdom.440"=MAP(out)[[1]][2],
                                 "anap.440"=MAP(out)[[1]][3])


#------------------------------------------------------------------------------------------
#Generate Rrs with the inversion retrieved params
forward.fit.optimized.sse <- Saber_forward(chl = Fit.optimized.ssobj$chl, 
                                           acdom440 = Fit.optimized.ssobj$acdom.440, 
                            anap440 =Fit.optimized.ssobj$anap.440,
                            #bbp.550 = HL.deep.iop$bbp550[j],
                            bbp.550 = Fit.input$bbp.550,
                            realdata = obsdata)
rrs.forward.fit.optimized.sse <- forward.fit.optimized.sse[[1]]$Rrs


forward.fit.optimized.mcmc <- Lee_forward(chl = Fit.optimized.mcmc$chl, 
                                            acdom440 = Fit.optimized.mcmc$acdom.440, 
                                           anap440 =Fit.optimized.mcmc$anap.440, 
                                           #bbp.550 = HL.deep.iop$bbp550[j],
                                           bbp.550 = Fit.input$bbp.550,
                                           realdata = obsdata)
rrs.forward.fit.optimized.mcmc <- forward.fit.optimized.mcmc[[1]]$Rrs


#Create Plot of actual Rrs and Rrs simulated with invrsion retrieved params
plot= FALSE
plotframe.rrs <- data.frame("wave"=wavelength, "rrs.est.sse"=rrs.forward.fit.optimized.sse,
                            "rrs.est.mcmc"=rrs.forward.fit.optimized.mcmc,
                            "rrs.obs"=insitu.data, "rrs.prefit"= rrs.prefit)
#Create labels
mcmc.label = paste0("theta[MAP]== {",signif(Fit.optimized.mcmc$chl,digits = 2),"*',",signif(Fit.optimized.mcmc$acdom.440,digits = 2),",",signif(Fit.optimized.mcmc$anap.440,digits = 2),"'}")
mle.label = paste0("theta[MLE]== {",signif(Fit.optimized.ssobj$chl,digits = 2),"*',",signif(Fit.optimized.ssobj$acdom.440,digits = 2),",",signif(Fit.optimized.ssobj$anap.440,digits = 2),"'}")
obs.label = paste0("theta[obs]== {",signif(Fit.input$chl,digits = 2),"*',",signif(Fit.input$acdom.440,digits = 2),",",signif(Fit.input$anap.440,digits = 2),"'}")
prefit.label = paste0("theta[prefit]== {",signif(prefit.best$C_ph,digits = 2),"*',",signif(prefit.best$a_cdom.440,digits = 2),",",signif(prefit.best$a.nap.440,digits = 2),"'}")

#Create AXIS
xmin = min(plotframe.rrs$wave); xmax= max(plotframe.rrs$wave); xstp=100
ymin= 0; ymax=max(plotframe.rrs$rrs.obs)+0.20*max(plotframe.rrs$rrs.obs);ystp= signif(ymax/5, digits = 1)
asp_rat <- (xmax-xmin)/(ymax-ymin)

#Create Plot
g1 <- ggplot()  + geom_line(data = plotframe.rrs,aes(y=rrs.est.sse,color="xx1",x = wave),
                            size=1.3,show.legend = TRUE) +
  geom_line(data = plotframe.rrs,aes(y=rrs.est.mcmc,color="xx2",x = wave),
            size=1.3,show.legend = TRUE) +
  geom_line(data = plotframe.rrs,aes(y=rrs.obs,x = wave,color="xx6"), 
            size=1.3,show.legend = TRUE)+
  geom_line(data = plotframe.rrs,aes(y=rrs.obs,x = wave,color="xx5"),linetype="dashed", 
            size=1.3,show.legend = TRUE)+
  
  scale_colour_manual(labels = c(expression(paste(italic("R")["rs,model,MLE"])),
                                 expression(paste(italic("R")["rs,model,MAP"])),
                                 expression(paste(italic("R")["rs,model,prefit"])),
                                 expression(paste(italic("R")["rs,actual"]))), 
                      values = c("blue","green","purple","red")) +
  
  scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, xstp))  +
  scale_y_continuous(name =expression(paste(italic("R"),{}[rs],"(",lambda,",", 0^"-",")[", sr^-1,"]")) , limits = c(ymin, ymax),
                     breaks = seq(ymin, ymax, ystp))+ 
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
  
  annotate("text",x=470, y= ymax*0.95, label = obs.label, parse=T, size=4, color="red")+
  annotate("text",x=470, y= ymax*0.90, label = prefit.label, parse=T, size=4, color="purple")+
  annotate("text",x=470, y= ymax*0.85, label = mle.label, parse=T, size=4, color="blue")+
  annotate("text",x=470, y= ymax*0.80, label = mcmc.label, parse=T, size=4, color="green")+
  
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
        axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.position=c(0.55, 0.9),
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 20, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_blank(),
        #legend.spacing.y = unit(2.0, 'cm'),
        plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
        legend.text.align = 0,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g1 

if (plot == "TRUE") {
  ggsave(paste0("./SABER.simulated.chl=",Fit.input$chl,"_acdom=",Fit.input$acdom.440,"_anap=",Fit.input$anap.440,"_bbp=",Fit.input$bbp.550,".png"), plot = g1,
         scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
}

#-------------------------------------------------------------------------------------------------
##Batch Inversion
Fit.optimized.ssobj.LUT <- matrix(nrow = length(Fit.input.LUT$C_ph), ncol = length(params),0)
colnames(Fit.optimized.ssobj.LUT) <- c("chl", "acdom440", "anap440")


for (i in 1:nrow(Fit.optimized.ssobj.LUT)) {
  
  inverse_output <- solve.objective.inverse(obj.fn = "log-LL", initial = par0, 
                                            obsdata =rrs.forward.LUT[i,] )
  cat(paste0("\033[0;42m",i," iterations over, ", (nrow(Fit.optimized.ssobj.LUT) - i), " remaining","\033[0m","\n"))
  Fit.optimized.ssobj.LUT[i,] <- as.numeric(inverse_output[[1]])
}
