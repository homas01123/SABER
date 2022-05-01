#===========================================================================
# main.R performs the inversion to retrieve the water components along with bathymetry and 
# bottom reflectance (for shallow water) given the wavelengths,

#Author: Mr. Soham Mukherjee - PhD Student, Aquatel, UQAR
#===========================================================================
rm(list=ls(all=TRUE))
source("./saber_forward.R")
source("./saber_inverse.R")
source("./snell_law.R")
source("./solve.objective.inverse.R")
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
setwd("S:/Work/UQAR/Analysis/WISE/R_inverse_wasi")
#--------------------------------------------------------------------------
## Inputs
#Water type specs
type_case_water = 2
type_Rrs_below = "deep"
type_Rrs_water = "below_surface"

# Desired Wavelength for the simulation
wavelength <- seq(400,800,10)

# Observed Rrs for following parameters to be input in RT
Rrs.cops <- read.csv("S:/Work/UQAR/Datasets/WISE/L2_SRIMGSAT/20190818_StationOUT-F18/IOP/hl_simulation/Cops_table_T-F18.csv", 
                    header = T)
Rrs.albert <- read.csv("S:/Work/UQAR/Datasets/WISE/L2_SRIMGSAT/20190818_StationOUT-F18/IOP/hl_simulation/Albert_table_T-F18.csv",
                       header = T)

Rrs_obs <- as.numeric(Rrs.albert[1,-1])
Rrs_obs_wl <- c(320, 330, 340, 380, 412, 443, 465, 490, 510, 
                532, 555, 589, 625, 665, 683, 694, 710, 780, 875)

Rrs_obs.interp = Hmisc::approxExtrap(Rrs_obs_wl, Rrs_obs, xout = wavelength,
                                     method = "linear")$y
# Viewing geometry in degrees
view = 0.98
sun  = 67

#bottom depth
zB=2

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
#--------------------------------------------------------------------------------------------------

##Single RUN
Fit.input <- data.frame("chl"= 3.66 , 
                        "acdom.440"= 1.31 ,
                        "anap.440"= 0.018 )

forward.op <- Saber_forward(chl = Fit.input$chl, acdom440 = Fit.input$acdom.440, 
                            anap440 =Fit.input$anap.440 )
rrs.forward <- forward.op[[1]]$Rrs

##Batch RUN
Fit <- data.frame("C_ph"=seq(1,10,0.5),
                  "a_cdom.440"=seq(0.5,5,0.25),
                  "a.nap.440"=seq(0.01,0.1,0.005))

Fit.input.LUT <- expand.grid(Fit) #Create Fit params LUT

rrs.forward.LUT <- matrix(nrow = length(Fit.input.LUT$C_ph), ncol = length(wavelength),0)

for (i in 1:length(Fit.input.LUT$C_ph)) { #Create Rrs LUT
  temp1 <- as.numeric(Fit.input.LUT[i,])
  temp2 <- Saber_forward(chl = temp1[1], acdom440 = temp1[2], 
                         anap440 =temp1[3] )
  rrs.forward.LUT[i,] <- temp2[[1]]$Rrs
  cat(paste0("\033[0;42m",i," iterations over, ", (nrow(rrs.forward.LUT) - i), " remaining","\033[0m","\n"))
}

#----------------------------------------------------------------------------------------------------
#Perform inverse modeling using optimization of inverse cost function
#-----------------------------------------------------------------------------------------------------
#Pre-FIT of initial values
pre.Fit <- data.frame("C_ph"=seq(1,10,0.5),
                  "a_cdom.440"=seq(0.5,5,0.25),
                  "a.nap.440"=seq(0.01,0.1,0.005))

pre.Fit.input.LUT <- expand.grid(pre.Fit) #Create pre-Fit params LUT

preFIT.rrs.forward.LUT <- matrix(nrow = length(pre.Fit.input.LUT$C_ph),
                                 ncol = length(wavelength),0)
reslist = vector()
for (i in 1:length(pre.Fit.input.LUT$C_ph)) { #Create Rrs LUT
  temp1 <- as.numeric(pre.Fit.input.LUT[i,])
  temp2 <- Saber_forward(chl = temp1[1], acdom440 = temp1[2], 
                         anap440 =temp1[3] )
  preFIT.rrs.forward.LUT[i,] <- temp2[[1]]$Rrs
  reslist[i] = temp2[[2]]
  cat(paste0("\033[0;42m",i," iterations over, ", (nrow(preFIT.rrs.forward.LUT) - i), " remaining","\033[0m","\n"))
}

prefit.best <- pre.Fit.input.LUT[which.min(reslist),] #retrieve best initial values using SSR

prefit.Rrs <- Saber_forward(chl = prefit.best$C_ph, acdom440 = prefit.best$a_cdom.440 , 
              anap440 =prefit.best$a.nap.440)[[1]]$Rrs

#plot(wavelength, Rrs_obs.interp, type="l", col="red")
#lines(wavelength, prefit.Rrs, col="green")

##Do the optimization 
#Initial values
par0 = c(chl = prefit.best$C_ph, acdom440 = prefit.best$a_cdom.440, anap440 = prefit.best$a.nap.440, 
         population.sd = 0.1) 

##Single RUN
obj = "log-LL"
inverse_output <- solve.objective.inverse(initial = par0, obsdata = Rrs_obs.interp,
                                          obj.fn = "log-LL")
if (obj == "log-LL") {
  Fit.optimized.ssobj <- data.frame("chl"=inverse_output[[1]]$estimates[1], 
                                    "acdom.440"=inverse_output[[1]]$estimates[2],
                                    "anap.440"=inverse_output[[1]]$estimates[3])
} else{
  if (obj == "SSR") {
    
    Fit.optimized.ssobj <- data.frame("chl"=inverse_output[[1]]$chl,
                                      "acdom.440"=inverse_output[[1]]$acdom.440,
                                "anap.440"=inverse_output[[1]]$anap.440)
  }
}

#Generate Rrs with the inversion retrieved params
forward.fit.optimized <- Saber_forward(chl = Fit.optimized.ssobj$chl, acdom440 = Fit.optimized.ssobj$acdom.440, 
                            anap440 =Fit.optimized.ssobj$anap.440 )
rrs.forward.fit.optimized <- forward.fit.optimized[[1]]$Rrs

#plot(wavelength, Rrs_obs.interp, type="l", col="red")
#lines(wavelength, rrs.forward.fit.optimized, col="blue")

#Create Plot of actual Rrs and Rrs simulated with invrsion retrieved params
plot= FALSE
plotframe.rrs <- data.frame("wave"=wavelength, "rrs.est"=rrs.forward.fit.optimized, 
                            "rrs.obs"=Rrs_obs.interp)

xmin = min(plotframe.rrs$wave); xmax= max(plotframe.rrs$wave); xstp=100
ymin= 0; ymax=max(plotframe.rrs$rrs.obs)+0.10*max(plotframe.rrs$rrs.obs);ystp= signif(ymax/5, digits = 1)
asp_rat <- (xmax-xmin)/(ymax-ymin)

g1 <- ggplot()  + geom_line(data = plotframe.rrs,aes(y=rrs.est,color="xx1",x = wave),
                            size=1.3,show.legend = TRUE) +
  geom_line(data = plotframe.rrs,aes(y=rrs.obs,x = wave,color="xx5"),linetype="dashed", 
            size=1.3,show.legend = TRUE)+
  scale_colour_manual(labels = c(expression(paste(italic("R")["rs,model"])),
                                 expression(paste(italic("R")["rs,actual"]))), 
                      values = c("blue","green")) +
  #ggtitle(paste0(stationlist[i])) +
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
        legend.position=c(0.70, 0.9),
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
  ggsave(paste0("./SABER.inv.Rrs.png"), plot = g1,
         scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
}

##Batch RUN
Fit.optimized.ssobj.LUT <- matrix(nrow = length(Fit.input.LUT$C_ph), ncol = length(params),0)
colnames(Fit.optimized.ssobj.LUT) <- c("chl", "acdom440", "anap440")


for (i in 1:nrow(Fit.optimized.ssobj.LUT)) {
  
  inverse_output <- solve.objective.inverse(obj.fn = "SSR", initial = par0, 
                                            obsdata =rrs.forward.LUT[i,] )
  cat(paste0("\033[0;42m",i," iterations over, ", (nrow(Fit.optimized.ssobj.LUT) - i), " remaining","\033[0m","\n"))
  Fit.optimized.ssobj.LUT[i,] <- as.numeric(inverse_output[[1]])
}
