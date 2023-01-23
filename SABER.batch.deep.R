rm(list=ls(all=TRUE))
require(dplyr)
require(readxl)
require(stats4)
require(MASS)
require(dglm)
library(fitdistrplus)
require(Riops)
require(Cops)

setwd("S:/Work/UQAR/Analysis/WISE/R_inverse_wasi")
#Declare the functions

read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#Load IOP and AOP data

#water IOP
insitu.data <- read_excel_allsheets(paste0(getwd(),"/IOP_AOP_Sun60.xls"))
water.data <- insitu.data$Basics[6:46, 1:3]                           
names(water.data) <- c("wavelength", "a_w", "bb_w") 

#chl, acdom440, anap440 & bbp550 
chldata <- insitu.data$a_ph[-1,1]
acdom440data <- insitu.data$a_g[which(names(insitu.data$a_g[1,]) == "440")]
anap440data <- insitu.data$a_dm[which(names(insitu.data$a_dm[1,]) == "440")]
bbp550data <- insitu.data$bb[which(names(insitu.data$bb[1,]) == "550")] -
  as.numeric(water.data$bb_w[water.data$wavelength == 550])

HL.deep.iop <- data.frame("chl"=chldata, "acdom440"= acdom440data$`440`,
                          "anap440"= anap440data$`440`, "bbp550"= bbp550data$`550`)
#sub-surface (0^-) Rrs
rrs.HL <- as.matrix(insitu.data$r_rs); rrs.HL.wl <- as.numeric(names(insitu.data$r_rs))

#Generate the AM03 simulated Rrs(0^-) using loaded IOPs
rrs.forward.SABER <- matrix(nrow = length(HL.deep.iop$chl), ncol = length(wavelength),0)

#Create the Progress Bar
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = length(HL.deep.iop$chl), # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")


for (i in 1:length(HL.deep.iop$chl)) { #Create forward SABER LUT
  temp1 <- as.numeric(HL.deep.iop[i,])
  temp2 <- Saber_forward(chl = temp1[1], acdom440 = temp1[2], 
                         anap440 =temp1[3], bbp.550 = temp1[4], verbose = F, 
                         realdata = rrs.HL[i,])
  rrs.forward.SABER[i,] <- temp2[[1]]$Rrs
  setTxtProgressBar(pb, i)
  if (i == length(HL.deep.iop$chl)) {
    cat(paste0("\033[0;32m","###############LUT generation for AM03 FINISHED################","\033[0m","\n"))
  }
  #cat(paste0("\033[0;42m",i," iterations over, ", (nrow(rrs.forward.SABER) - i), " remaining","\033[0m","\n"))
  
}

#Generate the LEE99 simulated Rrs(0^-) using loaded IOPs
rrs.forward.lee <- matrix(nrow = length(HL.deep.iop$chl), ncol = length(wavelength),0)

for (i in 1:length(HL.deep.iop$chl)) { #Create forward SABER LUT
  temp1 <- as.numeric(HL.deep.iop[i,])
  temp2 <- Lee_forward(chl = temp1[1], acdom440 = temp1[2], 
                         anap440 =temp1[3], bbp.550 = temp1[4], verbose = F, realdata = rrs.HL[i,] )
  rrs.forward.lee[i,] <- temp2[[1]]$Rrs
  #cat(paste0("\033[0;42m",i," iterations over, ", (nrow(rrs.forward.lee) - i), " remaining","\033[0m","\n"))
  setTxtProgressBar(pb, i)
  if (i == length(HL.deep.iop$chl)) {
    cat(paste0("\033[0;32m","###############LUT generation for LEE99 FINISHED################","\033[0m","\n"))
  }
}

#Test the 1:1 scatter b/w SABER forward and HL
plot(as.matrix(rrs.HL), as.numeric(rrs.forward.lee), xlim=c(0,0.05), ylim=c(0,0.05))
abline(0,1, col="green", lwd=3)

#Find the population parameters for noise (/epsilon ~= Rrs_HL - Rrs_SABER)
SABER.deep.residual <- rrs.HL - rrs.forward.SABER
hist(as.matrix(SABER.deep.residual))

SABER.deep.residual.norm <- fitdistrplus::fitdist(c(t(100*SABER.deep.residual)), "norm") #try to fit normal
plot(SABER.deep.residual.norm)
#---------------------------------------------------------------------------------------
# yx.lm <- lm(as.vector(rrs.HL) ~ as.vector(rrs.forward.SABER))
# summary(yx.lm)
# abline(yx.lm)
# 
# rrs.HL <- rrs.forward.SABER*yx.lm$coefficients[2] + yx.lm$coefficients[1]
#---------------------------------------------------------------------------------------
#PRE-FIT for inversion
pre.Fit <- data.frame("C_ph"=seq(min(chldata),max(chldata),0.5),
                      "a_cdom.440"=seq(min(HL.deep.iop$acdom440),max(HL.deep.iop$acdom440),
                                       (max(HL.deep.iop$acdom440) - min(HL.deep.iop$acdom440))/60)[-1],
                      "a.nap.440"=seq(min(HL.deep.iop$anap440),max(HL.deep.iop$anap440),
                                      (max(HL.deep.iop$anap440) - min(HL.deep.iop$anap440))/60)[-1])

#pre.Fit.input.LUT <- expand.grid(pre.Fit) #Create pre-Fit params LUT
pre.Fit.input.LUT <- pre.Fit #Create pre-Fit params LUT

reslist = vector()
#prefit.best <- data.frame(matrix(ncol = 3, nrow =dim(rrs.forward.SABER)[1], 0 ))
Fit.optimized.ssobj.batch <- data.frame("chl"=0, 
                                  "acdom.440"=0,
                                  "anap.440"=0,
                                  "chl.sd"=0, "acdom.440.sd"=0, "anap440.sd"=0)
for (j in 1:dim(rrs.HL)[1]) {
  
  for (i in 1:length(pre.Fit.input.LUT$C_ph)) { #Create Rrs LUT
    temp1 <- as.numeric(pre.Fit.input.LUT[i,])
    temp2 <- Saber_forward(chl = temp1[1], acdom440 = temp1[2], 
                           anap440 =temp1[3], verbose=F, bbp.550 = HL.deep.iop$bbp550[j], 
                           realdata = rrs.HL[j,] )
    #preFIT.rrs.forward.LUT[i,] <- temp2[[1]]$Rrs
    reslist[i] = temp2[[2]]
    cat(paste0("\033[0;33m",j, " no spectra: ",i," iterations over, ", (length(pre.Fit.input.LUT$C_ph) - i), " remaining","\033[0m","\n"))
    if (i == length(pre.Fit.input.LUT$C_ph)) {
      cat(paste0("\033[0;32m","###############PRE-FIT FINISHED for ", j, " th spectra################","\033[0m","\n"))
    }
  }
  
  prefit.best <- pre.Fit.input.LUT[which.min(reslist),] #retrieve best initial values using SSR
  names(prefit.best) <- c("chl", "acdom440", "anap440")
  
  ##Do the optimization 
  #Initial values from pre-fit
  par0 = c(chl = prefit.best$chl, acdom440 = prefit.best$acdom440, 
           anap440 = prefit.best$anap440)#, bb.550 = HL.deep.iop$bbp550[j])
  lower.bound <- par0 - 0.5*par0
  upper.bound <- par0 + par0
  
  ##BATCH RUN
  obj = "log-LL"
  methods.opt <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent")
  inverse_output <- solve.objective.inverse(initial = par0, obsdata = rrs.HL[j,],
                                            method.opt =  methods.opt[4], 
                                            obj.fn = "log-LL", sa.model = "am03",
                                            lower.b = lower.bound, upper.b = upper.bound, 
                                            batch = TRUE)
  
  Fit.optimized.ssobj.batch <- rbind(Fit.optimized.ssobj.batch, c(inverse_output[[1]]$estimates[1:3],
                                                      inverse_output[[1]]$`sd(+/-)`[1:3]))
  cat(paste0("\033[0;43m","############### INVERSION FINISHED for ", j, " th spectra ################","\033[0m","\n"))
}

Fit.optimized.ssobj.batch <- Fit.optimized.ssobj.batch[-1,]

rrs.HL.est.sse <- Saber_forward(chl = Fit.optimized.ssobj.batch$chl, 
                            acdom440 = Fit.optimized.ssobj.batch$acdom.440,
                              anap440 = Fit.optimized.ssobj.batch$anap.440, 
                            bbp.550 = HL.deep.iop$bbp550[j],
                            realdata = rrs.HL[j,])[[1]]$Rrs 

rrs.HL.est.mcmc <- Saber_forward(chl = Fit.optimized.mcmc$chl, 
                                            acdom440 = Fit.optimized.mcmc$acdom.440, 
                                            anap440 =Fit.optimized.mcmc$anap.440, 
                                            bbp.550 = HL.deep.iop$bbp550[j],
                                            #bbp.550 = Fit.input$bbp.550,
                                            realdata = Rrs_obs.interp)[[1]]$Rrs 

plotframe.rrs <- data.frame("wave"=wavelength, "rrs.est.sse"=rrs.HL.est.sse,
                            "rrs.est.mcmc"=rrs.HL.est.mcmc,
                            "rrs.obs"=rrs.HL[j,])

mcmc.label = "theta[MCMC]== {4.96*',1.50,0.026'}"
mle.label = "theta[MLE]== {3.73*',1.44,0.103'}"
obs.label = "theta[obs]== {5.00*',1.50,0.020'}"

xmin = min(plotframe.rrs$wave); xmax= max(plotframe.rrs$wave); xstp=100
ymin= 0; ymax=max(plotframe.rrs$rrs.obs)+0.20*max(plotframe.rrs$rrs.obs);ystp= signif(ymax/5, digits = 1)
asp_rat <- (xmax-xmin)/(ymax-ymin)

g1 <- ggplot()  + geom_line(data = plotframe.rrs,aes(y=rrs.est.sse,color="xx1",x = wave),
                            size=1.3,show.legend = TRUE) +
  geom_line(data = plotframe.rrs,aes(y=rrs.est.mcmc,color="xx2",x = wave),
            size=1.3,show.legend = TRUE) +
  geom_line(data = plotframe.rrs,aes(y=rrs.obs,x = wave,color="xx5"),linetype="dashed", 
            size=1.3,show.legend = TRUE)+
  scale_colour_manual(labels = c(expression(paste(italic("R")["rs,model,MLE"])),
                                 expression(paste(italic("R")["rs,model,MAP"])),
                                 expression(paste(italic("R")["rs,actual"]))), 
                      values = c("blue","green","red")) +
  
  scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, xstp))  +
  scale_y_continuous(name =expression(paste(italic("R"),{}[rs],"(",lambda,",", 0^"-",")[", sr^-1,"]")) , limits = c(ymin, ymax),
                     breaks = seq(ymin, ymax, ystp))+ 
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
  
  #annotate("text",x=470, y= 0.0021, label = obs.label, parse=T, size=4)+
  #annotate("text",x=470, y= 0.0019, label = mle.label, parse=T, size=4)+
  #annotate("text",x=470, y= 0.0017, label = mcmc.label, parse=T, size=4)+
  
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
  ggsave(paste0("./SABER.inv.Rrs_HL_AM_model.png"), plot = g1,
         scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
}
