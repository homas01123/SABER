#===========================================================================
# saber_forward.R simulates the remote sensing reflectance given the wavelengths,
# the water components along with bathymetry and bottom reflectance (for shallow water)

#Author: Mr. Soham Mukherjee - PhD Student, Aquatel, UQAR
#===========================================================================
#rm(list=ls(all=TRUE))
require(dplyr)
require(readxl)
require(stats4)
require(MASS)
require(dglm)
library(fitdistrplus)
require(Riops)
require(Cops)
require(ggplot2)
#--------------------------------------------------------------------------
setwd("S:/Work/UQAR/Analysis/WISE/R_inverse_wasi")
snell_law <- function(view,sun){#Function to convert above water to under water geometry
  
  # Index of refrations (real)
  n_air= 1;  # air index of refration (real part)
  n_w= 1.33; # water index of refration (real part)
  
  # Angles from the water
  
  # from deg to rad
  view=view*(180/pi); # rad
  sun=sun*(180/pi);   # rad
  
  # angles inside the water in rad
  view_w= asin((n_air/n_w)*sin(view));   # rad
  sun_w= asin((n_air/n_w)*sin(sun));     # rad
  
  # Fresnel Law
  
  rho_L= (1/2)*abs(((sin(view-view_w)^2)/(sin(view+view_w)^2))+((tan(view-view_w)^2)/(tan(view+view_w)^2)))
  return(data.frame("view_w"=view_w, "sun_w"=sun_w, "rho_L"=rho_L))
}

Saber_forward <-  function(chl=4.96, acdom440=0.9322314, anap440=0.07, bbp.550=0.00726002,z=2, 
                           rb.fraction=fA.set, plot=FALSE, verbose=FALSE,realdata){
  ## OAC and IOP initialization
  
  # Read input params for forward modeling
  
  base.CDOM <- acdom440 # CDOM absorption at 440nm [m^-1]
  C_ph <- chl           # Phytoplankton concentation [mg/m^-3]
  base.NAP <- anap440   # NAP absorption at 440nm [m^-1]
  base.bbp <- bbp.550   # particulate backscatter at 550nm [m^-1]
  lambda <- wavelength  #Desired spectral range for simulation
  fA <- rb.fraction     #Aerial fraction of bottom albedo
  zB <- z               #bottom depth
  Rrs_obs.interp <- realdata #Actual AOP data for which we try to simulate the forward model
  #lambda <- wavelength
  
  if (verbose == TRUE) {
    print("Function arguement initialized") 
  }
  ################################################################################################
  ##                                    Bio-Optical Development                                 ##
  ################################################################################################
  
  #--------------------------------------------------------------------------
  ## Absorption
  #--------------------------------------------------------------------------
  ## Pure water absorption (1/m)
  
  # wavelenght range [190;4000] [nm]
  abs.water <- read.table("./input-spectra/abs_W.A", header = F) 
  wavelength <- abs.water$V1
  absorpt_W <-  abs.water$V2
  
  a_W <- rep(0,length(lambda))# abs. of pure water [1/m]
  a_W <-  Hmisc::approxExtrap(wavelength, absorpt_W,xout = lambda,method = "linear")$y
  
  if (verbose == TRUE) {
    print(paste0("Water absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm.")) 
  }
  
  ## Plankthon absorption (1/m)
  
  # load plankton absorption data
  A0_A1_PhytoPlanc <- read.table("./input-spectra/A0_A1_PhytoPlanc.dat")
  # extract the values from the table
  lam_p <- A0_A1_PhytoPlanc$V1
  a0_p <- A0_A1_PhytoPlanc$V2
  a1_p <- A0_A1_PhytoPlanc$V3
  
  a0 <- rep(0,length(lambda))# [m^2/mg]
  a1 <- rep(0,length(lambda))# [m^2/mg]
  
  a0 <- Hmisc::approxExtrap(lam_p, a0_p,xout = lambda,method = "linear")$y
  a1 <- Hmisc::approxExtrap(lam_p, a1_p,xout = lambda,method = "linear")$y
  
  # Compute the value of plankton absorption as function of the concentration and wavelength
  aph_440 <- 0.06*(C_ph)^0.65# [mg/m^3]
  abs_ph <- rep(0,length(lambda))
  
  for (i in 1:length(lambda)){
    
    abs_ph[i] <- (a0[i] + a1[i]*log(aph_440))*aph_440
    
  }
  abs_ph[abs_ph < 0] <- 0
  
  if (verbose == TRUE) {
    print(paste0("Planktonic absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
  }
  
  ## CDOM absorption coefficient [1/m]
  
  Ga_CDOM <- 1# [m^2/mg]
  Oa_CDOM <- 0
  S_CDOM <- 0.014
  abs_CDOM_440 <-  (Ga_CDOM*base.CDOM)+Oa_CDOM# [1/m], CDOM abs. coeff. at 440 [nm]
  
  abs_CDOM <- rep(0,length(lambda))
  
  for (i in 1:length(lambda)){
    
    abs_CDOM[i]  <- abs_CDOM_440*exp(-S_CDOM*(lambda[i] - 440))
    
  }
  
  if (verbose == TRUE) {
    print(paste0("CDOM absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
  }
  
  ## SPM absorption coefficient [1/m]
  
  Ga_X <- 1# [m^2/mg]
  Oa_X <- 0
  S_X <- 0.01160
  abs_X_440 <-  (Ga_X*base.NAP)+Oa_X# [1/m], SPM abs. coeff. at 440 [nm]
  
  abs_X <- rep(0,length(lambda))
  
  for (i in 1:length(lambda)){
    
    abs_X[i]  <- abs_X_440*exp(-S_X*(lambda[i] - 440))
    
  }
  
  if (verbose == TRUE) {
    print(paste0("NAP absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
  }
  
  #Plot the componenet specific absoprtion
  plotframe.abs <- data.frame("wave"=lambda,"abs.ph"= abs_ph, "abs.cdom"=abs_CDOM,
                              "abs.nap"=abs_X)
  xmin = min(plotframe.abs$wave); xmax= max(plotframe.abs$wave); xstp=100
  ymin= 0; ymax=max(plotframe.abs$abs.cdom);ystp= signif(ymax/5, digits = 2)
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  g <- ggplot()  + geom_line(data = plotframe.abs,aes(y=abs.ph,color="xx1",x = wave),
                             size=1.3,show.legend = TRUE) +
    geom_line(data = plotframe.abs,aes(y=abs.cdom,x = wave,color="xx2"), 
              size=1.3,show.legend = TRUE)+
    geom_line(data = plotframe.abs,aes(y=abs.nap,x = wave,color="xx3"), 
              size=1.3,show.legend = TRUE)+
    scale_colour_manual(labels = c(expression(paste(italic("a")[phi])),
                                   expression(paste(italic("a")["CDOM"])),
                                   expression(paste(italic("a")["NAP"]))), 
                        values = c("green","darkgoldenrod2","navyblue")) +
    #ggtitle(paste0(stationlist[i])) +
    scale_x_continuous(name = expression(paste("Wavelength(", lambda, ")[nm]")), limits = c(xmin, xmax), 
                       breaks = seq(xmin, xmax, xstp))  +
    scale_y_continuous(name =expression(paste(italic("a"),{}["t-w"],"(",lambda,",", 0^"-",")[", m^-1,"]")) , limits = c(ymin, ymax),
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
  g 
  
  if (plot == "TRUE") {
    ggsave(paste0("./SABER.absorption.model_chl=",C_ph ,"_acdom440=", base.CDOM,"_anap440=", base.NAP, ".png"), plot = g,
           scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
  }
  
  #--------------------------------------------------------------------------
  ## Scattering and back-scattering
  #--------------------------------------------------------------------------
  
  ## Pure water backscattering (1/m)
  
  if (type_case_water == 1) {# case1 water
    b1 <-  0.00144#  [1/m]
    
  } else {# case2 water
    if (type_case_water == 2) {
      b1 <-  0.00111#  [1/m]
    } else{
      print("Wrong water case type choosen")
    }
    
  }
  
  lambda1 <-  500# [nm]
  bb_W <- b1*(lambda/lambda1)^(-4.32)# [1/m]
  
  if (verbose == TRUE) {
    print(paste0("Water backscatter calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
  }
  
  # ## SPM backscattering coeff
  
  # bbxS_W <- 0.0086# [m^2/g] specfic backscattering according to WASI
  # bbxS <- ((bbxS_W)*(33.57*(10^(-6))))/g_size# [m^2/g] specfic backscattering as function of grain size
  # # b_ratio = 0.019;                                 # ratio between specific scattering and specific backscatteirng
  # # bxS = b_ratio*bbxS;                              # [m^2/g] specific scattering
  
  refexponent <- 0.46
  bb_x <- rep(0,length(lambda))# Backscattering coefficient for suspended particles [1/m]
  
  for (i in 1:length(lambda)){
    
    #print(lambda[i])
    bb_x[i] = base.bbp*((lambda[i]/550)^-(refexponent)) #Implement power model
    
  }
  if (verbose == TRUE) {
    print(paste0("NAP backscatter calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
  }
  #--------------------------------------------------------------------------
  ## Water system IOPs
  #--------------------------------------------------------------------------
  
  ## Total Absorption Coefficient (1/m)
  
  abs <-  a_W + abs_ph + abs_CDOM + abs_X
  #abs= a_W + a_tw;
  
  if (verbose == TRUE) {
    print(paste0("Total absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
  }
  
  ## Total Backscattering Coefficient (1/m)
  
  bb <-  bb_W + bb_x
  #bb= bb_W + b_btw;
  
  if (verbose == TRUE) {
    print(paste0("Total backscatter calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
  }
  
  ## Extinction Coefficient (1/m) and Single Back Scattering Albedo
  
  ext <-  abs+bb# [1/m] extinction coeff.
  omega_b <-  bb/ext# single back scattering albedo
  
  if (verbose == TRUE) {
    print(paste0("omega albedo calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
  }
  
  ################################################################################################
  ##                                      RT Model                                              ##
  ################################################################################################
  
  #--------------------------------------------------------------------------
  ## Remote sensing reflectance below the surface
  #--------------------------------------------------------------------------
  geometry <- snell_law(view = view, sun = sun)
  sun_w <- geometry$sun_w; view_w <- geometry$view_w; rho_L <- geometry$rho_L
  
  if (verbose == TRUE) {
    print("Viewing Geometry and Fresnel reflectance below surface calculated")
  }
  
  ## Remote Sensing Reflectance below the water surface
  
  if(type_case_water == "1"){
    #case 1
    
    f_rs <- 0.095 # [1/sr] 
  } else {
    if (type_case_water == "2") {
      #case 2
      
      f_rs <- 0.0512*(1 + (4.6659*omega_b) +(-7.8387*(omega_b^2)) + (5.4571*(omega_b^3)) )*(1+(0.1098/cos(sun_w)))*(1+(0.4021/cos(view_w)))# [1/sr]
    }
  }
  
  
  Rrs_below_deep <- f_rs*omega_b# [1/sr]
  
  if (type_Rrs_below == "deep") {
    
    Rrs_below <- Rrs_below_deep
    
    if (verbose == TRUE) {
      print(paste0("Subsurface (0^-)Rrs for deep water calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
  } else{
    
    if (type_Rrs_below == "shallow") {
      #Bottom Contribution
      
      # Reflection factors of bottom surface [1/sr]
      B0 <- 1/pi; B1 <- 1/pi; B2 <- 1/pi; B3 <- 1/pi; B4 <- 1/pi; B5 <- 1/pi; 
      BOTTOM <- c(B0,B1,B2,B3,B4,B5)
      
      # Bottom Albedo (costant)
      # wavelength range [350;900] [nm]
      bott0<- read.table("./input-spectra/Bott0const.R")
      wavebottom <- bott0$V1
      Bott0 <-  bott0$V2
      abott0 <- rep(0,length(lambda))
      abott0 <-  Hmisc::approxExtrap(wavebottom, Bott0, xout = lambda, method = "linear")$y
      
      # Bottom Albedo Sand
      # wavelenght range [350;1000] [nm]
      bott1 <- read.table("./input-spectra/Bott1SAND.R")
      wavebottom <- bott1$V1
      Bott1 <-  bott1$V2
      abott1 <- rep(0,length(lambda))
      abott1 <-  Hmisc::approxExtrap(wavebottom, Bott1, xout = lambda, method = "linear")$y
      
      # Bottom Albedo of fine-grained sediment
      # wavelenght range [350;900] [nm]
      bott2 <- read.table("./input-spectra/Bott2silt.R")
      wavebottom <- bott2$V1
      Bott2 <-  bott2$V2
      abott2 <- rep(0,length(lambda))
      abott2 <-  Hmisc::approxExtrap(wavebottom, Bott2, xout = lambda, method = "linear")$y
      
      # Bottom Albedo of green makrophyte "Chara contraria"
      # wavelenght range [350;900] [nm]
      bott3 <- read.table("./input-spectra/Bott3chara.R")
      wavebottom <- bott3$V1
      Bott3 <-  bott3$V2
      abott3 <- rep(0,length(lambda))
      abott3 <-  Hmisc::approxExtrap(wavebottom, Bott3, xout = lambda, method = "linear")$y
      
      # Bottom Albedo of green makrophyte "Potamogeton perfoliatus"
      # wavelenght range [350;900] [nm]
      bott4 <- read.table("./input-spectra/Bott4perfol.R")
      wavebottom <- bott4$V1
      Bott4 <-  bott4$V2
      abott4 <- rep(0,length(lambda))
      abott4 <-  Hmisc::approxExtrap(wavebottom, Bott4, xout = lambda, method = "linear")$y
      
      # Bottom Albedo of green makrophyte "Potamogeton pectinatus"
      # wavelenght range [350;900] [nm]
      bott5 <- read.table("./input-spectra/Bott5pectin.R")
      wavebottom <- bott5$V1
      Bott5 <-  bott5$V2
      abott5 <- rep(0,length(lambda))
      abott5 <-  Hmisc::approxExtrap(wavebottom, Bott5, xout = lambda, method = "linear")$y
      
      abott <- rbind(abott0, abott1, abott2, abott3, abott4, abott5)
      
      Bottom <-  matrix(nrow = nrow(abott), ncol = ncol(abott), 0)
      Rrs_Bottom <- matrix(nrow = nrow(abott), ncol = ncol(abott), 0)# Bottom remote sensing reflectance [1/sr]
      
      for (i in 1:length(fA)){
        Bottom[i,] <-  fA[i]*abott[i,]
        Rrs_Bottom[i,] <-  BOTTOM[i]* Bottom[i,] #fA(i)*abott(:,i);
      }
      
      Bottom <- colSums(Bottom)
      Rrs_Bottom <- colSums(Rrs_Bottom)# [1/sr]
      
      if (verbose == TRUE) {
        print(paste0("Aerial fraction scaled bottom albedo generated from ", min(lambda), " nm to ", max(lambda), " nm."))
      }
      
      # Attenuation Coefficients
      if (type_case_water == 1) {
        
        k0 <- 1.0395 #case 1
      } else {
        if (type_case_water == 2) {
          
          k0 <- 1.0546 #case 2
        }
      }
      
      Kd <- k0*(ext/cos(sun_w))
      kuW <- (ext/cos(view_w))*((1+omega_b)^3.5421)*(1-(0.2786/cos(sun_w)))
      kuB <- (ext/cos(view_w))*((1+omega_b)^2.2658)*(1-(0.0577/cos(sun_w)))
      
      if (verbose == TRUE) {
        print(paste0("Attenuation coefficients calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
      }
      
      
      #Final calculation for shallow Rrs
      Ars1 <- 1.1576; Ars2 <- 1.0389 #Paramteric coeffs for shallow water
      Rrs_below_shallow <-  Rrs_below_deep*(1-(Ars1*exp(-zB*(Kd+kuW)))) + Ars2*Rrs_Bottom*exp(-zB*(Kd+kuB)) 
      Rrs_below <- Rrs_below_shallow
      
      if (verbose == TRUE) {
        print(paste0("Subsurface (0^-) Rrs in Shallow water calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
      }
    }
  }
  
  #--------------------------------------------------------------------------
  ## Remote sensing reflectance above the surface
  #--------------------------------------------------------------------------
  # Extraterrestrial solar irradiance [mW/m^2 nm]
  E0 <-  read.table("./input-spectra/E0.txt", header = F)
  E0.wavelength <-  E0$V1
  E0.ett <-  E0$V2
  E0 <- rep(length(lambda), 0)
  
  E0 <-  Hmisc::approxExtrap(x =E0.wavelength, y =E0.ett, xout = lambda, method = "linear")$y
  
  # Oxygen absorption [1/cm]
  absO2 <- read.table("./input-spectra/absO2.A", header = F)
  absO2.wavelength <-  absO2$V1
  absO2.oxy <-  absO2$V2
  abs_O2 <- rep(length(lambda), 0)
  
  abs_O2 <-  Hmisc::approxExtrap(x =absO2.wavelength, y =absO2.oxy, xout = lambda, method = "linear")$y
  
  # Ozone absorption [1/cm]
  absO3 <- read.table("./input-spectra/absO3.A", header = F)
  absO3.wavelength <-  absO3$V1
  absO3.oxy <-  absO3$V2
  abs_O3 <- rep(length(lambda), 0)
  
  abs_O3 <-  Hmisc::approxExtrap(x =absO3.wavelength, y =absO3.oxy, xout = lambda, method = "linear")$y
  
  # Water vapour absorption [1/cm]
  absWV <- read.table("./input-spectra/absWV.A", header = F)
  absWV.wavelength <-  absWV$V1
  absWV.wv <-  absWV$V2
  abs_WV <- rep(length(lambda), 0)
  
  abs_WV <-  Hmisc::approxExtrap(x =absWV.wavelength, y =absWV.wv, xout = lambda, method = "linear")$y
  
  # angles from deg to rad
  
  view_rad=view*(180/pi); # rad
  sun_rad=sun*(180/pi);   # rad
  
  # Downwelling Irradiance [mW/m^2 nm]
  M <- 1/(cos(sun_rad)+(0.50572*((90+ 6.079975-sun)^(-1.253))))
  M1 <-  (M*P)/1013.25
  Moz <-  1.0035/ (((cos(sun_rad)^2)+0.007)^0.5)
  
  
  # Air mass type
  if (type_case_water == "1") { 
    AM <- 1
  } else{
    
    AM <- 1
  }
  
  omega_a <- ((-0.0032*AM) + 0.972)*exp(RH*3.06*(10^-4))
  
  Ha <- 1; V <- 15# [km] aerosol scale height and horizontal visibility
  beta <- 3.91*(Ha/V)
  tau_a <-  beta*((lambda/550)^(-alpha))
  
  Tr <-   exp(-M1/((115.6406*(lambda^(4)))-(1.335*(lambda^(2)))))
  Taa <-  exp(-(1-omega_a)*tau_a*M)
  Tas <-  exp(-omega_a*tau_a*M)
  Toz <-  exp(-abs_O3*Hoz*Moz)
  To <-   exp((-1.41*abs_O2*M1)/((1+(118.3*abs_O2*M1))^0.45))
  Twv <-  exp((-0.2385*abs_WV*WV*M)/((1+(20.07*abs_WV*WV*M))^0.45))
  
  
  B3 <- 0.82 - (0.1417*alpha); B1 <-  B3*(1.459 +(B3*(0.1595+(0.4129*B3)))); B2 <-  B3*(0.0783 +(B3*(-0.3824-(0.5874*B3))))
  Fa <-  1-(0.5*exp((B1+(B2*cos(sun_rad)))*cos(sun_rad)))
  
  Edd <-  E0*Tr*Taa*Tas*Toz*To*Twv*cos(sun_rad)
  Edsr <- (1/2)*E0*(1-(Tr^(0.95)))*Taa*Tas*Toz*To*Twv*cos(sun_rad)
  Edsa <-  E0*(Tr^(1/2))*Taa*(1-Tas)*Toz*To*Twv*cos(sun_rad)*Fa
  Eds <-  Edsr+Edsa# diffuse downwelling irradiance (sum of Rayleigh and aerosol)
  
  Ed <-  (f_dd*Edd) + (f_ds*Eds)# downwelling irradiance
  
  # Sky Radiance
  Ls <-  (g_dd*Edd) + (g_dsr*Edsr) + (g_dsa*Edsa)# [mW/sr m^2 nm]
  
  # Remote sensing reflectance above the surface [1/sr]
  Rrs_above <-  rho_L*(Ls/Ed)
  
  if (verbose == TRUE) {
    print(paste0("Skyglint calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
  }
  #--------------------------------------------------------------------------
  ## Final Remote sensing reflectance computation
  #--------------------------------------------------------------------------
  sigma <- 0.03; nW <- 1.33; rho_U <- 0.54; Q <- 5#[sr]
  
  if (type_Rrs_water == "below_surface") {
    
    Rrs <-  Rrs_below
    
    if (verbose == TRUE) {
      print(paste0("Subsurface (0^-)Rrs calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
  } else {
    if (type_Rrs_water == "above_surface_with_glint") {
      
      Rrs <-  (((1-sigma)*(1-rho_L)/(nW^2))* (Rrs_below/(1-rho_U*Q*Rrs_below)))+ Rrs_above
      
      if (verbose == TRUE) {
        print(paste0("Above surface (0^+)Rrs with modelled glint calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
      }
      
    } else {
      if (type_Rrs_water == "above_surface_only") {
        
        Rrs <-  (((1-sigma)*(1-rho_L)/(nW^2))* (Rrs_below/(1-rho_U*Q*Rrs_below)))
        
        if (verbose == TRUE) {
          print(paste0("Above surface (0^+)Rrs calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
        
      }
    }
    }
  }
  
  #Compare actual vs modelled Rrs
  #Rrs_obs.interp <- Hmisc::approxExtrap(Rrs_obs_wl, Rrs_obs, xout = lambda, method = "linear")$y
  
  plotframe.rrs <- data.frame("wave"=lambda, "rrs.est"=Rrs, "rrs.obs"=Rrs_obs.interp)
  xmin = min(plotframe.rrs$wave); xmax= max(plotframe.rrs$wave); xstp=100
  ymin= 0; ymax=max(plotframe.rrs$rrs.obs);ystp= signif(ymax/5, digits = 1)
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
    ggsave(paste0("./SABER.forward.Rrs_chl=",C_ph ,"_acdom440=", base.CDOM,"_anap440=", base.NAP, ".png"), plot = g1,
           scale = 1.5, width = 4.5, height = 4.5, units = "in",dpi = 300)
  }
  
  #--------------------------------------------------------------------------
  ## Objective function for the inverse mode
  #--------------------------------------------------------------------------
  #Calculate minimizable residual fuction
  Res <- sum((Rrs_obs.interp-Rrs)^2) #eucledian space
  
  #Calculate residual as function of wavelength
  Res.spectral <- (Rrs - Rrs_obs.interp)*100/Rrs_obs.interp
  return(list(data.frame("wavelength"=lambda, "Rrs"=Rrs,
                         "p.bias"=Res.spectral),"ss.residual"=Res,"method"=c("SSR eucledian")))
  #return(Res)
}
