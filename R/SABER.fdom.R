



Rrs_fDOM <- function(use_qaa_adg = T, rrs_input_for_qaa = obsdata- obsdata_sicf,
                     use_qaa_spectral_slope = T,
                     use_manual_spectral_slope = F,
                     manual_spectral_slope_vals = c("s_g"=0.014, "s_d"=0.003, "gamma"=1),
                     
                     dg_443_val=bgc_params_rb$abs_invivo$a_cdom + bgc_params_rb$abs_invivo$a_nap,
                     #abs_cdom_443, abs_nap_443,
                     
                     Ed_path = "./data/input-spectra/Ed_HL.csv",
                     
                     use_analytic_Ed = TRUE,
                     sunzen_Ed = 60, lat_Ed = 49, lon_Ed = -68,
                     date_time_Ed = "2019-08-18 20:50 GMT",
                     
                     wavelength=seq(400,800,10), phi_f=0.01, use_fDOM_rad = F,
                     plot = T){
  
  print("Calculating the discretized wavelength redistribution function...")
  #Implement the wavelength redistribution functions with discretization 
  #for inelastic scattering 
  wave_redist_inel = wave.redist.inelastic(quant_y_phi = phi_f, include.chl = T,
                                           include.CDOM = T, wavelength = wavelength)
  
  #Plot the discretized wavelength redistribution functions for inelastic scattering 
  if (plot == TRUE) {
    plot(wavelength[-length(wavelength)], wave_redist_inel[["sicf"]][1,], pch=19, type="l", 
         lwd=2.5, col="green",
         ylab = expression(paste("Fluoroscence Efficieny (1/nm)")), xlab = "Wavelength" )
    
    lines(wavelength[-length(wavelength)], wave_redist_inel[["fDOM"]][1,], col="goldenrod2", lwd=2.5)
    
    legend(x = "topleft",          # Position
           legend = c("SICF","fDOM"),  # Legend texts
           lty = c(1, 1),           # Line types
           col = rev(c("goldenrod2", "green")),           # Line colors
           lwd = 2)
  }
  
  #=====================================================================
  #Obtain the fDOM scattering function components
  #=====================================================================
  
  waveb <-wavelength
  lambda <- wavelength
  Nwave <- length(waveb) - 1
  
  WRF_CDOM = wave_redist_inel$fDOM # extract the wave discretized quantum yield shape for fDOM
  
  ###### Scalar irradiance calculation #####
  
  #---------------------------------------------------------------------------
  #Get Ed0 from HL simulated Ed data
  #---------------------------------------------------------------------------
  if (use_analytic_Ed == FALSE) {
    
    #Ed_path = Ed_fDOM_path 
    
    print("NOTE: Pre-Defined Ed is used; results might be erroneous")
    
    Ed_sim = read.csv(file = Ed_path, header = TRUE, skip = 9)
    Ed_sim_interp = Hmisc::approxExtrap(x = Ed_sim$Wavelength, y = log(Ed_sim$Ed_total.W.m.2.nm.),
                                        xout = wavelength, method = "linear")$y
    Ed_sim_interp = exp(Ed_sim_interp)*100 #extrapolate and change units
    
    Ed_sim_interp_0m = 0.96*Ed_sim_interp #Convert to underwater E0
    
    sun_view = snell_law(sun = sunzen_below(sun_zen_aove = sun_above), view = view) #calculate sun zenith
    
    E0_sim_interp_scalar = Ed_sim_interp_0m/(cos(sun_view$sun_w)) #Convert to scalar irradiance
    
    E0_0m = E0_sim_interp_scalar
    Ed_0m = Ed_sim_interp_0m
    
  } else {
    
    #---------------------------------------------------------------------------
    #Calculate the Ed using Gergg & Carder model
    #---------------------------------------------------------------------------
    
    print("NOTE: Gregg & Carder 1990 is used to estimate Ed")
    
    Cops::GreggCarder.data()
    
    library(lubridate)
    library(dplyr)
    
    x = as_datetime(as.POSIXct(date_time_Ed))
    
    jday_no = yday(x) #Get Jullien Day no.
    time_no = format(x, "%T") #Get time in UTC
    
    time_dec = sapply(strsplit(as.character(time_no),":"), #Convert time in decimal format
                      function(x) {
                        x <- as.numeric(x)
                        x[1]+x[2]/60
                      }
    )
    
    if (sunzen_Ed < 0) {
      sunzen_Ed = Cops::GreggCarder.sunang(rad = 180/pi, iday = jday_no, 
                                           xlon = lon_Ed, ylat = lat_Ed, hr = time_dec)
      print(paste0("Sun Zenith was not provided,", sunzen_Ed," is calculated from Geometry"))
    }
    #browser()
    
    tryCatch({
      #Calculate the Ed following Gregg & Carder 1990
      #print("Entered try catch")
      test_Ed = GreggCarder.f.modified(the = sunzen_Ed, 
                                       lam.sel = lambda, hr = time_dec,
                                       jday = jday_no, rlon = lon_Ed, rlat = lat_Ed, debug = F)
      
      if (all(is.na(test_Ed))) {
        cat(paste0("Sun is below horizon with given sun-earth geometry, reseting to default"))
        
        sunzen_Ed = -999; lat_Ed = lat_Ed; lon_Ed = lon_Ed;
        date_time_Ed = date_time_Ed
        
        
        library(lubridate)
        library(dplyr)
        
        x = as_datetime(as.POSIXct(date_time_Ed))
        
        jday_no = yday(x) #Get Jullien Day no.
        time_no = format(x, "%T") #Get time in UTC
        
        time_dec = sapply(strsplit(as.character(time_no),":"), #Convert time in decimal format
                          function(x) {
                            x <- as.numeric(x)
                            x[1]+x[2]/60
                          }
        )
        
        test_Ed = GreggCarder.f.modified(the = sunzen_Ed, 
                                         lam.sel = lambda, hr = time_dec,
                                         jday = jday_no, rlon = lon_Ed, rlat = lat_Ed, debug = F)
      } }, 
      error = function(e){cat("ERROR in Irradiance model :",conditionMessage(e), "\n")}
    )
    
    
    Ed0.0p = test_Ed$Ed #Total Ed at 0+
    Ed0_dir.0p = test_Ed$Edir #Direct Ed at 0+
    Ed0_dif.0p = test_Ed$Edif #Diffused Ed at 0+
    
    rhoF <- GreggCarder.sfcrfl(rad = 180.0/pi, theta=sunzen_Ed, ws=5) #Calculate fresnel reflectance
    sun_view_op = snell_law(view = 0 , sun = sunzen_Ed) #calculate sun angle 0+ in radian
    
    #Translate Ed0+ to Ed0-
    Ed0.0m = (Ed0_dir.0p *(1 - rhoF$rod)) + # direct
      (Ed0_dif.0p * (1 - rhoF$ros))   # diffuse
    
    #Convert it to scalar irradiance
    sun_view_om = snell_law(view = 0 , sun = sunzen_below(sunzen_Ed))
    Ed0m_0 = (Ed0.0m*100) / cos(sun_view_om$sun_w)
    
    #Plot the simulated Ed profiles
    if (plot == TRUE) {
      plot(lambda, Ed0.0p*100, col="black", pch=19, ylab = "Ed", ylim=c(0,1.25*max(Ed0m_0)))
      lines(lambda, Ed0.0m*100, col="red", lwd=2.5, lty="dashed")
      lines(lambda, Ed0m_0, col="navyblue", lwd=2.5)
      
      legend(x = "bottom",          # Position
             legend = c("Ed0+","Ed0-", "E0"),  # Legend texts
             lty = c(1, 2, 1),           # Line types
             col = c("black", "red", "navyblue"), # Line colors
             lwd = 2, bty = "n")
      
    }
    
    E0_0m = Ed0m_0
    Ed_0m = Ed0.0m*100
    
  }
  
  ###### absorption of CDOM calculation #####
  
  if ((use_qaa_adg == TRUE & type_Rrs_below == "shallow") & (is.na(dg_443_val))) {
    stop("fDOM: The water type is shallow, QAA will obtain wrong a_CDM and s_CDM, please enter a_CDM(443) value manually")
  }
  
  if ((use_qaa_adg == TRUE & type_Rrs_below == "shallow") & (!is.na(dg_443_val))) {
    print("fDOM: The water type is shallow, QAA will obtain wrong a_CDM and s_CDM, thus, user given value of adg(443) will be used")
    abs_CDM_440 = dg_443_val
  }
  
  if ((use_qaa_adg == TRUE & type_Rrs_below == "deep") & (!is.na(dg_443_val))) {
    print("fDOM: The water type is deep,but user given value of adg(443) will be used")
    abs_CDM_440 = dg_443_val
  }
  
  if ((use_qaa_adg == FALSE & type_Rrs_below == "deep") & (!is.na(dg_443_val))) {
    print("fDOM: The water type is deep and user given value of adg(443) will be used")
    abs_CDM_440 = dg_443_val
  }
  
  if ((use_qaa_adg == FALSE & type_Rrs_below == "shallow") & (!is.na(dg_443_val))) {
    print("fDOM: The water type is shallow and user given value of adg(443) will be used")
    abs_CDM_440 = dg_443_val
  }
  
  if ((use_qaa_adg == TRUE & type_Rrs_below == "deep") & (is.na(dg_443_val))) {
    print("fDOM: The water type is deep and no adg(443) value is provided, hence QAA will be used")
    qaa_op = QAA.v5(waves = wavelength, Rrs = rrs_input_for_qaa)
    abs_CDM_440 =  qaa_op$a_dg_443
  }
  
  
  if (use_qaa_spectral_slope == TRUE & use_manual_spectral_slope == TRUE) {
    stop("Both QAA derived slope calculation and manual slope calculation cannot be set TRUE")
    
  }
  
  if (use_qaa_spectral_slope == FALSE & use_manual_spectral_slope == FALSE) {
    S_CDM <- 0.017 #<< Model Default >>
    print(paste0("The spectral slope for CDOM + NAP is kept constant as:", S_CDM))
    
  }
  
  if (use_qaa_spectral_slope == TRUE & use_manual_spectral_slope == FALSE) {
    #parametric formula to retrieve spectral slope of CDOM + NAP
    S_CDM = 0.015 + (0.002/(0.6 + (rrs_input_for_qaa[which.min(abs(wavelength - 443))]/rrs_input_for_qaa[which.min(abs(wavelength - 555))])))
    print(paste0("fDOM: The spectral slope for CDOM + NAP is calculated as:", S_CDM))
    
  }
  
  if (use_qaa_spectral_slope == FALSE & use_manual_spectral_slope == TRUE) {
    S_CDM = as.numeric(manual_spectral_slope_vals["s_g"] + manual_spectral_slope_vals["s_d"])
    print(paste0("The spectral slope for CDOM + NAP is supplied by user as:", S_CDM))
    
  }
  
  
  abs_CDM <- rep(0,length(wavelength))
  
  for (i in 1:length(wavelength)){
    
    abs_CDM[i]  <- abs_CDM_440*exp(-S_CDM*(wavelength[i] - 440))
    
  }
  
  
  # if ( use_true_IOPs == TRUE ) {
  #   qaa_op = QAA.v5(waves = wavelength, Rrs = Rrs_obs.interp)
  #   abs_CDM_440 =  qaa_op$a_dg_443
  #   
  #   if (type_Rrs_below == "deep") {
  #     if (slope.parametric == TRUE) {
  #       #parametric formula to retrieve spectral slope of CDOM + NAP
  #       S_CDM = 0.015 + (0.002/(0.6 + (Rrs_obs.interp[which.min(abs(wavelength - 443))]/Rrs_obs.interp[which.min(abs(wavelength - 555))])))
  #       print(paste0("fDOM: The spectral slope for CDOM + NAP is calculated as:", S_CDM))
  #     } else{
  #       
  #       if (use_manual_spectral_slope_vals  == TRUE) {
  #         S_CDM = as.numeric(manual_spectral_slope_vals["s_g"] + manual_spectral_slope_vals["s_d"])
  #         print(paste0("The spectral slope for CDOM + NAP is supplied by user as:", S_CDM))
  #       } else {
  #         S_CDM <- 0.017 #<< Model Default >>
  #         print(paste0("The spectral slope for CDOM + NAP is kept constant as:", S_CDM))
  #       }
  #       
  #     }
  #     
  #     
  #   } else {
  #     print("fDOM: The water type is shallow, QAA will obtain wrong s_CDOM, thus const. 0.017 is used")
  #     S_CDM <- 0.017 #<< USER INPUT >>
  #   }
  #   
  #   abs_CDM <- rep(0,length(wavelength))
  #   
  #   for (i in 1:length(wavelength)){
  #     
  #     abs_CDM[i]  <- abs_CDM_440*exp(-S_CDM*(wavelength[i] - 440))
  #     
  #   }
  #   
  # } else {
  #   #abs_CDM <- rep(0,length(wavelength))
  #   abs_CDM = abs_CDM_backup
  #   
  # }
  
  sum0 = matrix(0, ncol = Nwave, nrow = Nwave) #stores the wavelength specific integration of fDOM 
  
  #######Do the integration with respect to wavelength of the final source function of fDOM#####
  # the equation to be integrated is:a_CDOM*f_F(ex, em)
  
  for (jwave in 2:Nwave) { #loop for excitation (wavelength_ex)
    
    for (iwave in 1:(jwave-1)) { #loop for emission (wavelength_exm)
      
      Eozi <- E0_0m[iwave] # retrieve E0 at wavelength_ex
      
      absCDOM <- abs_CDM[iwave] # retrieve a_CDOM at wavelength_ex
      
      sum0[iwave,jwave] <- sum0[iwave,jwave] + (Eozi * absCDOM * WRF_CDOM[iwave, jwave]) #Integrate
      
    }
  }
  
  fDOM_rad = sum0[1,]/(4*pi) #multiply with fluoroscence anisotropic phase function
  print("fDOM equivalent radiance is calculated")
  
  #fDOM_rrs = fDOM_rad / Ed_interp[-length(Ed_interp)]
  #fDOM_rrs = fDOM_rad / Ed0_m_interp[-length(Ed0_m_interp)]
  
  fDOM_rrs = fDOM_rad / Ed_0m[-length(Ed_0m)]
  
  if (use_fDOM_rad == FALSE) {
    
    fdom = data_frame("wavelength" = wavelength[1:length(wavelength)-1], "fDOM" = fDOM_rrs)
    
  } else {
    fdom = data_frame("wavelength" = wavelength[1:length(wavelength)-1], "fDOM" = fDOM_rad)
  }
  
  print(paste0("Subsurface (0^-) Rrs equivalent to fDOM is calculated"))
  
  return(fdom)
  
  
}

