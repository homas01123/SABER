Rrs_sicf <- function(dg_comsposite=TRUE, dg_443,
                             c_chl, abs_cdom_443, abs_nap_443,
                             
                             Ed_path = "./data/input-spectra/Ed_HL.csv",
                             
                             use_analytic_Ed = TRUE,
                             sunzen_Ed = 60, lat_Ed = 49, lon_Ed = -68,
                             date_time_Ed = "2019-08-18 20:50 GMT",
                             
                             coeff_x=c(0.0992,0.40,0.078), 
                             use_sicf_rad = F,
                             wavelength=seq(400,800,10), phi_f=0.01) {
  lambda = wavelength
  
  ##approximated model :  A. Gilerson, J. Zhou, S. Hlaing, I. Ioannou, J. Schalles, 
  #B. Gross, F. Moshary, and S. Ahmed, "Fluorescence component in the reflectance spectra 
  #from coastal waters. Dependence on water composition," Opt. Express 15, 
  #15702-15721 (2007)
  
  
  if (use_analytic_Ed == FALSE) {
    
    print("SICF: NOTE: Pre-Defined Ed is used; results might be erroneous")
    Ed_sim = read.csv(file = Ed_path, header = TRUE, skip = 9)
    Ed_sim_interp = Hmisc::approxExtrap(x = Ed_sim$Wavelength, y = log(Ed_sim$Ed_total.W.m.2.nm.),
                                        xout = wavelength, method = "linear")$y
    Ed_sim_interp = exp(Ed_sim_interp)*100 #extrapolate and change units
    
    Ed_sim_interp_0m = 0.96*Ed_sim_interp #Convert to underwater E0
    Ed_0m_sicf = Ed_sim_interp_0m
    
  } else {
    print("SICF: NOTE: Gregg & Carder 1990 is used to estimate Ed")
    
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
      print("SICF: Sun Zenith was not provided, calculated from Geometry")
      sunzen_Ed = Cops::GreggCarder.sunang(rad = 180/pi, iday = jday_no, 
                                           xlon = lon_Ed, ylat = lat_Ed, hr = time_dec)
    }
    
    # test_Ed = GreggCarder.f.modified(the = sunzen_Ed, #Calculate the Ed following Gregg & Carder 1990
    #                                  lam.sel = lambda, hr = time_dec,
    #                                  jday = jday_no, rlon = lon_Ed, rlat = lat_Ed, debug = F)
    
    tryCatch({
      #Calculate the Ed following Gregg & Carder 1990
      test_Ed = GreggCarder.f.modified(the = sunzen_Ed, 
                                       lam.sel = lambda, hr = time_dec,
                                       jday = jday_no, rlon = lon_Ed, rlat = lat_Ed, debug = F)
      
      if (all(is.na(test_Ed))) {
        cat(paste0("Sun is below horizon with given sun-earth geometry, reseting to default"))
        
        sunzen_Ed = -999; lat_Ed = 49; lon_Ed = -68;
        date_time_Ed = "2019-08-18 20:50 GMT"
        
        
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
    Ed_0m_sicf = Ed0.0m*100
  }
  
  if (dg_comsposite == TRUE) {
    
    abs_cdom_nap_443 = dg_443
    print("SICF: Composite absorption for a_dg is obtained from SA subroutine, i.e. QAA")
    
  } else {
    abs_cdom_nap_443 = abs_cdom_443 + abs_nap_443
    print("SICF: Seperate absorption for a_dg is provided from user inputs")
  }
  
  
  Lf_685 <- phi_f * coeff_x[1] * c_chl / (1 + coeff_x[2] * abs_cdom_nap_443 + coeff_x[3] * c_chl)
  
  #Lf_685= 0.001 * Lf_685  #convert the unit 1/??m to 1/nm; Ed: (W/m^-2 nm-1)  v.s. 
  #Lf:  (W/m^-2 ??m-1) * (mg/m3)-1
  
  # The spectral shape of fluorescence, Lf, was modeled as a summation of two Gaussian 
  #spectral profile  centered at 685 nm, having a full width at half maximum (FWHM) of 25 nm
  
  Lf_distribute <- function(wv) {
    exp(-4 * log(2) * ((wv - 685)/25)^2) + 0.3 * exp(-4 * log(2) * ((wv - 730)/50)^2)
  }
  
  # Calculate Lf based on the spectral profile and the wavelength
  Lf <- Lf_685 * Lf_distribute(wavelength)
  Rrs_sicf = Lf/(Ed_0m_sicf)
  if (use_sicf_rad == FALSE) {
    
    sicf = data_frame("wavelength" = wavelength#[1:length(wavelength)-1]
                      , "sicf" = Rrs_sicf)
    
  } else {
    sicf = data_frame("wavelength" = wavelength#[1:length(wavelength)-1]
                      , "sicf" = Lf)
  }
  
  print(paste0("Subsurface (0^-) Rrs equivalent to SICF is calculated"))
  
  return(sicf)
  
}
