#=====================================================================================================
# Saber_forward_fast.R simulates the remote sensing reflectance given the wavelengths,
# the water components along with bathymetry and bottom reflectance (for shallow water).

# This code has all the fastest modes of bio-optical parametrizations for IOPs. For All B-OPT models,
# refer to Saber.forward.final().

# The SA algorithm is excluded for inelastic scattering equivalent Rrs. Refer to Saber.forward.final()
# for the inelastic scattering support.

#Author: Mr. Soham Mukherjee - PhD Student, Aquatel, UQAR
#=====================================================================================================
#rm(list=ls(all=TRUE))
library(dplyr)
library(readxl)
library(stats4)
library(MASS)
library(dglm)
library(fitdistrplus)
library(Riops)
library(Cops)
library(ggplot2)
#--------------------------------------------------------------------------

# #Function to convert above water to under water geometry
# snell_law <- function(view,sun){
#   
#   # Index of refrations (real)
#   n_air= 1;  # air index of refration (real part)
#   n_w= 1.33; # water index of refration (real part)
#   
#   # Angles from the water
#   
#   # from deg to rad
#   view=view*(180/pi); # rad
#   sun=sun*(180/pi);   # rad
#   
#   # angles inside the water in rad
#   view_w= asin((n_air/n_w)*sin(view));   # rad
#   sun_w= asin((n_air/n_w)*sin(sun));     # rad
#   
#   # Fresnel Law
#   
#   rho_L= (1/2)*abs(((sin(view-view_w)^2)/(sin(view+view_w)^2))+((tan(view-view_w)^2)/(tan(view+view_w)^2)))
#   return(data.frame("view_w"=view_w, "sun_w"=sun_w, "rho_L"=rho_L))
# }

#SABER FORWARD Model FAST version
Saber_forward_fast <-  function(use_true_IOPs = T, #Set TRUE if actual spectral IOPs exist
                                 a_non_water_path = "./data/rb_retrieve_demo_a.csv", #a path
                                 bb_non_water_path = "./data/rb_retrieve_demo_bb.csv", #bb path
                                 
                                 
                                 chl=4.96, #must be input if SICF=TRUE and use_true_IOPs = F
                                 
                                 
                                 a_dg = 1,  # provide value of acdom440/443 + anap440/443 (adg443)
                                
                                 
                                 bbp.550=0.00726002, #bbp value at 550/555 nm
                                 
                                 slope.parametric = TRUE, #Note, only possible if dg_composite = TRUE 
                                 
                                 use_manual_slope = FALSE, #use manual spectral slopes
                                 manual_slope = c("s_g"=0.015, "s_d"=0.01160, "gamma"=1), #Values of
                                 #manual spectral slopes. must be provided in a named vector as shown
                                 
                                 
                                 z=2, #bottom depth
                                 rb.fraction = fA.set, #aerial fraction of bottom types
                                 
                                 verbose = FALSE, #plot diagnostics and console output
                                 
                                wavelength = seq(400,800,10),
                                Rrs_input_for_slope, plot = T
                                 
){
  
  if (slope.parametric == TRUE & use_manual_slope == TRUE) {
    stop("both automatic and manual slope can't be set TRUE")
  }
  
  ## OAC and IOP initialization
  
  # Read input params for forward modeling
  
  base.CDM <- a_dg # CDM absorption at 440nm [m^-1]
  C_ph <- chl           # Phytoplankton concentration [mg/m^-3]
  base.bbp <- bbp.550   # particulate backscatter at 550nm [m^-1]
  lambda <- wavelength  #Desired spectral range for simulation
  Rrs_obs.interp <- Rrs_input_for_slope
  
  if (verbose == T) {
    
    cat(paste0("\033[0;36m","=====================================================================================================","\033[0m","\n"))
    cat(paste0("\033[0;36m","######################################### SIMULATION BEGINS #########################################","\033[0m","\n"))
    cat(paste0("\033[0;36m","=====================================================================================================","\033[0m","\n"))
    
  }
  
  ##If in situ IOPs are provided
  if (use_true_IOPs == TRUE) {
    
    if (!is.null(a_non_water_path) & !is.null(bb_non_water_path)) {
      
      a_non_water_df = read.csv(file = a_non_water_path, header = T)
      bb_non_water_df = read.csv(file = bb_non_water_path, header = T)
      
      if (is.data.frame(a_non_water_df) & is.data.frame(bb_non_water_df)) {
        
        a_non_water = a_non_water_df$at_w
        a_non_water_wave = a_non_water_df$wave
        
        bb_non_water = bb_non_water_df$bbp
        bb_non_water_wave = bb_non_water_df$wave
        
      } else {
        stop("The IOPs must be provided as data frames with <<wave>> as wavelength & <<a/bb>> as value")
      }
      
    }
  }
  
  #Shallow parameters
  fA <- rb.fraction     #Aerial fraction of bottom albedo
  zB <- z               #bottom depth
  
  if (verbose == TRUE) {
    print("Function arguement initialized") 
  }
  ################################################################################################
  ##                                    Bio-Optical Development                                 ##
  ################################################################################################
  
  #--------------------------------------------------------------------------
  ## Pure water IOPs
  #--------------------------------------------------------------------------
  ## absorption (1/m)
  
  # wavelength range [190;4000] [nm]
  abs.water <- read.table("./data/input-spectra/abs_W.A", header = F) 
  wavelength <- abs.water$V1
  absorpt_W <-  abs.water$V2
  
  a_W <- rep(0,length(lambda))# abs. of pure water [1/m]
  a_W <-  Hmisc::approxExtrap(wavelength, absorpt_W,xout = lambda,method = "linear")$y
  
  if (verbose == TRUE) {
    print(paste0("Water absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm.")) 
  }
  
  ## backscattering (1/m)
  
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
  
  #--------------------------------------------------------------------------
  ##Optically Active Constituents IOPs
  #--------------------------------------------------------------------------
  
  if (use_true_IOPs == FALSE) {
    
    if (verbose == TRUE) {
      print(paste0("Full spectral IOPs not provided, bio-optical models are used"))
    }
    
    ## Plankton absorption (1/m)
    
    # load plankton absorption data
    A0_A1_PhytoPlanc <- read.table("./data/input-spectra/A0_A1_PhytoPlanc.dat")
    
    #===============================================================================
    #Replace with the below for package environment for upper line
    #system.file("data", "input-spectra", "A0_A1_PhytoPlanc.dat", package = "SABER")
    #===============================================================================
    
    # extract the values from the table
    lam_p <- A0_A1_PhytoPlanc$V1
    a0_p <- A0_A1_PhytoPlanc$V2
    a1_p <- A0_A1_PhytoPlanc$V3
    
    a0 <- rep(0,length(lambda))# [m^2/mg]
    a1 <- rep(0,length(lambda))# [m^2/mg]
    
    a0 <- Hmisc::approxExtrap(lam_p, a0_p,xout = lambda,method = "linear")$y
    a1 <- Hmisc::approxExtrap(lam_p, a1_p,xout = lambda,method = "linear")$y
    
    aph_440 <- 0.06 * (C_ph)^0.65  # [mg/m^3] #Prieur & Satyendranath (1981)
    
    abs_ph <- sapply(1:length(lambda), function(i) (a0[i] + a1[i] * log(aph_440)) * aph_440) #Vectorization
    
    
    abs_ph[abs_ph < 0] <- 0
    
    
    if (verbose == TRUE) {
      print(paste0("Planktonic absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
    ## CDOM+NAP absorption coefficient [1/m]
    
    
    Ga_CDOM = 1
    Oa_CDOM = 0
    
    abs_CDM_440 <-  (Ga_CDOM*base.CDM)+Oa_CDOM# [1/m], CDOM abs. coeff. at 440 [nm]
    
    if (slope.parametric == TRUE) {
      
      if (type_Rrs_below == "deep") {
        #parametric formula to retrieve spectral slope of CDOM + NAP
        S_CDM = 0.015 + (0.002/(0.6 + (Rrs_obs.interp[which.min(abs(lambda - 443))]/Rrs_obs.interp[which.min(abs(lambda - 555))])))
        
        if (verbose == T) {
          print(paste0("The spectral slope for CDOM + NAP is calculated as:", S_CDM))
        }
        
        
      } else {
        
        if (verbose == T) {
          print("Shallow water type, QAA will derive wrong s_CDM, replaced with constant 0.017")
        }
        
        S_CDM <- 0.017 #<< HARD CODED>>
        
      }
      
      
    } else {
      if (use_manual_slope  == TRUE) {
        S_CDM = as.numeric(manual_slope["s_g"] + manual_slope["s_d"])
        
        if (verbose == T) {
          print(paste0("The spectral slope for CDOM + NAP is supplied by user as:", S_CDM))
        }
        
      } else {
        S_CDM <- 0.017 #<< Model Default >>
        
        if (verbose == T) {
          print(paste0("The spectral slope for CDOM + NAP is kept constant as:", S_CDM))
        }
        
      }
      
    }
    #Vectorization
    abs_CDM <- sapply(1:length(lambda), function(i) abs_CDM_440 * exp(-S_CDM * (lambda[i] - 440)))
    
    
    if (verbose == T) {
      print(paste0("CDOM+NAP absorption calculated using absolute values from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
    
    #--------------------------------------------------------------------------
    ## Scattering and back-scattering
    #--------------------------------------------------------------------------
    
    # Backscattering coefficient for suspended particles [1/m]
    if (slope.parametric == TRUE) {
      
      if (type_Rrs_below == "deep") {
        #parametric formula to retrieve spectral slope of bbp using QAAv5
        refexponent = 2*(1-(1.2*exp(-0.9 * (Rrs_obs.interp[which.min(abs(lambda - 443))]/Rrs_obs.interp[which.min(abs(lambda - 555))]))))
        
        if (verbose == T) {
          print(paste0("The spectral slope for bbp is calculated as:", refexponent))
        }
        
        
      } else {
        if (verbose == T) {
          print("Shallow water type, QAA will derive wrong spectral slope, manual slope will be used")
        }
        
        if (use_manual_slope  == TRUE) {
          refexponent = as.numeric(manual_slope["gamma"])
          
          if (verbose == T) {
            print(paste0("The spectral slope for bbp is supplied by user as:", refexponent))
          }
          
        } else {
          refexponent <- 0.46 #<< Model Default >>
          
          if (verbose == T) {
            print(paste0("The spectral slope for bbp is kept constant as:", refexponent))
          }
          
        }
        
      }
      
    } else {
      
      if (use_manual_slope  == TRUE) {
        refexponent = as.numeric(manual_slope["gamma"])
        
        if (verbose == T) {
          print(paste0("The spectral slope for bbp is supplied by user as:", refexponent))
        }
        
      } else {
        refexponent <- 0.46 #<< Model Default >>
        
        if (verbose == T) {
          print(paste0("The spectral slope for bbp is kept constant as:", refexponent))
        }
        
      }
      
    }
    #Vectorization
    bb_x <- sapply(1:length(lambda), function(i) base.bbp * ((lambda[i] / 550) ^ -refexponent))
    
    if (verbose == TRUE) {
      print(paste0("NAP backscatter calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
    #--------------------------------------------------------------------------
    ## Total Water system IOPs
    #--------------------------------------------------------------------------
    
    ## Total Absorption Coefficient (1/m)
    browser()
    abs <-  a_W + abs_ph + abs_CDM #+ abs_X
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
    
  } else {
    #--------------------------------------------------------------------------
    ## Water system total IOPs
    #--------------------------------------------------------------------------
    if (verbose == T) {
      print(paste0("Actual IOPs provided"))
    }
    
    a_non_water[a_non_water < 0] = 0
    
    ## Total Absorption Coefficient (1/m)
    if (!(all(length(a_non_water_wave) == length(lambda)) && all(a_non_water_wave == lambda))) {
      if (verbose == T) {
        print("Simulation wavelength and absorption wavelength are different, data will be interpolated")
      }
      
      
      a_non_water_interp = approx(x= a_non_water_wave, y = a_non_water,
                                  xout = lambda, method = "linear")$y
      
      a_non_water_interp[is.na(a_non_water_interp)] = 0
      abs <-  a_W + a_non_water_interp
      
    } else {
      
      abs <-  a_W + a_non_water
    }
    
    #Plot the fitted a_t-w
    if (plot == TRUE) {
      plot(a_non_water_df$wave, a_non_water_df$at_w, xlab="wavelength", 
           ylab="non-water absorption [m^-1]")
      lines(lambda, a_non_water_interp, col="red", lwd=3)
    }
    
    if (verbose == TRUE) {
      print(paste0("Total absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
    ## Total Backscattering Coefficient (1/m)
    
    if ( !(all(length(bb_non_water_wave) == length(lambda)) && all(bb_non_water_wave == lambda))) {
      if (verbose == T) {
        print("Simulation wavelength and backscatter wavelength are different")
      }
      
      #### Compute bbp and bb spectral slope
      if (length(bb_non_water_wave) == 6) {
        if (verbose == T) {
          print("HS-6 VSF is used")
        }
        
        HS6_wl = c(394, 420, 470, 532, 620, 700)
        x = 555/HS6_wl
        #nz=length(IOP.fitted.down$Depth)
        nz=1 #As we only are interested for surface bb
        bbP555.down =rep(0,nz)
        nuP.down    =rep(0,nz)
      }
      
      if (length(bb_non_water_wave) == 9) {
        if (verbose == T) {
          print("BB9 VSF is used")
        }
       
        HS6_wl = c(412, 440, 488, 510, 532, 595, 650, 676, 715)
        x = 555/HS6_wl
        #nz=length(IOP.fitted.down$Depth)
        nz=1 #As we only are interested for surface bb
        bbP555.down =rep(0,nz)
        nuP.down    =rep(0,nz)
      }
      
      for (i in 1:nz) {
        #y = IOP.fitted.down$HS6$bbP[i,]
        y = as.matrix(bb_non_water)[,i]
        y[y < 0] = NA
        y[y > 0.1] = NA
        if (any(is.na(y))) {
          bbP555.down[i]=NA
          nuP.down[i]=NA
        } else {
          model = nls(y~b*x^z, start = list(b = y[3], z = 1),data=data.frame(x,y), 
                      control=list(maxiter=100, warnOnly=T))
          bbP555.down[i]=coef(model)[1]
          nuP.down[i]=coef(model)[2]
        }
      }
      #### Apply the power-law model to get hyperspectral bb
      waveletngth_hs <- matrix(data=lambda)
      refbbp <- matrix(data=bbP555.down)
      refexponent <- matrix(nuP.down)
      
      bbp_hs <- matrix(nrow = length(refbbp), ncol =length(waveletngth_hs),data = 0 )
      i=1
      j=1
      
      while (j<=length(refbbp)) {
        #print(refbbp[j,])
        for (i in 1:length(waveletngth_hs)) {
          #print(waveletngth_hs[i,])
          
          #Implement power model
          bbp_hs[j,i] = refbbp[j,1]*((waveletngth_hs[i,1]/555)^-(refexponent[j,1])) 
          
        }
        j=j+1
      }
      if (verbose == T) {
        print("bbp for simulation wavelength is obtained using power-law fitting")
      }
      
      bb <-  bb_W + as.vector(bbp_hs)
      
    } else {
      
      bb <-  bb_W + bb_non_water
    }
    
    #Plot the fitted bbp
    if (plot == TRUE) {
      plot(lambda, as.vector(bbp_hs),xlab="wavelength", 
           ylab="non-water backscatter [m^-1]")
      lines(bb_non_water_df$wave, bb_non_water_df$bb, col="red", lwd=3)
    }
    
    
    if (verbose == TRUE) {
      print(paste0("Total backscatter calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
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
      
      
      if (verbose == T) {
        print("Bottom reflectance is assumed as LMM of user defined pure-spectra")
      }
      
      # Reflection factors of bottom surface [1/sr]
      #B0 <- 1/pi; 
      B1 <- 1/pi; B2 <- 1/pi; B3 <- 1/pi #;B4 <- 1/pi; B5 <- 1/pi; 
      BOTTOM <- c(B1,B2,B3)#,B4,B5)
      
      # Bottom Albedo Calculation
      
      abott1 <-  rb$class1 
      abott2 <-  rb$class2
      abott3 <-  rb$class3 
      
      abott <- rbind(abott1, abott2, abott3)#, abott4, abott5)
      
      Bottom <-  matrix(nrow = length(fA), ncol = ncol(abott), 0)
      Rrs_Bottom <- matrix(nrow = length(fA), ncol = ncol(abott), 0)# Bottom Rrs [1/sr]
      
      Bottom = fA * abott #BitWISE operation
      
      Rrs_Bottom = BOTTOM * Bottom #BitWISE operation

      Bottom <- colSums(Bottom) #[unitless]
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
      #browser()
      
      if (verbose == TRUE) {
        print(paste0("Subsurface (0^-) Rrs in Shallow water calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
      }
    }
  }
  
  #--------------------------------------------------------------------------
  ## Final Remote sensing reflectance computation
  #--------------------------------------------------------------------------
  Rrs <-  Rrs_below
  if (verbose == TRUE) {
    print(paste0("FINAL: Subsurface (0^-)Rrs calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    cat(paste0("\033[0;32m","=====================================================================================================","\033[0m","\n"))
    cat(paste0("\033[0;32m","######################################### SIMULATION ENDS #########################################","\033[0m","\n"))
    cat(paste0("\033[0;32m","=====================================================================================================","\033[0m","\n"))
    
  }
  #browser()
  return(list(data.frame("wavelength"=lambda, "Rrs"=Rrs)))
  #return(Rrs)
  
}

#=====================================================================================================
# Saber_forward_fast_sensitivity_test.R is EXACTLY same as Saber_forward_fast.R but the return data 
# structure is different (returns as a vector instead of a list of data-frame)

#Author: Mr. Soham Mukherjee - PhD Student, Aquatel, UQAR
#=====================================================================================================

#SABER FORWARD Model FAST version
Saber_forward_fast_sensitivity_test <-  function(use_true_IOPs = T, #Set TRUE if actual spectral IOPs exist
                                a_non_water_path = "./data/rb_retrieve_demo_a.csv", #a path
                                bb_non_water_path = "./data/rb_retrieve_demo_bb.csv", #bb path
                                
                                
                                chl=4.96, #must be input if SICF=TRUE and use_true_IOPs = F
                                
                                
                                a_dg = 1,  # provide value of acdom440/443 + anap440/443 (adg443)
                                
                                
                                bbp.550=0.00726002, #bbp value at 550/555 nm
                                
                                slope.parametric = TRUE, #Note, only possible if dg_composite = TRUE 
                                
                                use_manual_slope = FALSE, #use manual spectral slopes
                                manual_slope = c("s_g"=0.015, "s_d"=0.01160, "gamma"=1), #Values of
                                #manual spectral slopes. must be provided in a named vector as shown
                                
                                sun = 32.03672, #sun zenith angle
                                view = 30, #viewing zenith angle
                                
                                
                                z=2, #bottom depth
                                rb.fraction = fA.set, #aerial fraction of bottom types
                                
                                verbose = FALSE, #plot diagnostics and console output
                                
                                wavelength = seq(400,800,10),
                                Rrs_input_for_slope
                                
){
  
  if (slope.parametric == TRUE & use_manual_slope == TRUE) {
    stop("both automatic and manual slope can't be set TRUE")
  }
  
  ## OAC and IOP initialization
  
  # Read input params for forward modeling
  
  base.CDM <- a_dg # CDM absorption at 440nm [m^-1]
  C_ph <- chl           # Phytoplankton concentration [mg/m^-3]
  base.bbp <- bbp.550   # particulate backscatter at 550nm [m^-1]
  lambda <- wavelength  #Desired spectral range for simulation
  Rrs_obs.interp <- Rrs_input_for_slope
  
  if (verbose == T) {
    
    cat(paste0("\033[0;36m","=====================================================================================================","\033[0m","\n"))
    cat(paste0("\033[0;36m","######################################### SIMULATION BEGINS #########################################","\033[0m","\n"))
    cat(paste0("\033[0;36m","=====================================================================================================","\033[0m","\n"))
    
  }
  
  ##If in situ IOPs are provided
  if (use_true_IOPs == TRUE) {
    
    if (!is.null(a_non_water_path) & !is.null(bb_non_water_path)) {
      
      a_non_water_df = read.csv(file = a_non_water_path, header = T)
      bb_non_water_df = read.csv(file = bb_non_water_path, header = T)
      
      if (is.data.frame(a_non_water_df) & is.data.frame(bb_non_water_df)) {
        
        a_non_water = a_non_water_df$at_w
        a_non_water_wave = a_non_water_df$wave
        
        bb_non_water = bb_non_water_df$bbp
        bb_non_water_wave = bb_non_water_df$wave
        
      } else {
        stop("The IOPs must be provided as data frames with <<wave>> as wavelength & <<a/bb>> as value")
      }
      
    }
  }
  
  #Shallow parameters
  fA <- rb.fraction     #Aerial fraction of bottom albedo
  zB <- z               #bottom depth
  
  if (verbose == TRUE) {
    print("Function arguement initialized") 
  }
  ################################################################################################
  ##                                    Bio-Optical Development                                 ##
  ################################################################################################
  
  #--------------------------------------------------------------------------
  ## Pure water IOPs
  #--------------------------------------------------------------------------
  ## absorption (1/m)
  
  # wavelength range [190;4000] [nm]
  abs.water <- read.table("./data/input-spectra/abs_W.A", header = F) 
  wavelength <- abs.water$V1
  absorpt_W <-  abs.water$V2
  
  a_W <- rep(0,length(lambda))# abs. of pure water [1/m]
  a_W <-  Hmisc::approxExtrap(wavelength, absorpt_W,xout = lambda,method = "linear")$y
  
  if (verbose == TRUE) {
    print(paste0("Water absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm.")) 
  }
  
  ## backscattering (1/m)
  
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
  
  #--------------------------------------------------------------------------
  ##Optically Active Constituents IOPs
  #--------------------------------------------------------------------------
  
  if (use_true_IOPs == FALSE) {
    
    if (verbose == TRUE) {
      print(paste0("Full spectral IOPs not provided, bio-optical models are used"))
    }
    
    ## Plankton absorption (1/m)
    
    # load plankton absorption data
    A0_A1_PhytoPlanc <- read.table("./data/input-spectra/A0_A1_PhytoPlanc.dat")
    
    #===============================================================================
    #Replace with the below for package environment for upper line
    #system.file("data", "input-spectra", "A0_A1_PhytoPlanc.dat", package = "SABER")
    #===============================================================================
    
    # extract the values from the table
    lam_p <- A0_A1_PhytoPlanc$V1
    a0_p <- A0_A1_PhytoPlanc$V2
    a1_p <- A0_A1_PhytoPlanc$V3
    
    a0 <- rep(0,length(lambda))# [m^2/mg]
    a1 <- rep(0,length(lambda))# [m^2/mg]
    
    a0 <- Hmisc::approxExtrap(lam_p, a0_p,xout = lambda,method = "linear")$y
    a1 <- Hmisc::approxExtrap(lam_p, a1_p,xout = lambda,method = "linear")$y
    
    aph_440 <- 0.06 * (C_ph)^0.65  # [mg/m^3]
    abs_ph <- sapply(1:length(lambda), function(i) (a0[i] + a1[i] * log(aph_440)) * aph_440) #Vectorization
    
    
    abs_ph[abs_ph < 0] <- 0
    
    
    if (verbose == TRUE) {
      print(paste0("Planktonic absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
    ## CDOM+NAP absorption coefficient [1/m]
    
    
    Ga_CDOM = 1
    Oa_CDOM = 0
    
    abs_CDM_440 <-  (Ga_CDOM*base.CDM)+Oa_CDOM# [1/m], CDOM abs. coeff. at 440 [nm]
    
    if (slope.parametric == TRUE) {
      
      if (type_Rrs_below == "deep") {
        #parametric formula to retrieve spectral slope of CDOM + NAP
        S_CDM = 0.015 + (0.002/(0.6 + (Rrs_obs.interp[which.min(abs(lambda - 443))]/Rrs_obs.interp[which.min(abs(lambda - 555))])))
        
        if (verbose == T) {
          print(paste0("The spectral slope for CDOM + NAP is calculated as:", S_CDM))
        }
        
        
      } else {
        
        if (verbose == T) {
          print("Shallow water type, QAA will derive wrong s_CDM, replaced with constant 0.017")
        }
        
        S_CDM <- 0.017 #<< HARD CODED>>
        
      }
      
      
    } else {
      if (use_manual_slope  == TRUE) {
        S_CDM = as.numeric(manual_slope["s_g"] + manual_slope["s_d"])
        
        if (verbose == T) {
          print(paste0("The spectral slope for CDOM + NAP is supplied by user as:", S_CDM))
        }
        
      } else {
        S_CDM <- 0.017 #<< Model Default >>
        
        if (verbose == T) {
          print(paste0("The spectral slope for CDOM + NAP is kept constant as:", S_CDM))
        }
        
      }
      
    }
    #Vectorization
    abs_CDM <- sapply(1:length(lambda), function(i) abs_CDM_440 * exp(-S_CDM * (lambda[i] - 440)))
    
    
    if (verbose == T) {
      print(paste0("CDOM+NAP absorption calculated using absolute values from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
    
    #--------------------------------------------------------------------------
    ## Scattering and back-scattering
    #--------------------------------------------------------------------------
    
    # Backscattering coefficient for suspended particles [1/m]
    if (slope.parametric == TRUE) {
      
      if (type_Rrs_below == "deep") {
        #parametric formula to retrieve spectral slope of bbp using QAAv5
        refexponent = 2*(1-(1.2*exp(-0.9 * (Rrs_obs.interp[which.min(abs(lambda - 443))]/Rrs_obs.interp[which.min(abs(lambda - 555))]))))
        
        if (verbose == T) {
          print(paste0("The spectral slope for bbp is calculated as:", refexponent))
        }
        
        
      } else {
        if (verbose == T) {
          print("Shallow water type, QAA will derive wrong spectral slope, manual slope will be used")
        }
        
        if (use_manual_slope  == TRUE) {
          refexponent = as.numeric(manual_slope["gamma"])
          
          if (verbose == T) {
            print(paste0("The spectral slope for bbp is supplied by user as:", refexponent))
          }
          
        } else {
          refexponent <- 0.46 #<< Model Default >>
          
          if (verbose == T) {
            print(paste0("The spectral slope for bbp is kept constant as:", refexponent))
          }
          
        }
        
      }
      
    } else {
      
      if (use_manual_slope  == TRUE) {
        refexponent = as.numeric(manual_slope["gamma"])
        
        if (verbose == T) {
          print(paste0("The spectral slope for bbp is supplied by user as:", refexponent))
        }
        
      } else {
        refexponent <- 0.46 #<< Model Default >>
        
        if (verbose == T) {
          print(paste0("The spectral slope for bbp is kept constant as:", refexponent))
        }
        
      }
      
    }
    #Vectorization
    bb_x <- sapply(1:length(lambda), function(i) base.bbp * ((lambda[i] / 550) ^ -refexponent))
    
    if (verbose == TRUE) {
      print(paste0("NAP backscatter calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
    #--------------------------------------------------------------------------
    ## Total Water system IOPs
    #--------------------------------------------------------------------------
    
    ## Total Absorption Coefficient (1/m)
    browser()
    abs <-  a_W + array(unlist(abs_ph)) + array(unlist(abs_CDM)) #+ abs_X
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
    
  } else {
    #--------------------------------------------------------------------------
    ## Water system total IOPs
    #--------------------------------------------------------------------------
    if (verbose == T) {
      print(paste0("Actual IOPs provided"))
    }
    
    a_non_water[a_non_water < 0] = 0
    
    ## Total Absorption Coefficient (1/m)
    if (!(all(length(a_non_water_wave) == length(lambda)) && all(a_non_water_wave == lambda))) {
      if (verbose == T) {
        print("Simulation wavelength and absorption wavelength are different, data will be interpolated")
      }
      
      
      a_non_water_interp = approx(x= a_non_water_wave, y = a_non_water,
                                  xout = lambda, method = "linear")$y
      
      a_non_water_interp[is.na(a_non_water_interp)] = 0
      abs <-  a_W + a_non_water_interp
      
    } else {
      
      abs <-  a_W + a_non_water
    }
    
    #Plot the fitted a_t-w
    if (plot == TRUE) {
      plot(a_non_water_df$wave, a_non_water_df$a, xlab="wavelength", 
           ylab="non-water absorption [m^-1]")
      lines(lambda, a_non_water_interp, col="red", lwd=3)
    }
    
    if (verbose == TRUE) {
      print(paste0("Total absorption calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
    
    ## Total Backscattering Coefficient (1/m)
    
    if ( !(all(length(bb_non_water_wave) == length(lambda)) && all(bb_non_water_wave == lambda))) {
      if (verbose == T) {
        print("Simulation wavelength and backscatter wavelength are different")
      }
      
      #### Compute bbp and bb spectral slope
      if (length(bb_non_water_wave) == 6) {
        if (verbose == T) {
          print("HS-6 VSF is used")
        }
        
        HS6_wl = c(394, 420, 470, 532, 620, 700)
        x = 555/HS6_wl
        #nz=length(IOP.fitted.down$Depth)
        nz=1 #As we only are interested for surface bb
        bbP555.down =rep(0,nz)
        nuP.down    =rep(0,nz)
      }
      
      if (length(bb_non_water_wave) == 9) {
        if (verbose == T) {
          print("BB9 VSF is used")
        }
        
        HS6_wl = c(412, 440, 488, 510, 532, 595, 650, 676, 715)
        x = 555/HS6_wl
        #nz=length(IOP.fitted.down$Depth)
        nz=1 #As we only are interested for surface bb
        bbP555.down =rep(0,nz)
        nuP.down    =rep(0,nz)
      }
      
      for (i in 1:nz) {
        #y = IOP.fitted.down$HS6$bbP[i,]
        y = as.matrix(bb_non_water)[,i]
        y[y < 0] = NA
        y[y > 0.1] = NA
        if (any(is.na(y))) {
          bbP555.down[i]=NA
          nuP.down[i]=NA
        } else {
          model = nls(y~b*x^z, start = list(b = y[3], z = 1),data=data.frame(x,y), 
                      control=list(maxiter=100, warnOnly=T))
          bbP555.down[i]=coef(model)[1]
          nuP.down[i]=coef(model)[2]
        }
      }
      #### Apply the power-law model to get hyperspectral bb
      waveletngth_hs <- matrix(data=lambda)
      refbbp <- matrix(data=bbP555.down)
      refexponent <- matrix(nuP.down)
      
      bbp_hs <- matrix(nrow = length(refbbp), ncol =length(waveletngth_hs),data = 0 )
      i=1
      j=1
      
      while (j<=length(refbbp)) {
        #print(refbbp[j,])
        for (i in 1:length(waveletngth_hs)) {
          #print(waveletngth_hs[i,])
          
          #Implement power model
          bbp_hs[j,i] = refbbp[j,1]*((waveletngth_hs[i,1]/555)^-(refexponent[j,1])) 
          
        }
        j=j+1
      }
      if (verbose == T) {
        print("bbp for simulation wavelength is obtained using power-law fitting")
      }
      
      bb <-  bb_W + as.vector(bbp_hs)
      
    } else {
      
      bb <-  bb_W + bb_non_water
    }
    
    #Plot the fitted bbp
    if (plot == TRUE) {
      plot(lambda, as.vector(bbp_hs),xlab="wavelength", 
           ylab="non-water backscatter [m^-1]")
      lines(bb_non_water_df$wave, bb_non_water_df$bb, col="red", lwd=3)
    }
    
    
    if (verbose == TRUE) {
      print(paste0("Total backscatter calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    }
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
      
      
      if (verbose == T) {
        print("Bottom reflectance is assumed as LMM of user defined pure-spectra")
      }
      
      # Reflection factors of bottom surface [1/sr]
      #B0 <- 1/pi; 
      B1 <- 1/pi; B2 <- 1/pi; B3 <- 1/pi #;B4 <- 1/pi; B5 <- 1/pi; 
      BOTTOM <- c(B1,B2,B3)#,B4,B5)
      
      # Bottom Albedo for ALGAE-WISE
      
      abott1 <-  rb$class1 
      abott2 <-  rb$class2
      abott3 <-  rb$class3 
      
      abott <- rbind(abott1, abott2, abott3)#, abott4, abott5)
      
      Bottom <-  matrix(nrow = length(fA), ncol = ncol(abott), 0)
      Rrs_Bottom <- matrix(nrow = length(fA), ncol = ncol(abott), 0)# Bottom Rrs [1/sr]
      
      Bottom = fA * abott #BitWISE operation
      
      Rrs_Bottom = BOTTOM * Bottom #BitWISE operation
      
      Bottom <- colSums(Bottom) #[unitless]
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
      #browser()
      
      if (verbose == TRUE) {
        print(paste0("Subsurface (0^-) Rrs in Shallow water calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
      }
    }
  }
  
  #--------------------------------------------------------------------------
  ## Final Remote sensing reflectance computation
  #--------------------------------------------------------------------------
  Rrs <-  Rrs_below
  if (verbose == TRUE) {
    print(paste0("FINAL: Subsurface (0^-)Rrs calculated from ", min(lambda), " nm to ", max(lambda), " nm."))
    cat(paste0("\033[0;32m","=====================================================================================================","\033[0m","\n"))
    cat(paste0("\033[0;32m","######################################### SIMULATION ENDS #########################################","\033[0m","\n"))
    cat(paste0("\033[0;32m","=====================================================================================================","\033[0m","\n"))
    
  }
  #browser()
  #return(list(data.frame("wavelength"=lambda, "Rrs"=Rrs)))
  return(Rrs)
  
}


