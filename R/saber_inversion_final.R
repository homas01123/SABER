
source("./R/saber_forward_parametric_conc_wise.R") #FINAL VERSION
source("./R/SABER_forward_fast.R") #FASTEST VERSION

#Inelastic Scattering Models
source("./R/SABER.sicf.R")
source("./R/SABER.fdom.R")
source("./R/inelastic.wavelength.redistribution.R")

#QAAv5 for IOP and Spectral slopes
source("./R/QAA_v5.R")

#Inverse Am03 model with residual function
#source("./R/saber_inverse.R")

#Sun Earth Geometry and Downwelling Irradiance 
source("./R/snell_law.R")
source("./R/gregg_carder_ed.R")

#Spectral Bottom Reflectance Retrieval
source("./R/retrive.rb.spectral.wise.R")
source("./R/generate.rb.spectral.R")
source("./R/read.surface.IOPs.wise.R")
source("./R/retrieve.rb.spectral.minimal.R")

#Gradient Based Inverse Functions
source("./R/solve.objective.inverse.final.R") 
source("./R/solve.objective.inverse_fast.R") #FINAL VERSION

#MCMC Based Inverse Functions
source("./R/mcmc.functions.R")
source("./R/saber.inversion.vector.parallel.R")
source("./R/mcmc_bayes_parallel.R")

library(dplyr)
library(readxl)
library(stats4)
library(MASS)
library(dglm)
library(fitdistrplus)
library(Riops)
library(Cops)
library(ggplot2)
#library(rho)
library(marqLevAlg)
library(BayesianTools)
library(coda)
library(stats)
library(Deriv)
library(mcmc)
library(plot3D)
library(limSolve)
library(nloptr)
library(minpack.lm)
library(pracma)
library(alabama)
library(data.table)
library(ggalt)
library(ggExtra)
library(suncalc)
library(viridis)
library(doParallel)
library(sensitivity)
library(lubridate)
library(tidyr)
library(tidyverse)

# Function to trnaslate above water Sun Zenith angle to sub-surface Zenith angle ----                         
sunzen_below = function(sun_zen_aove=45){
  sun_zen_below_rad = asin(sin(sun_zen_aove*(pi/180))/1.333)
  sun_zen_below_deg = sun_zen_below_rad*(180/pi)
  return(sun_zen_below_deg)
}

# Function to translate Rrs0+ to Rrs0- ----
surface_rrs_translate <-  function(Rrs) { 
  rrs = Rrs / (0.52 + 1.7*Rrs)
  return(rrs)
}


# Global Declarations ----


## Water type and wavelength specifications ----
type_case_water = 2
type_Rrs_below = "shallow"
type_Rrs_below_rb = "shallow"
type_Rrs_water = "below_surface"
wavelength <- seq(400,750,10)


## Plotting specifications ----
plot=FALSE #Set TRUE when the ggplot2 outputs needed to be saved onto disk


## Inversion specifications ----
preFit = FALSE

use.lklhood= TRUE
use.wise.prior = FALSE 
use.nomad.prior = TRUE
use.ioccg.prior = FALSE

## Sun-sensor geometry specifications ----
view = 30
sun_above  = 45
sun = sunzen_below(sun_zen_aove = sun_above) #Make sure you pass "sun" for subsurface calculation

## Shallow parameters specifications ----

### Water depth ----
zB=2

### Areal fraction of bottom reflectance ----

# fA0=0; # constant (Excluded from Model) 
# fA1=0.2; # sand
# fA2=0.2; # sediment
# fA3=0; # Chara contraria
# fA4=0.6; # Potamogeton perfoliatus
# fA5=0; # Potamogeton pectinatus
#fA.set= c(fA1,fA2,fA3,fA4,fA5)

use_WASI_rb = FALSE
use_WISE_Man_rb = TRUE
use_algae_WISE_rb = FALSE

# if (use_WISE_Man_rb == TRUE) {
#   load("./data/WISE-Man.RData")
#   
# } else {
#   load("./data/algae-WISE.RData")
# }
# 
# if (!(all(length(rb$wavelength) == length(wavelength)) && all(rb$wavelength == wavelength))) {
#   rb_interp = data.frame("wavelength" = wavelength,
#                          "class1" = approx(rb$wavelength, rb$class1, xout = wavelength)$y,
#                          "class2" = approx(rb$wavelength, rb$class2, xout = wavelength)$y,
#                          "class3" = approx(rb$wavelength, rb$class3, xout = wavelength)$y
#   )
#   
# }
# rb = rb_interp
# rm(rb_interp)

fA1=0.5; # aerial fraction 1 
fA2=0.25; # aerial fraction 2 
fA3=0.25; # aerial fraction 3, i.e. fa1+fa2+fa3 = 1

fA.set= c(fA1,fA2,fA3)

source("./R/retrive.rb.spectral.R")
source("./R/generate.rb.spectral.R")

## Atmospheric specifications ----

### Irradiance intensities [1/sr] ----
g_dd=0.05; g_dsr=0; g_dsa=0;

### Intensities of light sources ----
f_dd= 1; f_ds= 1;

### Angstrom exponent ----
alpha = 1.317;

### Atmospheric pressure ----
P = 1013.25; # [mbar]

### Relative Humidity ----
RH = 0.60;

### Scale height for ozone ----
Hoz = 0.300; # [cm]

### Scale height of the precipitate water in the atmosphere ----
WV= 2.500; # [cm]


# Function to read the configuration file for inversion ----

read_inv_param <- function(inv_config_file){
  
  inv_config = read.table(file = inv_config_file, sep=",", header = T,colClasses = "character")
  rrs_type = as.character(trimws(inv_config$param_vals[1]))
  
  mode = as.character(trimws(inv_config$param_vals[2]))
  
  param_name = as.character(trimws((unlist((strsplit(inv_config$param_vals[3],";"))))))
  rb_type = as.character(trimws((unlist((strsplit(inv_config$param_vals[4],";"))))))
  min_par = as.numeric((unlist((strsplit(inv_config$param_vals[5],";")))))
  max_par = as.numeric((unlist((strsplit(inv_config$param_vals[6],";")))))
  init_par = as.numeric((unlist((strsplit(inv_config$param_vals[7],";")))))
  
  constrain_config = as.logical(trimws((unlist((strsplit(inv_config$param_vals[8],";"))))))
  names(constrain_config) = c("const_bbp", "const_bgc", "const_IOP")
  constrain_bbp_val = as.numeric(inv_config$param_vals[9])
  constrain_bgc_val = as.numeric(trimws((unlist((strsplit(inv_config$param_vals[10],";"))))))
  names(constrain_bgc_val) = c("chl", "adg443", "bbp555")
  constrain_iop = as.character(trimws((unlist((strsplit(inv_config$param_vals[11],";"))))))
  
  qaa_mode = as.logical(trimws(inv_config$param_vals[12]))
  qaa_prefit = as.logical(trimws(inv_config$param_vals[13]))
  
  auto_slope = as.logical(trimws(inv_config$param_vals[14]))
  manual_slope  =  as.logical(trimws(inv_config$param_vals[15]))
  manual_slope_vals = as.numeric((unlist((strsplit(inv_config$param_vals[16],";")))))
  names(manual_slope_vals) = c("s_g", "s_d", "gamma")
  
  sa_model = as.character(trimws((unlist((strsplit(inv_config$param_vals[17],";"))))))
  obj_fn = as.character(trimws((unlist((strsplit(inv_config$param_vals[18],";"))))))
  opt_method = as.character(trimws((unlist((strsplit(inv_config$param_vals[19],";"))))))
  
  iter_count = as.numeric(trimws(inv_config$param_vals[20]))
  sampler_mcmc = as.character(trimws(inv_config$param_vals[21]))
  hybrid_mode = as.logical(trimws(inv_config$param_vals[22]))
  
  config_list = list(rrs_type, mode, param_name, rb_type, min_par, max_par, init_par, constrain_config,
                     constrain_bbp_val, constrain_bgc_val, constrain_iop, qaa_mode,
                     qaa_prefit, auto_slope, manual_slope, manual_slope_vals, 
                     sa_model,obj_fn, opt_method, iter_count, sampler_mcmc, hybrid_mode )
  names(config_list) = inv_config$param_config
  return(config_list)

}

#test_param = read_inv_param(inv_config_file = "./data/inversion_config.txt")


# Function to create initial values and bounds for inversion ----

rb_count = length(fA.set)
create_init_bound <- function(pop.sd = "unknown", 
                              constrain.bbp= F,
                              constrain.shallow.bgc = F , constrain.shallow.iop = F,
                              manual_par0 = T, 
                              
                              init_par = c(chl = 2, adg440 = 0.8,
                                           bbp = 0.005,z = 2.5, 
                                           rb_c = rep(0.5, rb_count),
                                           population.sd = 0.1),
                              
                              upper_par = c(chl = 30, adg440 = 5,
                                            bbp = 0.01,z = 10, 
                                            rb_c = rep(1, rb_count),
                                            population.sd = 10),
                              
                              lower_par = c(chl = 1, adg440 = 0.2,
                                            bbp = 0.001,z = 0.5, 
                                            rb_c = rep(0, rb_count),
                                            population.sd = 1e-05),
                              
                              rrs_inv_type = "deep"){
  
  
  # Initial values for deep water 
  
  if (pop.sd == "unknown" & rrs_inv_type == "deep" & manual_par0 == TRUE & constrain.bbp == "TRUE") {# <<MANUAL INPUT>>
    par0 = c(chl = init_par[1],
             adg = init_par[2], 
             #bbp550 =0.005,
             population.sd = init_par[length(init_par)]) #@@  sigma UNKNOWN
  }
  
  if (pop.sd == "unknown" & rrs_inv_type == "deep" & manual_par0 == TRUE & constrain.bbp == "FALSE") {# <<MANUAL INPUT>>
    par0 = c(chl = init_par[1],
             adg = init_par[2], 
             bbp550 =init_par[3],
             population.sd = init_par[length(init_par)]) #@@  sigma UNKNOWN
  }
  
  #Set bounds (Instead of flat multiplier, implement range as a function of parameter sensitivity)
  if (rrs_inv_type == "deep" & constrain.bbp == TRUE & manual_par0 == FALSE) {
    upper.bound <- c(5*init_par[1:2], upper_par[length(upper_par)])  
    lower.bound <- c(0.2*init_par[1:2], lower_par[length(lower_par)])
    
  } else {
    
    if (rrs_inv_type == "deep" & manual_par0 == T) {
      if (constrain.bbp == FALSE) {
        upper.bound <- upper_par[c(1:3,length(upper_par))]  
        lower.bound <- lower_par[c(1:3,length(lower_par))]
      } else {
        upper.bound <- upper_par[c(1:2,length(upper_par))]  
        lower.bound <- lower_par[c(1:2,length(lower_par))]
      }
      
    }
    
  }  
  
  
  # Initial values for shallow water
  
  rb_count = length(fA.set)
  rep(0.5, rb_count)
  rb_c = paste0("rb", (1:rb_count))
  
  if (pop.sd == "unknown" & rrs_inv_type == "shallow" & constrain.shallow.bgc == "TRUE") { # <<MANUAL INPUT>>
    
    par0 = init_par[-(1:3)]
    
    if (manual_par0 == T) {
      
      lower.bound = lower_par[-(1:3)]
      upper.bound = upper_par[-(1:3)]
      
    } else {
      
      lower.bound <- c(z = 1,
                       rep(0, rb_count),
                       population.sd = 0.00001)
      
      upper.bound <- c(z = 10,
                       rep(1, rb_count),
                       population.sd = 10)
      
    }
    
  }
  
  
  if (pop.sd == "unknown" & rrs_inv_type == "shallow" & constrain.shallow.iop == "TRUE") {
    par0 = c(z = init_par[4], 
             rb_c = init_par[5:7],
             population.sd = init_par[length(init_par)])
    
    if (manual_par0 == T) {
      
      lower.bound = lower_par[4:length(lower_par)]
      upper.bound = upper_par[4:length(upper_par)]
      
    } else {
      
      lower.bound <- c(z = 1,
                       rep(0, rb_count),
                       population.sd = 0.00001)
      
      upper.bound <- c(z = 10,
                       rep(1, rb_count),
                       population.sd = 10)
      
    }
    
  }
  
  
  if (pop.sd == "unknown" & rrs_inv_type == "shallow" & constrain.shallow.bgc == "FALSE" & 
      constrain.shallow.iop == "FALSE") { 
    
    par0 = init_par
    
    if (manual_par0 == T) {
      
      lower.bound = lower_par
      upper.bound = upper_par
      
    } else {
      
      #Autoscale Intital values from pre-Ft
      increament.scale <- 1
      
      lower.bound <- c((par0[1:3] - 0.8*par0[1:3]),z = 1,
                       rep(0, rb_count),
                       population.sd = 0.00001)
      
      upper.bound <- c((par0[1:3] + 5*par0[1:3]),z = 10,
                       rep(1, rb_count),
                       population.sd = 10)
      
    }
  }
  
  return(list("par0"=par0, "upper.bound"=upper.bound, "lower.bound"=lower.bound))
}

# inv_bound = create_init_bound(rrs_inv_type = rrs_type, manual_par0 = T, 
#                               constrain.bbp = F, 
#                               constrain.shallow.bgc = T, 
#                               constrain.shallow.iop = F, 
#                               pop.sd =  "unknown",
#                               
#                               init_par = rowMeans(matrix(c(max_par,min_par), 
#                                                          ncol=2)),
#                               upper_par = max_par,
#                               lower_par = min_par
# )
# 
# par0 = inv_bound$par0 ; upper_bound = inv_bound$upper.bound  
# lower_bound= inv_bound$lower.bound


# A single Function to perform the inversion as per the configuration ----

run_inverse_mode <- function(rrs_input, 
                             wavelength_input, 
                             inv_config_file = "./data/inversion_config.txt"){
  
  config_list = read_inv_param(inv_config_file)
  
  print("READING CONFIGURATION FILE...")
  
  if (config_list$` rrs_type` == "deep") {
    type_Rrs_below = "deep"
    assign(x = "type_Rrs_below", value = type_Rrs_below, envir = .GlobalEnv)
    
    if (config_list$` mode` == "use_grad") {
    inv_output = inverse_runGrad(rrs_type = config_list$` rrs_type`, obsdata = rrs_input, 
                                 max_par = config_list$` max_par`, 
                                 min_par = config_list$` min_par`, 
                                 init_par = config_list$` init_par`, 
                                 param_lab = config_list$` param_name`, 
                                 constrain_config = config_list$` constrain_conf`, 
                                 bgc_const_val = config_list$` constrain_bgc_val`, 
                                 bbp_const_val = config_list$` constrain_bbp_val`,
                                 iop_const_path =  config_list$` constrain_iop`,
                                 qaa_prefit = config_list$` qaa_prefit`, 
                                 QAA_mode = config_list$` qaa_mode`, 
                                 qaa_slope = config_list$` auto_slope`, 
                                 manual_slope = config_list$` manual_slope`, 
                                 manual_slope_vals = config_list$` manual_slope_val`, 
                                 wavelength_sim = wavelength_input, 
                                 sa_model = config_list$` sa_model`, 
                                 obj_fn = config_list$` obj_fn`, 
                                 opt_method = config_list$` opt_method`)
    }
    
    if (config_list$` mode` == "use_mcmc") {
      inv_output = inverse_runBayes(obsdata = rrs_input, 
                                    rrs_type =  type_Rrs_below, 
                                    max_par = config_list$` max_par`, 
                                    min_par = config_list$` min_par`, 
                                    param_lab = config_list$` param_name`, 
                                    constrain_config = config_list$` constrain_conf`, 
                                    bbp_const_val = config_list$` constrain_bbp_val`, 
                                    bgc_const_val = config_list$` constrain_bgc_val`, 
                                    iop_const_path = config_list$` constrain_iop`,
                                    qaa_slope = config_list$` auto_slope`, 
                                    manual_slope = config_list$` manual_slope`, 
                                    manual_slope_val = config_list$` manual_slope_val`, 
                                    iter_count = config_list$` iter_count`, 
                                    sampler_mcmc = config_list$` sampler_mcmc`, 
                                    wavelngth_sim = wavelength_input, 
                                    sa_model = config_list$` sa_model`, 
                                    hybrid_mode = config_list$` hybrid_mode`, plot_rrs = T)
    }
    
    
  }
  
  if (config_list$` rrs_type` == "shallow") {
    type_Rrs_below = "shallow"
    assign(x = "type_Rrs_below", value = type_Rrs_below, envir = .GlobalEnv)
    
    load("./data/EGSL_bottom.RData")
    
    column_numbers <- which(names(rb) %in% config_list$` rb_types`)
    print(paste0("The set benthic classes are : ", config_list$` rb_types`))
    
    rb = rb[c(1,column_numbers)]
    names(rb) = c("wavelength", "class1", "class2", "class3")
    
    if (!(all(length(rb$wavelength) == length(wavelength_input)) && all(rb$wavelength == wavelength_input))) {
      
      print(paste0("Rb wavelength has " , length(rb$wavelength), " bands; Rrs wavelength has ", length(wavelength_input) ))
      
      rb_interp = data.frame("wavelength" = wavelength_input,
                             "class1" = approx(rb$wavelength, rb$class1, xout = wavelength_input)$y,
                             "class2" = approx(rb$wavelength, rb$class2, xout = wavelength_input)$y,
                             "class3" = approx(rb$wavelength, rb$class3, xout = wavelength_input)$y
      )
      
    }
    rb = rb_interp
    rm(rb_interp)
    
    assign(x = "rb", value = rb, envir = .GlobalEnv)
    
    #browser()
    
    if (config_list$` mode` == "use_grad") {
        
        inv_output = inverse_runGrad(rrs_type = config_list$` rrs_type`, 
                                     obsdata = rrs_input, max_par = config_list$` max_par`, 
                                     min_par = config_list$` min_par`, 
                                     init_par = config_list$` init_par`, 
                                     param_lab = config_list$` param_name`, 
                                     constrain_config = config_list$` constrain_conf`, 
                                     bgc_const_val = config_list$` constrain_bgc_val`,
                                     bbp_const_val = config_list$` constrain_bbp_val`,
                                     iop_const_path =  config_list$` constrain_iop`,
                                     qaa_prefit = config_list$` qaa_prefit`, 
                                     QAA_mode = config_list$` qaa_mode`, 
                                     qaa_slope = config_list$` auto_slope`, 
                                     manual_slope = config_list$` manual_slope`, 
                                     manual_slope_vals = config_list$` manual_slope_val`, 
                                     wavelength_sim = wavelength_input, 
                                     sa_model = config_list$` sa_model`, 
                                     obj_fn = config_list$` obj_fn`, 
                                     opt_method = config_list$` opt_method`)
      
    }
    
    if (config_list$` mode` == "use_mcmc") {
      type_Rrs_below = config_list$` rrs_type`
      assign(x = "type_Rrs_below", value = type_Rrs_below, envir = .GlobalEnv)
      
      inv_output = inverse_runBayes(obsdata = rrs_input, 
                                    rrs_type =  type_Rrs_below, 
                                    max_par = config_list$` max_par`, 
                                    min_par = config_list$` min_par`, 
                                    param_lab = config_list$` param_name`, 
                                    constrain_config = config_list$` constrain_conf`,
                                    bbp_const_val = config_list$` constrain_bbp_val`, 
                                    bgc_const_val = config_list$` constrain_bgc_val`, 
                                    iop_const_path = config_list$` constrain_iop`,
                                    qaa_slope = config_list$` auto_slope`, 
                                    manual_slope = config_list$` manual_slope`, 
                                    manual_slope_val = config_list$` manual_slope_val`, 
                                    iter_count = config_list$` iter_count`, 
                                    sampler_mcmc = config_list$` sampler_mcmc`, 
                                    wavelngth_sim = wavelength_input, 
                                    sa_model = config_list$` sa_model`, 
                                    hybrid_mode = config_list$` hybrid_mode`, plot_rrs = T)
      
    }
    
  }
  
  
  return(inv_output)
}


# Test the inversion function ----

# run_inverse_mode(rrs_input = obsdata #+ rnorm.trunc(n = length(obsdata), mean = 0, 
#                                              #sd = sd(obsdata), min = 0)*0.05 
#                       , wavelength_input = wavelength)
# 
# run_inverse_mode(rrs_input = as.numeric(rrs.forward.SABER[1,]), 
#                  wavelength_input = wavelength)
# param_vec[100,]
