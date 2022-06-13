#===========================================================================
# solve.objective.inverse.R solves the inverse objective function  to retrieve the water components 
#along with bathymetry and bottom reflectance (for shallow water) given the initial vals and Rrs.

#There are two optimization methods:

#1. Minimize the residual sum-square-error using non-linear multivariate optimization using 
#non-gradient based method
#2. Maximize the conditional log-likelihood (P(Data|Theta) of the model noise and obtain hessian
#derivative of the function, iff the function is differentiable over paramteric space)

#Author: Mr. Soham Mukherjee - PhD Student, Aquatel, UQAR
#===========================================================================
solve.objective.inverse <- function(obj.fn,initial,obsdata,sa.model, method.opt, 
                                    lower.b, upper.b, batch, pop.sd){
  #initial <- abs(initial)
  if (obj.fn == "log-LL") {
    
    ##Create log-likelihood function
    
    # NLL = function(pars, data) {
    #   
    #   # Values predicted by the forward model for single RUN
    #   
    #   Gpred = Saber_forward(chl = pars[1], acdom440 = pars[2],
    #                         anap440 =pars[3], bbp.550 =HL.deep.iop$bbp550[j],
    #                         verbose = F, realdata = data, plot = F)
    #   
    #   # # # Values predicted by the forward model for BATCH RUN
    #   #     Gpred = Saber_forward(chl = pars[1], acdom440 = pars[2],
    #   #                           anap440 =pars[3], bbp.550 = HL.deep.iop$bbp550[j],
    #   #                           verbose = F, realdata = data, plot = F)
    #   
    #   
    #   # Negative log-likelihood 
    #   
    #   ll= -sum(dnorm(x = data, mean = Gpred[[1]]$Rrs, sd = pars[4], log = TRUE)) #ON for unknown sigma
    #   #grad = numDeriv::grad(x=pars, func=NLL,data= insitu.data)
    #   #attr(ll, "gradient") <- grad
    #   return(ll)
    #   
    # }
    
    
    NLL = function(pars, data) {

      # Values predicted by the forward model for single RUN
        if (sa.model == "am03") {
          Gpred = Saber_forward(chl = pars[1], acdom440 = pars[2],
                                anap440 =pars[3], bbp.550 = Fit.input$bbp.550,
                                verbose = F, realdata = data, plot = F)

        } else {
          Gpred = Lee_forward(chl = pars[1], acdom440 = pars[2],
                              anap440 =pars[3], bbp.550 = Fit.input$bbp.550,
                              verbose = F, realdata = data, plot = F)
        }

      # Negative log-likelihood
       if(pop.sd == "TRUE") {
         if (batch == FALSE) {
      
           smull = -sum(dnorm(x = data, mean = Gpred[[1]]$Rrs, sd = 0.0001, log = TRUE)) #ON for Albert Rrs
      
         } else {
      
           smull = -sum(dnorm(x = data, mean = Gpred[[1]]$Rrs, sd = 0.13129354/100, log = TRUE)) #ON for IOCCG
      
         }
       } else {
         if (batch == TRUE | batch == FALSE) {

          smull = -sum(dnorm(x = 10000*data, mean = 10000*Gpred[[1]]$Rrs, sd = pars[4], log = TRUE)) #ON for unknown sigma

           }
      
       }
      return(smull)
    }
    
    ##Optimize the log-likelihood function
    cat(paste0("\033[0;44m","#################### INVERSION BEGINS #########################","\033[0m","\n"))
    
    #Turn ON known population sigma for likelihood
    if(pop.sd == TRUE) {
      par0 = c(chl = initial[1], acdom440 = initial[2], anap440 = initial[3])#, bbp550 = initial[4])
    } else {
      #Turn ON for uknown population sigma for likelihood
      par0 = c(chl = initial[1], acdom440 = initial[2], anap440 = initial[3], pop.sd = initial[4])
    }
    
    cat(paste0("\033[0;46m","Initial values are ",par0[1],"    ", par0[2],"    ", 
               par0[3],"\033[0m","\n"))
    
    cat(paste0("\033[0;41m","#################### OPTIMIZATION INITIALIZING #########################","\033[0m","\n"))
    Sys.sleep(2)
    start.time = Sys.time()
    
    if (sa.model == "am03") {
      print("Albert & Mobley 2003 SA model used for Likelihood calculation")
    } else {
      print("Lee et al. 1999 SA model used for Likelihood calculation")
    }
    
    if (method.opt == "L-BFGS-B") {
      
      MLE_estimates = optim(par = par0, fn = NLL, data = obsdata, 
                            lower = lower.b,     # Lower bound on parameters
                            upper = upper.b,  # Upper bound on parameters
                            #gr = ,
                            method = method.opt,
                            control = list(parscale = abs(par0)),
                            hessian = TRUE)
    } else {
      if (method.opt == "levenberg-marqardt") {
        
        LM = marqLevAlg::marqLevAlg(b = par0, fn = NLL,data = obsdata, print.info = F)
        MLE_estimates = data.frame("par"=LM_estimates$b)
        
      } else {
        
        MLE_estimates = optim(par = par0, fn = NLL, data = obsdata, 
                              #lower = c(0, 0, 0),     # Lower bound on parameters
                              #upper = c(10, 2, 0.3),  # Upper bound on parameters
                              method = method.opt,
                              control = list(parscale = abs(par0)),
                              hessian = FALSE)
      }
      
    }
    
    
    cat(paste0("\033[0;46m","#################### OPTIMIZATION ENDS #########################","\033[0m","\n"))
    
    #print(MLE_estimates$par)
    
    cat(paste0("\033[0;41m","#################### CALCULATE UNCERTAINITY: START #########################","\033[0m","\n"))
    Sys.sleep(2)
    
    #Calculate hessian matrix for var-covar matrix
    hessian.inverse <- numDeriv::hessian(x =MLE_estimates$par, func = NLL, data=obsdata)
    rownames(hessian.inverse) <- names(par0)
    colnames(hessian.inverse) <- names(par0)
    cat(paste0("\033[0;42m","#################### VAR-COV HESSIAN MATRIX #########################","\033[0m","\n"))
    prmatrix(hessian.inverse)
    cat(paste0("\033[0;46m","#################### CALCULATE UNCERTAINITY: END #########################","\033[0m","\n"))
    
    MLE_par <- MLE_estimates$par
    MLE_SE <- sqrt(diag(solve(hessian.inverse))) #solve for diagonal elements to get sd
    MLE <- data.table("param" = names(par0),
                      "estimates" = MLE_par,
                      "sd(+/-)" = MLE_SE)
    
    print("The retrieved parameters are:")
    prmatrix(MLE)
    cat(paste0("\033[0;42m","#################### INVERSION ENDS #########################","\033[0m","\n"))
    if (MLE_estimates$convergence == 0) {
      convergence <- "TRUE"
    } else {
      convergence = "FALSE"
    }
    end.time = Sys.time(); time_taken <- end.time - start.time
    return(list(MLE,"convergence"= convergence, "time.elapsed"= time_taken))
  } else {
    #@@@@@@@@@@@@@@ SSR Objective function @@@@@@@@@@@@@@@@@@@@@@@@
    if (obj.fn == "SSR") {
      cat(paste0("\033[0;44m","#################### INVERSION BEGINS #########################","\033[0m","\n"))
      Fit <- data.frame("C_ph"=initial[1],
                        "a_cdom.440"=initial[2],
                        "a.nap.440"=initial[3])
      
      params <- as.numeric(Fit)
      cat(paste0("\033[0;46m","Initial values are ",params[1],"    ", params[2],"    ", 
                 params[3],"\033[0m","\n"))
      cat(paste0("\033[0;41m","#################### OPTIMIZATION INITIALIZING #########################","\033[0m","\n"))
      Sys.sleep(5)
      start.time = Sys.time()
      
      ##Single RUN
      inverse_output <- pracma::fminunc(fn = Saber_inverse,x0 = params,valdata=obsdata,
                                        bbp.550= Fit.input$bbp.550)
      cat(paste0("\033[0;46m","#################### OPTIMIZATION ENDS #########################","\033[0m","\n"))
      Fit.optimized.ssobj <- data.frame("chl"=inverse_output$par[1], "acdom.440"=inverse_output$par[2],
                                        "anap.440"=inverse_output$par[3])
      prmatrix(Fit.optimized.ssobj)
      #Calculate hessian matrix for var-covar matrix
      hessian.inverse <- numDeriv::hessian(x =as.numeric(Fit.optimized.ssobj), 
                                           func = Saber_inverse, valdata=obsdata,
                                           bbp.550= Fit.input$bbp.550)
      rownames(hessian.inverse) <- names(par0)
      colnames(hessian.inverse) <- names(par0)
      cat(paste0("\033[0;42m","#################### VAR-COV HESSIAN MATRIX #########################","\033[0m","\n"))
      prmatrix(hessian.inverse)
      cat(paste0("\033[0;46m","#################### CALCULATE UNCERTAINITY: END #########################","\033[0m","\n"))
      
      MLE_par <- as.numeric(Fit.optimized.ssobj)
      MLE_SE <- sqrt(diag(solve(hessian.inverse))) #solve for diagonal elements to get sd
      MLE <- data.table("param" = names(par0),
                        "estimates" = MLE_par,
                        "sd(+/-)" = MLE_SE)
      
      print("The retrieved parameters are:")
      prmatrix(MLE)
      cat(paste0("\033[0;42m","#################### INVERSION ENDS #########################","\033[0m","\n"))
      end.time = Sys.time(); time_taken <- end.time - start.time
      return(list(Fit.optimized.ssobj,"time.elapsed"= time_taken))
    }
  }
  
  
}
