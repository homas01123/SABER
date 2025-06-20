# Load necessary libraries
library(ggradar)
library(dplyr)
library(scales)  
library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(fmsb)
library(colormap)
 
inv_param_wise_sens_optimal = list_of_2380_dfs[[15]]

# Function to simulate Rrs and calculate RMSE row-wise for Sobol analysis
inverse_model_function <- function(x, obsdata, inv_param, wavelength) {
  
  # Initialize a matrix to store the simulated Rrs for each row of input
  simulated_Rrs_all <- list()
  rmse_vec <- vector()
  rrs_diff_all <- list()
  
  
  for (i in 1:nrow(x)) {
    # Extract parameters for this row
    new_data <- data.frame(
      gamma = x[i, 1], 
      sg = x[i, 2], 
      sd = x[i, 3], 
      class1 = x[i, 4], 
      class2 = x[i, 5], 
      class3 = x[i, 6]
    )
    
    # Interpolate bottom reflectance for this row
    interp_rb(
      wavelength_input = wavelength, 
      bottom_type = c(new_data$class1, new_data$class2, new_data$class3)
    )
    
    # Simulate Rrs using inversion-retrieved model parameters
    simulated_Rrs = Saber_forward_fast_sensitivity_test(
      use_true_IOPs = F, 
      use_manual_slope = T, 
      chl = #rep(median(inv_param$chl),nrow(inv_param)),
      inv_param$chl,
      Rrs_input_for_slope = NA,
      a_dg = #rep(median(inv_param$adg443),nrow(inv_param)),
        inv_param$adg443,
      bbp.550 = #rep(median(inv_param$bbp555),nrow(inv_param)),
        inv_param$bbp555,
      slope.parametric = F, 
      manual_slope = c("s_g"=as.numeric(new_data$sg), 
                       "s_d"=as.numeric(new_data$sd), 
                       "gamma"=as.numeric(new_data$gamma)),
      z =  #rep(median(inv_param$H), nrow(inv_param)),  
      inv_param$H,
      # rb.fraction = as.data.frame(cbind(
      #   matrix(rep(apply(inv_param[c("fa1")], 2, median), 
      #              each = nrow(inv_param)), nrow = nrow(inv_param)),
      #   
      #   matrix(rep(apply(inv_param[c( "fa2")], 2, median), 
      #              each = nrow(inv_param)), nrow = nrow(inv_param)),
      #   
      #   matrix((inv_param["fa3"]$fa3), nrow = nrow(inv_param))
      # ))
      # ,
      inv_param[c("fa1", "fa2", "fa3")],
      wavelength = wavelength,
      verbose = F
    )
    
    rrs_diff = obsdata - simulated_Rrs
    
    RMSE <- apply(rrs_diff, 1, function(row) sqrt(mean(row^2)))
    
    rmse_vec[i] = mean(RMSE)
    simulated_Rrs_all[[i]] <- simulated_Rrs
    rrs_diff_all[[i]] <- rrs_diff
    
  }
  
  
  return(rmse_vec)
}


inverse_model_function(x = as.matrix(data.frame("gamma" = seq(0.1, 1, length.out = 20), 
                                      "sg" = seq(0.001, 0.0015, length.out = 20), 
                                      "sd" = seq(0.0001, 0.001, length.out = 20), 
                                      "class1" = combinations[,1],
                                      "class2" = combinations[,2], 
                                      "class3" = combinations[,3])), 
                       obsdata = (inverse_rrs_input_mod), 
                       inv_param = inv_param_wise_sens_optimal, wavelength = seq(400,750,1))


# Sampling from uniform distributions
set.seed(345)
n = 1820
x = (data.frame("gamma" = seq(0.1, 1, length.out = 20), 
                         "sg" = seq(0.001, 0.0015, length.out = 20), 
                         "sd" = seq(0.0001, 0.001, length.out = 20), 
                         "class1" = combinations[,1],
                         "class2" = combinations[,2], 
                         "class3" = combinations[,3]))
X1 <- data.frame(
  gamma = runif(n, min(x$gamma), max(x$gamma)),
  sg = runif(n, min(x$sg), max(x$sg)),
  sd = runif(n, min(x$sd), max(x$sd)),
  class1 = sample(x$class1, n, replace = TRUE),
  class2 = sample(x$class2, n, replace = TRUE),
  class3 = sample(x$class3, n, replace = TRUE)
)

X2 <- data.frame(
  gamma = runif(n, min(x$gamma), max(x$gamma)),
  sg = runif(n, min(x$sg), max(x$sg)),
  sd = runif(n, min(x$sd), max(x$sd)),
  class1 = sample(x$class1, n, replace = TRUE),
  class2 = sample(x$class2, n, replace = TRUE),
  class3 = sample(x$class3, n, replace = TRUE)
)

# Perform the Sobol sensitivity analysis
sobol_results <- sobolSalt(
  model = inverse_model_function,  # Use the inverse model function
  X1 = as.matrix(X1), 
  X2 = as.matrix(X2), 
  obsdata = (inverse_rrs_input_mod), 
  inv_param = (inv_param_wise_sens_optimal), 
  wavelength = seq(400, 750, 1),
  nboot = 1000, scheme = "A", conf = 0.95
)

# View results
print(sobol_results)



sobol_df <- data.frame(
  #group = "chl",
  "variance" = sobol_results[["T"]]$original*100
  
)

sobol_df = as.data.frame(t(sobol_df))
colnames(sobol_df) = c("gamma", "s_g", "s_d","benthic_class1", "benthic_class2", "benthic_class3")
sobol_df$s_dg <- sobol_df$s_g + sobol_df$s_d
sobol_df[c("s_g", "s_d")] = NULL
sobol_df <- sobol_df[, c(names(sobol_df)[1], "s_dg", 
                         names(sobol_df)[-c(1, which(names(sobol_df) == "s_dg"))])]
sobol_df$Group = "group"
sobol_df = sobol_df %>% dplyr::select(Group, everything())
sobol_df <- rbind(rep(100, nrow(sobol_df)), rep(0, nrow(sobol_df)), (sobol_df))
rownames(sobol_df) = seq(1, nrow(sobol_df),1)

sobol_df[3,-1] = c(7.57751, 64.7859, 7.40912, 4.20, 11.70) #for fa1

# Prepare color
colors_border = colormap(colormap=colormaps$viridis, nshades=6, alpha=1)
colors_in = colormap(colormap=colormaps$viridis, nshades=6, alpha=0.3)

par(mar = c(1, 2, 2, 2))
# Create radar chart
png("beautiful_radar_chart.png", width = 600, height = 600)  # Save plot as PNG
par(mar = c(1, 2, 2, 2))

input_df = lapply(sobol_df[,-1], function(col) if(is.numeric(col)) col else as.numeric(col))
input_df = as.data.frame(do.call(cbind, input_df))
input_df[3,] = (input_df[3,]/rowSums(input_df[3,]))*100

radarchart(input_df, axistype = 1, 
           pcol = c(colors_border), 
           pfcol = c(colors_in), 
           plwd=4, cglcol = "grey", cglty = 1, axislabcol = "grey", caxislabels = seq(0, 100, 25), 
           vlabels = c(expression(paste(gamma)),
                        expression(paste("s"["dg"])),
                       expression(paste("R"["B,1"])),
                        expression(paste("R"["B,2"])),
                        expression(paste("R"["B,3"]))
                       #expression(paste(epsilon))
                       ),
           title = expression(paste("[chl]")),
           cglwd = 1, vlcex = 1.3)
dev.off()





par(mar = c(5, 4, 4, 2) + 0.1)

labels <- c(
  bquote(s[dg]),                 # Spectral slope (sdg)
  bquote(R[B]),                   # Benthic reflectance (RB)
  bquote(R[rs,inel]),             # Remote sensing reflectance with inelastic scattering (Rrs,inel)
  bquote(epsilon)                 # Greek letter epsilon
)

labels = as.vector(unlist(labels))

radar_chl = ggradar(
  sobol_df[3,], 
  values.radar = c("0%", "50%", "100%"),
  grid.min = 0, grid.mid = 50, grid.max = 100,
  axis.labels = c(
    " ",           
    " ",         
    " ", 
    " ", " ") 
  #c("Spectral Slope", "Benthic Type", "Inelastic Scattering", "Model Noise")
  ,
  # Polygons
  group.line.width = 1.3, fill = T, fill.alpha = 0.3,
  
  group.point.size = 3, plot.legend = F, draw.points = T,
  group.colours =# "orange2" #"#00AFBB"
  colors_border[5] 
    , 
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "grey", grid.line.width = 1.1 
) 
ggsave(paste0("./outputs/uncertainity_contrib_chl.png"), plot = radar_chl, scale = 1, 
       width = 4.5, height = 4.5,
       units = "in",dpi = 300)
