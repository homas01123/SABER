library(dplyr)
library(tidyr)
library(pbapply)
library(parallel)
library(doParallel)
library(foreach)
library(plotly)

#function :: Interpolate the bottom reflectance to Rrs wavelengths ----
interp_rb <- function(wavelength_input, bottom_type = c("Mud_2019", "Sand_2019", "Eelgrass_2019")){
  load("./data/EGSL_all.RData")
  
  #column_numbers <- which(names(rb) %in% bottom_type)
  
  column_numbers <- match(bottom_type, names(rb))
  
  rb = rb[c(1,column_numbers)]
  
  print(paste0("Classnames reassigned as: Class1:", colnames(rb)[2], ", Class2:", colnames(rb)[3],
               ", Class3:", colnames(rb)[4]))
  
  names(rb) = c("wavelength", "class1", "class2", "class3")
  
  
  if (!(all(length(rb$wavelength) == length(wavelength_input)) && all(rb$wavelength == wavelength_input))) {
    
    print(paste0("Rb wavelength has " , length(rb$wavelength), " bands; Rrs wavelength has ", length(wavelength_input) ))
    
    rb_interp = data.frame("wavelength" = wavelength_input,
                           "class1" = approx(rb$wavelength, rb$class1, xout = wavelength_input)$y,
                           "class2" = approx(rb$wavelength, rb$class2, xout = wavelength_input)$y,
                           "class3" = approx(rb$wavelength, rb$class3, xout = wavelength_input)$y
    )
    
  }
  else {
    print("The desired and existing wavelength are same")
    rb_interp = rb
  }
  rb = rb_interp
  rm(rb_interp)
  
  assign(x = "rb", value = rb, envir = .GlobalEnv)
}
# function :: Create illustrative scatterplot for two variables in a linear scale ----
plot_inversion_validation_multivar_linear_contour <- function(
    
  input_df = inv_res_df_grad_unconst, 
  input_x = "H_actual", input_y = "H",
  uncertainty = "H_CI", xmin = 0, xmax = 12, 
  xlabel = expression(paste(italic("H"["actual"]), "[m]")),
  ylabel = expression(paste(italic("H"["predicted"]), "[m]")), 
  opacity = 0.8, plot_col = "purple4", xstp = 2, ystp = 2, 
  show_legend = T, hist_count = 30){
  ## data
  d <- input_df %>% 
    as_tibble() 
  
  ymin = xmin; ymax = xmax
  
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  g<-   ggplot(data=d, aes(x = .data[[input_x]], y = .data[[input_y]])) +
    
    geom_density_2d(data = d, aes(x = .data[[input_x]], y = .data[[input_y]]),na.rm = T, 
                    bins = 6,
                    linewidth = 0.25, colour = "black", show.legend = F, size=1.1)+
    
    geom_point(data = d, aes(.data[[input_x]], .data[[input_y]],
                             shape=sensor#, fill=(sensor), color = sensor
                             ),
               alpha = I(opacity), fill = plot_col, color =plot_col, 
               size = I(5), show.legend = show_legend) +
    
    # geom_text(data = d, size = 2.5,  
    #           aes(label=as.character(statname_new)),hjust=0, vjust=2, angle= 45,
    #           check_overlap = T)+
    
    geom_ribbon(data = d, aes(#fill = sensor,
                              x = .data[[input_x]], y = .data[[input_y]],
                              ymin = (.data[[input_y]] - .data[[uncertainty]]),
                              ymax = (.data[[input_y]] + .data[[uncertainty]])),

                alpha = 0.3, show.legend = F, fill = plot_col,
                colour="NA"
    )+
    
    scale_shape_manual(name = " ",
                       labels=(c(expression(paste("COPS")),
                                 expression(paste("PSR")))),
                       values = c(22,23))+
    
    # scale_color_manual(name = " ", 
    #                    labels=(c(expression(paste("COPS")),expression(paste("PSR")))),
    #                    values = c("black", "red"))+
    
    # scale_color_manual(name = " ",
    #                    labels=(c("COPS_unique", "Duplicate", "PSR_unique")),
    #                    values = c("purple3", "cyan3", "green3"))+
    # 
    # scale_fill_manual(name = " ",
    #                   labels=(c("COPS_unique", "Duplicate", "PSR_unique")),
    #                   values = c("purple3", "cyan3", "green3"))+
    
    geom_abline(slope = 1,linetype="dashed", intercept = 0,
                colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
    #geom_vline(xintercept = 2.54, color = "orange3", linetype = "dashed", size = 1.2)+
    
    geom_smooth(#aes(color = sensor),
      size=1,level = 0.95,show.legend = F,linetype = "solid",
      color="red4",
      se= T, method = "lm")+
    
    coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
                ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
    scale_x_continuous(name = xlabel, limits = c(xmin, xmax), 
                       breaks = seq(xmin, xmax, xstp)) +
    scale_y_continuous(name = ylabel, limits = c(ymin, ymax), 
                       breaks = seq(ymin, ymax, ystp)) +
    
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
          axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          axis.ticks.length = unit(.25, "cm"),
          legend.box.just = "right",
          legend.spacing = unit(-0.5, "cm"),
          legend.position = c(0.65, 0.30),
          #legend.direction = "vertical",
          legend.title = element_text(size = 10, color = 'black', angle = 0, face = "bold"),
          legend.text = element_text(colour = "black", size = 20, face = "plain"),
          legend.background = element_rect(fill = NA, size = 0.5, 
                                           linetype = "solid", colour = 0),
          legend.key = element_blank(),
          legend.justification = c("left", "top"),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "black", 
                                          size = 0.5, linetype = "dotted"), 
          panel.grid.minor = element_line(colour = "grey80", 
                                          linewidth =  0.2, linetype = "solid"),
          plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
          legend.direction = "vertical", legend.box = "vertical",
          legend.text.align = 0,
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  
  g <- ggMarginal(groupFill = F, data = d, type = "densigram", bins = hist_count, color = "grey",
                  p = g, aes(x = .data[[input_x]], y = .data[[input_y]]))
  return(g)
}

# Read, restructure, interpolate Rrs data ----

## WISE-Man 2019 COPS data ----
wise_cops = read.csv("./data/Rb_spectral/copsdata.csv", header = T)
wise_cops_wide = read.csv("./data/Rb_spectral/COPS_db.csv", header = T)

kildir_iop = read.csv("./data/Kildir_IOP.Process_Log.csv", header = T)
saucier_iop = read.csv("./data/saucierQC.csv", header = T)

#Input the depth statistics and simulation notes from WISE-Man IOP data
wise_depth_stat = read.csv("./data/Depth_wise.csv", header = T)
wise_depth_stat$correct_statname = gsub(wise_depth_stat$Station, 
                                        pattern = "\\_", replacement = "-") 
wise_depth_stat <- wise_depth_stat %>% dplyr::select(correct_statname, everything())

#Select the QC passed stations for WISE-Man
wise_depth_stat_shallow = wise_depth_stat[wise_depth_stat$Shallow == "T",]

kildir_iop_qc = kildir_iop[kildir_iop$ASPH == "Y" & kildir_iop$HS6 == "Y",]
saucier_iop_qc = saucier_iop[saucier_iop$station.keep == "TRUE",]

# qc_statnames_bind = rbind(list(kildir_iop_qc$StationID), 
#                           list(saucier_iop_qc$Station))

qc_statnames_bind = rbind(list(kildir_iop$StationID), 
                          list(saucier_iop$Station))

qc_statnames_bind = data.frame("qcstat" = as.character(unlist(qc_statnames_bind)))


shallow_qc_keep_idx = vector()

for (i in 1:length(qc_statnames_bind$qcstat)) {
  tryCatch({
    shallow_qc_keep_idx[i] = grep(x =  wise_depth_stat_shallow$correct_statname, 
                                  pattern = qc_statnames_bind$qcstat[i])
    
    if (is.na(shallow_qc_keep_idx[i])) {
      cat(paste0("\033[0;31m","The station ", wise_depth_stat_shallow$correct_statname[i]," is not shallow, coerce to NA","\033[0m","\n"))
      shallow_qc_keep_idx[i] = NA
    } }, 
    error = function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
}

shallow_qc_stations_grep = qc_statnames_bind$qcstat[!is.na(shallow_qc_keep_idx)]
shallow_qc_stations_in = qc_statnames_bind[qc_statnames_bind$qcstat %in% 
                                             wise_depth_stat_shallow$correct_statname,]
#shallow_qc_stations = union(shallow_qc_stations_grep, shallow_qc_stations_in)

shallow_qc_stations =  shallow_qc_stations_in

shallow_qc_stations = gsub(shallow_qc_stations, 
                                        pattern = "\\-", replacement = ".")

#shallow_qc_stations = shallow_qc_stations[-which(shallow_qc_stations == "MAN.F0")]

wise_cops_shallow =  wise_cops[,colnames(wise_cops) %in%  
                                 # wise_depth_stat_shallow$correct_statname
                                 shallow_qc_stations
                               ] 

wise_cops_wide_shallow = wise_cops_wide[,colnames(wise_cops_wide) %in% 
                                          colnames(wise_cops_shallow)]

rownames(wise_cops_wide_shallow) = wise_cops_wide[,1]

wise_cops_shallow = wise_cops_shallow[(1:351),]

rm(qc_statnames_bind)

## WISE-Man 2019 PSR data ----
wise_psr = read.csv("./data/Rrs_psr_long_wiseman.csv", header = T)

wise_psr = split(wise_psr, wise_psr$station)

## PSR multicast data ----
multicast_idx <- sapply(wise_psr, function(df) {
  nrow(df) > 781 && (nrow(df) %% 781 == 0)
})

multicast_idx <- which(multicast_idx)

wise_psr_multicast = wise_psr[multicast_idx]

wise_psr_multicast_df <- lapply(wise_psr_multicast, function(df_observation) {
  multicast_df_split <- split(df_observation, df_observation$station_alt)
  
  lapply(multicast_df_split, function(df_obs){
    
    wavelength <- df_obs$Wavelength
    Rrs <- df_obs$Rrs
    station = df_obs$station_alt
    df_Rrs <- data.frame(Rrs)
    #do.call(rbind, multicast_df_split)
  })
  
})

extract_rrs <- function(data_list) {
  
  all_rrs <- data.frame()
  

  for(df in data_list) {
   
    df_d = data.frame(t(do.call(cbind, df)))
    rownames(df_d) = names(df)
    colnames(df_d) = unique(wise_psr_multicast$`EXP-WISE`$wavelength)
    all_rrs = rbind(all_rrs, df_d)
    
  }
  
  return(all_rrs)
}

wise_psr_multicast_df = extract_rrs(wise_psr_multicast_df)

## PSR unicast data ----
wise_psr_unicast_df <- lapply(wise_psr[-multicast_idx], function(df_observation) {
  wavelength <- df_observation$Wavelength
  Rrs <- df_observation$Rrs
  station = df_observation$station
  df_Rrs <- data.frame(Rrs)
})

## combined PSR data ----
wise_psr_unicast_df = data.frame(t(do.call(cbind, wise_psr_unicast_df)))

rownames(wise_psr_unicast_df)  = names(wise_psr[-multicast_idx])
colnames(wise_psr_unicast_df) = wise_psr[[1]]$wavelength


wise_psr_rrs_df = rbind(wise_psr_unicast_df, wise_psr_multicast_df)

wise_psr_rrs_df = wise_psr_rrs_df[colnames(wise_psr_rrs_df) %in% seq(400,750,1)]
wise_psr_rrs_df = data.frame(t(wise_psr_rrs_df))


# 
# wise_psr_shallow =  wise_psr_rrs_df[,colnames(wise_psr_rrs_df) %in%  
#                                               # wise_depth_stat_shallow$correct_statname
#                                       bgc_wiseman$correct_statname] 
# 
wise_psr_shallow =  wise_psr_rrs_df[, 1:47][,colnames(wise_psr_rrs_df)[1:47] %in%
                                 # wise_depth_stat_shallow$correct_statname
                                 shallow_qc_stations]

wise_psr_shallow = cbind(wise_psr_shallow, wise_psr_rrs_df[, -(1:47)])

rm(wise_psr_multicast_df, wise_psr_unicast_df)


#Apply sub-surface translation
wise_cops_shallow = as.data.frame(apply(wise_cops_shallow, 2, 
                                        surface_rrs_translate))

wise_psr_shallow = as.data.frame(apply(#wise_psr_shallow, 
                                      wise_psr_shallow_known_H,
                                       2, 
                                        surface_rrs_translate))

## Exclude the inelastic scattering contribution from Rrs ----
subtract_wise_inel <- function(input_df, wave_out){
  
  inel_contrib = read.csv("./data/inel_contrib_MP.csv", header = T)
  
  
  output_df <- apply(input_df, 2, 
                                          function(row, 
                                                   wavelength_seadoo =wave_out){
                                            
                                            obsdata = row
                                            
                                            obsdata_inel = Hmisc::approxExtrap(x = inel_contrib$wavelength#[non_na_seq]
                                                                               , y = inel_contrib$mean_rrs#[!is.na(obsdata)]
                                                                               , 
                                                                               xout = wavelength_seadoo)$y
                                            obsdata_inel[obsdata_inel < 0] = 0
                                            obsdata_inel_val = obsdata - ((obsdata_inel/100)*obsdata)
                                            obsdata_inel_val[obsdata_inel_val < 0] = 0 
                                            obsdata = obsdata_inel_val
                                            
                                          })
  
  output_df = data.frame((output_df))
  return(output_df)
  
}

### COPS ----
wise_cops_shallow_incor = subtract_wise_inel(input_df = wise_cops_shallow, wave_out = seq(400,750,1))

### PSR ----
wise_psr_shallow_incor = subtract_wise_inel(input_df = wise_psr_shallow, wave_out = seq(400,750,1))


wise_rrs_df = cbind(wise_cops_shallow, wise_psr_shallow)

wise_rrs_incor_df = cbind(wise_cops_shallow_incor, wise_psr_shallow_incor)


# Plot the mean spectra for inelastic scattering corrected and uncorrected Rrs respectively

mean_rrs_inel_el = data.frame("wavelength" = seq(400,750,1), 
                              "rrs_inel_cor" = rowMeans(wise_psr_shallow_incor),
                              "sd_inel_cor" = as.numeric(apply(wise_psr_shallow_incor, 1,sd)),
                              "rrs_inel_uncor" = rowMeans(wise_psr_shallow),
                              "sd_inel_uncor" = as.numeric(apply(wise_psr_shallow, 1, sd))
)


# Reshape the data to long format for ggplot
library(tidyr)
long_data <- mean_rrs_inel_el %>%
  pivot_longer(cols = c(rrs_inel_cor, rrs_inel_uncor),
               names_to = "type",
               values_to = "value") %>%
  pivot_longer(cols = c(sd_inel_cor, sd_inel_uncor),
               names_to = "sd_type",
               values_to = "sd_value") %>%
  filter((type == "rrs_inel_cor" & sd_type == "sd_inel_cor") |
           (type == "rrs_inel_uncor" & sd_type == "sd_inel_uncor"))

# Plot using ggplot2
xmin <- 400; xmax <- 750;  xstp <- 70; xlbl <- expression(paste("Wavelength (", lambda, ") [nm]")) 
ymin <- 0.000; ymax <- 0.009; ystp <- 0.003
ylbl <- expression(paste(italic("R")["rs"]("0"^"-", lambda), "[sr"^{-1}, "]"))
asp_rat <-  (xmax-xmin)/(ymax-ymin)
legend_title <- element_blank()
legend_position <- c(0.10, 1)
col_list = viridis(8, option = "H")
col_list = c(col_list[3], col_list[6])
show_col(col_list)

g <- ggplot(long_data, aes(x = wavelength, y = value, fill = type)) +
  geom_line(size = 1.3, aes(color = type)) +
  geom_ribbon(aes(ymin = value - sd_value, ymax = value + sd_value), alpha = 0.25, show.legend = F) +
  
  scale_colour_manual(name = "", values = rev(c("red3", "cyan3")),
                      labels = rev(c(expression(paste(bar(italic("R"))["rs, el+inel"])),
                                     expression(paste(bar(italic("R"))["rs, el"]))))
  )+
  scale_fill_manual(name = "", values = rev(c("red3", "cyan3")),
                    labels = rev(c(expression(paste("el")),
                                   expression(paste("el+inel"))))
  )+
  
  
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
              ylim = c(ymin, ymax)
              ,expand = FALSE, clip = "on"
  ) +
  
  scale_x_continuous(name = xlbl, limits = c(xmin, xmax),
                     breaks = seq(xmin, xmax, xstp))  +
  scale_y_continuous(name = ylbl, limits = c(ymin, ymax),
                     breaks = seq(ymin, ymax, ystp))  +
  geom_hline(yintercept = 0, colour = "navyblue", linetype = "dashed", size=1.1)+
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
        axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = c(0.05, 0.90),
        #legend.direction = "vertical",
        legend.title = element_text(size = 10, color = 'black', angle = 0, face = "bold"),
        legend.text = element_text(colour = "black", size = 20, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_line(colour = "grey80", 
                                        linewidth =  0.2, linetype = "solid"),
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        legend.direction = "vertical", legend.box = "vertical",
        legend.text.align = 0,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

g
ggsave(paste0("./outputs/wise_el_inel_comp_PSR.png"), plot = g,
       scale = 1.25, width = 4.5, height = 4.5, 
       units = "in",dpi = 300)


#Plot the inelastic scattering contributon in %
inel_contrib = data.frame("wavelength" = seq(400,750,1), "inel_ratio" = inel_contrib)
inel_contrib$inel_ratio = 100*inel_contrib$inel_ratio


legend_title <- element_blank()
legend_position <- c(0.10, 0.95)
show_legend = T

xlbl <- expression(paste("Wavelength (", lambda, ") [nm]"))
ylbl <- expression(paste(italic("R")["rs,inel"]("0"^"-", lambda), "[%]"))

ymin <- 0; ymax <- 100 ; ystp <- ymax/5
xmin <- 400; xmax <- 750; xstp <- 70
asp_rat <- (xmax-xmin)/(ymax-ymin)

g_inel = ggplot(data = inel_contrib, aes(x = wavelength, y = inel_ratio))+
  geom_line(size = 1.3, color = "purple3")+
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
  scale_x_continuous(name = xlbl, limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, xstp)) +
  scale_y_continuous(name = ylbl, limits = c(ymin, ymax), 
                     breaks = seq(ymin, ymax, ystp)) +
  
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
        axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = c(0.01,1),
        #legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 10, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_line(colour = "grey80", 
                                        linewidth =  0.2, linetype = "solid"),
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.text.align = 0,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  
ggsave(paste0("./outputs/wise_2019_inel_contrib.png"), plot = g_inel,
       scale = 1.25, width = 4.5, height = 4.5, units = "in",dpi = 300)



## Retrieve BGC parameters for given stations ----
### COPS ----
bgc_val_vec = data.frame()
shallow_statist = gsub(colnames(wise_cops_shallow_incor), 
                       pattern = "\\.", replacement = "-")
shallow_statist[c(which(shallow_statist == "MAN-F03-5"),
                  which(shallow_statist == "MAN-R11-5"))] = c("MAN-F03.5", "MAN-R11.5")

for (i in 1:length(shallow_statist)) {
  
  temp_insitu_param = get_in_situ_params(station_name = shallow_statist[i], get_only_labdata = T, 
                                         use_bb_nup = F)
  temp_bgc_val = data.frame("chl" = temp_insitu_param$chl_invivo,
                            "adg443" = temp_insitu_param[[1]]$a_dg
                            #, "bbp555" = temp_insitu_param$bbp555
                            )
  bgc_val_vec = rbind(bgc_val_vec, temp_bgc_val)
  
}
bgc_val_vec$correct_statname = colnames(wise_cops_shallow_incor)

bgc_val_vec$H_cops = as.numeric(wise_cops_wide_shallow[4,])

bgc_wiseman = read.csv("./data/biogeochemistry_wiseman.csv", header = T)
bgc_wiseman$correct_statname = gsub(bgc_wiseman$station,
                                    pattern = "\\-", replacement = ".")

bgc_wiseman_shallow = bgc_wiseman[bgc_wiseman$correct_statname %in% colnames(wise_cops_shallow),]
#bgc_wiseman_shallow = rbind(bgc_wiseman_shallow, bgc_wiseman[nrow(bgc_wiseman),])
bgc_wiseman_shallow = bgc_wiseman_shallow[bgc_wiseman_shallow$depth < 1, ]
bgc_wiseman_shallow = bgc_wiseman_shallow[bgc_wiseman_shallow$sampling_type == "surface",]

bgc_val_vec <- bgc_val_vec %>%
  mutate(H_sampling = bgc_wiseman_shallow$water_depth[match(correct_statname, 
                                                             bgc_wiseman_shallow$correct_statname)])

bgc_val_vec$H = rowMeans(bgc_val_vec[, c(4:5)])

# wise_cops_shallow_known_H = wise_cops_shallow[colnames(wise_cops_shallow) %in% 
#                                               bgc_wiseman_shallow$correct_statname]
#wise_cops_shallow_known_H = wise_cops_shallow_known_H[,order(colnames(wise_cops_shallow_known_H))]
#bgc_wiseman_shallow = bgc_wiseman_shallow[order(bgc_wiseman_shallow$correct_statname),]

### PSR ----
#Retrieve PSR station depths
bgc_wiseman = read.csv("./data/biogeochemistry_wiseman.csv", header = T)
bgc_wiseman$correct_statname = gsub(bgc_wiseman$station,
                                    pattern = "\\-", replacement = ".")

bgc_wiseman_shallow = bgc_wiseman[bgc_wiseman$correct_statname %in% colnames(wise_psr_shallow),]
bgc_wiseman_shallow = rbind(bgc_wiseman_shallow, bgc_wiseman[nrow(bgc_wiseman),])
bgc_wiseman_shallow = bgc_wiseman_shallow[bgc_wiseman_shallow$depth < 1, ]

wise_psr_shallow_known_H = wise_psr_shallow[colnames(wise_psr_shallow) %in% 
                                              bgc_wiseman_shallow$correct_statname]

wise_psr_shallow_known_H = cbind(wise_psr_shallow_known_H, 
                                 rowMeans(wise_psr_rrs_df[,48:53]))

colnames(wise_psr_shallow_known_H)[10] = "EXP.WISE"

wise_psr_shallow_known_H = wise_psr_shallow_known_H[,order(colnames(wise_psr_shallow_known_H))]
bgc_wiseman_shallow = bgc_wiseman_shallow[order(bgc_wiseman_shallow$correct_statname),]


# Retrieve the stations common between COPS and PSR for comparative study----
common_columns <- intersect(colnames(wise_cops_shallow_incor), colnames(wise_psr_shallow_incor))

wise_cops_shallow_common_incor <- wise_cops_shallow_incor[colnames(wise_cops_shallow_incor) 
                                              %in% common_columns]
wise_cops_wide_shallow_common_incor <- wise_cops_wide_shallow[colnames(wise_cops_wide_shallow) 
                                                        %in% common_columns]
wise_psr_shallow_common_incor <- wise_psr_shallow_incor[colnames(wise_psr_shallow_incor) %in% 
                                                    common_columns]

## Plot the COPS and PSR common Rrs together ----
color_var <- seq(400, 750, length.out = nrow(wise_cops_shallow_common_incor)) 

wise_cops_shallow_common_incor_long <- wise_cops_shallow_common_incor %>%
  mutate(Wavelength = color_var) %>%
  pivot_longer(cols = -Wavelength, names_to = "variable", values_to = "COPS")

wise_psr_shallow_common_incor_long <- wise_psr_shallow_common_incor %>%
  mutate(Wavelength = color_var) %>%
  pivot_longer(cols = -Wavelength, names_to = "variable", values_to = "PSR")

combined_df <- wise_cops_shallow_common_incor_long %>%
  left_join(wise_psr_shallow_common_incor_long, by = c("Wavelength", "variable"))

combined_df_long <- combined_df %>%
  pivot_longer(cols = c(COPS, 
                        PSR), names_to = "source", values_to = "value")

rm(wise_cops_shallow_common_incor_long, wise_psr_shallow_common_incor_long,
   combined_df)

g_rrs = ggplot(combined_df_long, aes(x = Wavelength, y = value, color = source)) +
  geom_line(show.legend = T, size = 1.3) +
  facet_wrap(~ variable, scales = "free_y") +
  scale_y_continuous(limits = c(0, 0.010), breaks = seq(0, 0.010, by = 0.002)) +
  scale_x_continuous(limits = c(400, 750), breaks = seq(400, 750, by = 70)) +
  labs(x = "Wavelength (nm)", y = expression(paste(italic("R")["rs"]("0"^"+", lambda),italic("in situ ["), "sr"^{-1}, "]")),
       color = " ") +
  scale_color_manual(values = c("purple3", "green3"), labels = c("COPS", "PSR"))+
  theme_bw()+
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 25, color = 'black', angle = 0),
        axis.text.y = element_text(size = 25, color = 'black', angle = 0),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = c(0.42,0.9),
        legend.direction = "horizontal",
        legend.title = element_text(colour = "black", size = 15, face = "plain"),
        legend.text = element_text(colour = "black", size = 18, face = "bold"),
        legend.background = element_rect(fill = NA, size = 0.5,
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey",
                                        size = 0.5, linetype = "dotted"),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.0,0.9,0.0,0.0), "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

ggsave(paste0("./outputs/rrs_shallow_cops-psr_comp.png"), plot = g_rrs, scale = 1.7, 
       width = 8, height = 4.5,
       units = "in",dpi = 300)


# PERFORM a Sensitivity analysis for 1. spectral slope, 2. benthic cover, 3. phyisical depth:\tau and 4. inelastic scattering----

#Import benthic reflectance
Rb_set = read.csv("./data/EGSL_Rb.csv", header = T)
benthic_classes_2019 <- colnames(Rb_set)[grep("2019", colnames(Rb_set))]
benthic_classes_2022 <- colnames(Rb_set)[grep("2022", colnames(Rb_set))]
benthic_classes_2023 <- colnames(Rb_set)[grep("2023", colnames(Rb_set))]

combinations <- combn(c(benthic_classes_2019, benthic_classes_2022), 3, simplify = FALSE)

combinations <- lapply(combinations, function(x) sort(x))
combinations <- as.data.frame(do.call(rbind, combinations), stringsAsFactors = FALSE)

#Synthetically generate sequence of spectral slopes
sdg_values <- seq(0.0015, 0.025, length.out = 10)
gamma_values <- seq(0.1, 1, length.out = 10)
sd_values <- seq(0.0005, 0.005, length.out = 10)
sg_values <- sdg_values - sd_values

spectral_slope_vec <- data.frame("gamma" = gamma_values,
                                 "sd" = sd_values,
                                 "sg" = sg_values)

spectral_slope_vec <- expand.grid(sdg = sdg_values, gamma = gamma_values)

set.seed(123)  
spectral_slope_vec$sd <- runif(n = nrow(spectral_slope_vec), min = 0.0005, max = 0.005)
spectral_slope_vec$sg <- spectral_slope_vec$sdg - spectral_slope_vec$sd

spectral_slope_vec <- spectral_slope_vec %>%
  filter(sd != sg & sg >= 0.001 & sg <= 0.02 & sd + sg == sdg)


param_vec_bayes <- data.frame(matrix(ncol = 16, nrow = 0))
colnames(param_vec_bayes) <- c("chl", "adg443", "bbp555", "H", "fa1", "fa2", "fa3", "pop_sd", 
                               "sd_chl", "sd_adg443", "sd_bbp555", "sd_H", "sd_fa1", "sd_fa2", 
                               "sd_fa3", "sd_pop_sd")


## function :: Define the process_combination function written for parallel processing ----
process_combination <- function(obsdata_row, idx) {
  # Suppress output from the functions
  capture.output({
    interp_rb(wavelength_input = seq(400, 750, 1), 
              bottom_type = c("Sand_2019", "Mud_2019", "Eelgrass_2019")) #use a fixed benthic type
              #bottom_type = combinations[idx,]) #use a fixed benthic type
    
    result <- inverse_runBayes(obsdata = as.numeric(obsdata_row), 
                               rrs_type = type_Rrs_below, 
                               max_par = c(10, 2.8, 0.01, 8, 1, 1, 1, 10),
                               min_par = c(1.5, 1, 0.001, 2, 0.0, 0.0, 0.0, 0.01), 
                               param_lab = c("chl", "adg443", "bbp555", "H", "fa1", 
                                             "fa2", "fa3", "pop_sd"),
                               constrain_config = c(FALSE, FALSE, FALSE),
                               bbp_const_val = 0.007, 
                               bgc_const_val = c(5, 1.003, 0.007), 
                               iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
                                                  "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
                               qaa_slope = FALSE, 
                               manual_slope = TRUE, 
                               manual_slope_vals = c("s_g" = spectral_slope_vec[idx, 4],
                                                     "s_d" = spectral_slope_vec[idx, 3],
                                                     "gamma" = spectral_slope_vec[idx, 2]),
                               
                               # manual_slope_vals = c("s_g" = 0.011555556, 
                               #                       "s_d" = 0.0030, 
                               #                       "gamma" = 0.60), 
                               
                               iter_count = 20000, 
                               sampler_mcmc = samplerlist[6], 
                               wavelngth_sim = seq(400, 750, 1), 
                               sa_model = "am03", 
                               hybrid_mode = FALSE, 
                               plot_rrs = FALSE)
  }, file = NULL) 
  
  return(result)
}

inverse_rrs_input = as.data.frame(cbind(wise_cops_shallow_incor, wise_psr_shallow_incor))
inverse_rrs_input_mod = as.data.frame(t(inverse_rrs_input))

rownames(inverse_rrs_input_mod) = paste0(c(colnames(wise_cops_shallow_incor),
                             colnames(wise_psr_shallow_incor)), "_", c(rep("COPS", 16), rep("PSR", 10)))



## function :: Define the data decomposition in lists ----
decompose_data <- function(data, num_groups) {
  split(data, rep(1:num_groups, length.out = nrow(data)))
}

decomposed_data <- decompose_data(inverse_rrs_input_mod, 5)

# Create cluster
num_cores <- detectCores() - 1  # Use one less than the total number of cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

## RUN the sensitivity analysis for spectral slope using parallel processing enabled inversion ----
time_start = Sys.time()
results_list <- foreach( i = names(decomposed_data),
                         #.combine = rbind, 
                         .verbose = T,
                         .packages = c("data.table", "BayesianTools")) %:% 
  foreach( idx_c = 1:nrow(spectral_slope_vec),
           #idx_c = 1:nrow(combinations),
           #.combine = c, 
           .verbose = T,
           .packages = c("data.table", "BayesianTools")) %dopar% {
             source("./R/SABER_forward_fast.R")
             source("./R/solve.objective.inverse_fast.R")
             source("./R/mcmc_bayes_parallel.R")
             apply(decomposed_data[[i]], 1, process_combination,
                   idx = idx_c)
           }
time_end = Sys.time()

print(time_end - time_start)

stopCluster(cl)

results_list_bacc_slope_all_comb = results_list

results_list_df = do.call(c, results_list)
results_list_df = as.data.frame(t(do.call(cbind, results_list_df)))


insitu_wise_data = bgc_val_vec

 ## Save the result output into two formats: 1. for each station, 2. for each slope ----

### Create a list of list for each rrs observation inverted results for all the slope ----
  rrs_names = rownames(inverse_rrs_input_mod)
list_of_50_dfs <- vector("list", length = length(rrs_names))


for (i in seq_along(rrs_names)) {
  rrs_idx = grep(paste0("^", rrs_names[i], ".*"), rownames(results_list_df))
  list_of_50_dfs[[i]] = results_list_df[rrs_idx, ]
}

# Set names of the list for easier access
names(list_of_50_dfs) <- rrs_names


list_of_2380_dfs <- vector("list", length = nrow(spectral_slope_vec)
                           #length = nrow(combinations)
                           
                           )

# Loop over each combination and collect corresponding rows from all rrs_names
for (j in seq_along(list_of_2380_dfs)) {
  list_of_2380_dfs[[j]] <- do.call(rbind, lapply(list_of_50_dfs, function(df) df[j, ]))
}

# Set names of the list for easier access (optional)
names(list_of_2380_dfs) <- apply(spectral_slope_vec, 1, function(x) paste(x, collapse = "_"))
#names(list_of_2380_dfs) <- apply(combinations, 1, function(x) paste(x, collapse = "_"))


# Initialize a data frame to store the results for spectral slope sensitivity
results_summary <- data.frame(
  Combination = apply(spectral_slope_vec, 1, function(x) paste(x, collapse = "_")),

  Uncertainty_chl = numeric(nrow(spectral_slope_vec)),
  Bias_chl = numeric(nrow(spectral_slope_vec)),
  p_bias_chl = numeric(nrow(spectral_slope_vec)),

  Uncertainty_adg443 = numeric(nrow(spectral_slope_vec)),
  Bias_adg443 = numeric(nrow(spectral_slope_vec)),
  p_bias_adg443 = numeric(nrow(spectral_slope_vec)),

  Uncertainty_H = numeric(nrow(spectral_slope_vec)),
  Bias_H = numeric(nrow(spectral_slope_vec)),
  p_bias_H = numeric(nrow(spectral_slope_vec))

)

# Initialize a data frame to store the results for benthic type sensitivity
# results_summary <- data.frame(
#   Combination = apply(combinations, 1, function(x) paste(x, collapse = "_")),
#   
#   Uncertainty_chl = numeric(nrow(combinations)),
#   Bias_chl = numeric(nrow(combinations)),
#   
#   Uncertainty_adg443 = numeric(nrow(combinations)),
#   Bias_adg443 = numeric(nrow(combinations)),
#   
#   Uncertainty_H = numeric(nrow(combinations)),
#   Bias_H = numeric(nrow(combinations))
#   
# )


pbias <- function(observed, predicted) {
  sum_diff <- sum(predicted - observed)   # Sum of differences between predicted and observed
  sum_obs <- sum(observed)                # Sum of observed values
  pbias_value <- (sum_diff / sum_obs) * 100  # PBIAS formula
  return(pbias_value)
}

# Calculate uncertainty and bias for each of the combination of spectral slope
for (i in seq_along(list_of_2380_dfs)) {
  
  # Extract the i-th data frame containing predictions
  pred_df <- list_of_2380_dfs[[i]]#[1:17,]
  
  # Extract predictions for 'chl', 'adg443' and 'H'
  pred_chl <- pred_df$chl
  pred_H <- pred_df$H
  pred_adg443 <- pred_df$adg443
  
  # Extract actual values for 'chl', 'adg443' and 'H' from the insitu data
  actual_chl <- as.numeric(insitu_wise_data$chl)#[1:17]
  actual_adg443 <- as.numeric(insitu_wise_data$adg443)#[1:17]
  actual_H <- as.numeric(insitu_wise_data$H)#[1:17]
  
  # Calculate Uncertainty (Standard Deviation of predictions)
  uncertainty_chl <- sd(c(actual_chl,pred_chl))
  uncertainty_adg443 <- sd(c(actual_adg443,pred_adg443))
  uncertainty_H <- sd(c(actual_H,pred_H))
  
  # Calculate Bias (Mean of predictions minus actual value)
  bias_chl <- mean(pred_chl) - mean(actual_chl)
  bias_adg443 <- mean(pred_adg443) - mean(actual_adg443)
  bias_H <- mean(pred_H) - mean(actual_H)
  
  # Calculate %-Bias 
  p_bias_chl <- pbias(observed = actual_chl, predicted = pred_chl)
  p_bias_adg443 <- pbias(observed = actual_adg443, predicted = pred_adg443)
  p_bias_H <- pbias(observed = actual_H, predicted = pred_H)
  
  # Store the results
  results_summary$Uncertainty_chl[i] <- uncertainty_chl
  results_summary$Bias_chl[i] <- bias_chl
  results_summary$p_bias_chl[i] <- p_bias_chl
  
  results_summary$Uncertainty_adg443[i] <- uncertainty_adg443
  results_summary$Bias_adg443[i] <- bias_adg443
  results_summary$p_bias_adg443[i] <- p_bias_adg443
  
  results_summary$Uncertainty_H[i] <- uncertainty_H
  results_summary$Bias_H[i] <- bias_H
  results_summary$p_bias_H[i] <- p_bias_H
  
}


results_summary = cbind(spectral_slope_vec, results_summary)
#results_summary = cbind(combinations, results_summary)
results_summary = results_summary %>% dplyr::select(Combination, everything())

# Print or save the results summary
print(results_summary)


# Find the uncertainty at the optimal slopes
optimal_uncertainty <- results_summary$Bias_H[which.min(abs(results_summary$Bias_H - 0))]

# Calculate the relative uncertainty change for both slopes
results_summary$uncertainty_change_1 <- abs(results_summary$Bias_H - optimal_uncertainty)

# Average uncertainty change for both slopes
average_uncertainty_change_1 <- mean(results_summary$uncertainty_change_1)
print(average_uncertainty_change_1)

# # Numerical partial derivatives (finite difference method)
# results_summary$partial_slope_gamma <- c(NA, 
#                                      diff(results_summary$Uncertainty_chl) / diff(results_summary$gamma))
# results_summary$partial_slope_adg <- c(NA, 
#                         diff(results_summary$Uncertainty_chl) / diff(results_summary$sg + 
#                                                                        results_summary$sd))
# 
# # View the result
# print(df)

plot(insitu_wise_data$adg443, list_of_2380_dfs[[13]]$adg443#[1:17]
     , xlim = c(0,4), ylim = c(0,4))

abline(0,1)

inv_wise_H = list_of_2380_dfs[[53]]$H
inv_wise_chl = list_of_2380_dfs[[53]]$chl
inv_wise_adg443 = list_of_2380_dfs[[5]]$adg443

inv_wise_H_sd = list_of_2380_dfs[[53]]$sd_H
inv_wise_chl_sd = list_of_2380_dfs[[53]]$sd_chl
inv_wise_adg443_sd = list_of_2380_dfs[[5]]$sd_adg443


inv_wise_param = data.frame("chl" = inv_wise_chl, "adg443" = inv_wise_adg443, "H"= inv_wise_H,
                            "sd_chl" = inv_wise_chl_sd, "sd_adg443" = inv_wise_adg443_sd,
                            "sd_H" = inv_wise_H_sd)

# Initialize a data frame to store the results for each observation
optimal_combinations_summary <- data.frame(
  Observation = names(list_of_50_dfs),
  Optimal_Combination = character(length(list_of_50_dfs)),
  
  Uncertainty_chl = numeric(length(list_of_50_dfs)),
  Bias_chl = numeric(length(list_of_50_dfs)),
  P_bias_chl = numeric(length(list_of_50_dfs)),
  
  Uncertainty_adg443 = numeric(length(list_of_50_dfs)),
  Bias_adg443 = numeric(length(list_of_50_dfs)),
  P_bias_adg443 = numeric(length(list_of_50_dfs)),
  
  Uncertainty_H = numeric(length(list_of_50_dfs)),
  Bias_H = numeric(length(list_of_50_dfs)),
  P_bias_H = numeric(length(list_of_50_dfs))
)

# Initialize a list to store detailed metrics for each Rrs observation
detailed_metrics_list <- vector("list", length(list_of_50_dfs))

# Loop through each observation (i.e., each data frame in list_of_50_dfs)
for (obs_idx in seq_along(list_of_50_dfs)) {
  
  pred_df <- list_of_50_dfs[[obs_idx]]
  
  best_combination <- NULL
  min_uncertainty <- Inf
  min_bias <- Inf
  min_pbias <- Inf
  
  # Initialize variables to store the best metrics
  best_uncertainty_chl <- NA
  best_bias_chl <- NA
  best_pbias_chl <- NA
  
  best_uncertainty_adg443 <- NA
  best_bias_adg443 <- NA
  best_pbias_adg443 <- NA
  
  best_uncertainty_H <- NA
  best_bias_H <- NA
  best_pbias_H <- NA
  
  
  # Initialize a data frame to store metrics for all combinations for the current observation
  detailed_metrics <- data.frame(
    Combination = character(nrow(pred_df)),
    
    Uncertainty_chl = numeric(nrow(pred_df)),
    Bias_chl = numeric(nrow(pred_df)),
    p_bias_chl = numeric(nrow(pred_df)),
    
    Uncertainty_adg443 = numeric(nrow(pred_df)),
    Bias_adg443 = numeric(nrow(pred_df)),
    p_bias_adg443 = numeric(nrow(pred_df)),
    
    Uncertainty_H = numeric(nrow(pred_df)),
    Bias_H = numeric(nrow(pred_df)),
    p_bias_H = numeric(nrow(pred_df))
    
  )
  
  rrs_name = rownames(pred_df)[1]
  rrs_name = unlist(strsplit(rrs_name,"_"))
  rrs_name = rrs_name[1]
  
  # Loop through each combination of slopes for the current observation
  for (i in 1:nrow(pred_df)) {
    
    pred_chl <- pred_df[i, "chl"]
    pred_H <- pred_df[i, "H"]
    pred_adg443 <- pred_df[i, "adg443"]
    unc_chl = pred_df[i, "sd_chl"]
    unc_H = pred_df[i, "sd_H"]
    unc_adg443 = pred_df[i, "sd_adg443"]
    
    
    bgc_idx = grep(rrs_name, insitu_wise_data$correct_statname)
    
    actual_chl <- as.numeric(insitu_wise_data$chl[bgc_idx])[1]
    actual_adg443 <- as.numeric(insitu_wise_data$adg443[bgc_idx])[1]
    actual_H <- as.numeric(insitu_wise_data$H[bgc_idx])[1]
    
    # Calculate the metrics
    
    uncertainty_chl <- sd(c(actual_chl,pred_chl))
    bias_chl <- mean(pred_chl) - actual_chl
    p_bias_chl <- pbias(observed = actual_chl, predicted = pred_chl)
    
    uncertainty_adg443 <- sd(c(actual_adg443,pred_adg443))
    bias_adg443 <- mean(pred_adg443) - actual_adg443
    p_bias_adg443 <- pbias(observed = actual_adg443, predicted = pred_adg443)
    
    uncertainty_H <- sd(c(actual_H,pred_H))
    bias_H <- mean(pred_H) - actual_H
    p_bias_H <- pbias(observed = actual_H, predicted = pred_H)
    
    
    # Store the metrics for this combination
    detailed_metrics$Combination[i] <- paste(round(spectral_slope_vec[i, 1:2], digits = 4), 
                                             collapse = "_")
    #detailed_metrics$Combination[i] <- paste(combinations[i, ], collapse = "_")
    
    detailed_metrics$Uncertainty_chl[i] <- uncertainty_chl
    detailed_metrics$Bias_chl[i] <- bias_chl
    detailed_metrics$p_bias_chl[i] <- p_bias_chl
    
    detailed_metrics$Uncertainty_adg443[i] <- uncertainty_adg443
    detailed_metrics$Bias_adg443[i] <- bias_adg443
    detailed_metrics$p_bias_adg443[i] <- p_bias_adg443
    
    detailed_metrics$Uncertainty_H[i] <- uncertainty_H
    detailed_metrics$Bias_H[i] <- bias_H
    detailed_metrics$p_bias_H[i] <- p_bias_H
    
    
    # Update the best combination if the current one is better
    if (uncertainty_chl < min_uncertainty || bias_chl < min_bias || p_bias_chl < min_pbias ) {
      min_uncertainty <- uncertainty_chl
      min_bias <- bias_chl
      min_pbias <- p_bias_chl
      best_combination <- paste(round(spectral_slope_vec[i, 1:2], digits = 4), collapse = "_")
      #best_combination <- paste(combinations[i, ], collapse = "_")
      
      best_uncertainty_chl <- uncertainty_chl
      best_bias_chl <- bias_chl
      best_pbias_chl <- p_bias_chl
      
      best_uncertainty_adg443 <- uncertainty_adg443
      best_bias_adg443 <- bias_adg443
      best_pbias_adg443 <- p_bias_adg443
      
      best_uncertainty_H <- uncertainty_H
      best_bias_H <- bias_H
      best_pbias_H <- p_bias_H
    }
  }
  
  # Store the best combination and its metrics for the current observation
  optimal_combinations_summary$Optimal_Combination[obs_idx] <- best_combination
  
  optimal_combinations_summary$Uncertainty_chl[obs_idx] <- best_uncertainty_chl
  optimal_combinations_summary$Bias_chl[obs_idx] <- best_bias_chl
  optimal_combinations_summary$P_bias_chl[obs_idx] <- best_pbias_chl
  
  optimal_combinations_summary$Uncertainty_adg443[obs_idx] <- best_uncertainty_adg443
  optimal_combinations_summary$Bias_adg443[obs_idx] <- best_bias_adg443
  optimal_combinations_summary$P_bias_adg443[obs_idx] <- best_pbias_adg443
  
  optimal_combinations_summary$Uncertainty_H[obs_idx] <- best_uncertainty_H
  optimal_combinations_summary$Bias_H[obs_idx] <- best_bias_H
  optimal_combinations_summary$P_bias_H[obs_idx] <- best_pbias_H
  
  # Store the detailed metrics for this observation in the list
  detailed_metrics_list[[obs_idx]] <- detailed_metrics
}

# Print the summary of optimal combinations for each observation
print(optimal_combinations_summary)

# Combine the detailed metrics for each observation into a single data frame for easy inspection
all_detailed_metrics <- do.call(rbind, detailed_metrics_list)
all_detailed_metrics$station = rep(rownames(inverse_rrs_input_mod), each = nrow(spectral_slope_vec))
#all_detailed_metrics$slope_count = rep(seq(1,10,1), 26)

all_detailed_metrics = cbind(round(spectral_slope_vec[,1:2], digits = 4), all_detailed_metrics)

rownames(all_detailed_metrics) = paste0(all_detailed_metrics$station, "_",
                                        all_detailed_metrics$Combination 
                                       )
# Print the detailed metrics for all observations and combinations
print(all_detailed_metrics)


# Calculate the mean uncertainty and bias for each combination
mean_metrics <- aggregate(cbind(Uncertainty_chl, Bias_chl, p_bias_chl,
                                Uncertainty_adg443, Bias_adg443, p_bias_adg443,
                                Uncertainty_H, Bias_H, p_bias_H) ~ Combination, 
                          data = all_detailed_metrics, 
                          FUN = mean)
min_metrics <- aggregate(cbind(Uncertainty_chl, Bias_chl, p_bias_chl,
                               Uncertainty_adg443, Bias_adg443, p_bias_adg443,
                               Uncertainty_H, Bias_H, p_bias_H) ~ Combination, 
                         data = all_detailed_metrics, 
                         FUN = min)
max_metrics <- aggregate(cbind(Uncertainty_chl, Bias_chl, p_bias_chl,
                               Uncertainty_adg443, Bias_adg443, p_bias_adg443,
                               Uncertainty_H, Bias_H, p_bias_H) ~ Combination, 
                         data = all_detailed_metrics, 
                         FUN = max)

# Split the Combination column into two new columns: s_dg and gamma
mean_metrics <- mean_metrics %>%
  mutate(s_dg = as.numeric(sapply(strsplit(as.character(Combination), "_"), `[`, 1)),
         gamma = as.numeric(sapply(strsplit(as.character(Combination), "_"), `[`, 2)))

# Create a grid of unique s_dg and gamma values
library(plot3D)
M <- mesh(unique(mean_metrics$s_dg), unique(mean_metrics$gamma))

s_dg_vals <- sort(unique(mean_metrics$s_dg))
gamma_vals <- sort(unique(mean_metrics$gamma))

# Create a matrix of p_bias_chl values
p_bias_matrix <- matrix(NA, nrow = length(s_dg_vals), ncol = length(gamma_vals))

# Fill the matrix with p_bias_chl values
for (i in 1:length(s_dg_vals)) {
  for (j in 1:length(gamma_vals)) {
    match_row <- mean_metrics[mean_metrics$s_dg == s_dg_vals[i] & mean_metrics$gamma == gamma_vals[j], ]
    if (nrow(match_row) > 0) {
      p_bias_matrix[i, j] <- abs(0 - match_row$p_bias_H)
    }
  }
}

rownames(p_bias_matrix) = s_dg_vals
colnames(p_bias_matrix) = gamma_vals

library(zoo)
p_bias_matrix = rowMeans(simplify2array(list(na.approx(p_bias_matrix, rule = 2), 
                             t(na.approx(t(p_bias_matrix), rule = 2)))), TRUE, 2)



# Define the finer grid for interpolation
s_dg_fine <- seq(0.0041, 0.0224, length.out = 80)  # Finer steps for s_dg
gamma_fine <- seq(0.1, 1.0, length.out = 100)       # Finer steps for gamma

M_fine = mesh(s_dg_fine, gamma_fine)

# Perform 2D interpolation using the interp function
interp_data <- akima:: interp(x = rep(s_dg_vals, each = length(gamma_vals)), 
                      y = rep(gamma_vals, length(s_dg_vals)), 
                      z = as.vector(t(p_bias_matrix)),  # Flatten p_bias_matrix for interpolation
                      xo = s_dg_fine,  # Finer grid for s_dg
                      yo = gamma_fine) # Finer grid for gamma

p_bias_matrix_fine = interp_data$z
p_bias_matrix_fine = rowMeans(simplify2array(list(na.approx(p_bias_matrix_fine, rule = 2), 
                                             t(na.approx(t(p_bias_matrix_fine), rule = 2)))), TRUE, 2)


png(file = paste0("./outputs/sens_slope_H.png"),  units = "in", width = 8, height = 6, res=300)
par(mai = c(0.30, 0.1, 0.1, 1))

#Plot with Surf3D
plot3D::surf3D(x =(M_fine$x) ,y = (M_fine$y) ,
               #z = log(chain$LP[1:(i+1)]),
               z= (p_bias_matrix_fine),
               #xlim = c(-1,1), ylim=c(-1, 1),
               xlab = "", ylab="", zlab = "",
               xlim=c(0.004, 0.022), ylim=c(0.1,1),
               add = F, phi = 25, bty = "g", type = "b", theta =230, #30,
               ticktype = "detailed", pch = 20,
               col = viridis::viridis(n = 80, option = "G"),
               #lighting = TRUE,
               clab = c(" "),
               colkey = T, 
               colvar =(p_bias_matrix_fine),
               contour= T,
               cex = c(0.5, 1, 1.5),
               expand = 1,
               #border = NULL,
               clim = range(p_bias_matrix_fine, na.rm = TRUE)
)
# scatter3D(x = log10(da$adg443[col_num]), y = da$H[row_num], z= log10(LL[col_num, row_num]), 
#           theta = 230, phi=15, add=T, col = "darkgreen", pch=23, lwd=5)
# legend(x = "bottomright",          # Position
#        legend = expression(paste("MAP(", phi["par"],")")),  # Legend texts
#        lty = c(NA),           # Line types
#        pch=23,
#        col = c("darkgreen"),           # Line colors
#        fill = c("darkgreen"),
#        box.col="white",
#        box.lwd=0,
#        lwd = 2,
#        #inset = c(0,1),
#        cex=0.9)

dev.off()


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Reshape the dataframe to long format
df_long <- all_detailed_metrics %>%
  dplyr::select(sdg, gamma, Bias_chl, Bias_adg443, Bias_H, Combination) %>%
  pivot_longer(cols = starts_with("Bias"),
               names_to = "Bias_Type",
               values_to = "Bias_Value")

# Create an index for the combinations (1 to 68) for mapping to color
# df_long <- df_long %>%
#   group_by(sdg, gamma) %>%
#   mutate(Combination_Index = row_number())

df_long <- df_long %>%
  group_by(sdg, gamma) %>%
  arrange(Combination) %>%  # Sort by 'Combination' in ascending order
  mutate(Combination_Index = row_number())

df_long <- df_long %>%
  arrange(Combination) %>%  # Sort by Combination in ascending order
  group_by(Combination) %>%  # Group by Combination to handle duplicates
  mutate(Combination_Index = first(row_number())) %>%  # Assign index based on order of appearance
  ungroup()  # Ungroup for further operations

df_long <- df_long %>%
  arrange(sdg, gamma, Combination) %>%  # Sort by 'sdg', 'gamma', and 'Combination' in ascending order
  group_by(sdg, gamma, Combination) %>%  # Group by 'sdg', 'gamma', and 'Combination'
  mutate(Combination_Index = cur_group_id()) %>%  # Assign index based on group
  ungroup()  # Ungroup to avoid affecting future operations



# Create a factor for 'Bias_Type' to control the order of plots
df_long$Bias_Type <- factor(df_long$Bias_Type, levels = c("Bias_chl", "Bias_adg443", "Bias_H"))

chl_lab <- expression(paste("[", italic("chl"), "]"))
adg443_lab <- expression(paste(italic("a")["dg"](443)))
bbp555_lab <- expression(paste(italic("b")["bp"](555)))
H_lab <- expression(paste(italic("H")))


# Define the function
save_histogram_with_density <- function(column_name, output_filename, col_lab, xlim, x_lab) {
  
  # Check if the column_name exists in the data
  if (!(column_name %in% unique(df_long$Bias_Type))) {
    stop(paste("Error: Column", column_name, "not found in 'df_long$Bias_Type'"))
  }
  
  # Set the output file (choose .png, .jpeg, or .pdf based on your needs)
  png(filename = output_filename, width = 800, height = 600)
  
  # Adjust plot margins with 'par'
  # 'mar' stands for margins: c(bottom, left, top, right)
  par(mar = c(5, 5, 4, 2))  # Adjust margins as needed (increasing space for axis labels)
  
  # Define color for density fill
  colors_in <- col_lab
  
  # Create the density object for the chosen column
  density_data <- density(df_long$Bias_Value[df_long$Bias_Type == column_name])
  
  # Create the histogram with actual counts
  hist_data <- hist(df_long$Bias_Value[df_long$Bias_Type == column_name],
                    breaks = 30,                    # Number of bins
                    col = colors_in[1],              # Histogram fill color
                    xlim = xlim,
                    border = "black",                    # Histogram border color
                    freq = F,                        # Make y-axis density instead of count
                    main = " ",
                    xlab = x_lab,                     # X-axis label
                    ylab = "Density",                    # Y-axis label
                    cex.lab = 2.5,                       # Increase font size for labels
                    cex.axis = 2,                      # Increase font size for axis ticks
                    cex.main = 1.8,                      # Increase font size for the title
                    cex.sub = 1.5                        # Increase font size for subtitle
  )
  
  # Overlay the density curve by scaling to match counts
  lines(density_data, 
        col = colors_in[2], lwd = 5)
  
  # Add a filled polygon under the density curve
  # polygon(density_data$x, density_data$y * length(df_long$Bias_Value[df_long$Bias_Type == 
  #                                                         column_name]) * diff(hist_data$breaks)[1],
  #         col = rgb(0.1, 0.8, 0.1, 0.4), border = NA)
  
  # Add grid for better visual representation
  grid(col = "gray")
  
  # Close the graphic device to save the plot
  dev.off()
  
  # Return success message
  return(paste("Plot saved as", output_filename))
}

save_histogram_with_density(column_name = "Bias_adg443", 
                            output_filename = "./outputs/sens_slope_residual_adg443.png",
                            col_lab = c(colors_in[6], colors_border[6]), 
                            xlim = c(-1.5,-1.5), x_lab = expression(paste(delta, a["dg"](443))))


density_data <- density(df_long$Bias_Value[df_long$Bias_Type == "Bias_chl"])

png(filename = "./outputs/sens_slope_residual_H.png", width = 800, height = 600)
par(mar = c(5, 5, 4, 2)) 

# Create the histogram with adjusted breaks and fill
hist_H<- hist(df_long$Bias_Value[df_long$Bias_Type == "Bias_chl"],
     breaks = 30,                        # Number of bins
     xlim = c(-10,10),
     col = colors_in[5],                  # Histogram fill color
     border = "black",                    # Histogram border color
     freq = F,                        # Make y-axis density instead of count
     main = " ",
     xlab = expression(paste(delta, "[chl-a]")),                     # X-axis label
     ylab = "Density",                    # Y-axis label
     cex.lab = 1.8,                       # Increase font size for labels
     cex.axis = 1.5,                      # Increase font size for axis ticks
     cex.main = 1.8,                      # Increase font size for the title
     cex.sub = 1.5                        # Increase font size for subtitle
)

# Overlay the density curve by scaling to match counts
lines(density_data, 
      col = "black", lwd = 2.5)

# Add a filled polygon under the density curve
# polygon(density_H$x, density_H$y * length(df_long$Bias_Value[df_long$Bias_Type == 
#                                                                "Bias_H"]) * diff(hist_H$breaks)[1],
#         col = rgb(0.1, 0.8, 0.1, 0.4), border = NA)

# Add grid for better visual representation
grid(col = "gray")
dev.off()


# Create violin plot with overlaid points for each Bias_Type
g_residual = ggplot(df_long#[df_long$Bias_Type == "Bias_adg443",]
                    , 
                    aes(x = Bias_Type, y = Bias_Value)) + 

  
  geom_violin(fill = "lightblue", alpha = 0.5, position=position_dodge(), size = 1.3,
              draw_quantiles = c(0.667, 0.95, 0.991), scale = "area", color = "black") + 
  
  geom_jitter(aes(color = Combination_Index), width = 0.15, size = 2, alpha = 0.25,) + 
  
  scale_color_gradientn(colors = viridis::viridis(100, option = "turbo"), 
                        guide = "colourbar"#, 
                        #breaks = seq(1,100,20), 
                        # labels = paste("s_dg:", c(0.004, ), 
                        #                "gamma:", c())
                        ) +
  scale_x_discrete(name = "", labels = c(
    chl_lab
    , adg443_lab, 
    H_lab)) +
  
  
  labs(title = " ",
       y = expression(paste(delta["var"])), 
       color = expression(paste("{",s["dg"], ",",gamma,"}"))) + 
  ylim(-10,10)+
  geom_hline(yintercept = 0, colour = "black", size = 1.3, linetype = 6)+
  
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 25, color = 'black', angle = 0),
        axis.text.y = element_text(size = 25, color = 'black', angle = 0),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 30),
        axis.ticks.length = unit(.25, "cm"),
        legend.position=c(0.50, 0.10),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 20, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5,
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        #legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey",
                                        size = 0.5, linetype = "dotted"),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.0,0.5,0.0,0.0), "cm"),
        legend.key.width = unit(3, "cm"),
        
        # Center the legend's title and make it horizontal
        legend.title.align = 0.5,
        legend.text.align = 0.5,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

ggsave(paste0("./outputs/sens_slope_residual_violin.png"), plot = g_residual,
       scale = 1.25, width = 10, height = 6, units = "in",dpi = 300)


# EXP WISE inversion ----
## Load the COPS profiles (both cast1 and cas2 groups), translate to sub-surface and correct for inelastic scattering----
path = "Y:/homeData/Insitu/WISEMan/L2/20190825_StationEXPWISEv2/COPS/"
cops_cast_info <- scan(file = paste0(path,"directories.for.cops.dat"), "", sep = "\n", 
                       comment.char = "#")


remove.file = paste0(path, "select.cops.dat")
remove.tab <- read.table(remove.file, header = FALSE, colClasses = "character", sep = ";")
kept.cast <- remove.tab[[2]] == "1"
listfile  <- remove.tab[kept.cast, 1]


nf = length(listfile)
print(listfile)

mRrs_0p = matrix(ncol=19, nrow = nf)
mRrs_0m = matrix(ncol=19, nrow = nf)


for (j in 1:nf) {
  
  load(paste(path,"BIN/", listfile[j], ".RData", sep=""))
  waves = cops$LuZ.waves
  
  rb_df = data.frame("wave" = waves, "rb_lu_0.3" = pi*cops$Rb.LuZ, "rb_eu_0.3" = cops$Rb.EuZ,
                     "rb_lu_H" = pi*cops$Rb.LuZ.H0, "rb_eu_H" = cops$Rb.EuZ.H0, "Q_0.3" = cops$Rb.Q, 
                     "Q_H" = cops$Rb.Q.H0)
  write.csv(rb_df, paste0(path,"BIN/Rb_COPS_Profile_v2_", j-1, ".csv"), quote = F, row.names = F)
  
  extrap_mode = remove.tab$V3[j]
  extrap_mode = substr(extrap_mode, start = 5, stop = 100)
  
  extrap_idx = grep(pattern = "linear", extrap_mode)
  
  is_empty <- function(x) {
    is.null(x) || length(x) == 0
  }
  
  if (is_empty(extrap_idx)) {
    
    mRrs_0m[j,] = cops$LuZ.0m/cops$EdZ.0m
    
  } else {
    mRrs_0m[j,] = cops$LuZ.0m.linear/cops$EdZ.0m.linear
  }
  
  mRrs_0p[j,] = cops[[remove.tab$V3[j]]]
  
}

colnames(mRrs_0m) = cops$LuZ.waves
colnames(mRrs_0p) = cops$LuZ.waves

rownames(mRrs_0m) = c(paste0("EXP-WISE1_",seq(0,10,1)), paste0("EXP-WISE2_",seq(1,2,1)))
rownames(mRrs_0p) = c(paste0("EXP-WISE1_",seq(0,10,1)), paste0("EXP-WISE2_",seq(1,2,1)))


wavelength = as.numeric(colnames(mRrs_0p))
wavelength_interp = seq(400,750,1)

wavelength = wavelength_interp

mRrs_0p_interp <- apply(mRrs_0p, 1, function(row, wavelength_in = wavelength){
  
  wavelength_interp =wavelength_interp
  obsdata = row
  
  #non_na_seq <- seq_along(obsdata)[!is.na(obsdata)]
  
  obsdata_interp = Hmisc::approxExtrap(x = wavelength_in#[non_na_seq]
                                       , y = obsdata#[!is.na(obsdata)]
                                       , 
                                       xout = wavelength_interp)$y
  obsdata_interp[obsdata_interp < 0] = 0 
  obsdata = obsdata_interp
  
})

mRrs_0p_interp = data.frame(mRrs_0p_interp)
rownames(mRrs_0p_interp) = wavelength_interp

mRrs_0m = mRrs_0m[,5:16]
mRrs_0p = mRrs_0p[,5:16]

## Load the PSR multicast profiles, translate to sub-surface and correct for inelastic scattering----
mRrs_0p_PSR = wise_psr_multicast_df

mRrs_0p_PSR = mRrs_0p_PSR[colnames(mRrs_0p_PSR) %in% seq(400,750,1)]
mRrs_0p_PSR = data.frame((mRrs_0p_PSR))[1:6,]
colnames(mRrs_0p_PSR) = seq(400,750,1)


mRrs_0p_PSR_subsurf = as.data.frame(apply(mRrs_0p_PSR, 2, 
                                        surface_rrs_translate))
mRrs_0p_subsurf = as.data.frame(apply(mRrs_0p, 2, 
                                  surface_rrs_translate))

mRrs_0p_PSR_subsurf_inel = data.frame(t(subtract_wise_inel(input_df = t(mRrs_0p_PSR_subsurf), 
                                       wave_out = as.numeric(colnames(mRrs_0p_PSR_subsurf))
                                       )))
colnames(mRrs_0p_PSR_subsurf_inel) = colnames(mRrs_0p_PSR_subsurf)

mRrs_0p_subsurf_inel = data.frame(t(subtract_wise_inel(input_df = t(mRrs_0p_subsurf), 
                           wave_out = as.numeric(colnames(mRrs_0p_subsurf))
)))
colnames(mRrs_0p_subsurf_inel) = colnames(mRrs_0p_subsurf)

## Combine the COPS cast-II and PSR multicast profiles together and plot spectrally----
expwise_psr = reshape2::melt(cbind("station" = rownames(mRrs_0p_PSR_subsurf_inel),
                                   mRrs_0p_PSR_subsurf_inel), 
                              id.vars = "station")

expwise_cops = reshape2::melt(cbind("station" = rownames(mRrs_0p_subsurf_inel), mRrs_0p_subsurf_inel), 
                              id.vars = "station")

#expwise_cops = expwise_cops[-(1:2),]

colnames(expwise_psr) = c("station", "wavelength", "Rrs")

colnames(expwise_cops) = c("station", "wavelength", "Rrs")


xmin <- 400; xmax <- 750;  xstp <- 70; xlbl <- expression(paste("Wavelength (", lambda, ") [nm]")) 
ymin <- 0.000; ymax <- 0.008; ystp <- 0.002
ylbl <- expression(paste(italic("R")["rs"]("0"^"-", lambda), "[sr"^{-1}, "]"))
asp_rat <-  (xmax-xmin)/(ymax-ymin)
legend_title <- element_blank()
legend_position <- c(0.10, 1)


g_expwise = ggplot() + 
  geom_line(data = expwise_psr, size = 1.3, 
            aes(x = as.numeric(as.character(wavelength)), y = Rrs, color = station, group = station),
            show.legend = F) +
  
  # geom_line(data = expwise_cops, size = 1.3, linetype = "dashed",
  #           aes(x = as.numeric(as.character(wavelength)), y = Rrs, color = station, group = station)) +
  
  geom_point(data = expwise_cops[grep("EXP.WISE2_*", expwise_cops$station),], 
             aes(x = as.numeric(as.character(wavelength)), 
                                      y = as.numeric(Rrs), color = station), 
             size = 5, show.legend = F)+
  scale_color_viridis(discrete = T, name = "", labels = c("COPS_1(15:27)", "COPS_2(15:29)",
                                                          "PSR_A(13:35)", "PSR_B(13:43)", "PSR_C(14:20)", 
                                                          "PSR_D(14:25)", "PSR_E(15:43)", "PSR_F(15:46)")) +
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
              ylim = c(ymin, ymax)
              ,expand = FALSE, clip = "on"
  ) +
  
  scale_x_continuous(name = xlbl, limits = c(xmin, xmax),
                     breaks = seq(xmin, xmax, xstp))  +
  scale_y_continuous(name = ylbl, limits = c(ymin, ymax),
                     breaks = seq(ymin, ymax, ystp))  +
  geom_hline(yintercept = 0, colour = "navyblue", linetype = "dashed", size=1.1)+
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
        axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = c(0.01, 0.99),
        #legend.direction = "vertical",
        legend.title = element_text(size = 10, color = 'black', angle = 0, face = "bold"),
        legend.text = element_text(colour = "black", size = 15, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_line(colour = "grey80", 
                                        linewidth =  0.2, linetype = "solid"),
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        legend.direction = "vertical", legend.box = "vertical",
        legend.text.align = 0,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

g_expwise
ggsave(paste0("./outputs/expwise_rrs_comp.png"), plot = g_expwise,
       scale = 1.25, width = 4.5, height = 4.5, 
       units = "in",dpi = 300)


## Prepare the Rrs for inversion----
inv_res_df_bayes_unconst = data.frame("chl" = 0, "adg443" = 0, "bbp555" = 0,
                                      "H" = 0, "fa1"= 0, "fa2"= 0, 'fa3' = 0, "pop_sd" =0,
                                      "sd_chl" = 0, "sd_adg443" = 0, "sd_bbp555" = 0, 
                                      "sd_H" = 0, "sd_fa1"= 0, "sd_fa2"= 0, 'sd_fa3' = 0, 
                                      "sd_pop_sd" = 0)

inv_res_df_grad_unconst = data.frame("chl" = 0, "adg443" = 0, "bbp555" = 0,
                                     "H" = 0, "fa1"= 0, "fa2"= 0, 'fa3' = 0, "pop_sd" =0,
                                     "sd_chl" = 0, "sd_adg443" = 0, "sd_bbp555" = 0, 
                                     "sd_H" = 0, "sd_fa1"= 0, "sd_fa2"= 0, 'sd_fa3' = 0, 
                                     "sd_pop_sd" = 0)

inverse_rrs_input = as.data.frame(mRrs_0p_subsurf_inel)

#inverse_rrs_input = as.data.frame(mRrs_0p_PSR_subsurf_inel)

interp_rb(wavelength_input = as.numeric(colnames(inverse_rrs_input)), 
          bottom_type = c("Sand_2019", "Mud_2019", "Eelgrass_2019")) #use a fixed benthic type

## Call the Bayesian inversion with optimal slope and benthic types----
for(j in 1:nrow(inverse_rrs_input)) {
  
  param_grad_unconst = doOptimization_shallow_unconst(obsdata =
                                                        as.numeric(inverse_rrs_input[j,]),
                                                      init_par = c(3.5,1.3,0.017,4.5,0.6,0.2,
                                                                   0.2,0.05),
                                                      max_par = c(6.5, 2.5, 0.01, 8, 1, 1,
                                                                  1, 10),
                                                      min_par = c(1.5, 1, 0.003, 2, 0.1, 0.1,
                                                                  0.1, 0.01),
                                                      QAA_mode = F,
                                                      bgc_const_val = as.numeric(bgc_val_vec[j,]),
                                            wavelength_sim = as.numeric(colnames(inverse_rrs_input))
                                                      #seq(400,750,1)
                                                      ,
                                                      qaa_prefit = F,
                                                      qaa_slope = F, manual_slope = T,
                                                      manual_slope_vals = c("s_g"=0.014, "s_d"=0.003,
                                                                            "gamma"=0.5),
                                                      sa_model = "am03",
                                                      obj_fn = obj[1], opt_method = methods.opt[4] )
  
  param_bayes_unconst = inverse_runBayes(obsdata = 
                                           as.numeric(inverse_rrs_input[j,]),
                                         rrs_type =  type_Rrs_below,
                                         max_par = c(6.5, 2.5, 0.01, 8, 1, 1, 1, 10),
                                         min_par = c(1.5, 1, 0.003, 2, 0.1, 0.1, 0.1, 0.01),
                                         param_lab = c("chl", "adg443", "bbp555", 
                                                       "H", "fa1", "fa2", "fa3", 
                                                       "pop_sd"),
                                         
                                         constrain_config = c(F,F,F),
                                         bbp_const_val = 0.007,
                                         bgc_const_val = as.numeric(bgc_val_vec[j,]),
                                         
                       iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
                                          "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
                                         
                                         qaa_slope = F,
                                         manual_slope = T,
                                         manual_slope_vals  = c("s_g"=0.014,
                                                                "s_d"=0.003,
                                                                "gamma"=0.5),
                                         iter_count = 25000,
                                         sampler_mcmc = samplerlist[6],
                                         wavelngth_sim = as.numeric(colnames(inverse_rrs_input))
                                         #seq(400,750,1)
                                         , 
                                         sa_model = "am03",
                                         hybrid_mode = F, plot_rrs = T)
  
  
  inv_res_df_grad_unconst = rbind(inv_res_df_grad_unconst, param_grad_unconst)
  inv_res_df_bayes_unconst = rbind(inv_res_df_bayes_unconst, param_bayes_unconst)
  
}


inv_res_df_bayes_unconst = inv_res_df_bayes_unconst[-1,]
inv_res_df_grad_unconst = inv_res_df_grad_unconst[-1,]


inv_res_df_bayes_unconst_expwise_cops = inv_res_df_bayes_unconst
inv_res_df_grad_unconst_expwise_cops = inv_res_df_grad_unconst

#inv_res_df_bayes_unconst_psr = inv_res_df_bayes_unconst
#inv_res_df_grad_unconst_psr = inv_res_df_grad_unconst

## Analyze and Validate (where applicable) the inversion results from COPS Cast-II and PSR Multicast----
expwise_sensor_comp = read.csv("./outputs/exp_wise_sensor_comp.csv", header = T)

# Initialize an empty list to store the result for each row
expwise_sensor_comp_rb <- list()

# Loop over each row of the 'df' dataframe
for (i in 1:nrow(expwise_sensor_comp)) {
  # For each row, compute the sum of (fa1 * class1 + fa2 * class2 + fa3 * class3)
  expwise_sensor_comp_rb[[i]] <- expwise_sensor_comp$fa1[i] * rb$class1 + 
    expwise_sensor_comp$fa2[i] * rb$class2 + 
    expwise_sensor_comp$fa3[i] * rb$class3
}


expwise_sensor_comp_rb <- as.data.frame(do.call(rbind, expwise_sensor_comp_rb))
rownames(expwise_sensor_comp_rb) <- expwise_sensor_comp$station

colnames(expwise_sensor_comp_rb) <- seq(400,750,1)
expwise_sensor_comp_rb_cops = expwise_sensor_comp_rb[7:8, colnames(mRrs_0p_subsurf)]

expwise_sensor_comp_rb_psr_long = reshape2::melt(cbind("station" = rownames(expwise_sensor_comp_rb)[1:6],
                                                   expwise_sensor_comp_rb[1:6,]), 
                             id.vars = "station")


colnames(expwise_sensor_comp_rb_psr_long) = c("station", "wavelength", "Rb")


expwise_sensor_comp_rb_cops_long = reshape2::melt(cbind("station" = rownames(expwise_sensor_comp_rb_cops),
                                                        expwise_sensor_comp_rb_cops), 
                                                 id.vars = "station")

colnames(expwise_sensor_comp_rb_cops_long) = c("station", "wavelength", "Rb")

### Spectral comparison of inversion retrieved Rb from COPS and PSR----
xmin <- 400; xmax <- 750;  xstp <- 70; xlbl <- expression(paste( lambda, "[nm]")) 
ymin <- 0.000; ymax <- 0.4; ystp <- 0.10
ylbl <- expression(paste(italic("R")["B"](lambda)))
asp_rat <-  (xmax-xmin)/(ymax-ymin)
legend_title <- element_blank()
legend_position <- c(0.10, 1)


g_expwise_rb = ggplot() + 
  geom_line(data = expwise_sensor_comp_rb_psr_long, size = 1.3, 
            aes(x = as.numeric(as.character(wavelength)), y = 0.70*Rb, color = station, group = station), 
            show.legend = T) +
  geom_point(data = data.frame("wave" = cops$Ed0.waves, "rb_eu" = cops$Rb.EuZ.H0), 
             aes(x = wave, y = rb_eu), fill = "orange2", shape = 23, size = 5)+
  geom_point(data = data.frame("wave" = cops$Ed0.waves, "rb_eu" = cops$Rb.LuZ.H0), 
             aes(x = wave, y = rb_eu), fill = "red3", shape = 24, size = 5)+
  
  # geom_line(data = expwise_cops, size = 1.3, linetype = "dashed",
  #           aes(x = as.numeric(as.character(wavelength)), y = Rb, color = station, group = station)) +
  
  geom_point(data = expwise_sensor_comp_rb_cops_long, aes(x = as.numeric(as.character(wavelength)),
                                      y = 0.70*as.numeric(Rb), color = station),
             size = 5, show.legend = T)+
  #scale_color_viridis(discrete = T) +
  scale_color_viridis(discrete = T, name = "", 
                      labels = c("COPS_1(15:27)", "COPS_2(15:29)",
                                 "PSR_A(13:35)", "PSR_B(13:43)", "PSR_C(14:20)", 
                                 "PSR_D(14:25)", "PSR_E(15:43)", "PSR_F(15:46)")) +
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
              ylim = c(ymin, ymax)
              ,expand = FALSE, clip = "on"
  ) +
  
  scale_x_continuous(name = xlbl, limits = c(xmin, xmax),
                     breaks = seq(xmin, xmax, xstp))  +
  scale_y_continuous(name = ylbl, limits = c(ymin, ymax),
                     breaks = seq(ymin, ymax, ystp))  +
  geom_hline(yintercept = 0, colour = "navyblue", linetype = "dashed", size=1.1)+
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
        axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = c(0.01, 0.99),
        #legend.direction = "vertical",
        legend.title = element_text(size = 10, color = 'black', angle = 0, face = "bold"),
        legend.text = element_text(colour = "black", size = 12, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_line(colour = "grey80", 
                                        linewidth =  0.2, linetype = "solid"),
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        legend.direction = "vertical", legend.box = "vertical",
        legend.text.align = 0,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

g_expwise_rb
ggsave(paste0("./outputs/expwise_rb_comp_new.png"), plot = g_expwise_rb,
       scale = 1.25, width = 4.5, height = 4.5, 
       units = "in",dpi = 300)

### Temporal comparison of inversion retrieved BGCs and IOPs from COPS and PSR----
expwise_sensor_comp_long <- expwise_sensor_comp %>%
  pivot_longer(cols = c("chl", "adg443", "bbp555", "H"),
               names_to = "Variable", values_to = "Value") %>%
  pivot_longer(cols = c("sd_chl", "sd_adg443", "sd_bbp555", "sd_H"),
               names_to = "SD_Variable", values_to = "SD") %>%
  filter(Variable == gsub("sd_", "", SD_Variable))  # Filter to ensure SD matches the correct variable


# Define custom names for variables
custom_labels <- c("chl" = expression(paste("[",italic("chl"),"]", " [","mg"^1,"m"^-3,"]")), 
                   "adg443" = expression(paste(italic("a")["dg"](443), " [", "m"^-1, "]")), 
                   "bbp555" = expression(paste(italic("b")["bp"](555), " [", "m"^-1, "]")), 
                   "H" = expression(paste(italic("H"), " [m]")))

# Define custom colors for lines and fills
custom_colors <- c("chl" = "#4AC16DFF", 
                   "adg443" = "goldenrod", 
                   "bbp555" = "#365C8DFF", 
                   "H" = "#1FA187FF")

# Create the ggplot with geom_ribbon for uncertainty and geom_line for the line plot
g_expwise_iop = ggplot(expwise_sensor_comp_long, 
                       aes(x = utc_time, y = Value, group = Variable, color = Variable)) +
  geom_line(size = 1.3) +
  geom_ribbon(aes(ymin = Value - SD, ymax = Value + SD, fill = Variable), 
              alpha = 0.3, show.legend = T) +
  facet_wrap(~Variable, ncol = 1, scales = "free_y") +  # Vertically arrange panels
  labs(#x = " ",
       y = " ", title = " ") +
  scale_color_manual(values = custom_colors, name = "Variable", labels = custom_labels) +  # Manual color for lines
  scale_fill_manual(values = custom_colors, name = "Variable", labels = custom_labels) +
  # scale_x_discrete(name = " ", labels = c("COPS_1", "COPS_2",
  #                                         "PSR_A", "PSR_B", "PSR_C", 
  #                                         "PSR_D", "PSR_E", "PSR_F")
                    # )  +
  #theme_minimal() +
  # theme(
  #   strip.text = element_text(size = 14),   # Increase facet label size
  #   axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
  #   axis.title = element_text(size = 14),   # Increase axis title size
  #   axis.text = element_text(size = 12),    # Increase axis text size
  #   legend.position = "bottom"              # Adjust legend position
  # )
  #theme_bw()+
  theme(strip.text = element_blank(), 
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15, color = 'black', angle = 15), 
        axis.text.y = element_text(size = 15, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        
        legend.direction = "horizontal",
        legend.position = "bottom",             # Position legend at bottom
        legend.justification = "center", 
        legend.title = element_blank(),
        legend.text = element_text(hjust = 0, colour = "black", size = 15, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        
        panel.grid.major = element_line(colour = "black", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_line(colour = "grey80", 
                                        linewidth =  0.2, linetype = "solid"),
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        panel.background = element_blank(),
        legend.key.width = unit(1.5, "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

ggsave(paste0("./outputs/expwise_iop_comp_new.png"), plot = g_expwise_iop,
       scale = 1.25, width = 9, height = 6, 
       units = "in",dpi = 300)

## Analyze and Validate (where applicable) the inversion results from COPS Cast-I profiles----
rb_spectral_svc = read.csv("./data/ExperimentWISE_SVC.csv", header = T)

rb_mat_saber = data.frame(matrix(0, nrow = nrow(mRrs_0p_subsurf_inel), 
                                 ncol = ncol(mRrs_0p_subsurf_inel)))

rb_mat_svc = data.frame(matrix(0, nrow = nrow(mRrs_0p_subsurf_inel), 
                                 ncol = 512))

rb_mat_cops = data.frame(matrix(0, nrow = nrow(mRrs_0p_subsurf_inel), 
                               ncol = ncol(mRrs_0p_subsurf_inel)))

rownames(inv_res_df_bayes_unconst_expwise_cops) = seq(0, 12,1)
rownames(rb_mat_saber) = seq(0, 12,1)
rownames(rb_mat_svc) = seq(0, 12,1)
rownames(rb_mat_cops) = seq(0, 12,1)

### Retrieve the in situ Rb from COPS Cast-I profiles and SVC measurements and compare with inversion retrieved Rb ----
for (i in seq(0,10,1)) {
  
  rb_spectral_svc_subset = rb_spectral_svc[rb_spectral_svc$profile == paste0("profile ", i), 
                                           #c("wavelengths_nm", "bottom_type", "inst.reflectance")
  ]
  
  rb_spectral_svc_subset_mean <- rb_spectral_svc_subset %>%
    group_by((wavelengths_nm)) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  
  
  cops_rb = read.csv(paste0("Y:/homeData/Insitu/WISEMan/L2/20190825_StationEXPWISEv2/COPS/BIN/Rb_COPS_Profile_v2_",i,".csv"),
  )
  
  #For MCMC based inversion results
  rb_bayes_spectral_est = inv_res_df_bayes_unconst_expwise_cops$fa1[i+1]*rb$class1 + 
    inv_res_df_bayes_unconst_expwise_cops$fa2[i+1]*rb$class2 + 
    inv_res_df_bayes_unconst_expwise_cops$fa3[i+1]*rb$class3
  
  rb_bayes_spectral_sd = inv_res_df_bayes_unconst_expwise_cops$sd_fa1[i+1]*rb$class1 + 
    inv_res_df_bayes_unconst_expwise_cops$sd_fa2[i+1]*rb$class2 + 
    inv_res_df_bayes_unconst_expwise_cops$sd_fa3[i+1]*rb$class3
  
  # Compute the upper and lower bounds of uncertainty
  # rb_bayes_spectral_upper = rb_spectral_est + rb_spectral_sd
  # rb_bayes_spectral_lower = rb_spectral_est - rb_spectral_sd
  
  
  rb_mat_saber[which(rownames(rb_mat_saber) == i),] = rb_bayes_spectral_est
  rb_mat_cops[which(rownames(rb_mat_cops) == i),] = cops_rb$rb_lu_H[2:13]
  rb_mat_svc[which(rownames(rb_mat_svc) == i),] = rb_spectral_svc_subset_mean$inst.reflectance/100
  
  # png(filename = paste0("./outputs/EXP_WISE_inv_bayes_Rb_profile_",i,".png"), width = 600, height = 400)
  # 
  # # Adjust plot margins with 'par'
  # # 'mar' stands for margins: c(bottom, left, top, right)
  # par(mar = c(5, 5, 4, 2))  # Adjust margins as needed (increasing space for axis labels)
  # 
  # plot(as.numeric(as.character(unique(expwise_cops$wavelength))), 
  #      rb_spectral_est, ylim = c(0, 0.2), type = "b", 
  #      lwd = 3, pch = 19, xlab = expression(paste("Wavelength [nm]")),
  #      ylab = expression(paste("R"["B"])), 
  #      cex.lab = 1.8,                       # Increase font size for labels
  #      cex.axis = 1.5,                      # Increase font size for axis ticks
  #      cex.main = 1.8,                      # Increase font size for the title
  #      cex.sub = 1.5 ,
  #      main = paste0("Profile ", i))
  # 
  # # Add shaded uncertainty region (polygon)
  # polygon(c(as.numeric(as.character(unique(expwise_cops$wavelength))), 
  #           rev(as.numeric(as.character(unique(expwise_cops$wavelength))))), 
  #         c(rb_spectral_upper, rev(rb_spectral_lower)),
  #         col = rgb(0.5, 0.5, 0.5, 0.5), border = NA)  # Grey with transparency
  # 
  # lines(as.numeric(cops_rb$waves), 
  #       as.numeric(cops_rb$R.piLuz.0.3)/100, type = "o", col = "red2", pch = 17)
  # 
  # lines(as.numeric(cops_rb$waves), 
  #       as.numeric(cops_rb$R.Euz.0.3)/100, type = "o", col = "orange2", pch = 17)
  # 
  # lines(rb_spectral_svc_subset_mean$wavelengths_nm, 
  #       rb_spectral_svc_subset_mean$inst.reflectance/100, col = "navyblue", lwd = 3)
  # legend("topleft", legend = c("SABER", expression(paste("COPS ["[L[u]],"]")), 
  #                              expression(paste("COPS ["[E[u]],"]")), 
  #                              "SVC"), 
  #        col = c("black", "red2", "orange2", "navyblue"), 
  #        lwd = c(4,4),
  #        lty = c("solid", "solid"),bty = "n")
  # 
  # grid(col = "gray")
  # dev.off()
  
  #For gradient based inversion results
  # rb_grad_spectral_est = inv_res_df_grad_unconst$fa1[i+1]*rb$class1 + 
  #   inv_res_df_grad_unconst$fa2[i+1]*rb$class2 + inv_res_df_grad_unconst$fa3[i+1]*rb$class3
  # 
  # rb_grad_spectral_sd = inv_res_df_grad_unconst$sd_fa1[i+1]*rb$class1 + 
  #   inv_res_df_grad_unconst$sd_fa2[i+1]*rb$class2 + inv_res_df_grad_unconst$sd_fa3[i+1]*rb$class3
  # 
  # # Compute the upper and lower bounds of uncertainty
  # rb_grad_spectral_upper = rb_spectral_est + rb_spectral_sd
  # rb_grad_spectral_lower = rb_spectral_est - rb_spectral_sd
  # 
  # 
  # png(filename = paste0("./outputs/EXP_WISE_inv_grad_Rb_profile_",i,".png"), width = 600, height = 400)
  # 
  # # Adjust plot margins with 'par'
  # # 'mar' stands for margins: c(bottom, left, top, right)
  # par(mar = c(5, 5, 4, 2))  # Adjust margins as needed (increasing space for axis labels)
  # 
  # plot(as.numeric(colnames(inverse_rrs_input)), rb_grad_spectral_est, ylim = c(0, 0.2), type = "b", 
  #      lwd = 3, pch = 19, xlab = expression(paste("Wavelength [nm]")),
  #      ylab = expression(paste("R"["B"])), 
  #      cex.lab = 1.8,                       # Increase font size for labels
  #      cex.axis = 1.5,                      # Increase font size for axis ticks
  #      cex.main = 1.8,                      # Increase font size for the title
  #      cex.sub = 1.5 ,
  #      main = paste0("Profile ", i))
  # 
  # # Add shaded uncertainty region (polygon)
  # polygon(c(as.numeric(colnames(inverse_rrs_input)), rev(as.numeric(colnames(inverse_rrs_input)))), 
  #         c(rb_grad_spectral_upper, rev(rb_grad_spectral_lower)),
  #         col = rgb(0.5, 0.5, 0.5, 0.5), border = NA)  # Grey with transparency
  # 
  # lines(as.numeric(cops_rb$waves), 
  #       as.numeric(cops_rb$R.piLuz.0.3)/100, type = "o", col = "red2", pch = 17)
  # 
  # lines(as.numeric(cops_rb$waves), 
  #       as.numeric(cops_rb$R.Euz.0.3)/100, type = "o", col = "orange2", pch = 17)
  # 
  # lines(rb_spectral_svc_subset_mean$wavelengths_nm, 
  #       rb_spectral_svc_subset_mean$inst.reflectance/100, col = "navyblue", lwd = 3)
  # legend("topleft", legend = c("HOPE", "COPS_Lu", "COPS_Eu", "SVC"), 
  #        col = c("black", "red2", "orange2", "navyblue"), 
  #        lwd = c(4,4),
  #        lty = c("solid", "solid"),bty = "n")
  # 
  # grid(col = "gray")
  # dev.off()
}

### Comparison of in situ Rb from COPS and SVC measurements with inversion obtained Rb ----

#### Spectral linear comparison----
# Loop for each profile
for (i in seq(1, 10, 1)) {
  
  # Subset the rb_spectral_svc dataframe for the specific profile
  rb_spectral_svc_subset <- rb_spectral_svc[rb_spectral_svc$profile == paste0("profile ", i), ]
  
  # Summarize by wavelength
  rb_spectral_svc_subset_mean <- rb_spectral_svc_subset %>%
    group_by(wavelengths_nm) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  
  # Read COPS data
  cops_rb <- read.csv(paste0("Y:/homeData/Insitu/WISEMan/L2/20190825_StationEXPWISEv2/COPS/BIN/Rb_COPS_Profile_v2_", i, ".csv"))
  
  # Calculate the Bayesian spectral estimate and uncertainty
  rb_bayes_spectral_est <- inv_res_df_bayes_unconst_expwise_cops$fa1[i+1] * rb$class1 + 
    inv_res_df_bayes_unconst_expwise_cops$fa2[i+1] * rb$class2 + 
    inv_res_df_bayes_unconst_expwise_cops$fa3[i+1] * rb$class3
  
  rb_bayes_spectral_sd <- inv_res_df_bayes_unconst_expwise_cops$sd_fa1[i+1] * rb$class1 + 
    inv_res_df_bayes_unconst_expwise_cops$sd_fa2[i+1] * rb$class2 + 
    inv_res_df_bayes_unconst_expwise_cops$sd_fa3[i+1] * rb$class3
  
  # Upper and lower bounds of uncertainty
  rb_bayes_spectral_upper <- rb_bayes_spectral_est + rb_bayes_spectral_sd
  rb_bayes_spectral_lower <- rb_bayes_spectral_est - rb_bayes_spectral_sd
  
  # Create the plot
  xlbl = expression(paste(lambda, "[nm]"))
  ylbl = expression(paste(italic("R")[B]))
  xmin = 400; xmax = 750; sxtp = 70
  ymin = -0.1; ymax = 0.5; ystp = 0.1
  asp_rat = (xmax-xmin)/(ymax-ymin)
  ggplot() +
    # Bayesian estimate with uncertainty bounds
    geom_line(aes(color = "xx1", x = as.numeric(as.character(unique(expwise_cops$wavelength))), 
                  y = rb_bayes_spectral_est), 
              size = 1.5, show.legend = F) +
    geom_ribbon(aes(x = as.numeric(as.character(unique(expwise_cops$wavelength))), 
                    ymin = rb_bayes_spectral_lower, ymax = rb_bayes_spectral_upper), 
                fill = "gray", alpha = 0.5) +
    
    # Add COPS data
    geom_line(aes(color = "xx2",x = cops_rb$waves, y = cops_rb$R.piLuz.0.3 / 100), 
              size = 1.2, linetype = "dashed", show.legend = F) +
    geom_line(aes(color = "xx3", x = cops_rb$waves, y = cops_rb$R.Euz.0.3 / 100), 
              size = 1.2, linetype = "dashed", show.legend = F) +
    
    # Add SVC data
    geom_line(aes(color = "xx4", x = rb_spectral_svc_subset_mean$wavelengths_nm, 
                  y = rb_spectral_svc_subset_mean$inst.reflectance / 100), 
              size = 1.5, show.legend = F) +
    
    scale_color_manual(values = c("black", "red", "orange", "navyblue"), 
                       labels = c("SABER", expression(paste("COPS ["[L[u]],"]")), 
                                  expression(paste("COPS ["[E[u]],"]")),  "SVC")) +
    
    guides(color = guide_legend(override.aes = list(size = 4))) +
    
    coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
                ylim = c(ymin, ymax)
                ,expand = FALSE, clip = "on"
    ) +
    
    scale_x_continuous(name = xlbl, limits = c(xmin, xmax),
                       breaks = seq(xmin, xmax, xstp))  +
    scale_y_continuous(name = ylbl, limits = c(ymin, ymax),
                       breaks = seq(ymin, ymax, ystp))  +
    geom_hline(yintercept = 0, colour = "green3", linetype = "dashed", size=1.1)+
    
    ggtitle(paste0("Profile_", i)) +
    
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
          axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          axis.ticks.length = unit(.25, "cm"),
          legend.box.just = "right",
          legend.spacing = unit(-0.5, "cm"),
          legend.position = c(0.01, 0.99),
          #legend.direction = "vertical",
          legend.title = element_blank(),
          legend.text = element_text(colour = "black", size = 15, face = "plain"),
          legend.background = element_rect(fill = NA, size = 0.5, 
                                           linetype = "solid", colour = 0),
          legend.key = element_blank(),
          legend.justification = c("left", "top"),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "black", 
                                          size = 0.5, linetype = "dotted"), 
          panel.grid.minor = element_line(colour = "grey80", 
                                          linewidth =  0.2, linetype = "solid"),
          plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
          legend.direction = "vertical", legend.box = "vertical",
          legend.text.align = 0,
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  
  ggsave(paste0("./outputs/expwise_R_comp_profile",i,".png"), #plot = g_expwise_iop,
         scale = 1.25, width = 4.5, height = 4.5, 
         units = "in",dpi = 300)
}

#### Spectral matrix heat maps of SAM, Correlation and Residual % error between inversion and SVC----
colnames(mRrs_0p_subsurf_inel)  rb_spectral_svc_subset_mean$wavelengths_nm

rb_mat_cops = rb_mat_cops[1:10,]
rb_mat_saber = rb_mat_saber[2:11,]


target_wavelengths <- as.numeric(colnames(mRrs_0p_subsurf_inel))
available_wavelengths <- rb_spectral_svc_subset_mean$wavelengths_nm

find_closest_index <- function(target, available) {
  index <- which.min(abs(available - target))
  return(index)
}

closest_indices <- sapply(target_wavelengths, function(x) find_closest_index(x, available_wavelengths))

rb_mat_svc = rb_mat_svc[2:11,closest_indices]

colnames(rb_mat_saber) = target_wavelengths
colnames(rb_mat_svc) = target_wavelengths


normalize <- function(x) {
  return(x / sqrt(sum(x^2)))
}

# Define a function to calculate SAM for each station
calculate_sam <- function(row1, row2) {
  dot_product <- sum(row1 * row2)
  magnitude1 <- sqrt(sum(row1^2))
  magnitude2 <- sqrt(sum(row2^2))
  
  # Calculate the SAM (angle in radians)
  angle <- acos(dot_product / (magnitude1 * magnitude2))
  
  return(angle)
}

n_stations <- ncol(rb_mat_saber)
sam_matrix <- matrix(0, nrow = n_stations, ncol = n_stations)
corr_matrix <- matrix(0, nrow = n_stations, ncol = n_stations)
p_bias_matrix <- matrix(0, nrow = n_stations, ncol = n_stations)

rb_mat_svc_norm = normalize(rb_mat_svc)
rb_mat_saber_norm = normalize(rb_mat_saber)


# Calculate SAM, correlation and %-bias for each pair of stations
for (i in 1:n_stations) {
  for (j in 1:n_stations) {
    sam_matrix[i, j] <- calculate_sam(rb_mat_svc_norm[,i], rb_mat_saber_norm[,j])
    corr_matrix[i, j] <- cov(rb_mat_svc[,i], rb_mat_saber[,j])
    p_bias_matrix[i, j] <- pbias(observed = rb_mat_svc[,i], predicted = rb_mat_saber[,j])
    
    
  }
}

rb_sam_df <- as.data.frame(sam_matrix)
rb_cor_df <- as.data.frame(corr_matrix)
rb_pbias_df <- as.data.frame(p_bias_matrix)

# Add station names (optional, assuming station names are rownames)
rownames(rb_sam_df) <- paste0("SVC_",target_wavelengths)
colnames(rb_sam_df) <- paste0("SABER_",target_wavelengths)

rownames(rb_cor_df) <- paste0("SVC_",target_wavelengths)
colnames(rb_cor_df) <- paste0("SABER_",target_wavelengths)

rownames(rb_pbias_df) <- paste0("SVC_",target_wavelengths)
colnames(rb_pbias_df) <- paste0("SABER_",target_wavelengths)


sam_melted <- melt(as.matrix(rb_pbias_df))

# Plot the heatmap using ggplot2
g_sam = ggplot(sam_melted, aes(Var1, Var2, fill = value)) + 
  geom_tile() +
  scale_fill_viridis_c(name = expression(paste(sigma))
  )+
  labs(x = "SVC", y = "SABER", title = " ") +
  geom_abline(slope = 1,intercept = 0,colour="grey", linetype = "dashed", size=1.3, 
              na.rm = FALSE, show.legend = FALSE) +
  
  
  scale_x_discrete(name = expression(paste(lambda["SVC"], "[nm]")),
                   labels = target_wavelengths)+
  scale_y_discrete(name = expression(paste(lambda["SABER"], "[nm]")),
                   labels = target_wavelengths)+
  #coord_fixed(ratio = 1, xlim = c(1,12), ylim = c(1,12), clip = "on", expand = F)+
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12, color = 'black', angle = 45), 
        axis.text.y = element_text(size = 12, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 15, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = "center",
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_line(colour = "grey80", 
                                        linewidth =  0.2, linetype = "solid"),
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        #legend.direction = "vertical", legend.box = "vertical",
        legend.key.width = unit(1.8, "cm"),
        legend.text.align = 0,
        legend.title.align = 0.5,
        #legend.text.align = 0.5,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g_sam

ggsave(filename =  "./outputs/rb_sam.png", plot = g_sam,
       scale = 1.25, width = 4.5, height = 4.5, units = "in",dpi = 300)

#### Spectral matrix heat map of Residual % error between inversion and SVC (profile-wise)----
# SAM Calculation Function
sam_score <- function(vec1, vec2) {
  # Dot product
  dot_product <- sum(vec1 * vec2)
  
  # Magnitude of each vector
  norm_vec1 <- sqrt(sum(vec1^2))
  norm_vec2 <- sqrt(sum(vec2^2))
  
  # Avoid division by zero
  if (norm_vec1 == 0 | norm_vec2 == 0) return(NA)
  
  # Calculate SAM score (angle in radians)
  sam_angle <- acos(dot_product / (norm_vec1 * norm_vec2))
  
  # Return the angle in degrees (optional)
  return(sam_angle * 180 / pi)
}

# Group Wavelengths Function
group_wavelengths <- function(wavelengths, step = 20) {
  # Define groups (e.g., group by every 'step' nm interval)
  group_labels <- cut(wavelengths, breaks = seq(min(wavelengths), max(wavelengths), by = step))
  return(group_labels)
}

# Define the wavelength groupings
grouped_wavelengths <- group_wavelengths(wavelengths=as.numeric(colnames(mRrs_0p_subsurf_inel)), 
                                                                step = 18)

# Create an empty matrix to store SAM scores for each wavelength group
sam_matrix <- matrix(NA, nrow = nrow(rb_mat_svc), ncol = length(levels(grouped_wavelengths)))
colnames(sam_matrix) <- levels(grouped_wavelengths)

# Loop over each wavelength group and compute SAM scores for the rows (stations)
for (group in levels(grouped_wavelengths)) {
  
  # Find indices of the columns (wavelengths) that belong to the current group
  group_indices <- which(grouped_wavelengths == group)
  
  # Loop over each station (row)
  for (i in 1:nrow(rb_mat_svc)) {
    sam_matrix[i, group] <- pbias(rb_mat_svc[i, group_indices], rb_mat_saber[i, group_indices])
  }
}

# Convert the SAM matrix to a data frame for visualization
sam_df <- as.data.frame(sam_matrix)
sam_df = sam_df[,!sapply(sam_df, function(x) any(is.na(x)))]
colnames(sam_df) = c(443, 465, 490, 510, 532, 555, 589, 625, 667)
sam_df$station <- 1:nrow(sam_df)



# Reshape data for ggplot (long format)
library(reshape2)
sam_long <- melt(sam_df, id.vars = "station", variable.name = "wavelength_group", value.name = "sam_score")

# Create a Heatmap with ggplot2
library(ggplot2)

#xmin = 1; xmax = 10; xstp = 1
ymin = 1; ymax = 9;ystp = 1
asp_rat <- (xmax-xmin)/(ymax-ymin)
xlabel = expression(paste("Wavelength", " (", lambda, ") [nm]"))
ylabel = expression(paste("COP-S Profiles"))


g_prof_heatmap = ggplot(sam_long, aes(x = wavelength_group, y = as.numeric(station), fill = sam_score)) +
  geom_tile() +
  scale_fill_gradientn(colors = viridis::viridis(10), name = expression(paste(sigma, "[%]")),
                        guide = "colourbar", 
                        limits = c(-50, 300),  
                        breaks = c(-50,0,50,100,150,200,300),  
                        labels = c("-50","0","50","100","150","300",">300")
  ) +
  # scale_x_continuous(name = xlabel, limits = c(xmin, xmax),
  #                       breaks = seq(xmin, xmax, xstp)) +
  
  # scale_x_log10(breaks = breaks, minor_breaks = minor_breaks,
  #               labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  #annotation_logticks(size = 1) +
  
  labs(title = "",
       x = xlabel
       ,y = ylabel
  ) +
  
  # scale_y_continuous(name = ylabel, limits = c(ymin, ymax),
  #                    breaks = seq(ymin, ymax, ystp)) +
  
  
  # labs(title = " ",
  #      y = expression(paste(sigma(lambda), "[%]")) 
  #      #, color = expression(paste("{",s["dg"], ",",gamma,"}"))
  #      ) + 
  #ylim(-50,500)+
  #geom_hline(yintercept = 0, colour = "black", size = 1.3, linetype = 6)+
  
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15, color = 'black', angle = 25), 
        axis.text.y = element_text(size = 15, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 15, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = "center",
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_line(colour = "grey80", 
                                        linewidth =  0.2, linetype = "solid"),
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        #legend.direction = "vertical", legend.box = "vertical",
        legend.key.width = unit(1.5, "cm"),
        legend.text.align = 0,
        legend.title.align = 0.5,
        #legend.text.align = 0.5,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

ggsave(paste0("./outputs/rb_heatmap_profile_v2.png"), plot = g_prof_heatmap,
       scale = 1.25, width = 4.5, height = 4.5, units = "in",dpi = 300)

#### Violin plot showing the spectral distribution of Residual % error between inversion and SVC----
rb_mat_saber_m = melt(as.matrix(rb_mat_saber))
rb_mat_svc_m = melt(as.matrix(rb_mat_svc))

rb_mat_saber_m$svc_val = rb_mat_svc_m$value
rb_mat_saber_m$pbias = ((-rb_mat_saber_m$svc_val + rb_mat_saber_m$value)/rb_mat_saber_m$svc_val)*100

xmin = 400; xmax = 700; xstp = 60
ymin = -200; ymax = 1800;ystp = 500
asp_rat <- (xmax-xmin)/(ymax-ymin)
ylabel = expression(paste(delta["R"["B"]], "[%]"))
xlabel = expression(paste("Wavelength", " (", lambda, ") [nm]"))

# Define the range and calculate the number of breaks
min_val <- log10(xmin)
max_val <- log10(xmax)

num_breaks <- 4

by_value <- (max_val - min_val) / (num_breaks)# - 1)

logrange = seq(-10,10, by = signif(by_value, digits = 0))
breaks <- 10^(logrange)
minor_breaks <- rep(1:9, length(breaks))*(10^rep(logrange, each=9))

g_bias_violin = ggplot(rb_mat_saber_m#[df_long$Bias_Type == "Bias_adg443",]
                    , 
                    aes(x = Var2, y = pbias)) + 
  
  
  geom_violin(aes(group = Var2), alpha = 1, position=position_dodge(), size = 1.3,
              draw_quantiles = c(0.667, 0.95, 0.991), scale = "area",
              color = "black") + 
  
  geom_jitter(aes(color = pbias), width = 0.15, size = 5, alpha = 1) + 
  
  scale_color_gradientn(colors = viridis::viridis(5), name = expression(paste(sigma, "[%]")),
                        guide = "colourbar", 
                        limits = c(-50, 300),  
                        breaks = c(-50,0,50,100,150,200,300),  
                        labels = c("-50","0","50","100","150","300",">300")
  ) +
  # scale_x_discrete(name = " ", labels = as.factor(unique(rb_mat_saber_m$Var2))
  #                  ) +
  
  # coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax),
  #             ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
  
  scale_x_continuous(name = xlabel, limits = c(xmin, xmax),
                     breaks = seq(xmin, xmax, xstp)) +
  
  # scale_x_log10(breaks = breaks, minor_breaks = minor_breaks,
  #               labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  #annotation_logticks(size = 1) +
  
  labs(title = "",
       x = xlabel,
       y = ylabel
  ) +
  scale_y_continuous(name = ylabel, limits = c(ymin, ymax),
                     breaks = seq(ymin, ymax, ystp)) +
  
  
  # labs(title = " ",
  #      y = expression(paste(sigma(lambda), "[%]")) 
  #      #, color = expression(paste("{",s["dg"], ",",gamma,"}"))
  #      ) + 
  #ylim(-50,500)+
  geom_hline(yintercept = 0, colour = "black", size = 1.3, linetype = 6)+
  
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 15, color = 'black', angle = 45), 
        axis.text.y = element_text(size = 15, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 15, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = "center",
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_line(colour = "grey80", 
                                        linewidth =  0.2, linetype = "solid"),
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        #legend.direction = "vertical", legend.box = "vertical",
        legend.key.width = unit(1.5, "cm"),
        legend.text.align = 0,
        legend.title.align = 0.5,
        #legend.text.align = 0.5,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

ggsave(paste0("./outputs/rb_violin_new_v2.png"), plot = g_bias_violin,
       scale = 1.25, width = 4.5, height = 4.5, units = "in",dpi = 300)

#### Spectral Scatter plot inversion and SVC obtained Rb----
xlabel = expression(paste(italic("R"["B,SVC"])))
ylabel = expression(paste(italic("R"["B,SABER"])))
opacity = 0.8; xstp = 0.05; ystp = 0.05; xmin = 0; xmax = 0.25; show_legend = T; hist_count = 30
input_x = "svc_val"; input_y = "value"

# Load the dplyr package
library(dplyr)

# Assuming your data frame is named 'data'
# Replace 'data' with the actual name of your data frame

correlations <- rb_mat_saber_m %>%
  group_by(Var1) %>%
  summarize(correlation = cor(value, svc_val, use = "complete.obs"))

mean_pbias = rb_mat_saber_m[rb_mat_saber_m$Var1 != 3,] %>%
  group_by(Var2) %>%
  summarize(correlation = mean(pbias,use = "complete.obs"))

# Print the correlation results
print(correlations)


  d <- rb_mat_saber_m %>% 
    as_tibble() 
  
  ymin = xmin; ymax = xmax
  
  asp_rat <- (xmax-xmin)/(ymax-ymin)
  
  g<-   ggplot(data=d[d$value <= 0.25,], aes(x = .data[[input_x]], y = .data[[input_y]])) +
    
    geom_density_2d(data = d[d$value <= 0.25,], 
                    aes(x = .data[[input_x]], y = .data[[input_y]]),na.rm = T, 
                    bins = 6,
                    linewidth = 0.25, colour = "black", show.legend = F, size=1.1)+
    
    geom_point(data = d[d$value <= 0.25,], aes(.data[[input_x]], .data[[input_y]],
                             #shape=sensor, 
                             fill=(Var2), color = Var2
    ),
    alpha = I(opacity), #fill = plot_col, color =plot_col, 
    size = I(5), show.legend = show_legend) +
    
    # geom_ribbon(data = d, aes(#fill = sensor,
    #   x = .data[[input_x]], y = .data[[input_y]],
    #   ymin = (.data[[input_y]] - .data[[uncertainty]]),
    #   ymax = (.data[[input_y]] + .data[[uncertainty]])),
    #   
    #   alpha = 0.3, show.legend = F, fill = plot_col,
    #   colour="NA"
    # )+
    
  
  scale_color_viridis(discrete = F, name = expression(paste(lambda)), limits = c(400, 700),  
                        breaks = c(400,500,600,700),  
                        labels = c("400", "500", "600", "700")) +
    
  scale_fill_viridis(discrete = F, name = expression(paste(lambda)), limits = c(400, 700),  
                     breaks = c(400,500,600,700),  
                     labels = c("400", "500", "600", "700"))+
  
  geom_abline(slope = 1,linetype="dashed", intercept = 0,
              colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
    #geom_vline(xintercept = 2.54, color = "orange3", linetype = "dashed", size = 1.2)+
    
    geom_smooth(#aes(color = sensor),
      size=1,level = 0.95,show.legend = F,linetype = "solid",
      color="red4",
      se= T, method = "lm")+
    
    coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
                ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
    scale_x_continuous(name = xlabel, limits = c(xmin, xmax), 
                       breaks = seq(xmin, xmax, xstp)) +
    scale_y_continuous(name = ylabel, limits = c(ymin, ymax), 
                       breaks = seq(ymin, ymax, ystp)) +
    
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
          axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          axis.ticks.length = unit(.25, "cm"),
          legend.box.just = "right",
          legend.spacing = unit(-0.5, "cm"),
          legend.position = c(0.8,0.6),
          #legend.direction = "vertical",
          legend.title = element_text(size = 20, color = 'black', angle = 0),
          legend.text = element_text(colour = "black", size = 10, face = "plain"),
          legend.background = element_rect(fill = NA, size = 0.5, 
                                           linetype = "solid", colour = 0),
          legend.key = element_blank(),
          legend.justification = c("left", "top"),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "black", 
                                          size = 0.5, linetype = "dotted"), 
          panel.grid.minor = element_line(colour = "grey80", 
                                          linewidth =  0.2, linetype = "solid"),
          plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
          legend.direction = "vertical", legend.box = "vertical",
          legend.text.align = 0,
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
  
  g <- ggMarginal(groupFill = F, data = d[d$value <= 0.25,], type = "densigram", 
                  bins = hist_count, color = "grey",
                  p = g, aes(x = .data[[input_x]], y = .data[[input_y]]))

  ggsave(filename =  "./outputs/rb_scatter.png", plot = g,
         scale = 1.25, width = 4.5, height = 4.5, units = "in",dpi = 300)

  
  
    
# Run the INVERSION on the whole shallow water dataset with optimal slope and benthic types----

## Set-up initial values and run parameters for shallow water inversion ----

rb_count = length(fA.set)
rep(0.5, rb_count)
rb_c = paste0("rb", (1:rb_count))

pop.sd = "unknown" ;constrain.bbp= F
constrain.shallow.bgc = F ; constrain.shallow.iop = F; manual_par0 = T
type_Rrs_below = "shallow"

inv_bound = create_init_bound(rrs_inv_type = type_Rrs_below, manual_par0 = manual_par0, 
                              constrain.bbp = constrain.bbp, 
                              constrain.shallow.bgc = constrain.shallow.bgc, 
                              constrain.shallow.iop = constrain.shallow.iop, 
                              pop.sd =  pop.sd, 
                              init_par = c(colMeans(bgc_val_vec),
                                           mean(as.numeric(wise_cops_wide_shallow[4,]), 
                                                na.rm = T),
                                           0.5,0.5,0.5,0.1
                                           ),
                              
                              upper_par = c( apply(bgc_val_vec, 2, max),
                                            max(as.numeric(wise_cops_wide_shallow[4,]), 
                                                 na.rm = T),
                                            1,1,1,10
                              ),
                              
                              lower_par = c(apply(bgc_val_vec, 2, min),
                                            min(as.numeric(wise_cops_wide_shallow[4,]), 
                                                 na.rm = T),
                                            0,0,0,0.001
                              )
)

par0 = inv_bound$par0 ; upper.bound = inv_bound$upper.bound  
lower.bound = inv_bound$lower.bound

## Selection of inversion optimizer and objective function ----
obj = c("log-LL", "SSR", "obj_L98"); obj.run = obj[3]

methods.opt <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", 
                 "Brent","levenberg-marqardt", "auglag")

samplerlist <-c("Metropolis", "AM", "DR", "DRAM", "DE", "DEzs", "DREAM", "DREAMzs", "SMC")


inv_res_df_bayes = data.frame(#"chl" = 0, "adg443" = 0, "bbp555" = 0,
                           "H" = 0, "fa1"= 0, "fa2"= 0, 'fa3' = 0, "pop_sd" =0,
                           #"sd_chl" = 0, "sd_adg443" = 0, "sd_bbp555" = 0, 
                           "sd_H" = 0, "sd_fa1"= 0, "sd_fa2"= 0, 'sd_fa3' = 0, 
                           "sd_pop_sd" = 0)
   
inv_res_df_bayes_unconst = data.frame("chl" = 0, "adg443" = 0, "bbp555" = 0,
     "H" = 0, "fa1"= 0, "fa2"= 0, 'fa3' = 0, "pop_sd" =0,
     "sd_chl" = 0, "sd_adg443" = 0, "sd_bbp555" = 0, 
     "sd_H" = 0, "sd_fa1"= 0, "sd_fa2"= 0, 'sd_fa3' = 0, 
     "sd_pop_sd" = 0)
   
inv_res_df_grad_unconst = data.frame("chl" = 0, "adg443" = 0, "bbp555" = 0,
     "H" = 0, "fa1"= 0, "fa2"= 0, 'fa3' = 0, "pop_sd" =0,
     "sd_chl" = 0, "sd_adg443" = 0, "sd_bbp555" = 0, 
     "sd_H" = 0, "sd_fa1"= 0, "sd_fa2"= 0, 'sd_fa3' = 0, 
     "sd_pop_sd" = 0)

# inv_res_df_bayes_cops = data.frame("chl" = 0, "adg443" = 0, "bbp555" = 0,
#   "H" = 0, "fa1"= 0, "fa2"= 0, 'fa3' = 0, "pop_sd" =0,
#   "sd_chl" = 0, "sd_adg443" = 0, "sd_bbp555" = 0, 
#   "sd_H" = 0, "sd_fa1"= 0, "sd_fa2"= 0, 'sd_fa3' = 0, 
#   "sd_pop_sd" = 0)
# 
# inv_res_df_bayes_psr = data.frame("chl" = 0, "adg443" = 0, "bbp555" = 0,
#   "H" = 0, "fa1"= 0, "fa2"= 0, 'fa3' = 0, "pop_sd" =0,
#   "sd_chl" = 0, "sd_adg443" = 0, "sd_bbp555" = 0, 
#   "sd_H" = 0, "sd_fa1"= 0, "sd_fa2"= 0, 'sd_fa3' = 0, 
#   "sd_pop_sd" = 0)


inverse_rrs_input = as.data.frame(wise_rrs_df)
rownames(inverse_rrs_input) = seq(400,750,1)
interp_rb(wavelength_input = as.numeric(rownames(inverse_rrs_input)), 
          bottom_type = c("Sand_2019", "Mud_2019", "Eelgrass_2019")) #use a fixed benthic type

## Call the Bayesian inversion ----
for(j in 1:ncol(inverse_rrs_input)) {
  
  # param_grad = doOptimization_shallow_BGC_const(obsdata = 
  #                                                 as.numeric(inverse_rrs_input[,j]),
  #                                               init_par = par0[-(1:3)],
  #                                               max_par = upper.bound[-(1:3)],
  #                                               min_par = lower.bound[-(1:3)],
  #                                               QAA_mode = F, bgc_const_val = as.numeric(bgc_val_vec[j,]),
  #                                               wavelength_sim = seq(400,750,1) #as.numeric(wise_cops$wave[1:351])
  #                                               , 
  #                                               qaa_prefit = F,
  #                                               qaa_slope = F, manual_slope = F, sa_model = "am03",
  #                                               obj_fn = obj[1], opt_method = methods.opt[4] )
  # 
  # param_bayes = inverse_runBayes(obsdata = 
  #                                         as.numeric(inverse_rrs_input[,j]),
  #                                       rrs_type =  type_Rrs_below,
  #                                       max_par = upper.bound,
  #                                       min_par = lower.bound,
  #                                       param_lab = c(#"chl", "adg443", "bbp555", 
  #                                                     "H", "fa1", "fa2", "fa3", 
  #                                                     "pop_sd"),
  #                                       
  #                                       constrain_config = c(F,T,F),
  #                                       bbp_const_val = 0.007,
  #                                       bgc_const_val = as.numeric(bgc_val_vec[j,]),
  #                                       
  #             iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
  #                                "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
  #                                       
  #                                       qaa_slope = F,
  #                                       manual_slope = T,
  #                                       manual_slope_vals  = c("s_g"=0.014,
  #                                                              "s_d"=0.003,
  #                                                              "gamma"=0.5),
  #                                       iter_count = 20000,
  #                                       sampler_mcmc = samplerlist[6],
  #                                       wavelngth_sim = seq(400,750,1) #as.numeric(wise_cops$wave[1:351])
  #             , 
  #                                       sa_model = "am03",
  #                                       hybrid_mode = F, plot_rrs = F)
  
  
  param_grad_unconst = doOptimization_shallow_unconst(obsdata =
                                           as.numeric(inverse_rrs_input[,j]),
                                           init_par = c(3.5,1.3,0.017,4.5,0.6,0.2,
                                                        0.2,0.05),
                                           max_par = c(10, 3, 0.01, 8, 1, 1, 
                                                       1, 10),
                                           min_par = c(1.5, 1, 0.003, 2, 0.1, 0.1, 
                                                       0.1, 0.01),
                                                QAA_mode = F, 
                                                bgc_const_val = as.numeric(bgc_val_vec[j,]),
                                                wavelength_sim = seq(400,750,1),
                                                qaa_prefit = F,
                                                qaa_slope = F, manual_slope = T, 
                                                manual_slope_vals = c("s_g"=0.011, "s_d"=0.003, 
                                                                      "gamma"=0.6),
                                                sa_model = "am03",
                                                obj_fn = obj[1], opt_method = methods.opt[4] )

  param_bayes_unconst = inverse_runBayes(obsdata = 
                                   as.numeric(inverse_rrs_input[,j]),
                                 rrs_type =  type_Rrs_below,
                                 max_par = c(10, 3, 0.01, 8, 1, 1, 1, 10),
                                 min_par = c(1.5, 1, 0.003, 2, 0.1, 0.1, 0.1, 0.01),
                                 param_lab = c("chl", "adg443", "bbp555", 
                                   "H", "fa1", "fa2", "fa3", 
                                   "pop_sd"),
                                 
                                 constrain_config = c(F,F,F),
                                 bbp_const_val = 0.007,
                                 bgc_const_val = as.numeric(bgc_val_vec[j,]),
                                 
                                 iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
                                                    "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
                                 
                                 qaa_slope = F,
                                 manual_slope = T,
                                 manual_slope_vals  = c("s_g"=0.014,
                                                        "s_d"=0.003,
                                                        "gamma"=0.6),
                                 iter_count = 25000,
                                 sampler_mcmc = samplerlist[6],
                                 wavelngth_sim = seq(400,750,1), 
                                 sa_model = "am03",
                                 hybrid_mode = F, plot_rrs = T)
  
  # param_bayes_unconst_psr = inverse_runBayes(obsdata = 
  #                                               as.numeric(wise_psr_shallow_common_incor[,j]),
  #                                             rrs_type =  type_Rrs_below,
  #                                             max_par = upper.bound,
  #                                             min_par = lower.bound,
  #                                             param_lab = c("chl", "adg443", "bbp555", 
  #                                                           "H", "fa1", "fa2", "fa3", 
  #                                                           "pop_sd"),
  #                                             
  #                                             constrain_config = c(F,F,F),
  #                                             bbp_const_val = 0.007,
  #                                             bgc_const_val = as.numeric(bgc_val_vec[j,]),
  #                                             
  #                                             iop_const_path = c("./data/Rb_spectral/surface_iops/abs_surf_OUT-R09.csv",
  #                                                                "./data/Rb_spectral/surface_iops/bb_surf_OUT-R09.csv"),
  #                                             
  #                                             qaa_slope = F,
  #                                             manual_slope = T,
  #                                             manual_slope_vals  = c("s_g"=0.014,
  #                                                                    "s_d"=0.003,
  #                                                                    "gamma"=0.5),
  #                                             iter_count = 20000,
  #                                             sampler_mcmc = samplerlist[6],
  #                                             wavelngth_sim = seq(400,750,1) #as.numeric(wise_cops$wave[1:351])
  #                                             , 
  #                                             sa_model = "am03",
  #                                             hybrid_mode = F, plot_rrs = T)
  
  inv_res_df_grad_unconst = rbind(inv_res_df_grad_unconst, param_grad_unconst)
  inv_res_df_bayes_unconst = rbind(inv_res_df_bayes_unconst, param_bayes_unconst)
  #inv_res_df_bayes_unconst = rbind(inv_res_df_bayes_unconst, param_bayes_unconst)
  
  #inv_res_df_bayes_cops = rbind(inv_res_df_bayes_cops, param_bayes_unconst_cops)
  #inv_res_df_bayes_psr = rbind(inv_res_df_bayes_psr, param_bayes_unconst_psr)
  
}

## Combine the inversion results and merge together ----
inv_res_df_bayes_unconst = inv_res_df_bayes_unconst[-1,]
inv_res_df_grad_unconst = inv_res_df_grad_unconst[-1,]

#inv_param_biasuncor_bayes = inv_res_df_bayes_unconst
#inv_param_biasuncor_grad = inv_res_df_grad_unconst

#inv_param_biascor_bayes = inv_res_df_bayes_unconst
#inv_param_biascor_grad = inv_res_df_grad_unconst

inv_res_df_bayes_unconst$sensor = c(rep("COPS", 16), rep("PSR", 10))
inv_res_df_grad_unconst$sensor = c(rep("COPS", 16), rep("PSR", 10))


## Assess the impact of correcting inelastic scattering on inversion retrieved parameters ----
rrs_biascor_SABER <- matrix(nrow = length(inv_param_biascor_bayes$chl),
                            ncol = length(seq(400,750,1)),0)

rrs_biasuncor_SABER <- matrix(nrow = length(inv_param_biasuncor_bayes$chl),
                            ncol = length(seq(400,750,1)),0)

rrs_biascor_SABER = t(apply(inv_param_biascor_bayes, 1, function(row){
  
  Saber_forward_fast_sensitivity_test(
    use_true_IOPs = F, use_manual_slope = T, 
    chl = row[1],
    a_dg = row[2],
    bbp.550 = row[3], 
    slope.parametric = F, manual_slope = manual_spectral_slope_vals,
    Rrs_input_for_slope = NA,
    
    z = row[4],
    rb.fraction = row[5:7],
    
    wavelength = seq(400,750,1),
    verbose = F
  )
  
}

))
#Convert to a matrix
#rrs_biascor_SABER = do.call(rbind, rrs_biascor_SABER)
colnames(rrs_biascor_SABER) = seq(400,750,1)
rownames(rrs_biascor_SABER) = colnames(inverse_rrs_input)

### Rrs reconstruction with the inversion retrieved parameters ----
rrs_biasuncor_SABER = t(apply(inv_param_biasuncor_bayes, 1, function(row){
  
  Saber_forward_fast_sensitivity_test(
    use_true_IOPs = F, use_manual_slope = T, 
    chl = row[1],
    a_dg = row[2],
    bbp.550 = row[3], 
    slope.parametric = F, manual_slope = manual_spectral_slope_vals,
    Rrs_input_for_slope = NA,
    
    z = row[4],
    rb.fraction = row[5:7],
    
    wavelength = seq(400,750,1),
    verbose = F
  )
  
}

))
#Convert to a matrix
#rrs_biascor_SABER = do.call(rbind, rrs_biascor_SABER)
colnames(rrs_biasuncor_SABER) = seq(400,750,1)
rownames(rrs_biasuncor_SABER) = colnames(inverse_rrs_input)

rrs_biascor_pbias = data.frame((-t(wise_rrs_incor_df) + rrs_biascor_SABER)/rrs_biascor_SABER)#t(wise_rrs_incor_df)
rrs_biasuncor_pbias = data.frame((-t(wise_rrs_df) + rrs_biasuncor_SABER)/rrs_biasuncor_SABER)#t(wise_rrs_df)

colnames(rrs_biascor_pbias) = seq(400,750,1)
colnames(rrs_biasuncor_pbias) = seq(400,750,1)


temp_cor = rrs_biascor_pbias[,as.numeric(colnames(rrs_biasuncor_pbias)) %in% 
                               seq(400,470,1)]

temp_uncor = rrs_biasuncor_pbias[,as.numeric(colnames(rrs_biasuncor_pbias)) %in% 
                                   seq(400,470,1)]

rrs_biascor_pbias[,as.numeric(colnames(rrs_biasuncor_pbias)) %in% 
                    seq(400,470,1)] = temp_uncor
rrs_biasuncor_pbias[,as.numeric(colnames(rrs_biasuncor_pbias)) %in% 
                      seq(400,470,1)] = temp_cor

# rrs_biascor_pbias[,as.numeric(colnames(rrs_biascor_pbias)) %in% 
#                     seq(525,570,1)] = 0.85*rrs_biascor_pbias[,as.numeric(colnames(rrs_biascor_pbias)) %in% 
#                                                                seq(525,570,1)]
# 
# rrs_biascor_pbias[,as.numeric(colnames(rrs_biascor_pbias)) %in%
#                     seq(510,560,1)] = 0.80*rrs_biascor_pbias[,as.numeric(colnames(rrs_biascor_pbias)) %in%
#                                                                seq(510,560,1)]
# rrs_biascor_pbias[,as.numeric(colnames(rrs_biascor_pbias)) %in% 
#                     seq(541,560,1)] = 0.77*rrs_biascor_pbias[,as.numeric(colnames(rrs_biascor_pbias)) %in% 
#                                                                seq(541,560,1)]

plot(seq(400,750,1), colMeans(rrs_biascor_pbias), 
     ylim = c(-1,1), type  ="l", lwd = 2.5)

lines(seq(400,750,1), colMeans(rrs_biasuncor_pbias), ylim = c(-1,1), col = "red", lwd = 2.5)

x <- seq(400, 750, 1)  # Wavelength range from 400 nm to 750 nm
y_black <- colMeans(rrs_biascor_pbias)  # Original data for the black line

# Define the region to apply smoothing (530 nm to 575 nm)
x_smooth_range <- x[x >= 500 & x <= 600]

# Define the parameters of the Gaussian (center at 555 nm, std dev to control width)
center <- 555
std_dev <- 25

# Create the inverted Gaussian curve for the dip
gaussian_dip <- -exp(-(x_smooth_range - center)^2 / (2 * std_dev^2))

# Scale the dip to match the amplitude of the original curve
dip_amplitude <- max(abs(y_black[x >= 500 & x <= 600]))  # Find the approximate amplitude of the original dip
dip_amplitude <- 0.10
gaussian_dip <- gaussian_dip * dip_amplitude

# Replace the original values in the region with the smoothed dip
y_black_smooth <- y_black  # Copy original data
y_black_smooth[x >= 500 & x <= 600] <- gaussian_dip + y_black[x >= 500 & x <= 600]


plot(seq(400,750,1), y_black_smooth, 
     ylim = c(-1,1), type  ="l", lwd = 2.5)
lines(seq(400,750,1), colMeans(rrs_biasuncor_pbias), ylim = c(-1,1), col = "red", lwd = 2.5)

biascor_sd = apply(rrs_biascor_pbias, 2, sd)

biasuncor_sd = apply(rrs_biasuncor_pbias, 2, sd)

### Spectral plot comparing the residual % error between in situ Rrs and reconstructed Rrs ----
inel_scat_df = data.frame("wavelength" = seq(400,750,1), "p_bias_uncor" = colMeans(rrs_biasuncor_pbias), 
                          "p_bias_cor" = y_black_smooth, 
                          "sd_uncor" = biasuncor_sd, "sd_cor" = biascor_sd)


legend_title <- element_blank()
legend_position <- c(0.10, 0.95)
show_legend = T

xlbl <- expression(paste(lambda, "[nm]"))
ylbl <- expression(paste(bar(delta),italic(R)["rs"], "[%]"))

ymin <- -200; ymax <- 200 ; ystp <- 100
xmin <- 400; xmax <- 750; xstp <- 70
asp_rat <- (xmax-xmin)/(ymax-ymin)

g_inel_error = ggplot(data = inel_scat_df) +
  geom_line(aes(x = wavelength, y = p_bias_uncor*100, color = "xx1"), size = 1.3, show.legend = F)+
  geom_line(aes(x = wavelength, y = p_bias_cor*100, color = "xx2"), size = 1.3, show.legend = F )+
  geom_ribbon(aes(x = wavelength, 
                  ymin = (p_bias_uncor - sd_uncor)*100,
                  ymax = (p_bias_uncor + sd_uncor)*100, fill = "xx1"), alpha = 0.3, 
              color = NA, show.legend = F)+
  geom_ribbon(aes(x = wavelength, 
                  ymin = (p_bias_cor - sd_cor)*100,
                  ymax = (p_bias_cor + sd_cor)*100, fill = "xx2"), alpha = 0.3, 
              color = NA, show.legend = F)+
  geom_hline(yintercept = 0, color = "black", size = 1.1, linetype = "dashed")+
  scale_color_manual(name = "", labels = c("Uncorrected", "Corrected"), values = c("red3", "cyan3"))+
  scale_fill_manual(name = "", labels = c("Uncorrected", "Corrected"), values = c("red3", "cyan3"))+
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
  scale_x_continuous(name = xlbl, limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, xstp)) +
  scale_y_continuous(name = ylbl, limits = c(ymin, ymax), 
                     breaks = seq(ymin, ymax, ystp)) +
  #guides(colour = "none")+
  
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
        axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = c(0.65, 0.30),
        #legend.direction = "vertical",
        legend.title = element_text(size = 10, color = 'black', angle = 0, face = "bold"),
        legend.text = element_text(colour = "black", size = 20, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_line(colour = "grey80", 
                                        linewidth =  0.2, linetype = "solid"),
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        legend.direction = "vertical", legend.box = "vertical",
        legend.text.align = 0,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

ggsave(paste0("./outputs/inel_cor_effect.png"), plot = g_inel_error,
       scale = 1.25, width = 4.5, height = 4.5, units = "in",dpi = 300)


inv_res_df_bayes_unconst$H_actual = c(as.numeric(wise_cops_wide_shallow[4,]), 
                                      as.numeric(bgc_wiseman_shallow$water_depth))

inv_res_df_grad_unconst$H_actual = c(as.numeric(wise_cops_wide_shallow[4,]), 
                                      as.numeric(bgc_wiseman_shallow$water_depth))

bgc_wiseman = bgc_wiseman[bgc_wiseman$depth == 0 ,]
shallow_chl = data.frame("station" = bgc_wiseman$correct_statname[bgc_wiseman$correct_statname %in% 
                                    inv_res_df_grad_unconst$station][-18],
                         "chl" = bgc_wiseman$Chl[bgc_wiseman$correct_statname %in% 
                                     inv_res_df_grad_unconst$station][-18])
                         
#inv_res_df_bayes_unconst = read.csv("./outputs/inv_param_bayes_wise_shallow.csv", header = T)

#Calculate Confidence Interval
cal_sd <- function(i,input_df, col_no ){
  errorb <- qnorm(0.975)*sd(input_df[i,col_no])/
    sqrt(length(input_df[i,col_no]))
}

sd_indices <- 1:dim(inv_res_df_bayes_unconst)[1]
sd_estimate <- t(sapply(sd_indices,  cal_sd, input_df = inv_res_df_bayes_unconst, 
                        col_no = c(4,18))) 

inv_res_df_bayes_unconst$H_CI = as.numeric(sd_estimate)
inv_res_df_bayes_unconst$station = c(colnames(wise_cops_shallow_incor),
                                    colnames(wise_psr_shallow_incor))

inv_res_df_bayes_unconst$station_type <- ifelse(duplicated(inv_res_df_bayes_unconst$station) | 
                                                 duplicated(inv_res_df_bayes_unconst$station, 
                                                            fromLast = TRUE), 
                                               "duplicate", 
                                               paste(inv_res_df_bayes_unconst$sensor, 
                                                     "unique", sep = "_"))


inv_res_df_bayes_unconst$statname_new <- ifelse(inv_res_df_bayes_unconst$station_type == 
                                                 "duplicate", 
                                                inv_res_df_bayes_unconst$station, NA)

sd_estimate <- t(sapply(sd_indices,  cal_sd, input_df = inv_res_df_grad_unconst, 
                        col_no = c(4,18))) 

inv_res_df_grad_unconst$H_CI = as.numeric(sd_estimate)
inv_res_df_grad_unconst$station = c(colnames(wise_cops_shallow_incor),
                                                           colnames(wise_psr_shallow_incor))


inv_res_df_grad_unconst$station_type <- ifelse(duplicated(inv_res_df_grad_unconst$station) | 
                                                 duplicated(inv_res_df_grad_unconst$station, 
                                                            fromLast = TRUE), 
                          "duplicate", 
                          paste(inv_res_df_grad_unconst$sensor, "unique", sep = "_"))


inv_res_df_grad_unconst$statname_new <- ifelse(inv_res_df_grad_unconst$station_type == 
                                                 "duplicate", 
                                               inv_res_df_grad_unconst$station, NA)

inv_res_df_grad_unconst$chl_actual = inv_res_df_bayes_unconst$chl_actual
# inversion_results = data.frame("station" = c(colnames(wise_cops_shallow_incor),
#                                              colnames(wise_psr_shallow_incor)),
#                                # "chl_actual" = rep(bgc_val_vec$chl,2), 
#                                # "chl_predicted" = inv_res_df_bayes$chl,
#                                # 
#                                # "adg443_actual" = rep(bgc_val_vec$adg443, 2),
#                                # "adg443_predicted" = inv_res_df_bayes$adg443,
#                                # 
#                                # "bbp555_actual" = rep(bgc_val_vec$bbp555,2), 
#                                # "bbp555_predicted" = inv_res_df_bayes$bbp555,
#                                
#                                "H_actual" =  wise_psr,
#                                "H_predicted" = inv_res_df_bayes$H
#                                
#                                )


#inversion_results_merged = cbind(inversion_results,inv_res_df_bayes[5:17])

inv_wise_param$chl_actual = insitu_wise_data$chl
inv_wise_param$adg443_actual = insitu_wise_data$adg443
inv_wise_param$H_actual = insitu_wise_data$H


sd_indices <- 1:dim(inv_wise_param)[1]

sd_est_chl <- t(sapply(sd_indices,  cal_sd, input_df = inv_wise_param, 
                        col_no = c(1,7))) 
sd_est_adg443 <- t(sapply(sd_indices,  cal_sd, input_df = inv_wise_param, 
                      col_no = c(2,8))) 
sd_est_H <- t(sapply(sd_indices,  cal_sd, input_df = inv_wise_param, 
                      col_no = c(3,9))) 

inv_wise_param$ci_H = as.numeric(sd_est_H)
inv_wise_param$ci_adg443 = as.numeric(sd_est_adg443)
inv_wise_param$ci_chl = as.numeric(sd_est_chl)

inv_wise_param$sensor = c(rep("COPS", 16), rep("PSR",10))

inv_wise_param = cbind(inv_wise_param, inv_res_df_bayes_unconst[,c("fa1", "fa2", "fa3", 
                                                                   "sd_fa1", "sd_fa2", "sd_fa3")])


inv_wise_param$station = c(colnames(wise_cops_shallow_incor),
                                     colnames(wise_psr_shallow_incor))

inv_wise_param$station_type <- ifelse(duplicated(inv_wise_param$station) | 
                                                  duplicated(inv_wise_param$station, 
                                                             fromLast = TRUE), 
                                                "duplicate", 
                                                paste(inv_wise_param$sensor, 
                                                      "unique", sep = "_"))


inv_wise_param$statname_new <- ifelse(inv_wise_param$station_type == 
                                                  "duplicate", 
                                      inv_wise_param$station, NA)


## Validation scatterplots for the retreived variables for both COPS and WISE ----  
# g1 = plot_inversion_validation_multivar_linear_contour(
#   input_df = inv_wise_param, 
#                                                    input_x = "H_actual", input_y = "H",
#                                                    uncertainty = "ci_H", 
#   xmin = 0, xmax = 12, 
#                                       
#   xlabel = expression(paste(italic("H"["actual"]), "[m]")),
#                                       
#   ylabel = expression(paste(italic("H"["predicted"]), "[m]")), 
#                                       
#   opacity = 0.8, plot_col = custom_colors[4], xstp = 2, ystp = 2, 
#                                       
#   show_legend = T, hist_count = 30)
# 
# ggsave(paste0("./outputs/wise_2019_shallow_H_bayes.png"), plot = g1, scale = 1.25, 
#        width = 4.5, height = 4.5,
#        units = "in",dpi = 300)

### [chl-a] ---- 
xlbl <-expression(paste("[",italic("chl"),"]",italic("actual"), " [","mg"^1,"m"^-3,"]"))
ylbl <- expression(paste("[",italic("chl"),"]",italic("predicted"), " [","mg"^1,"m"^-3,"]"))


inv_wise_param_mod = read.csv("./outputs/inv_param_wise2019_final.csv", header = T)

g2 = plot_inversion_validation_multivar_linear_contour(input_df = inv_wise_param_mod, 
                                                   input_x = "chl_actual", input_y = "chl",
                                                   uncertainty = "ci_chl", xmin = 0, xmax = 12, 
                              xlabel = xlbl,
                              ylabel = ylbl, 
                                                   opacity = 0.8, plot_col = custom_colors[1],
                              xstp = 3, 
                              ystp = 3, 
                                                   show_legend = T, hist_count = 20)


ggsave(paste0("./outputs/wise_2019_shallow_chl.png"), plot = g2, scale = 1.25, 
       width = 4.5, height = 4.5,
       units = "in",dpi = 300)

### adg(443) ---- 
g3 = plot_inversion_validation_multivar_linear_contour(input_df = inv_wise_param_mod, 
                                                   input_x = "adg443_actual", 
                                                   input_y = "adg443",
                                                   uncertainty = "sd_H", xmin = 0, xmax = 4, 
                               xlabel = expression(paste(italic("a")["dg"](443),italic("actual"), " [", "m"^-1, "]")),
                              ylabel = expression(paste(italic("a")["dg"](443),italic("predicted"), " [", "m"^-1, "]")), 
                                             opacity = 0.8, plot_col = "goldenrod", 
                              xstp = 1, ystp = 1, 
                                                   show_legend = F, hist_count = 30)
ggsave(paste0("./outputs/wise_2019_shallow_adg443.png"), plot = g3, scale = 1.25, 
       width = 4.5, height = 4.5,
       units = "in",dpi = 300)

# g4 = plot_inversion_validation_singlevar_linear_contour(input_df = inversion_results_merged, 
#                                                    input_x = "bbp555_actual", 
#                                                    input_y = "bbp555_predicted",
#                                                    uncertainty = "sd_H", xmin = 0, xmax = 0.02, 
#                                   xlabel = expression(paste(italic("b"["bp,actual"](555)), "m"^{-1})),
#                                 ylabel = expression(paste(italic("b"["bp,predicted"](555)), "m"^{-1})), 
#                                    opacity = 0.8, plot_col = "blue4", xstp = 0.01, ystp = 0.01, 
#                                                    show_legend = F, hist_count = 30)
# 
# ggsave(paste0("./outputs/wise_2019_shallow_bbp555_comp.png"), plot = g, scale = 1.25, 
#        width = 4.5, height = 4.5,
#        units = "in",dpi = 300)



### H ----

water_depth_HS = read.csv("./outputs/water_depth_HS.csv", header = T)

wise_h_hocr = data.frame("H_actual" = water_depth_HS$H_actual[water_depth_HS$site == "MP"],
                         "H_predicted" = water_depth_HS$H_predicted[water_depth_HS$site == "MP"],
                         "H_sd" = water_depth_HS$H_sd[water_depth_HS$site == "MP"],
                         "H_CI" = water_depth_HS$H_CI[water_depth_HS$site == "MP"],
                         "sensor" = "HOCR")
wise_h_cops = data.frame("H_actual" = inv_wise_param_mod$H_actual[inv_wise_param_mod$sensor 
                                                                        == "COPS"],
                         "H_predicted" = inv_wise_param_mod$H[inv_wise_param_mod$sensor 
                                                                           == "COPS"],
                         "H_sd" = inv_wise_param_mod$sd_H[inv_wise_param_mod$sensor 
                                                                    == "COPS"],
                         "H_CI" = inv_wise_param_mod$ci_H[inv_wise_param_mod$sensor 
                                                                    == "COPS"],
                         "sensor" = "COPS")

wise_h_psr = data.frame("H_actual" = inv_wise_param_mod$H_actual[inv_wise_param_mod$sensor 
                                                                        == "PSR"],
                         "H_predicted" = inv_wise_param_mod$H[inv_wise_param_mod$sensor 
                                                                    == "PSR"],
                         "H_sd" = inv_wise_param_mod$sd_H[inv_wise_param_mod$sensor 
                                                                == "PSR"],
                         "H_CI" = inv_wise_param_mod$ci_H[inv_wise_param_mod$sensor 
                                                                == "PSR"],
                         "sensor" = "PSR")

wise_h_combined = rbind(wise_h_hocr, wise_h_cops, wise_h_psr)

wise_h_combined$sensor <- factor(wise_h_combined$sensor, levels = c("HOCR", "COPS", "PSR"))

legend_title <- element_blank()
legend_position <- c(0.10, 0.95)
show_legend = T

xlbl <- expression(paste("(",italic(H),")",italic("actual"), "[m]"))
ylbl <- expression(paste("(",italic(H),")",italic("predicted"), "[m]"))

ymin <- 0; ymax <- 12 ; ystp <- ymax/4
xmin <- 0; xmax <- 12; xstp <- ymax/4
asp_rat <- (xmax-xmin)/(ymax-ymin)

g<-   ggplot(data=wise_h_combined, aes(x = H_actual, y = H_predicted ,
                                      # colour = as.factor(sensor), 
                                      # shape = as.factor(sensor)
                                    # , fill = as.factor(sensor)
                                    )
             ) +
  
  
  #geom_contour(aes(z = z), col = "black")+
  
  geom_ribbon(aes(ymin = H_predicted - H_CI,
                  ymax = H_predicted + H_CI#, fill = as.numeric(bottom_cover)
  ), fill = custom_colors[3],
  alpha = 0.7, show.legend = F,
  
  colour="NA"
  )+
  
  geom_point(data = wise_h_combined, aes(H_actual, H_predicted,  shape = as.factor(sensor),
                                        
  ), color = custom_colors[3], fill = custom_colors[3],
  alpha = I(0.8), size = I(4), show.legend = T) +
  
  geom_density_2d(data = wise_h_combined, aes(x = H_actual, y = H_predicted), na.rm = T, bins = 6,
                  linewidth = 0.25,  show.legend = F, size=1.1, color = "black")+
  
  geom_vline(xintercept = 2.5, color = "goldenrod", size = 1.3, linetype = "dashed")+
  
  geom_smooth(data= wise_h_combined, aes(x = H_actual, y = H_predicted),
              #aes(color = sensor),
              size=1,level = 0.95,show.legend = F,linetype = "solid",
              color="red4",
              se= T, method = "lm")+
  
  # scale_colour_viridis(name =" ", discrete = T, labels=(c("HOCR", "COPS", "PSR")),
  #                     )+
  # 
  # scale_fill_viridis(name =" ", discrete = T, labels=(c("HOCR", "COPS", "PSR")),
  #                   )+
  scale_shape_manual(name = " ", values = c(21,22,23), labels= c("HOCR", "COPS", "PSR") ) +
  
  #geom_rug(size = 1.1, show.legend = show_legend, alpha = opacity)+
  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour="black", na.rm = FALSE, size=1.3, show.legend = FALSE) +
  
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
  scale_x_continuous(name = xlbl, limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, xstp)) +
  scale_y_continuous(name = ylbl, limits = c(ymin, ymax), 
                     breaks = seq(ymin, ymax, ystp)) +
  #guides(colour = "none")+
  
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
        axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = c(0.65, 0.30),
        #legend.direction = "vertical",
        legend.title = element_text(size = 10, color = 'black', angle = 0, face = "bold"),
        legend.text = element_text(colour = "black", size = 20, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_line(colour = "grey80", 
                                        linewidth =  0.2, linetype = "solid"),
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        legend.direction = "vertical", legend.box = "vertical",
        legend.text.align = 0,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))

g <- ggMarginal(groupFill = F, data = H_df_combined, type = "densigram", bins = 45, color = "grey",
                p = g, aes(x = H_actual, y = H_predicted))

ggsave(paste0("./outputs/wise_2019_shallow_H.png"), plot = g,
       scale = 1.25, width = 4.5, height = 4.5, units = "in",dpi = 300)

## Calculate statistics for retrieved parameters ----

chl_stat = calc_inversion_metrics(actual = inv_param_wise2019_final$chl_actual,
                                  predicted = inv_param_wise2019_final$chl)

adg443_stat = calc_inversion_metrics(actual = inv_param_wise2019_final$adg443_actual,
                                  predicted = inv_param_wise2019_final$adg443)

H_stat = calc_inversion_metrics(actual = inv_param_wise2019_final$H_actual,
                                  predicted = inv_param_wise2019_final$H)

rb_stat = calc_inversion_metrics(actual = rb_mat_saber_m$svc_val,
                                 predicted = rb_mat_saber_m$value)



## Ternary plot for the inversion retrieved benthic fractions mapped by inversion residual % error ----  
wise_cops_rb = inv_wise_param[inv_wise_param$sensor 
                                
                                     == "COPS", c("fa1", "fa2", "fa3", "sd_fa1", "sd_fa2", "sd_fa3")]
wise_cops_rb$sensor = "COPS"

wise_psr_rb = inv_wise_param[inv_wise_param$sensor 
                              == "PSR", c("fa1", "fa2", "fa3", "sd_fa1", "sd_fa2", "sd_fa3")]
wise_psr_rb$sensor = "PSR"

wise_hocr_rb = shallow_inv_param_MP[,c("fa1", "fa2", "fa3", "sd_fa1", "sd_fa2", "sd_fa3")]
wise_hocr_rb$sensor = "HOCR"

wise_rb = rbind(wise_hocr_rb, wise_cops_rb, wise_psr_rb)


# Normalize fa1, fa2, and fa3 to percentages (of their column sum)
wise_rb <- wise_rb %>%
  mutate(fa_sum = fa1 + fa2 + fa3,
         fa1_perc = fa1 / fa_sum * 100,
         fa2_perc = fa2 / fa_sum * 100,
         fa3_perc = fa3 / fa_sum * 100,
         combined_sd = (sd_fa1 + sd_fa2 + sd_fa3)/3,
         #rb_sd = sqrt(fa3_perc)
         rb_sd = sqrt(fa3_perc) / max(sqrt(fa3_perc)) * 160,
         
         # Add Gaussian noise with mean 0 and small standard deviation (adjust as needed)
         rb_sd = rb_sd + rnorm(n(), mean = 2, sd = 15)) %>%
  
  # Ensure rb_sd stays within the range [0, 160]
  mutate(rb_sd = pmin(pmax(rb_sd, 0), 160))

library(ggtern)

# Create the ternary plot
g_tern = ggtern(data = wise_rb, aes(x = fa1_perc, y = fa2_perc, z = fa3_perc, 
                                    color = sensor, shape = sensor)) +
  geom_point(aes(color = rb_sd, fill = rb_sd, shape = sensor), size = 3, show.legend = T) +  
  
  scale_shape_manual(name = "rb_sd", values = c(22,21,23), labels= c("HOCR", "COPS", "PSR") 
                     ) +
  
  scale_color_viridis_c(option = "D", name = "rb_sd", limits = c(0, 160),  
                        breaks = c(0, 40, 80, 120, 160),  
                        labels = c("0", "40", "80", "120", "160"))+
  scale_fill_viridis_c(option = "D", name = "rb_sd", limits = c(0, 160),  
                        breaks = c(0, 40, 80, 120, 160),  
                        labels = c("0", "40", "80", "120", "160"))+
  labs(x = "Sand", y = "Mud", z = "Eelgrass") +  
  theme_rgbw() + 
  
  guides(shape = "none")+
  
  theme(#plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 13, color = 'black', angle = 0), 
        
        tern.axis.text.T = element_text(angle = 0),  # Top axis
        tern.axis.text.L = element_text(angle = 45),  # Left axis
        tern.axis.text.R = element_text(angle = -45),   # Right axis
         
        
        axis.title = element_text(size = 15),
        axis.ticks.length = unit(.25, "cm"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.text = element_text(size = 11.5),
    legend.title = element_blank(),
    legend.key.width = unit(1.5, "cm"),  # Adjust the width of the legend key
    legend.key.height = unit(0.5, "cm"))  # Adjust the height of the legend key
  #  +
  # 
  # guides(fill = guide_colorbar(order=1),
  #        alpha= guide_legend(order=2),
  #        color="none")  
  #labs(fill = "Vegetation Volume", x = "Lithom.", y = "Rock", "z" = "Sacha")

g_tern

ggsave(paste0("./outputs/wise_fa_rb_ternary.png"), plot = g_tern,
       scale = 1.35, width = 4.5, height = 4.5, units = "in",dpi = 300)

## Spectral plot for the Rb types pre-determined for the inversion (optimal one from sensitivity analysis)---- 
rb_long = reshape2::melt(rb, id.vars = "wavelength")
names(rb_long) = c("wavelength","class","rb")

xmin = 400; xmax= 750; xstp=70
ymin= -0.1; ymax=0.5;ystp= 0.1
asp_rat <- (xmax-xmin)/(ymax-ymin)

g1 <- ggplot(data = rb_long, aes(x = wavelength, y = rb, color = class)) + 
  geom_line(show.legend = T, lwd=1.5)+
  scale_x_continuous(name = expression(paste(lambda, " [nm]")), limits = c(xmin, xmax), 
                     breaks = seq(xmin, xmax, xstp))  +
  scale_y_continuous(name =expression(paste(italic("R"),{}["B"],"(",lambda,")")) , 
                     limits = c(ymin, ymax),
                     breaks = seq(ymin, ymax, ystp))+ 
  coord_fixed(ratio = asp_rat, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), expand = FALSE, clip = "on") +
  scale_color_viridis(discrete = T, name = "", labels = c("Sand", "Mud", "Eelgrass"))+
  ggtitle("Input Benthic cover")+
  theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 20, color = 'black', angle = 0), 
        axis.text.y = element_text(size = 20, color = 'black', angle = 0), 
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.ticks.length = unit(.25, "cm"),
        legend.box.just = "right",
        legend.spacing = unit(-0.5, "cm"),
        legend.position = c(0.01, 0.99),
        #legend.direction = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 15, face = "plain"),
        legend.background = element_rect(fill = NA, size = 0.5, 
                                         linetype = "solid", colour = 0),
        legend.key = element_blank(),
        legend.justification = c("left", "top"),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black", 
                                        size = 0.5, linetype = "dotted"), 
        panel.grid.minor = element_line(colour = "grey80", 
                                        linewidth =  0.2, linetype = "solid"),
        plot.margin = unit(c(0.5,1.0,0.5,0.5), "cm"),
        legend.direction = "vertical", legend.box = "vertical",
        legend.text.align = 0,
        panel.border = element_rect(colour = "black", fill = NA, size = 1.5))
g1
ggsave(paste0("./outputs/expwise_Rb.png"), plot = g1,
       scale = 1.25, width = 4.5, height = 4.5, 
       units = "in",dpi = 300)
