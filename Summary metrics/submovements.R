## ---------------------------
##
## Script name: 
##
## Purpose of script:
##
## Author: Dr. Cory Holdom
##
## Date Created: 2025-06-23
##
## Copyright (c) Cory Holdom, 2025
## Email: cory.holdom@ndcn.ox.ac.uk
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)

require(dplyr)
require(FactoMineR)

## ---------------------------

## load up our functions into memory

# source("functions/summarise_data.R") 

## ---------------------------

clean_acc = function(raw_signal, filter_order = 4, fs = 100, freq_cutoff = c(0.1, 20)){
  
  # Estimate vector magnitude and direction 
  # Filter data with a band-pass Butterworth filter to smooth motion
  # and remove gravity component
  
  w_cutoff = freq_cutoff / (fs / 2)
  
  butt_filt = gsignal::butter(filter_order, w = freq_cutoff / (fs / 2), type = "pass")
  
  clean_signal = gsignal::filtfilt(butt_filt, raw_signal)
  
  return(clean_signal)
  
}

calculate_velocity = function(filtered_accelerations, fs = 100){
  
  # Expects signal to have been (bandpass) filtered
  
  col_names = names(filtered_accelerations)
  
  if(!("time"         %in% col_names)) stop("Missing time column in mapped bouts")
  if(!("x"            %in% col_names)) stop("Missing x column in mapped bouts")
  if(!("y"            %in% col_names)) stop("Missing y column in mapped bouts")
  if(!("z"            %in% col_names)) stop("Missing z column in mapped bouts")
  if(!("bout_status"  %in% col_names)) stop("Missing bout_status column in mapped bouts")
  if(!("bout_index"   %in% col_names)) stop("Missing bout_index column in mapped bouts")
  if(!("x_f"          %in% col_names)) stop("Missing x_f column in mapped bouts - make sure to run clean_acc")
  if(!("y_f"          %in% col_names)) stop("Missing y_f column in mapped bouts - make sure to run clean_acc")
  if(!("z_f"          %in% col_names)) stop("Missing z_f column in mapped bouts - make sure to run clean_acc")
  
  
  estimated_velocities = filtered_accelerations |>
    select(time, bout_index, x_f:z_f)
  
  estimated_velocities = estimated_velocities |>
    group_by(bout_index) |>
    mutate(vel_x = if_else(time - min(time) < 0.0001, 0, cumsum(x_f * 0.001)),
           vel_y = if_else(time - min(time) < 0.0001, 0, cumsum(y_f * 0.001)),
           vel_z = if_else(time - min(time) < 0.0001, 0, cumsum(z_f * 0.001))) |>
    ungroup() |>
    select(vel_x, vel_y, vel_z)
    
  return(estimated_velocities)
  
}


project_velocity = function(data_with_velocities){
  
  # Identify principle axes of movement using PCA over the triaxial velocities
  # Returns velocities within each of the three principle axes of movement *during each bout
  
  col_names = names(data_with_velocities)
  
  if(!("time"         %in% col_names)) stop("Missing time column in mapped bouts")
  if(!("x"            %in% col_names)) stop("Missing x column in mapped bouts")
  if(!("y"            %in% col_names)) stop("Missing y column in mapped bouts")
  if(!("z"            %in% col_names)) stop("Missing z column in mapped bouts")
  if(!("bout_status"  %in% col_names)) stop("Missing bout_status column in mapped bouts")
  if(!("bout_index"   %in% col_names)) stop("Missing bout_index column in mapped bouts")
  if(!("x_f"          %in% col_names)) stop("Missing x_f column in mapped bouts - make sure to run clean_acc")
  if(!("y_f"          %in% col_names)) stop("Missing y_f column in mapped bouts - make sure to run clean_acc")
  if(!("z_f"          %in% col_names)) stop("Missing z_f column in mapped bouts - make sure to run clean_acc")
  if(!("vel_x"        %in% col_names)) stop("Missing vel_x column in mapped bouts - make sure to run calculate_velocity")
  if(!("vel_y"        %in% col_names)) stop("Missing vel_y column in mapped bouts - make sure to run calculate_velocity")
  if(!("vel_z"        %in% col_names)) stop("Missing vel_z column in mapped bouts - make sure to run calculate_velocity")
  
  pcs = data_with_velocities |>
    group_by(bout_index) |>
    select(bout_index, vel_x:vel_z) |>
    group_map(~ FactoMineR::PCA(.x, graph = F)$ind$coord)
  
  est_pcs = data.table::rbindlist(lapply(pcs, as.data.frame))
  
  return(est_pcs)

}



identify_submovement_boundaries = function(data_with_velocities) {
  
  # Submovements are periods of movement (+ or -ve velocities) flanked by zero-velocity crossings
  # Thus, a submovement boundary can be found where the sign of the velocity changes.
  # As with the activity bouts, the submovements themselves can be identified using run-length encodings of the sign of the velocity
  
  col_names = names(data_with_velocities)
  
  if(!("time"         %in% col_names)) stop("Missing time column in mapped bouts")
  if(!("x"            %in% col_names)) stop("Missing x column in mapped bouts")
  if(!("y"            %in% col_names)) stop("Missing y column in mapped bouts")
  if(!("z"            %in% col_names)) stop("Missing z column in mapped bouts")
  if(!("bout_status"  %in% col_names)) stop("Missing bout_status column in mapped bouts")
  if(!("bout_index"   %in% col_names)) stop("Missing bout_index column in mapped bouts")
  if(!("x_f"          %in% col_names)) stop("Missing x_f column in mapped bouts - make sure to run clean_acc")
  if(!("y_f"          %in% col_names)) stop("Missing y_f column in mapped bouts - make sure to run clean_acc")
  if(!("z_f"          %in% col_names)) stop("Missing z_f column in mapped bouts - make sure to run clean_acc")
  if(!("vel_x"        %in% col_names)) stop("Missing vel_x column in mapped bouts - make sure to run calculate_velocity")
  if(!("vel_y"        %in% col_names)) stop("Missing vel_y column in mapped bouts - make sure to run calculate_velocity")
  if(!("vel_z"        %in% col_names)) stop("Missing vel_z column in mapped bouts - make sure to run calculate_velocity")
  if(!("vel_pc1"      %in% col_names)) stop("Missing vel_pc1 column in mapped bouts - make sure to run project_velocity")
  if(!("vel_pc2"      %in% col_names)) stop("Missing vel_pc2 column in mapped bouts - make sure to run project_velocity")
  if(!("vel_pc3"      %in% col_names)) stop("Missing vel_pc3 column in mapped bouts - make sure to run project_velocity")
  
  velocity_signs = sign(data_with_velocities[, c("vel_pc1", "vel_pc2", "vel_pc3")])
  
  flanked_vels = apply(sign(data_with_velocities[, c("vel_pc1", "vel_pc2", "vel_pc3")]), 2, rle)
  
  #print(flanked_vels)
  
  submovement_bouts = list(
    vel_pc1 = data.frame(
      submovement_length = flanked_vels$vel_pc1$lengths,
      bout_status = flanked_vels$vel_pc1$values
    ),
    
    vel_pc2 = data.frame(
      submovement_length = flanked_vels$vel_pc2$lengths,
      bout_status = flanked_vels$vel_pc2$values
    ),
    
    vel_pc3 = data.frame(
      submovement_length = flanked_vels$vel_pc3$lengths,
      bout_status = flanked_vels$vel_pc3$values
    )
  )
  
  print(head(submovement_bouts$vel_pc1))
  
  submovement_bouts$vel_pc1$index = cumsum(dplyr::lag(submovement_bouts$vel_pc1$submovement_length, default = 1))
  submovement_bouts$vel_pc2$index = cumsum(dplyr::lag(submovement_bouts$vel_pc2$submovement_length, default = 1))
  submovement_bouts$vel_pc3$index = cumsum(dplyr::lag(submovement_bouts$vel_pc3$submovement_length, default = 1))
  
  submovement_bouts$vel_pc1 = submovement_bouts$vel_pc1[1:(nrow(submovement_bouts$vel_pc1)-1),]
  submovement_bouts$vel_pc2 = submovement_bouts$vel_pc2[1:(nrow(submovement_bouts$vel_pc2)-1),]
  submovement_bouts$vel_pc3 = submovement_bouts$vel_pc3[1:(nrow(submovement_bouts$vel_pc3)-1),]
  
  rm(velocity_signs, flanked_vels)
  
  return(submovement_bouts)
  
}



map_submovements = function(){
  
  # Goal is to map submovements for each PC back onto original dataset to
  # have labelled submovements within labelled activity bouts
  
  
  bout_status = rep(0, nrow(raw_acc))
  
  origin = as.POSIXct(raw_acc[[time_col]][1])
  
  print(paste("Recording origin: ", origin))
  
  b_indices = bout_indices[[bout_index_col]]
  
  print(paste("Length of index array :", length(b_indices)))
  
  for(i in 1:length(b_indices)){
    
    if(i < length(b_indices)){
      
      t_start = (b_indices[i] - 1) * 100 + 1
      t_end = (b_indices[i+1] - 1) * 100
      
      bout_status[t_start:t_end] = bout_indices[[bout_status_col]][i]
      
    }
    
  }
  
  return(bout_status)
  
}












# 
# dist_xy = dtw::dtw(x = flanked_vels$vel_x$lengths, y = flanked_vels$vel_y$lengths)
# dist_xz = dtw::dtw(x = flanked_vels$vel_x$lengths, y = flanked_vels$vel_z$lengths)
# dist_yz = dtw::dtw(x = flanked_vels$vel_y$lengths, y = flanked_vels$vel_z$lengths)
# 
# dist_xy$distance / max(length(flanked_vels$vel_x$lengths), length(flanked_vels$vel_y$lengths))
# dist_xz$distance / max(length(flanked_vels$vel_x$lengths), length(flanked_vels$vel_z$lengths))
# dist_yz$distance / max(length(flanked_vels$vel_y$lengths), length(flanked_vels$vel_z$lengths))
# 
# dtw_pc1_x = dtw::dtw(x = bout_pca_res$Dim.1, y = bout_pca_res$vel_x)
# 
# dtw_pc1_x$normalizedDistance
# 
# dist_yz$normalizedDistance



# gsignal::filtfilt(butt_filt, v_data$x)
# 
# 
# 
# clean_acc(v_data$x)
# 
# Data$data$x_f = 0
# Data$data$y_f = 0
# Data$data$z_f = 0
# 
# rbind(Data$data, clean_acc(Data$data))
# 
# 
# 
# v_data = calculate_velocity(Data$data)
# 
# v_data |>
#   filter(bout_index %in% seq(1,17,2)) |>
#   mutate(p_time = (int_time - first(int_time)) / max(int_time - first(int_time))) |>
#   ungroup() |>
#   select(bout_index, p_time, x, vel_x) |>
#   ggplot(aes(x = p_time)) +
#   theme_bw() +
#   geom_line(aes(y = x), colour = "#228b22") +
#   geom_line(aes(y = vel_x)) +
#   facet_wrap(~ bout_index)


# 
# head(submovement_bouts)
# 
# hist(submovement_bouts$bout_length[submovement_bouts$bout_length > 2])
# 
# 
# ggplot(submovement_bouts) +
#   geom_histogram(aes(x = bout_length)) +
#   scale_x_log10() +
#   scale_y_log10()
# 
# 
# 






