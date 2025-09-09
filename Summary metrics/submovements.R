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

# clean_acc = function(raw_signal, filter_order = 4, fs = 100, freq_cutoff = c(0.1, 20)){
#   
#   # Estimate vector magnitude and direction 
#   # Filter data with a band-pass Butterworth filter to smooth motion
#   # and remove gravity component
#   
#   w_cutoff = freq_cutoff / (fs / 2)
#   
#   butt_filt = gsignal::butter(filter_order, w = freq_cutoff / (fs / 2), type = "pass")
#   
#   clean_signal = gsignal::filtfilt(butt_filt, raw_signal)
#   
#   return(clean_signal)
#   
# }

clean_acc = function(raw_signals, filter_order = 4, fs = 100, freq_cutoff = c(0.1, 20)){
  
  # Estimate vector magnitude and direction 
  # Filter data with a band-pass Butterworth filter to smooth motion
  # and remove gravity component
  
  w_cutoff = freq_cutoff / (fs / 2)
  
  butt_filt = gsignal::butter(filter_order, w = freq_cutoff / (fs / 2), type = "pass")
  
  clean_signal = data.table::as.data.table(sapply(raw_signals, gsignal::filtfilt, filt = butt_filt))
  
  return(clean_signal)
  
}

# freq_cutoff = c(0.1, 10)
# fs = 50
# filter_order = 4
# butt_filt = gsignal::butter(filter_order, w = freq_cutoff / (fs / 2), type = "pass")
# clean_signal = gsignal::filtfilt(butt_filt, select(walk_data, userAcceleration.x:userAcceleration.z))
# 
# sapply(select(walk_data, userAcceleration.x), gsignal::filtfilt, filt = butt_filt)

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
    select(bout_index, vel_x:vel_z) |>
    nest(.by = bout_index) |>
    mutate(PComp = lapply(data, function(df) as.data.table(FactoMineR::PCA(df, graph = F)$ind$coord))) |>
    select(bout_index, PComp) |>
    unnest(PComp, .drop = T) |>
    select(-bout_index)
    
  return(pcs)

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
  
  #velocity_signs = sign(data_with_velocities[, c("vel_pc1", "vel_pc2", "vel_pc3")])
  
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
  
  #print(head(submovement_bouts$vel_pc1))
  
  submovement_bouts$vel_pc1$index = cumsum(dplyr::lag(submovement_bouts$vel_pc1$submovement_length, default = 1))
  submovement_bouts$vel_pc2$index = cumsum(dplyr::lag(submovement_bouts$vel_pc2$submovement_length, default = 1))
  submovement_bouts$vel_pc3$index = cumsum(dplyr::lag(submovement_bouts$vel_pc3$submovement_length, default = 1))
  
  # submovement_bouts$vel_pc1 = submovement_bouts$vel_pc1[1:(nrow(submovement_bouts$vel_pc1)-1),]
  # submovement_bouts$vel_pc2 = submovement_bouts$vel_pc2[1:(nrow(submovement_bouts$vel_pc2)-1),]
  # submovement_bouts$vel_pc3 = submovement_bouts$vel_pc3[1:(nrow(submovement_bouts$vel_pc3)-1),]
  
  submovement_bouts$vel_pc1 = submovement_bouts$vel_pc1[c(3,1,2)]
  submovement_bouts$vel_pc2 = submovement_bouts$vel_pc2[c(3,1,2)]
  submovement_bouts$vel_pc3 = submovement_bouts$vel_pc3[c(3,1,2)]
  
  return(submovement_bouts)
  
}


## Logic to map bouts back: for each index in submovement_bouts, assign the corresponding row in data$Data the 
# 
# submovement_bout_indices$vel_pc1$bout_status = 1:length(submovement_bout_indices$vel_pc1$bout_status)
# 
# rep(submovement_bout_indices$vel_pc1$bout_status, submovement_bout_indices$vel_pc1$submovement_length)
# 


map_submovement_bouts = function(data_with_velocities, submovement_bout_indices){
  
  submovement_bout_indices$vel_pc1$bout_status = 1:length(submovement_bout_indices$vel_pc1$bout_status)
  submovement_bout_indices$vel_pc2$bout_status = 1:length(submovement_bout_indices$vel_pc2$bout_status)
  submovement_bout_indices$vel_pc3$bout_status = 1:length(submovement_bout_indices$vel_pc3$bout_status)
  
  mapped_bouts = list(
    "vel_pc1_sm" = rep(submovement_bout_indices$vel_pc1$bout_status,
                       submovement_bout_indices$vel_pc1$submovement_length),
    
    "vel_pc2_sm" = rep(submovement_bout_indices$vel_pc2$bout_status,
                       submovement_bout_indices$vel_pc2$submovement_length),
    
    "vel_pc3_sm" = rep(submovement_bout_indices$vel_pc3$bout_status,
                       submovement_bout_indices$vel_pc3$submovement_length)
  )
  
  return(mapped_bouts)
  
}



label_submovements = function(data_with_submovements, short_sm_limits = c(0.05, 0.6), long_sm_limits = c(0.6, 5.0), fs = 100){
  
  col_names = names(data_with_submovements)
  
  if(!("time"            %in% col_names)) stop("Missing time column in mapped bouts")
  if(!("x"               %in% col_names)) stop("Missing x column in mapped bouts")
  if(!("y"               %in% col_names)) stop("Missing y column in mapped bouts")
  if(!("z"               %in% col_names)) stop("Missing z column in mapped bouts")
  if(!("bout_status"     %in% col_names)) stop("Missing bout_status column in mapped bouts")
  if(!("bout_index"      %in% col_names)) stop("Missing bout_index column in mapped bouts")
  if(!("x_f"             %in% col_names)) stop("Missing x_f column in mapped bouts - make sure to run clean_acc")
  if(!("y_f"             %in% col_names)) stop("Missing y_f column in mapped bouts - make sure to run clean_acc")
  if(!("z_f"             %in% col_names)) stop("Missing z_f column in mapped bouts - make sure to run clean_acc")
  if(!("vel_x"           %in% col_names)) stop("Missing vel_x column in mapped bouts - make sure to run calculate_velocity")
  if(!("vel_y"           %in% col_names)) stop("Missing vel_y column in mapped bouts - make sure to run calculate_velocity")
  if(!("vel_z"           %in% col_names)) stop("Missing vel_z column in mapped bouts - make sure to run calculate_velocity")
  if(!("vel_pc1"         %in% col_names)) stop("Missing vel_pc1 column in mapped bouts - make sure to run project_velocity")
  if(!("vel_pc2"         %in% col_names)) stop("Missing vel_pc2 column in mapped bouts - make sure to run project_velocity")
  if(!("vel_pc3"         %in% col_names)) stop("Missing vel_pc3 column in mapped bouts - make sure to run project_velocity")
  if(!("vel_pc1_sm"      %in% col_names)) stop("Missing vel_pc1_sm column in mapped bouts - make sure to run map_submovements")
  if(!("vel_pc2_sm"      %in% col_names)) stop("Missing vel_pc2_sm column in mapped bouts - make sure to run map_submovements")
  if(!("vel_pc3_sm"      %in% col_names)) stop("Missing vel_pc3_sm column in mapped bouts - make sure to run map_submovements")
  
  # PC1 submovement categorisation
  res = data_with_submovements |>
    mutate(vel_pc1_sm_type = case_when((n() > (short_sm_limits[1] * fs)) & (n() < (short_sm_limits[2] * fs)) ~ "Short",
                                       (n() > (long_sm_limits[1]  * fs)) & (n() < (long_sm_limits[2]  * fs)) ~ "Long"), .by = vel_pc1_sm) |>
    mutate(vel_pc2_sm_type = case_when((n() > (short_sm_limits[1] * fs)) & (n() < (short_sm_limits[2] * fs)) ~ "Short",
                                       (n() > (long_sm_limits[1]  * fs)) & (n() < (long_sm_limits[2]  * fs)) ~ "Long"), .by = vel_pc2_sm) |>
    mutate(vel_pc3_sm_type = case_when((n() > (short_sm_limits[1] * fs)) & (n() < (short_sm_limits[2] * fs)) ~ "Short",
                                       (n() > (long_sm_limits[1]  * fs)) & (n() < (long_sm_limits[2]  * fs)) ~ "Long"), .by = vel_pc3_sm)
  
  return(res)
  
}



summarise_submovement_features = function(data_with_submovements, filter_invalid_submovements = T, limits = list(), fs = 100){
  
  # Within each submovement (PC1 and PC2), calculate duration, distance, peak Vel
  #    - Duration = length (/fs)
  #    - Distance = cumsum of velocity
  #    - Peak Vel = max
  
  col_names = names(data_with_submovements)
  
  if(!("time"            %in% col_names)) stop("Missing time column in mapped bouts")
  if(!("x"               %in% col_names)) stop("Missing x column in mapped bouts")
  if(!("y"               %in% col_names)) stop("Missing y column in mapped bouts")
  if(!("z"               %in% col_names)) stop("Missing z column in mapped bouts")
  if(!("bout_status"     %in% col_names)) stop("Missing bout_status column in mapped bouts")
  if(!("bout_index"      %in% col_names)) stop("Missing bout_index column in mapped bouts")
  if(!("x_f"             %in% col_names)) stop("Missing x_f column in mapped bouts - make sure to run clean_acc")
  if(!("y_f"             %in% col_names)) stop("Missing y_f column in mapped bouts - make sure to run clean_acc")
  if(!("z_f"             %in% col_names)) stop("Missing z_f column in mapped bouts - make sure to run clean_acc")
  if(!("vel_x"           %in% col_names)) stop("Missing vel_x column in mapped bouts - make sure to run calculate_velocity")
  if(!("vel_y"           %in% col_names)) stop("Missing vel_y column in mapped bouts - make sure to run calculate_velocity")
  if(!("vel_z"           %in% col_names)) stop("Missing vel_z column in mapped bouts - make sure to run calculate_velocity")
  if(!("vel_pc1"         %in% col_names)) stop("Missing vel_pc1 column in mapped bouts - make sure to run project_velocity")
  if(!("vel_pc2"         %in% col_names)) stop("Missing vel_pc2 column in mapped bouts - make sure to run project_velocity")
  if(!("vel_pc3"         %in% col_names)) stop("Missing vel_pc3 column in mapped bouts - make sure to run project_velocity")
  if(!("vel_pc1_sm"      %in% col_names)) stop("Missing vel_pc1_sm column in mapped bouts - make sure to run map_submovements")
  if(!("vel_pc2_sm"      %in% col_names)) stop("Missing vel_pc2_sm column in mapped bouts - make sure to run map_submovements")
  if(!("vel_pc3_sm"      %in% col_names)) stop("Missing vel_pc3_sm column in mapped bouts - make sure to run map_submovements")
  if(!("vel_pc1_sm_type" %in% col_names)) stop("Missing vel_pc1_sm_type column in mapped bouts - make sure to run label_submovements")
  if(!("vel_pc2_sm_type" %in% col_names)) stop("Missing vel_pc2_sm_type column in mapped bouts - make sure to run label_submovements")
  if(!("vel_pc3_sm_type" %in% col_names)) stop("Missing vel_pc3_sm_type column in mapped bouts - make sure to run label_submovements")
  
  # if(length(limits) == 0){
  #   
  #   
  #   
  # }
  
  
  # PC1
  data_with_submovements |>
    group_by(vel_pc1_sm) |> # Need to group into short/long SMs per dimension
    summarise(duration  = n()/fs,
              distance  = sum(abs(vel_pc1))/fs,
              peak_vel  = max(abs(vel_pc1)),
              peak_acc  = max(abs(vel_pc1 - lag(vel_pc1, default = 0)))/fs,
              peak_jerk = max(abs(vel_pc1 - 2 * lag(vel_pc1, n = 1, default = 0) - lag(vel_pc1, n = 2, default = 0)))
              )
              
    
  #TODO: Add PC2/PC3
  
  
  
}


resample_submovement = function(submovement_velocities, length_out = 40, fs = 100){
  
  t = seq(1, length(submovement_velocities))/fs
  
  if(length(t) < 2){return(as.double(NA))}
  
  res = approx(t, submovement_velocities, n = length_out)
  
  return(res$y)
  
}


# Normalise submovement from 0-1
normalise_submovement_velocity = function(submovement_velocities){
  
  submovement_velocities = abs(submovement_velocities)
  
  return(submovement_velocities / max(submovement_velocities))
  
}


# Remove time-velocity curves that don't return to (near) 0
remove_extreme_velocity = function(velocity_curves, col_name, submovement_group_col_name, vel_threshold = 0.1){
  
  col_name = enquo(col_name)
  
  submovement_group_col_name = enquo(submovement_group_col_name)
  
  (velocity_curves |>
    mutate(exclude = (abs(first(!!col_name) - last(!!col_name)) > vel_threshold), .by = !!submovement_group_col_name) |>
    mutate(col_name = if_else(exclude, NA, !!col_name)))[["col_name"]]
  
  
}






