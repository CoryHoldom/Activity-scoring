## ---------------------------
##
## Script name: 
##
## Purpose of script: A compilation of functions used for estimating activity index based on Bai et al., 2016
##
## Author: Dr. Cory Holdom
##
## Date Created: 2025-06-10
## Date Updated: 2025-06-11
## 
## Copyright (c) Cory Holdom, 2025
## Email: cory.holdom@ndcn.ox.ac.uk
##
## ---------------------------
##
## Notes: 
##   Bai et al., 2016:
##   doi.org/10.1371/journal.pone.0160644  
##   
##
## ---------------------------

require(dplyr)
require(tibble)

# Activity Index: the average of the variances over the 3 axes, normalised by
# the device noise. For this measure, I've assumed the recordings during sitting
# are stationary (maybe an overly optimistic assumption, but the exploration of
# the data look alright). The AI is further transformed to original units (g),
# to be on the same scale as the original measure.
# Based on Bai et al., 2016: doi.org/10.1371/journal.pone.0160644

calculate_stationary_variance = function(stationary_signal){
  
  # The variance of the acceleration (in SD units, g^2) is measured by taking
  # the sum of the variance of the three axes at rest
  
  stationary_variance = sd(unlist(stationary_signal[,1], use.names = F)) +
    sd(unlist(stationary_signal[,2], use.names = F)) +
    sd(unlist(stationary_signal[,3], use.names = F))
  
  return(stationary_variance) # Sigma_hat^2 in Bai, J. et al., 2016
  
}


calculate_windowed_sd = function(accelerations, window_secs = 10, fs = 100){
  
  # Estimates standard deviation of accelerations over a series of chunks (default: 10s @ 100Hz)
  # This is NOT a rolling SD, it is a serial measure
  
  accelerations = as_tibble(accelerations)
  
  accelerations = accelerations |>
    mutate(rn = row_number()) |>
    group_by(grp = ceiling(rn / (window_secs * fs))) |>
    mutate(stdev = case_when(
      length(grp) == window_secs * fs ~ sd(value),
      TRUE ~ NA
    )) |>
    summarise(block_10s = first(grp), stdev = first(stdev)) |>
    ungroup()
  
  return(accelerations)

}

# test_data = sapply(Data$data[,2:4], FUN = calculate_windowed_sd, window_secs = 10, fs = 100)["stdev",] |>
#   bind_rows() |>
#   mutate(max_sd = x + y + z)
# 
# ggplot(test_data, aes(x = max_sd)) +
#   theme_bw() +
#   geom_histogram() +
#   scale_x_log10()
# 
# table(test_data$max_sd < 0.013) |> prop.table()


aggregate_signal = function(raw_signal, window_length, func = "sd"){
  
  # Group data into "bins" to aggregate over. E.g., if data are recorded at 50Hz
  # and desired output is 1-second summaries, window_length = 50
  
  raw_signal$bin = floor(row_number(raw_signal) / window_length)
  
  # Estimate summary metric (e.g., SD) over each bin
  
  aggregated_signal = aggregate(raw_signal[,2:4],
                                by = raw_signal["bin"],
                                FUN = func)
  
  return(aggregated_signal)
  
}

calculate_activity_index = function(aggregated_signal_variance, stationary_variance = 0, relative_activity_index = FALSE){
  # The sum of the SD in each axis is normalised by the device variance at rest
  # and then divided by 3 to get the average SD over the three axes. The result
  # is transformed with a square root to obtain an activity index on the same
  # units as the original signal.
  # 
  # If relative_activity_index is TRUE, the AI is normalised by the systemic
  # variance of the device
  
  
  
  if(!relative_activity_index){
    
    # print("Estimated AI is not scaled for systemic noise")
    
    mean_of_normalised_variances = (aggregated_signal_variance[2] +
                                      aggregated_signal_variance[3] +
                                      aggregated_signal_variance[4] -
                                      3 * stationary_variance) / 3 
    
    activity_index = max(mean_of_normalised_variances, 0) ^ 0.5
    
    
  }
  
  else{
    
    # print("Estimated AI is scaled for systemic noise")
    
    mean_of_normalised_variances = (aggregated_signal_variance[2] +
                                      aggregated_signal_variance[3] +
                                      aggregated_signal_variance[4] -
                                      3 * stationary_variance) / stationary_variance / 3
    
    activity_index = max(mean_of_normalised_variances, 0) ^ 0.5
    
  }
  
  return(activity_index)
  
}