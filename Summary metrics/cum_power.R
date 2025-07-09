## ---------------------------
##
## Script name: 
##
## Purpose of script: Estimate cumulative spectral power based on Welch's method
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

require(gsignal, exclude = "filter")

## ---------------------------

## load up our functions into memory

# source("functions/summarise_data.R") 

## ---------------------------

cum_power = function(filtered_signal, freq_min = 0.1, freq_max = 5, fs = 50, window_length = 1000, show_plot = F, proportion = T) {
  
  ### Estimates cumulative power based on Welch's power spectral desnity
  ### Recommend smoothing with bandpass to remove gravity and high-freq
  ### fs needs to be equal to the sampling frequency
  ### *** need to work out an ideal window length for the 24hr measures ***
  
  if(freq_max > (fs/2)){
    stop("Maximum frequency (freq_max) must be below half of sampling freq (fs)")
  }
  
  psd = gsignal::pwelch(x = filtered_signal, fs = 50, 
                        window = 1000, detrend = "short-mean")
  
  spec_data = data.frame(
    freq = psd$freq,
    spec = psd$spec
  )
  
  spec_data$power = cumsum(spec_data$spec)
  
  if(show_plot){
    plot(spec_data$freq, spec_data$power)
  }
  
  if(proportion){  
    spec_data$power = spec_data$power / max(spec_data$power)
  }
    
    
  # Identify the start and end range of interest using a fuzzy search - should make a function for this
  spec_data = spec_data[min(which((spec_data$freq - freq_min) >= 0)) : max(which((spec_data$freq - freq_max) <= 0)), ]
  
  return(max(spec_data$power))
  
}













