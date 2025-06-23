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

memory.limit(30000000)     # this is needed on some PCs to increase memory allowance, but has no impact on macs.

## ---------------------------

## load up the packages we will need:  (uncomment as required)

require(dplyr)

## ---------------------------

## load up our functions into memory

#source("Summary metrics/activity_index.R") 

## ---------------------------


identify_activity_bouts = function(activity_indices){
  
  # Reduce variance to binary activity/no activity
  activity_status = as.numeric(activity_indices > 0) # Return a seq like 000111110000111011100001
  
  # Identify runs of activity
  activity_rle = rle(activity_status)
  
  # Collate outcomes
  activity_bouts = data.frame(
    bout_length = activity_rle$lengths,
    bout_status = activity_rle$values
  )
  
  # Flag indices of acticity bouts in original data 
  activity_bouts$index = cumsum(lag(activity_bouts$bout_length, default = 1))
  
  return(activity_bouts)
  
}
