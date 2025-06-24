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

## ---------------------------

## load up our functions into memory

#source("Summary metrics/activity_index.R") 

## ---------------------------


identify_activity_bouts = function(activity_indices){ #NEED TO FILTER OUT BAD BOUTS AFTER (not 4-18s)
  
  # Reduce variance to binary activity/no activity
  activity_status = as.numeric(activity_indices > 0) # Return a seq like 000111110000111011100001
  
  # Identify runs of activity
  activity_rle = rle(activity_status)
  
  # Collate outcomes
  activity_bouts = data.frame(
    bout_length = activity_rle$lengths,
    bout_status = activity_rle$values
  )
  
  # Flag indices of activity bouts in original data 
  activity_bouts$index = cumsum(dplyr::lag(activity_bouts$bout_length, default = 1))
  
  return(activity_bouts[c(3,1,2)])
  
  # ----------- # ----------- # ----------- #
  # bout_length # bout_status # index       #
  # ----------- # ----------- # ----------- #
  #     1       #     0       #     1       #
  #     5       #     1       #     2       #
  #     7       #     0       #     7       #
  
}


# maybe need an aggregate function idk



calculate_bout_features = function(activity_bout, fs = 50){ 
  
  # Daily measures involve mean and SD of bout acceleration and jerk
  
  b_acceleration = max(activity_bout)
  
  b_jerk = mean((activity_bout - dplyr::lag(activity_bout)) / fs, na.rm = T)
  
  bout_features = list(
    bout_acceleration = b_acceleration,
    bout_jerk = b_jerk
  )
  
  return(bout_features)
  
}

test = c(
  rep(0,8),
  runif(17),
  rep(0,3),
  runif(6)
)

abouts = identify_activity_bouts(test)



f = rep(abouts$index, abouts$bout_length)

splote = split(test, f)

cat = sapply(splote, calculate_bout_features)

as.data.frame(cat)
