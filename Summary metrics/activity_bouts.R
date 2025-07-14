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

source(here::here("find_subvec.R"))

## ---------------------------

identify_activity_bouts = function(activity_indices, report_short_bouts = F){ #NEED TO FILTER OUT BAD BOUTS AFTER (not 4-18s)
  
  # Reduce variance to binary activity/no activity
  activity_status = as.numeric(activity_indices > 0) # Return a seq like 000111110000111011100001
  
  # Replace "bouts" of 1s with inactivity to mirror previous works (doi.org/10.1007/s12311-022-01385-5)
  sub_vec = c(0,1,0)
  repl_vec = c(0,0,0)
  
  sub_indices = find_subvec(activity_status, sub_vec)
  
  for (i in sub_indices) activity_status[i + seq_along(sub_vec) - 1] = repl_vec
  
  # Identify runs of activity
  activity_rle = rle(activity_status)
  
  # Collate outcomes
  activity_bouts = data.frame(
    bout_length = activity_rle$lengths,
    bout_status = activity_rle$values
  )
  
  # Flag indices of activity bouts in original data 
  activity_bouts$index = cumsum(dplyr::lag(activity_bouts$bout_length, default = 1))
  
  activity_bouts = activity_bouts[1:(nrow(activity_bouts)-1),]
  
  if(report_short_bouts){
    results = list(
      activity_bouts = activity_bouts[c(3,1,2)],
      replaced_indices = sub_indices
    )
  }
  
  else{
    results = list(
      activity_bouts = activity_bouts[c(3,1,2)]
    )
  }
  
  return(results)
  
  # ----------- # ----------- # ----------- #
  #    index    # bout_length # bout_status #
  # ----------- # ----------- # ----------- #
  #     1       #      1      #      0      #
  #     2       #      5      #      1      #
  #     7       #      7      #      0      #
  
}


# maybe need an aggregate function idk


# Map activity bouts back onto original data

map_activity_bouts = function(raw_acc, bout_indices, time_col = "time", bout_index_col = "index", bout_length_col = "bout_length", bout_status_col = "bout_status") {
  
  # 
  
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

calculate_bout_features = function(mapped_bouts, fs = 100){ 
  
  # Daily measures involve mean and SD of bout acceleration and jerk
  
  col_names = names(mapped_bouts)
  
  if(!("time"        %in% col_names)) stop("Missing time column in mapped bouts")
  if(!("x"           %in% col_names)) stop("Missing x column in mapped bouts")
  if(!("y"           %in% col_names)) stop("Missing y column in mapped bouts")
  if(!("z"           %in% col_names)) stop("Missing z column in mapped bouts")
  if(!("bout_status" %in% col_names)) stop("Missing bout_status column in mapped bouts")
  if(!("bout_index"  %in% col_names)) stop("Missing bout_index column in mapped bouts")
  
  if(!("vm"          %in% col_names)) mapped_bouts$vm = (mapped_bouts$x^2 + mapped_bouts$y^2 + mapped_bouts$z^2)^0.5
  
  bout_features = mapped_bouts |>
    group_by(bout_index) |>
    summarise(bout_acceleration = max(vm),
              bout_jerk = mean((vm - dplyr::lag(vm)) / (1/100), na.rm = T))
  
  return(bout_features)
  
}

# Data$data |>
#   group_by(bout_index) |>
#   summarise(bout_acceleration = max(vm),
#             bout_jerk = mean((vm - dplyr::lag(vm)) / (1/100), na.rm = T))
# 
# calculate_bout_features(Data$data)
# 
# 
# ggplot(agd) +
#   theme_bw() +
#   geom_point(aes(x = bout_acceleration, y = bout_jerk)) +
#   scale_x_log10()
# 
# ggplot(agd) +
#   theme_bw() +
#   geom_histogram(aes(x = bout_jerk)) +
#   scale_y_log10()


# test = c(
#   rep(0,8),
#   runif(17),
#   rep(0,3),
#   1,
#   rep(0,3),
#   1,
#   0,
#   runif(6)
# )
# 
# abouts = identify_activity_bouts(test, report_short_bouts = T)
# 
# abouts
# 
# f = rep(abouts$activity_bouts$index, abouts$activity_bouts$bout_length)
# 
# splote = split(test, f)
# 
# cat = sapply(splote, calculate_bout_features)
# 
# as.data.frame(cat)
