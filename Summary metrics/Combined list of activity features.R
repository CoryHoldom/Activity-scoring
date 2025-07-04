## ---------------------------
##
## Script name: Structure of Activity Outputs
##
## Purpose of script: An example structure of (week) activity outputs
##
## Author: Cory Holdom
##
## Date Created: 2025-06-13
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

activity_features = list(
  
  activity_index = list( T
    ai_mean = ,
    ai_median = ,
    ai_mode = ,
    
    ai_entropy = ,
    
    ai_percent_low = ,
    ai_percent_moderate = ,
    ai_percent_high = ,
    
    ai_percent_single_direction = list(
      
      ai_low = ,
      ai_moderate = ,
      ai_high = ,
      
    )
    
  ),
  
  spectral = list(
    
    spectral_total_power = T
    
  ),
  
  activity_bouts = list( T # Need to use mean and sd of bout accs and jerks over a 24hr period
                         # I think this should be okay now
    
    bout_acceleration = list( ~
      
      acceleration_mean = ,
      acceleration_sd = 
        
    ),
    
    bout_jerk = list( ~
      
      jerk_mean = ,
      jerk_sd = 
      
    )
    
  ),
  
  submovement = list(
    
    sm_distance = ,
    sm_velocity = ,
    sm_acceleration = 
    
  )
  
)

# Test addition from dnanexus
