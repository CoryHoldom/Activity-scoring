# Actigraphy Metrics
Cory Holdom

## Overview

Pipeline for the derivation of activity summary metrics for wrist-worn
wearables outcomes.

Examples are sourced from the [MotionSense
dataset](https://www.kaggle.com/datasets/malekzadeh/motionsense-dataset).

``` r
library(tidyverse)
```

    ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ✔ ggplot2   3.5.2     ✔ tibble    3.2.1
    ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ✔ purrr     1.0.4     
    ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ✖ dplyr::filter() masks stats::filter()
    ✖ dplyr::lag()    masks stats::lag()
    ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(gsignal)
```


    Attaching package: 'gsignal'

    The following object is masked from 'package:lubridate':

        dst

    The following object is masked from 'package:dplyr':

        filter

    The following objects are masked from 'package:stats':

        filter, gaussian, poly

``` r
filter = dplyr::filter
```

``` r
Data_raw_walk <- readr::read_csv(
  "../Motionsense Data/A_DeviceMotion_data/wlk_7/sub_1.csv",
  col_types = cols_only(`...1` = col_number(),
                        userAcceleration.x = col_double(),
                        userAcceleration.y = col_double(),
                        userAcceleration.z = col_double())
)

names(Data_raw_walk) <- c("ts",
                          "userAcceleration_X",
                          "userAcceleration_Y",
                          "userAcceleration_Z")

Data_raw_walk$ts = Data_raw_walk$ts / 50 # Converting index to a relative time in seconds

head(Data_raw_walk)
```

    # A tibble: 6 × 4
         ts userAcceleration_X userAcceleration_Y userAcceleration_Z
      <dbl>              <dbl>              <dbl>              <dbl>
    1  0               0.0917             0.416               0.0937
    2  0.02            0.367              0.00457            -0.106 
    3  0.04            0.172             -0.217              -0.163 
    4  0.06            0.00496           -0.238              -0.0191
    5  0.08           -0.0403            -0.241               0.0152
    6  0.1            -0.0908            -0.286              -0.0346

``` r
ggplot(Data_raw_walk[1000:1500,], aes(x = ts, y = userAcceleration_Z)) +
  theme_bw() +
  geom_line() +
  labs(title = "Example 10s of walking data", x = "Duration (s)", y = "Z-Acceleration (g)")
```

![](Actigraphy-Metrics_files/figure-commonmark/Importing%20Data-1.png)

``` r
clean_acc = function(raw_signal, filter_order = 6, fs = 50, freq_cutoff = c(0.1, 20)){
  
  # Estimate vector magnitude and direction 
  # Filter data with a band-pass Butterworth filter to smooth motion
  # and remove gravity component
  
  w_cutoff = freq_cutoff / (fs / 2)
  
  butt_filt = gsignal::butter(n = filter_order, w = w_cutoff, type = "pass")
  
  clean_signal = gsignal::filter(butt_filt, raw_signal)
  
  return(clean_signal)
  
}

combine_acc = function(raw_signal, clean = TRUE){
  
  if(clean){
    raw_signal = sapply(raw_signal, clean_acc)
  }
  
  vm = sqrt(raw_signal[,1] ** 2 + raw_signal[,2] ** 2 + raw_signal[,3] ** 2)
  
  return(vm)
  
}
```

``` r
Data_clean_walk = data.frame(
  ts = Data_raw_walk$ts,
  sapply(Data_raw_walk[,2:4], clean_acc)
)

Data_clean_walk$userAcceleration_VM = combine_acc(Data_clean_walk[,2:4], clean = F)

ggplot(Data_clean_walk[1000:1500,], aes(x = ts, y = userAcceleration_VM)) +
  theme_bw() +
  geom_line() +
  labs(title = "Example 10s of (cleaned) walking data", x = "Duration (s)", y = "Acceleration magnitude (g)")
```

![](Actigraphy-Metrics_files/figure-commonmark/unnamed-chunk-1-1.png)
