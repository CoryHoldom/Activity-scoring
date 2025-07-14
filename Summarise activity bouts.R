summarise_avtivity_bouts = function(desired_eid){
  
  ## Data import
  
  header = GGIRread::readAxivity(paste0(here::here("CWAs/"), desired_eid, "_90001_0_0.cwa"))$header
  
  Data = GGIRread::readAxivity(paste0(here::here("CWAs/"), desired_eid, "_90001_0_0.cwa"), start = 1, end = header$blocks, progressBar = T)
  
  # Time is encoded as number of seconds since the epoch - multiply by 100 to get number of centiseconds since epoch
  Data$data$int_time = Data$data$time * 100
  
  Data$data$int_time = as.integer(round(Data$data$int_time) - round(Data$data$int_time[1]))
  
  # Truncate last fraction of a second on end of recording
  Data$data = Data$data[1:(floor(nrow(Data$data) / 100) * 100), ]
  
  cat(paste0("\nData successfully imported: eid ", desired_eid, "\n"))
  
  gc(reset = T)
  
  
  
  # Estimate Activity Indices
  
  Data_10s = Data$data |>
    mutate(int_time10 = floor(int_time/100)) |>
    group_by(int_time10) |>
    summarise(s = max(sd(x), sd(y), sd(z))) |>
    filter(s < 0.013)
  
  gc(reset = T)
  
  stat_var = (Data$data |>
                mutate(int_time10 = floor(int_time/100)) |>
                group_by(int_time10) |>
                summarise(sx = sd(x), sy = sd(y), sz = sd(z)) |>
                mutate(s = pmax(sx, sy, sz)) |>
                filter(s < 0.013))$s |>
    mean()
  
  print("Stationary variance estimated")
  
  gc(reset = T)
  
  agg_data = aggregate_signal(Data$data,
                              window_length = 100)
  
  agg_data$AI = apply(agg_data, MARGIN = 1,
                      FUN = calculate_activity_index,
                      stationary_variance = stat_var,
                      relative_activity_index = F)
  
  
  gc(reset = T)
  
  print("AI Estimates done")
  
  ## Estimate activity bouts
  
  a_bouts = identify_activity_bouts(agg_data$AI,
                                    report_short_bouts = F)
  
  
  print("Activity bouts done")
  
  gc(reset = T)
  
  
  ## Map activity bouts back onto raw acc data
  
  Data$data$bout_status = map_activity_bouts(Data$data, a_bouts$activity_bouts)
  
  bout_rl_encoding = rle(Data$data$bout_status)
  
  bout_rl_encoding = tibble(
    lengths = bout_rl_encoding$lengths,
    values = 1:length(bout_rl_encoding$lengths)
  )
  
  Data$data$bout_index = inverse.rle(bout_rl_encoding)
  
  print("AI bouts mapped")
  
  gc(reset = T)
  
  hist_data = Data$data |>
    group_by(bout_index) |>
    summarise(Bout_Status = first(bout_status), Time = last(time) - first(time)) |>
    mutate(eid = desired_eid) |>
    ungroup()
  
  bout_file_name = paste0("AI bouts/", desired_eid, "_AI_bouts.csv")
  
  readr::write_csv(x = hist_data, file = here::here(bout_file_name))
  
}
