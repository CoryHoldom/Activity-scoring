library(GGIRread)
library(ggplot2)

# system('gzip --decompress --keep "Example CWA/sample.cwa.gz')

source("Summary metrics/activity_bouts.R")
source("Summary metrics/activity_index.R")
source("Summary metrics/cum_power.R")
source("Summary metrics/submovements.R")

header = GGIRread::readAxivity("Example cwa/sample.cwa")$header

Data = GGIRread::readAxivity("Example cwa/sample.cwa", start = 1, end = header$blocks, frequency_tol = 0, progressBar = T)

GGIRread::

t_start = 2700000
t_diff = 1000

t = Data$data$time[t_start:(t_start + t_diff)]/100
xx = Data$data$z[t_start:(t_start + t_diff)]

butt = gsignal::butter(n = 4, w = c(0.1/50, 15/50), type = "pass", output = "Sos")
xf = gsignal::filtfilt(butt, x = xx)

#plot(t, xx, type = "l")
#lines(t, xf, col = "red")

t1 = as.POSIXct.numeric(Data$data$time[1])
t2 = as.POSIXlt.numeric(Data$data$time[2])

t1$sec
t2$sec

sort(sapply(ls(), function(x) format(object.size(get(x)), unit = 'auto')))

as.POSIXct.POSIXlt(t2)

ggplot(Data$data[c(12000:13000),]) +
  theme_bw() +
  geom_line(aes(x = time, y = x), colour = "red") +
  geom_line(aes(x = time, y = y), colour = "green") +
  geom_line(aes(x = time, y = z), colour = "blue")

stable_recoring = Data$data[c(12000:13000),]

stat_var = calculate_stationary_variance(stationary_signal = stable_recoring[, 2:4])

agg_data = aggregate_signal(Data$data, window_length = 100)

agg_data$AI = apply(agg_data, MARGIN = 1,
                    FUN = calculate_activity_index,
                    stationary_variance = stat_var,
                    relative_activity_index = F)


ggplot(dplyr::filter(agg_data, AI > 0), aes(x = AI)) +
  theme_bw() +
  geom_histogram(aes(x = x), fill = "red", alpha = 0.4) +
  geom_histogram(aes(x = y), fill = "green", alpha = 0.4) +
  geom_histogram(aes(x = z), fill = "blue", alpha = 0.4)

agg_long = tidyr::pivot_longer(agg_data, cols = x:AI)

bouts = identify_activity_bouts(agg_data$AI)

ggplot(bouts[bouts$bout_length>1,], aes(x = bout_length, group = as.factor(bout_status), colour = as.factor(bout_status))) +
  theme_bw() +
  geom_density() +
  lims(x = c(0,40))

## NEED TO WORK OUT HOW TO REMOVE BOUTS < 2s EFFICIENTLY AND RECALCULATE THE INDICES AND DURATIONS

active_bouts = bouts[bouts$bout_length>1,]

head(active_bouts, 10)
head(bouts, 10)

active_bouts$index = 






