library(GGIR)
library(ggplot2)

# GGIR::GGIR(mode = c(1, 2, 3, 4, 5),
#            datadir = "Example CWA",
#            outputdir = "CWA Outs")#, configfile = "CWA Outs/output_Example CWA/config.csv")
# 
# load("CWA Outs/output_Example CWA/meta/ms2.out/sample.cwa.RData")

inspect_file = g.inspectfile("ALS TDI Data/Patient_1039_Feb2018_MOS2D38156710 (2018-02-20).gt3x")

calibrate_params = g.calibrate(datafile = "ALS TDI Data/Patient_1039_Feb2018_MOS2D38156710 (2018-02-20).gt3x", inspectfileobject = inspect_file)

inspect_file = g.inspectfile("Example CWA/sample.cwa")

calibrate_params = g.calibrate(datafile = "Example CWA/sample.cwa", inspectfileobject = inspect_file)


Data_cwa = GGIRread::readAxivity("Example CWA/sample.cwa", end = inspect_file$header$value$blocks)

data_length = nrow(Data_cwa$data)

Data_cal = list(
  header = Data_cwa$header,
  data   = data.frame(
    time = Data_cwa$data$time,
       x = rep(0, data_length),
       y = rep(0, data_length),
       z = rep(0, data_length),
    temp = Data_cwa$data$temp
  ),
  QClog  = Data_cwa$QClog 
)

Data_cal$data$x = calibrate_params$scale[1] * (Data_cwa$data$x +
                                               calibrate_params$offset[1] +
                                               calibrate_params$tempoffset[1])

Data_cal$data$y = calibrate_params$scale[2] * (Data_cwa$data$y +
                                                 calibrate_params$offset[2] +
                                                 calibrate_params$tempoffset[2])

Data_cal$data$z = calibrate_params$scale[3] * (Data_cwa$data$z +
                                                 calibrate_params$offset[3] +
                                                 calibrate_params$tempoffset[3])



Data_cal$data$x.sub = calibrate_params$scale[1] * (Data_cwa$data$x -
                                                 calibrate_params$offset[1] -
                                                 calibrate_params$tempoffset[1])

Data_cal$data$y.sub = calibrate_params$scale[2] * (Data_cwa$data$y -
                                                 calibrate_params$offset[2] -
                                                 calibrate_params$tempoffset[2])

Data_cal$data$z.sub = calibrate_params$scale[3] * (Data_cwa$data$z -
                                                 calibrate_params$offset[3] -
                                                 calibrate_params$tempoffset[3])



Data_cal$data$x.div = (Data_cwa$data$x -
                     calibrate_params$offset[1] -
                     calibrate_params$tempoffset[1])/
  calibrate_params$scale[1]

Data_cal$data$y.div = (Data_cwa$data$y -
                     calibrate_params$offset[2] -
                     calibrate_params$tempoffset[2])/
  calibrate_params$scale[2]

Data_cal$data$z.div = (Data_cwa$data$z -
                     calibrate_params$offset[3] -
                     calibrate_params$tempoffset[3])/
  calibrate_params$scale[3]



data_ba = data.frame(
  M_x = (Data_cal$data$x[1:10000] + Data_cwa$data$x[1:10000])/2,
  D_x = (Data_cal$data$x[1:10000] - Data_cwa$data$x[1:10000]),
  
  M_y = (Data_cal$data$y[1:10000] + Data_cwa$data$y[1:10000])/2,
  D_y = (Data_cal$data$y[1:10000] - Data_cwa$data$y[1:10000]),
  
  M_z = (Data_cal$data$z[1:10000] + Data_cwa$data$z[1:10000])/2,
  D_z = (Data_cal$data$z[1:10000] - Data_cwa$data$z[1:10000])
)

ggplot(data_ba) +
  theme_bw() +
  geom_point(aes(x = M_x, y = D_x))

lm(D_z ~ M_z, data_ba)

plot(Data_cwa$data[1:10000, 2:4])
plot(Data_cal$data[1:10000, 2:4])

ggpubr::ggarrange(
  plotlist = list(
    ggplot(Data_cwa$data[c(12000:20000), ]) +
      theme_bw() +
      theme(axis.text.x = element_blank()) +
      geom_line(aes(x = time, y = x), colour = "#d95f02") +
      geom_line(aes(x = time, y = y), colour = "#1b9e77") +
      geom_line(aes(x = time, y = z), colour = "#7570b3") +
      geom_line(aes(x = time, y = sqrt(x^2+y^2+z^2)), colour = "#4d4d4d") +
      lims(y = c(-1, 1)),
    
    ggplot(Data_cal$data[c(12000:20000), ]) +
      theme_bw() +
      theme(axis.text.x = element_blank()) +
      geom_line(aes(x = time, y = x), colour = "#d95f02") +
      geom_line(aes(x = time, y = y), colour = "#1b9e77") +
      geom_line(aes(x = time, y = z), colour = "#7570b3") +
      geom_line(aes(x = time, y = sqrt(x^2+y^2+z^2)), colour = "#4d4d4d") +
      lims(y = c(-1, 1))
    
  ),
  align = "hv",
  nrow = 1
)

Data_cwa$data[c(12000:13000), ]$x |> sd()

Data_cal$data[c(12000:13000), ]$x |> sd()

lm(sqrt(x^2+y^2+z^2) ~ time, data = Data_cal$data[c(12000:20000), ]) |> summary()

Data_cal$data[c(12000:20000), ] |>
  mutate(vm = sqrt(x^2+y^2+z^2)) |>
  select(c(x,y,z,vm)) |>
  summarise_all(.funs = c("M" = mean, "S" = sd))

Data_cwa$data[c(12000:20000), ] |>
  mutate(vm = sqrt(x^2+y^2+z^2)) |>
  select(c(x,y,z,vm)) |>
  summarise_all(.funs = c("M" = mean, "S" = sd))

Data_dif = list()

Data_dif$data = mutate(Data_cal$data[c(12000:20000), 2:4], vm = sqrt(x^2+y^2+z^2)) - mutate(Data_cwa$data[c(12000:20000), 2:4], vm = sqrt(x^2+y^2+z^2))

ggplot(abs(Data_dif$data)) +
  theme_bw() +
  geom_violin(aes(x = 1, y = x), colour = "#d95f02") +
  geom_violin(aes(x = 2, y = y), colour = "#1b9e77") +
  geom_violin(aes(x = 3, y = z), colour = "#7570b3") +
  geom_violin(aes(x = 4, y = vm), colour = "#4d4d4d")



1-mean((mutate(Data_cwa$data[c(12000:20000), 2:4], vm = sqrt(x^2+y^2+z^2)))$vm)

plot_data_cwa = Data_cwa$data[c(19000:20000), ] |>
  mutate(Time = time - first(time),
         vm = (x^2 + y^2 + z^2)^0.5)

plot_data_cal = Data_cal$data[c(19000:20000), ] |>
  mutate(Time = time - first(time),
         vm = (x^2 + y^2 + z^2)^0.5)

plot_data_cwa$phi = acos(plot_data_cwa$z / plot_data_cwa$vm)

plot_data_cwa$theta = asin(plot_data_cwa$y / (plot_data_cwa$vm * sin(plot_data_cwa$phi)))

plot_data_cal$phi = acos(plot_data_cal$z / plot_data_cal$vm)

plot_data_cal$theta = asin(plot_data_cal$y / (plot_data_cal$vm * sin(plot_data_cal$phi)))

ggpubr::ggarrange(
  plotlist = list(
    ggplot(plot_data_cwa) +
      theme_bw() +
      theme(#axis.text.x = element_blank(),
            axis.ticks = element_blank(), 
            text = element_text(size = 16)) +
      geom_line(aes(x = Time, y = x), colour = "#d95f02") +
      geom_line(aes(x = Time, y = y), colour = "#1b9e77") +
      geom_line(aes(x = Time, y = z), colour = "#7570b3") +
      geom_line(aes(x = Time, y = sqrt(x^2+y^2+z^2)), colour = "#4d4d4d") +
      lims(y = c(-1, 1.1)) +
      labs(x = "Time (s)", y = "Acceleration (g)"),
    
    ggplot(plot_data_cal) +
      theme_bw() +
      theme(#axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_blank(), 
            text = element_text(size = 16)) +
      geom_line(aes(x = Time, y = x), colour = "#d95f02") +
      geom_line(aes(x = Time, y = y), colour = "#1b9e77") +
      geom_line(aes(x = Time, y = z), colour = "#7570b3") +
      geom_line(aes(x = Time, y = sqrt(x^2+y^2+z^2)), colour = "#4d4d4d") +
      lims(y = c(-1, 1.1)) +
      labs(x = "Time (s)")
    
  ),
  align = "hv",
  nrow = 1, labels = c("Raw","Calibrated"), label.x = c(0.24,0.15), label.y = 0.98
)

sd(plot_data_cal$vm)
sd(plot_data_cwa$vm)


ggpubr::ggarrange(plotlist = list(
ggplot(plot_data_cwa) +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    #axis.text.y = element_blank(),
    #axis.title.y = element_blank(),
    axis.ticks = element_blank(), 
    text = element_text(size = 16)) +
  geom_line(aes(x = Time, y = phi / pi * 180), colour = "#d95f02") +
  geom_line(aes(x = Time, y = theta / pi * 180), colour = "#1b9e77") +
  # geom_line(aes(x = Time, y = z), colour = "#7570b3") +
  # geom_line(aes(x = Time, y = sqrt(x^2+y^2+z^2)), colour = "#4d4d4d") +
  scale_y_continuous(breaks = c(-180, -90, 0, 90, 180), limits = c(-180, 180)) +
  labs(x = "Time (s)", y = "Orientation"),

ggplot(plot_data_cal) +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(), 
    text = element_text(size = 16)) +
  geom_line(aes(x = Time, y = phi / pi * 180), colour = "#d95f02") +
  geom_line(aes(x = Time, y = theta / pi * 180), colour = "#1b9e77") +
  # geom_line(aes(x = Time, y = z), colour = "#7570b3") +
  # geom_line(aes(x = Time, y = sqrt(x^2+y^2+z^2)), colour = "#4d4d4d") +
  scale_y_continuous(breaks = c(-180, -90, 0, 90, 180), limits = c(-180, 180)) +
  labs(x = "Time (s)", y = "Orientation")),
align = "hv",
nrow = 1, labels = c("Raw","Calibrated"), label.x = c(0.24,0.15), label.y = 0.98)

sd(plot_data_cal$phi / pi * 180)

