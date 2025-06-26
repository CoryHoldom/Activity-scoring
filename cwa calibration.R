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
      lims(y = c(-1, 0.5)),
    
    ggplot(Data_cal$data[c(12000:20000), ]) +
      theme_bw() +
      theme(axis.text.x = element_blank()) +
      geom_line(aes(x = time, y = x), colour = "#d95f02") +
      geom_line(aes(x = time, y = y), colour = "#1b9e77") +
      geom_line(aes(x = time, y = z), colour = "#7570b3") +
      lims(y = c(-1, 0.5))
  ),
  align = "hv",
  nrow = 1
)

Data_cwa$data[c(12000:13000), ]$x |> sd()

Data_cal$data[c(12000:13000), ]$x |> sd()






