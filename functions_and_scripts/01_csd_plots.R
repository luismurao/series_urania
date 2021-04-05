#' Project: Migratory patterns of Urania boisduvalii (Lepidoptera: Uraniidae)
#' as a function of biotic interactions.
#' CSD plots
#' Code to estimate the connectivity suitability dispersal plot (CSD-plot)
#' Date 28-01-2021
#' Paper authors: Claudia Nunez-Penichet
#'                Jorge Soberon
#                 Luis Osorio-Olvera
#' Code author: Luis Osorio-Olvera

library(bam)
library(raster)
library(purrr)
library(magrittr)
library(ggplot2)
library(ggthemes)
library(dplyr)
set.seed(111)
rm(list = ls())
source("functions/00_time_series_csd_functions.R")
# Read suitability models
omp_path <- "modelo_omphalea.tif"
omp <- raster::raster(omp_path)
ura_path <- "modelo_urania.tif"
ura <- raster::raster(ura_path)
# Suitability parches for both species
ura_omp <- ura*omp
# Convert suitability map to a sparse model
spar_mod <- bam::model2sparse(model = ura_omp)
# Connectivity matrices, for 1,2,..,8,9,10,15, and 22 neighbors
dsteps <- c(seq(1,10,2),c(15,20,22))
csd_vals <- csd_estimate(model = spar_mod,dispersal_steps = dsteps)
csd_vals$plot_data
#write.csv(csd_vals$plot_data,"paper/csd_plot_data.csv",row.names = F)
df <-read.csv("paper/csd_plot_data.csv")
#df <- csd_vals$plot_data

# Code to plot the CSD diagram

labs <- c(0,50,100,150,200,"","","",600)

p2 <- ggplot(df, aes(d, Clusters)) +
  #geom_point(shape = 16, size = 5, show.legend = FALSE) +
  geom_point(size= log10(df$mean_area)+2.5,fill=NA,shape=21, stroke = 1.2) +
  theme_classic() +
  #geom_point() +
  #scale_color_gradient(low = "#32aeff", high = "#f2aeff") +
  scale_alpha(range = c(.25, .6)) +
  scale_y_continuous(trans = squish_trans(200, 600, 300),
                     breaks= c(seq(0,300,50),c(400,660)),
                     labels = labs) +
  theme( axis.line = element_line(colour = "black",
                                  size = 1, linetype = "solid"),
         text = element_text(size=22,family = "Arial"),
         axis.title = element_text(size=20,family = "Arial"),
         axis.text = element_text(size=20,family = "Arial"),
         legend.position="none") +
  xlab("Dispersal distance (km)") +
  ylab("Number of clusters") +
  geom_line(size=1.2,alpha=0.75,col="black")
p2
