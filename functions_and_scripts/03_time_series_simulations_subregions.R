#' Porject: Migratory patterns of Urania boisduvalii (Lepidoptera: Uraniidae)
#' as a function of biotic interactions.
#' Time series simulations for different sub-regions,toxicity and connectivity
#' scenarios. The dispersal process start from a random point in a sub-region.
#' The period for Omphalea to became toxic vary across
#' the connectivity scenarios (5,10,15 and 20 neighbors).
#' Paper authors: Claudia Nunez-Penichet
#'                Jorge Soberon
#                 Luis Osorio-Olvera
#' Code author: Luis Osorio-Olvera


# Packages
library(bam)
library(raster)
library(purrr)
library(magrittr)
library(furrr)
set.seed(111)
rm(list=ls())
# Functions
source("functions/00_time_series_csd_functions.R")
# Read suitability rasters
ura <- raster::raster("modelo_urania.tif")
omp <- raster::raster("modelo_omphalea.tif")
# Suitability parches for both species
ura_omp <- ura*omp
sparse_mod <- bam::model2sparse(ura_omp)
# Paths to subregions
rutas_path <- list.files("paper/final_rasters/",full.names = TRUE)


# Periods to became toxic
periods_negativeL = c(2,3,3,3,3,6,6,6,9,9,9)
# Periods to become suitable
periods_positiveL = c(3,2,3,6,9,3,6,9,3,6,9)
ngbsL <- c(5,10,15,20)
# List of adjacency  matrices
plan(multisession(workers = length(ngbsL)+1))
adj_matrixL <- ngbsL %>% furrr::future_map(function(x){
  adj_matrix <- bam::adj_mat(sparse_mod,ngbs=x)

},.progress = TRUE,.options = furrr_options(seed=NULL))
plan(sequential)
gc()

# Number of scenarios
nprocesss <- length(rutas_path)*length(adj_matrixL)*
  length(periods_negativeL)
# Regions indexes
seqs_regions <- rep(seq_along(rutas_path),
                    times=length(adj_matrixL)*
                      length(periods_negativeL))
# Adjacency matrices indexes
seqs_matrices <- rep(seq_along(adj_matrixL),
                     each=length(rutas_path)*
                       length(periods_negativeL))
# Periods indexes
seqs_periods <- rep(rep(seq_along(periods_negativeL),
                        each=length(rutas_path)),
                    times=length(adj_matrixL))

sim_scenarios <- data.frame(seqs_matrices,seqs_regions,seqs_periods)

x=1
plan(multisession(workers = 14))
options(future.globals.maxSize= 8500*1024^2)

df_all_tseries <- 1:nprocesss %>% furrr::future_map_dfr(function(x){
  set.seed(111)
  matrix_index <- seqs_matrices[x]
  adj_matrix <- adj_matrixL[[matrix_index]]
  periods_negative <- periods_negativeL[seqs_periods[x]]
  periods_positive <- periods_positiveL[seqs_periods[x]]

  path_focal <- rutas_path[seqs_regions[x]]
  focal_region <- raster::raster(path_focal)

  ts_r1 <- ts_simulte(sp1 = ura,
                      sp2 = omp,
                      adj_matrix = adj_matrix,
                      periods_negative = periods_negative,
                      periods_positive = periods_positive,
                      nsteps = 220,
                      npops = 1,
                      focal_region = focal_region)


  rcoords <- t(sapply(7:ncol(ts_r1), function(y){
    ids <- which(ts_r1[,y]>0 )
    if(length(ids)==0L) return(c(x=NA,y=NA))
    colMeans(ts_r1[ids,c("x","y")])
  }))

  rasname <- rutas_path[seqs_regions[x]] %>%
    stringr::str_split(.,"/",simplify = T)


  df_summary <- data.frame(time = colnames(ts_r1[,-(1:6)]),
                           region_id = seqs_regions[x],
                           rasname,
                           ts_r1[1,(1:3)],rcoords,
                           mean_occ=colMeans(ts_r1[,-(1:6)]),
                           row.names = NULL)


  return(df_summary)

},.progress = TRUE,.options = future_options(seed = NULL))

plan(sequential)

df_all <- data.frame(df_all_tseries)
rio::export(df_all,"time_series_subregions_coordinates.csv")
