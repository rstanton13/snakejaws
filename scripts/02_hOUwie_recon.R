## set wd and load packages
{
  rm(list=ls())
  library(phytools)
  library(ape)
  library(OUwie)
  library(phylolm)
  library(tidyverse)
  `%nin%` <- Negate(`%in%`)
  options(scipen=999)
}

## load best fitting RQL model (CD OUMV)
model = readRDS("hOUwie/ratio/MA/MA_9cats_cd_oum.RDS")

## make the reconstruction (this takes a long time)
hrecon = hOUwie.recon(model, nodes = "internal"); write_rds(hrecon, "hOUwie/hrecon.RDS")
