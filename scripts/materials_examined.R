## set wd and load packages
{
  rm(list=ls())
  library(readxl)
  library(stringr)
  library(tidyverse)
  library(viridis)
  library(jsonlite)
  library(reptinames)
  library(sqldf)
  library(phytools)
  library(ape)
  library(ggplot2)
  library(patchwork)
  library(OUwie)
  `%nin%` <- Negate(`%in%`)
  options(scipen=999)
}

sheet = read_xlsx("input_sheets/serpentes_MA_data_Title2024_MANUAL.xlsx")

sh = sheet[!is.na(sheet$round2) & !grepl("skip",sheet$Notes),]
spec = sh$Specimen
inst = toupper(word(spec, 1, 1, sep=":"))
spec = paste(inst, (str_extract_all(word(spec, 3, 3, sep = ":"),"\\d+")), sep = " ")

spec = spec[-c(which(grepl("character",spec)))]
spec = sort(spec)
spec = paste(spec, collapse = ", ")

write(spec, file = "specimens_examined.txt")
