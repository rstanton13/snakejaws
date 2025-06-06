## Set wd and load packages
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
  library(patchwork)
  library(corHMM)
  library(OUwie)
  library(phylolm)
  library(ggplot2)
  library(reshape2)
  library(grid)
  library(gridExtra)
  library(ggbeeswarm)
  library(cowplot)
  `%nin%` <- Negate(`%in%`)
  options(scipen=999)
}

## Some info on what the non-obvious constraint group names mean
{
  ## "9cats" = all nine categories; 
  ## "bigbulky" = fish+crusts+otherverts; 
  ## "ingestionratio" = fish+crusts+verts & bugs+worms & eels+lizards & mollusks & snakes
  ## "nolimbs" = snakes+worms+eels
  ## "clades" = snakes+lizards & fish+eels & bugs+crusts & mollusks+worms
  ## "shapesize" = fish+crusts+otherverts & snakes+worms+eels
  ## "RQLgroup1" = crusts+fish+mollusks + bugs+worms+lizards+eels
  ## "RQLgroup2" = crusts+fish+mollusks+otherverts + bugs+worms+lizards+eels
  ## "RQLgroup3" = crusts+fish+mollusks+otherverts + bugs+worms+lizards+eels+snakes
  ## "SEWB" = snakes+eels+worms+bugs
}

## Load the data and the tree
{
  ## morphological and phylogenetic data
  dat <- read.csv("hOUwie/skull_dat.csv",row.names=1)
  rdat = read.csv("hOUwie/resid_dat.csv",row.names=1)
  dat = cbind(dat[,c("species","diet","MA","RQL","ma_var","rql_var")], rdat[,c("inlever","quad","inlever_error","quad_error")])
  phy <- read.tree("trees/Title2024_macrostomata_names_replaced.tre")
  phy <- keep.tip(phy, dat$species)
  Tmax <- max(branching.times(phy))
  rdat = read.csv("hOUwie/resid_dat.csv",row.names=1)
}

## Set up the analysis function
hOUwie_big_run = function(data, phy, output_dir, scaling=NULL, n_sim=100,
                          trait_to_analyze=c("MA","RQL"),
                          overwrite_old=FALSE,show_diagn_msg=FALSE,
                          grouping=c("9cats","bigbulky","ingestionratio",
                                     "nolimbs","clades","shapesize",
                                     "RQLgroup1","RQLgroup2","RQLgroup3",
                                     "tet-crusts","tet-fish","fish-crusts",
                                     "snakes-eels","snakes-worms","eels-worms",
                                     "SEWB", "bigbulkybugs",
                                     "FC-SW","FC-SE","FC-nolimbs",
                                     "bigbulky-SW","bigbulky-SE"
                          ), 
                          model_types=c("null","cd","hyb"), 
                          model_classes=c("bm1","bmv","ou1",
                                          "oua","ouv","ouva",
                                          "oum","oumv","oumva")
){
  ## Create output folders as necessary
  sapply(trait_to_analyze, function(name) {
    folder_path <- file.path(output_dir, name)
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
  })
  
  ## Set up internal functions 
  {
    renumber_matrix = function(mat){
      flattened = as.vector(t(mat))
      namevec = setNames(1:length(unique(flattened)),
                         flattened[!duplicated(flattened)])
      if(any(is.na(names(namevec)))){
        namevec[(which(is.na(names(namevec)))+1):length(namevec)]=namevec[(which(is.na(names(namevec)))+1):length(namevec)]-1
        namevec = namevec[-which(is.na(names(namevec)))]
      }
      newmat = matrix(namevec[as.character(flattened)], nrow = nrow(mat), byrow = TRUE)
      rownames(newmat)=rownames(mat)
      colnames(newmat)=colnames(mat)
      return(newmat)
    }
    
    make_cont_cd_model = function(g, j, data){
      ## Set up the model
      cont_model = getOUParamStructure(toupper(word(j,2,2,sep="_")), 9)
      colnames(cont_model) = getStateMat4Dat(data[,c("species","diet")])$legend
      ## Homogenize parameters depending on specified grouping
      if(g=="bigbulky"){
        cont_model[,c("crusts","fish")]=cont_model[,c("other verts")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="ingestionratio"){
        cont_model[,c("crusts","fish")]=cont_model[,c("other verts")]
        cont_model[,c("bugs")]=cont_model[,c("worms")]
        cont_model[,c("eels")]=cont_model[,c("lizards")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="nolimbs"){
        cont_model[,c("snakes","worms")]=cont_model[,c("eels")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="clades"){
        cont_model[,c("snakes")]=cont_model[,c("lizards")]
        cont_model[,c("fish")]=cont_model[,c("eels")]
        cont_model[,c("bugs")]=cont_model[,c("crusts")]
        cont_model[,c("mollusks")]=cont_model[,c("worms")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="shapesize"){
        cont_model[,c("crusts","fish")]=cont_model[,c("other verts")]
        cont_model[,c("snakes","worms")]=cont_model[,c("eels")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="RQLgroup1"){
        cont_model[,c("crusts","fish")]=cont_model[,c("mollusks")]
        cont_model[,c("bugs","worms","lizards")]=cont_model[,c("eels")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="RQLgroup2"){
        cont_model[,c("crusts","fish","other verts")]=cont_model[,c("mollusks")]
        cont_model[,c("bugs","worms","lizards")]=cont_model[,c("eels")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="RQLgroup3"){
        cont_model[,c("crusts","fish","other verts")]=cont_model[,c("mollusks")]
        cont_model[,c("bugs","worms","lizards","snakes")]=cont_model[,c("eels")]
        cont_model = renumber_matrix(cont_model)
      }else  if(g=="tet-crusts"){
        cont_model[,c("crusts")]=cont_model[,c("other verts")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="tet-fish"){
        cont_model[,c("fish")]=cont_model[,c("other verts")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="fish-crusts"){
        cont_model[,c("fish")]=cont_model[,c("crusts")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="snakes-eels"){
        cont_model[,c("snakes")]=cont_model[,c("eels")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="snakes-worms"){
        cont_model[,c("snakes")]=cont_model[,c("worms")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="eels-worms"){
        cont_model[,c("eels")]=cont_model[,c("worms")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="SEWB"){
        cont_model[,c("snakes","eels","worms")]=cont_model[,c("bugs")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="bigbulkybugs"){
        cont_model[,c("crusts","fish","bugs")]=cont_model[,c("other verts")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="FC-SW"){
        cont_model[,c("fish")]=cont_model[,c("crusts")]
        cont_model[,c("snakes")]=cont_model[,c("worms")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="FC-SE"){
        cont_model[,c("fish")]=cont_model[,c("crusts")]
        cont_model[,c("snakes")]=cont_model[,c("eels")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="FC-nolimbs"){
        cont_model[,c("fish")]=cont_model[,c("crusts")]
        cont_model[,c("snakes","eels")]=cont_model[,c("eels")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="bigbulky-SW"){
        cont_model[,c("crusts","fish")]=cont_model[,c("other verts")]
        cont_model[,c("snakes")]=cont_model[,c("worms")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="bigbulky-SE"){
        cont_model[,c("crusts","fish")]=cont_model[,c("other verts")]
        cont_model[,c("snakes")]=cont_model[,c("eels")]
        cont_model = renumber_matrix(cont_model)
      }else if(g!="9cats") stop(cat(paste("The specified grouping",g,"is not valid; stopping."),fill=TRUE))
      if(j=="cd_bmv"){
        cont_model["alpha",] = NA
      }
      colnames(cont_model)=colnames(getOUParamStructure(toupper(word(j,2,2,sep="_")), 9))
      return(cont_model)
    }
    
    make_cont_hyb_model = function(g, j, data){
      ## Set up the model
      cont_model = getOUParamStructure(toupper(word(j,2,2,sep="_")), 9, rate.cat = 2)
      colnames(cont_model) = as.vector(outer(getStateMat4Dat(data[,c("species","diet")])$legend,c(1,2),"paste"))
      if(word(j,1,1,sep="_")=="null"){
        cont_model[,1:9]=cont_model[,1]
        cont_model[,10:18]=cont_model[,10]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="bigbulky"){
        cont_model[,c("crusts 1","fish 1")]=cont_model[,c("other verts 1")]
        cont_model[,c("crusts 2","fish 2")]=cont_model[,c("other verts 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="ingestionratio"){
        cont_model[,c("crusts 1","fish 1")]=cont_model[,c("other verts 1")]
        cont_model[,c("bugs 1")]=cont_model[,c("worms 1")]
        cont_model[,c("eels 1")]=cont_model[,c("lizards 1")]
        cont_model[,c("crusts 2","fish 2")]=cont_model[,c("other verts 2")]
        cont_model[,c("bugs 2")]=cont_model[,c("worms 2")]
        cont_model[,c("eels 2")]=cont_model[,c("lizards 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="nolimbs"){
        cont_model[,c("snakes 1","worms 1")]=cont_model[,c("eels 1")]
        cont_model[,c("snakes 2","worms 2")]=cont_model[,c("eels 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="clades"){
        cont_model[,c("snakes 1")]=cont_model[,c("lizards 1")]
        cont_model[,c("fish 1")]=cont_model[,c("eels 1")]
        cont_model[,c("bugs 1")]=cont_model[,c("crusts 1")]
        cont_model[,c("mollusks 1")]=cont_model[,c("worms 1")]
        cont_model[,c("snakes 2")]=cont_model[,c("lizards 2")]
        cont_model[,c("fish 2")]=cont_model[,c("eels 2")]
        cont_model[,c("bugs 2")]=cont_model[,c("crusts 2")]
        cont_model[,c("mollusks 2")]=cont_model[,c("worms 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="shapesize"){
        cont_model[,c("crusts 1","fish 1")]=cont_model[,c("other verts 1")]
        cont_model[,c("snakes 1","worms 1")]=cont_model[,c("eels 1")]
        cont_model[,c("crusts 2","fish 2")]=cont_model[,c("other verts 2")]
        cont_model[,c("snakes 2","worms 2")]=cont_model[,c("eels 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="RQLgroup1"){
        cont_model[,c("crusts 1","fish 1")]=cont_model[,c("mollusks 1")]
        cont_model[,c("bugs 1","worms 1","lizards 1")]=cont_model[,c("eels 1")]
        cont_model[,c("crusts 2","fish 2")]=cont_model[,c("mollusks 2")]
        cont_model[,c("bugs 2","worms 2","lizards 2")]=cont_model[,c("eels 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="RQLgroup2"){
        cont_model[,c("crusts 1","fish 1","other verts 1")]=cont_model[,c("mollusks 1")]
        cont_model[,c("bugs 1","worms 1","lizards 1")]=cont_model[,c("eels 1")]
        cont_model[,c("crusts 2","fish 2","other verts 2")]=cont_model[,c("mollusks 2")]
        cont_model[,c("bugs 2","worms 2","lizards 2")]=cont_model[,c("eels 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="RQLgroup3"){
        cont_model[,c("crusts 1","fish 1","other verts 1")]=cont_model[,c("mollusks 1")]
        cont_model[,c("bugs 1","worms 1","lizards 1","snakes 1")]=cont_model[,c("eels 1")]
        cont_model[,c("crusts 2","fish 2","other verts 2")]=cont_model[,c("mollusks 2")]
        cont_model[,c("bugs 2","worms 2","lizards 2","snakes 2")]=cont_model[,c("eels 2")]
        cont_model = renumber_matrix(cont_model)
      }else  if(g=="tet-crusts"){
        cont_model[,c("crusts 1")]=cont_model[,c("other verts 1")]
        cont_model[,c("crusts 2")]=cont_model[,c("other verts 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="tet-fish"){
        cont_model[,c("fish 1")]=cont_model[,c("other verts 2")]
        cont_model[,c("fish 1")]=cont_model[,c("other verts 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="fish-crusts"){
        cont_model[,c("fish 1")]=cont_model[,c("crusts 2")]
        cont_model[,c("fish 1")]=cont_model[,c("crusts 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="snakes-eels"){
        cont_model[,c("snakes 1")]=cont_model[,c("eels 1")]
        cont_model[,c("snakes 2")]=cont_model[,c("eels 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="snakes-worms"){
        cont_model[,c("snakes 1")]=cont_model[,c("worms 1")]
        cont_model[,c("snakes 2")]=cont_model[,c("worms 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="eels-worms"){
        cont_model[,c("eels 1")]=cont_model[,c("worms 1")]
        cont_model[,c("eels 2")]=cont_model[,c("worms 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="SEWB"){
        cont_model[,c("snakes 1","eels 1","worms 1")]=cont_model[,c("bugs 1")]
        cont_model[,c("snakes 2","eels 2","worms 2")]=cont_model[,c("bugs 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="bigbulkybugs"){
        cont_model[,c("crusts 1","fish 1","bugs 1")]=cont_model[,c("other verts 1")]
        cont_model[,c("crusts 2","fish 2","bugs 2")]=cont_model[,c("other verts 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="FC-SW"){
        cont_model[,c("fish 1")]=cont_model[,c("crusts 1")]
        cont_model[,c("snakes 1")]=cont_model[,c("worms 1")]
        cont_model[,c("fish 2")]=cont_model[,c("crusts 2")]
        cont_model[,c("snakes 2")]=cont_model[,c("worms 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="FC-SE"){
        cont_model[,c("fish 1")]=cont_model[,c("crusts 1")]
        cont_model[,c("snakes 1")]=cont_model[,c("eels 1")]
        cont_model[,c("fish 2")]=cont_model[,c("crusts 2")]
        cont_model[,c("snakes 2")]=cont_model[,c("eels 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="FC-nolimbs"){
        cont_model[,c("fish 1")]=cont_model[,c("crusts 1")]
        cont_model[,c("snakes 1","eels 1")]=cont_model[,c("eels 1")]
        cont_model[,c("fish 2")]=cont_model[,c("crusts 2")]
        cont_model[,c("snakes 2","eels 2")]=cont_model[,c("eels 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="bigbulky-SW"){
        cont_model[,c("crusts 1","fish 1")]=cont_model[,c("other verts 1")]
        cont_model[,c("snakes 1")]=cont_model[,c("worms 1")]
        cont_model[,c("crusts 2","fish 2")]=cont_model[,c("other verts 2")]
        cont_model[,c("snakes 2")]=cont_model[,c("worms 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g=="bigbulky-SE"){
        cont_model[,c("crusts 1","fish 1")]=cont_model[,c("other verts 1")]
        cont_model[,c("snakes 1")]=cont_model[,c("eels 1")]
        cont_model[,c("crusts 2","fish 2")]=cont_model[,c("other verts 2")]
        cont_model[,c("snakes 2")]=cont_model[,c("eels 2")]
        cont_model = renumber_matrix(cont_model)
      }else if(g!="9cats") stop(cat(paste("The specified grouping",g,"is not valid; stopping."),fill=TRUE))
      if(j=="hyb_bmv"|j=="null_bmv"){
        cont_model["alpha",] = NA
      }
      colnames(cont_model)=colnames(getOUParamStructure(toupper(word(j,2,2,sep="_")), 9, rate.cat = 2))
      return(cont_model)
    }
    
    make_Qmat_1rcat = function(n_cont){
      ##                    b  c  e  f  l  m  v  s  w
      Q_mat = matrix(data=c(1, 0, 0, 0, 0, 0, 0, 0, 1, # b
                            0, 1, 0, 1, 0, 0, 0, 0, 0, # c
                            0, 0, 1, 1, 0, 0, 0, 0, 1, # e
                            0, 1, 1, 1, 0, 0, 1, 0, 1, # f
                            1, 0, 1, 0, 1, 1, 1, 1, 1, # l
                            0, 0, 0, 0, 0, 1, 0, 0, 1, # m
                            1, 0, 1, 1, 1, 1, 1, 1, 1, # v
                            1, 0, 0, 0, 1, 0, 1, 1, 0, # s
                            1, 0, 0, 0, 0, 1, 0, 0, 1  # w
      ))
      Q_mat = matrix(data=1,nrow=9,ncol=9)
      Q_mat[Q_mat!=0] = 1:length(Q_mat[Q_mat!=0])
      Qmat = matrix(Q_mat, nrow = 9, byrow = TRUE); rm(Q_mat)
      diag(Qmat)=NA
      Qmat = renumber_matrix(Qmat)#-1
      colnames(Qmat)=rownames(Qmat)=colnames(getOUParamStructure("BMV",9,1))
      
      Qmat = apply(Qmat, MARGIN = c(1,2), FUN = function(x){
        if(is.na(x)){
          x=NA
        }else if(x==0){
          x=NA
        }else if(x!=0){
          x=x+n_cont
        }
      })
      
      return(Qmat)
    }
    
    make_Qmat_2rcat = function(n_cont){
      ##                    b  c  e  f  l  m  v  s  w
      t_mat = matrix(data=c(1, 0, 0, 0, 0, 0, 0, 0, 1, # b
                            0, 1, 0, 1, 0, 0, 0, 0, 0, # c
                            0, 0, 1, 1, 0, 0, 0, 0, 1, # e
                            0, 1, 1, 1, 0, 0, 1, 0, 1, # f
                            1, 0, 1, 0, 1, 1, 1, 1, 1, # l
                            0, 0, 0, 0, 0, 1, 0, 0, 1, # m
                            1, 0, 1, 1, 1, 1, 1, 1, 1, # v
                            1, 0, 0, 0, 1, 0, 1, 1, 0, # s
                            1, 0, 0, 0, 0, 1, 0, 0, 1  # w
      ))
      t_mat = matrix(data=1,nrow=9,ncol=9)
      diag9 = matrix(data=NA,nrow=9,ncol=9)
      diag(diag9) = 1
      
      Q_mat = matrix(data=NA,nrow=18,ncol=18)
      Q_mat[1:9,1:9] = t_mat
      Q_mat[10:18,10:18] = t_mat
      Q_mat[1:9,10:18] = diag9
      Q_mat[10:18,1:9] = diag9
      diag(Q_mat) = NA
      
      vQ = as.vector(Q_mat)
      for(i in 1:length(vQ)){
        qi = vQ[i]
        if(is.na(qi)){
          next
        }else if(qi==0){
          next
        }else if(qi==1){
          vQ[i] = i
        }
      }
      Qmat = renumber_matrix(matrix(vQ,nrow=18))#-1
      colnames(Qmat)=rownames(Qmat)=colnames(getOUParamStructure("BMV",9,2))
      
      Qmat = apply(Qmat, MARGIN = c(1,2), FUN = function(x){
        if(is.na(x)){
          x=NA
        }else if(x==0){
          x=0
        }else if(x!=0){
          x=x+n_cont
        }
      })
      
      return(Qmat)
    }
  }
  
  ## Setting up models to run and some parameters
  {
    all_models = c(as.vector(outer(model_types, model_classes[!grepl("bm1|ou1",model_classes)], paste, sep = "_")))
    if("null" %in% model_types){
      if("ou1" %in% model_classes) all_models=c("null_ou1",all_models)
      if("bm1" %in% model_classes) all_models=c("null_bm1",all_models)
    }
    null_models = all_models[grepl("null_",all_models)]
    cd_models = all_models[grepl("cd_",all_models)]
    hyb_models = all_models[grepl("hyb_",all_models)]
    Tmax = max(branching.times(phy))
    if(is.null(scaling)){
      scaling = word(output_dir,2,sep="/")
    }
  }
  
  ## Analysis loop
  for(tr in trait_to_analyze){
    if(tr =="MA"){
      if(scaling == "ratio"){
        trait = data[,c("species","diet","MA","ma_var")]
      }else if(scaling == "resid"){
        trait = data[,c("species","diet","inlever","inlever_error")]
      }else stop(cat("The scaling is improperly specified; must be either 'ratio' or 'resid'\n"))
      plus_factor = abs(min(trait[,3])) + 2*sd(trait[,3]) + 1
      trait[,3] = trait[,3] + plus_factor
      lower.bounds = c(1e-10, 1e-10, 0)
      upper.bounds = c(log(2)/(0.01*Tmax), log(2)/(0.01*Tmax), max(trait[,3])+sd(trait[,3]))
    }else if(tr=="RQL"){
      if(scaling == "ratio"){
        trait = data[,c("species","diet","RQL","rql_var")]
      }else if(scaling == "resid"){
        trait = data[,c("species","diet","quad","quad_error")]
      }else stop(cat("The scaling is improperly specified; must be either 'ratio' or 'resid'\n"))
      plus_factor = abs(min(trait[,3])) + 2*sd(trait[,3])
      trait[,3] = trait[,3] + plus_factor
      lower.bounds = c(1e-10, 1e-10, 0)
      upper.bounds = c(log(2)/(0.01*Tmax), log(2)/(0.01*Tmax), max(trait[,3])+sd(trait[,3]))
    }else{
      stop(cat("The trait is improperly specified; must be MA or RQL",fill=TRUE))
    }
    
    ## Run null models (since null, grouping does not matter, so only once per trait)
    if(!is_empty(null_models)){
      for(j in null_models){
        if(grepl("ou1|bm1",j)){
          ## For null CID models
          if(!file.exists(paste(output_dir,tr,"/",tr,"_",j,".RDS",sep="")) | overwrite_old){
            
            ## Set up the continuous and discrete character models given the model class and the grouping
            cont_model = make_cont_cd_model(g="9cats",j,dat)
            disc_model = make_Qmat_1rcat(length(unique(na.omit(as.vector(cont_model)))))
            
            ## Fit the model
            cat(paste0("Fitting ",tr," data to model ",toupper(gsub("_"," ",j)),"; started at ", Sys.time(), sep=""),fill=TRUE)
            inp = data.frame(col1="In progress...")
            write_rds(inp,paste(output_dir,tr,"/",tr,"_",j,".RDS",sep=""))
            mod = hOUwie(phy, trait, rate.cat = 1, discrete_model = "ARD",
                         continuous_model = cont_model, nSim = n_sim, 
                         sample_nodes = TRUE, adaptive_sampling = TRUE, 
                         tip.fog = "known", quiet = TRUE, diagn_msg = show_diagn_msg,
                         lb_continuous_model = lower.bounds, ub_continuous_model = upper.bounds, 
                         lb_discrete_model = 0, ub_discrete_model = 1
            )
            write_rds(mod,paste(output_dir,tr,"/",tr,"_",j,".RDS",sep=""))
            cat(paste0("Finished at ",Sys.time(),sep=""),fill=TRUE)
            
          }
        }else{
          ## For null CID+ models
          if(!file.exists(paste(output_dir,tr,"/",tr,"_",j,".RDS",sep="")) | overwrite_old){
            ## Set up the continuous and discrete character models given the model class and the grouping
            cont_model = make_cont_hyb_model(g="9cats",j,dat)
            disc_model = make_Qmat_2rcat(length(unique(na.omit(as.vector(cont_model)))))
            
            ## Fit the model
            cat(paste0("Fitting ",tr," data to model ",toupper(gsub("_"," ",j)),"; started at ", Sys.time(), sep=""),fill=TRUE)
            inp = data.frame(col1="In progress...")
            write_rds(inp,paste(output_dir,tr,"/",tr,"_",j,".RDS",sep=""))
            mod = hOUwie(phy, trait, rate.cat = 2, discrete_model = "ARD",
                         continuous_model = cont_model, 
                         sample_nodes = TRUE, adaptive_sampling = TRUE, 
                         null.model = TRUE, tip.fog = "known", quiet = TRUE, diagn_msg = show_diagn_msg,
                         lb_continuous_model = lower.bounds, ub_continuous_model = upper.bounds, 
                         lb_discrete_model = 0, ub_discrete_model = 1
            )
            write_rds(mod,paste(output_dir,tr,"/",tr,"_",j,".RDS",sep=""))
            cat(paste0("Finished at ",Sys.time(),sep=""),fill=TRUE)
          }
        }
      }
    }
    
    ## Run the CD models
    if(!is_empty(cd_models)){
      for(g in grouping){
        for(j in cd_models){
          if(!file.exists(paste(output_dir,tr,"/",tr,"_",g,"_",j,".RDS",sep="")) | overwrite_old){
            
            ## Set up the continuous and discrete character models given the model class and the grouping
            cont_model = make_cont_cd_model(g,j,dat)
            disc_model = make_Qmat_1rcat(length(unique(na.omit(as.vector(cont_model)))))
            
            
            ## Run the model
            cat(paste0("Fitting ",tr," data to model ",toupper(gsub("_"," ",j))," for ",g,"; started at ", Sys.time(), sep=""),fill=TRUE)
            inp = data.frame(col1="In progress...")
            write_rds(inp,paste(output_dir,tr,"/",tr,"_",g,"_",j,".RDS",sep=""))
            mod = hOUwie(phy, trait, rate.cat = 1, discrete_model = "ARD",
                         continuous_model = cont_model, nSim = n_sim,
                         tip.fog = "known", quiet = TRUE, diagn_msg = show_diagn_msg,
                         lb_continuous_model = lower.bounds, ub_continuous_model = upper.bounds, 
                         lb_discrete_model = 0, ub_discrete_model = 1
            )
            write_rds(mod,paste(output_dir,tr,"/",tr,"_",g,"_",j,".RDS",sep=""))
            cat(paste0("Finished at ",Sys.time(),sep=""),fill=TRUE)
            
          }
        } ## end of j loop
      } ## end of g loop
    }
    
    ## Run the hybrid models
    if(!is_empty(hyb_models)){
      for(g in grouping){
        for(j in hyb_models){
          if(!file.exists(paste(output_dir,tr,"/",tr,"_",g,"_",j,".RDS",sep="")) | overwrite_old){
            
            ## Set up the continuous and discrete character models given the model class and the grouping
            cont_model = make_cont_hyb_model(g,j,dat)
            disc_model = make_Qmat_2rcat(length(unique(na.omit(as.vector(cont_model)))))
            
            ## Run the model
            cat(paste0("Fitting ",tr," data to model ",toupper(gsub("_"," ",j))," for ",g,"; started at ", Sys.time(), sep=""),fill=TRUE)
            inp = data.frame(col1="In progress...")
            write_rds(inp,paste(output_dir,tr,"/",tr,"_",g,"_",j,".RDS",sep=""))
            mod = hOUwie(phy, trait, rate.cat = 2, discrete_model = "ARD",
                         continuous_model = cont_model, nSim = n_sim,
                         sample_nodes = TRUE, adaptive_sampling = TRUE,
                         tip.fog = "known", quiet = TRUE, diagn_msg = show_diagn_msg,
                         lb_continuous_model = lower.bounds, ub_continuous_model = upper.bounds, 
                         lb_discrete_model = 0, ub_discrete_model = 1
            )
            write_rds(mod,paste(output_dir,tr,"/",tr,"_",g,"_",j,".RDS",sep=""))
            cat(paste0("Finished at ",Sys.time(),sep=""),fill=TRUE)
            
          }
        } ## end of j loop
      } ## end of g loop
    }
    
  } ## end of tr loop
}

#ggplot(dat) + geom_violin(aes(y = inlever, x = diet, color = diet)) + ggtitle("Mechanical advantage")
#ggplot(dat) + geom_violin(aes(y = quad, x = diet, color = diet)) + ggtitle("Relative quadrate length")

## Run the analysis
hOUwie_big_run(dat, phy, output_dir = "hOUwie/ratio/", n_sim = 100,
               model_types=c("cd"), model_classes = c("bmv","ouv", "oum","oumv"), trait_to_analyze = c("MA","RQL"), 
               overwrite_old = FALSE
)





## Investigate the results
output_dir = "hOUwie/ratio/" ## change dir to ratio or resid to change which one you analyze/pull results from
ror = ifelse(grepl("ratio",output_dir),"ratio","residual")
rql_avg_pars = readRDS(paste("hOUwie/",ror,"_RQL_results.RDS",sep=""))
ma_avg_pars = readRDS(paste("hOUwie/",ror,"_MA_results.RDS",sep=""))
reduce_runtime = TRUE; n_goodmod = 10

## investigate the results for MA
{
  model_set = list()
  file.names = setNames(list.files(paste(output_dir,"MA",sep="/"), pattern = "\\.RDS"),word(list.files(paste(output_dir,"MA",sep="/"), pattern = "\\.RDS"),1,1,sep="\\."))
  if(reduce_runtime){
    mtab = read.csv(paste(output_dir,"model_table_MA.csv",sep=""),row.names=1)
    goodmods = paste(rownames(mtab)[1:n_goodmod],".RDS",sep="")
    goodmods = file.names[file.names %in% goodmods]
    nin.files = file.names[file.names %nin% paste(rownames(mtab),".RDS",sep="")]
    file.names = c(nin.files,goodmods)
  }
  for(i in 1:length(file.names)){
    filetoread = paste(output_dir,"MA/",file.names[i],sep="")
    if(filetoread %in% paste(output_dir,"MA/",names(model_set),".RDS",sep="")) next
    else model_set[[names(file.names)[i]]] <- read_rds(paste(output_dir,"MA",file.names[i],sep="/"))
  }
  model_set = model_set[sapply(model_set, function(x) class(x)) == "houwie"]
  
  ma_plus = ifelse(grepl("resid",output_dir), (model_set[[1]]$data$inlever[1])-(rdat$inlever[1]), (model_set[[1]]$data$MA[1])-(dat$MA[1]))
  
  model_table = getModelTable(model_set, type="AICc")[order(getModelTable(model_set, type="AICc")$dAICc),]
  head(round(model_table,4),10)
  if(reduce_runtime){
    model_table = rbind(model_table, mtab[rownames(mtab) %nin% rownames(model_table),])
    model_table = model_table[order(model_table$dAICc),]
  }
  write.csv(model_table, paste(output_dir,"model_table_MA.csv",sep=""))
  
  model_avg_pars <- getModelAvgParams(model_set,type = "AICc")
  model_avg_pars$theta = model_avg_pars$theta - ma_plus
  model_avg_pars$expected_mean = model_avg_pars$expected_mean - ma_plus
  ma_avg_pars <- model_avg_pars
  ma_exp <- setNames(model_avg_pars$expected_mean, rownames(model_avg_pars))
  print(head(model_avg_pars))
  
  plot_data <- melt(model_avg_pars)
  MA_plot_data <- plot_data
  plot_data$tip_state <- factor(plot_data$tip_state)
  plot_data$tip_state <- fct_relevel(plot_data$tip_state, c("other verts","lizards","snakes","fish","eels","mollusks","worms","bugs","crusts"))
  pma = ggplot(plot_data, aes(x = tip_state, y = value, color = tip_state)) +
    geom_point(size = 3, shape = 21) +
    stat_summary(fun=mean,geom="point",aes(group=1, size = 2)) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group=1), width = 0.15, color = "black") +
    theme_classic() +
    facet_wrap(~variable, scales = "free") + scale_x_discrete(labels=c("Nonelongate
tetrapods", "Lizards", "Elongate
squamates","Nonelongate
fishes", 
                                                                       "Anguilliform
vertebrates","Gastropods","Earthworms","Arhtropods","Crustaceans"), name = NULL) +
    scale_color_manual(breaks=c("other verts","lizards","snakes","fish","eels","mollusks","worms","bugs","crusts"), guide = "none",
                       values = c("gray50","forestgreen","mediumpurple2","royalblue2","turquoise2","saddlebrown","lightpink2", "goldenrod2","tomato2"))+
    ggtitle("Mechanical advantage")
  print(pma)
  
  ggsave(plot = pma,
         filename = paste("figs/",ror,"MA_averaged_results.png",sep=""),
         height = 6,
         width = 9)
  
  write_rds(ma_avg_pars, paste("hOUwie/",ror,"_MA_results.RDS",sep=""))
  
}

## investigate the results for RQL
{
  model_set = list()
  file.names = setNames(list.files(paste(output_dir,"RQL",sep="/"), pattern = "\\.RDS"),word(list.files(paste(output_dir,"RQL",sep="/"), pattern = "\\.RDS"),1,1,sep="\\."))
  if(reduce_runtime){
    mtab = read.csv(paste(output_dir,"model_table_RQL.csv",sep=""),row.names=1)
    goodmods = paste(rownames(mtab)[1:n_goodmod],".RDS",sep="")
    goodmods = file.names[file.names %in% goodmods]
    nin.files = file.names[file.names %nin% paste(rownames(mtab),".RDS",sep="")]
    file.names = c(nin.files,goodmods)
  }
  for(i in 1:length(file.names)){
    filetoread = paste(output_dir,"RQL/",file.names[i],sep="")
    if(filetoread %in% paste(output_dir,"RQL/",names(model_set),".RDS",sep="")) next
    else model_set[[names(file.names)[i]]] <- read_rds(paste(output_dir,"RQL",file.names[i],sep="/"))
  }
  model_set = model_set[sapply(model_set, function(x) class(x)) == "houwie"]
  
  rql_plus = ifelse(grepl("resid",output_dir), (model_set[[1]]$data$quad[1])-(rdat$quad[1]), (model_set[[1]]$data$RQL[1])-(dat$RQL[1]))
  
  model_table = getModelTable(model_set, type="AICc")[order(getModelTable(model_set, type="AICc")$dAICc),]
  head(round(model_table,4),10)
  if(reduce_runtime){
    model_table = rbind(model_table, mtab[rownames(mtab) %nin% rownames(model_table),])
    model_table = model_table[order(model_table$dAICc),]
  }
  write.csv(model_table, paste(output_dir,"model_table_RQL.csv",sep=""))
  
  model_avg_pars <- getModelAvgParams(model_set, type = "AICc")
  model_avg_pars$theta = model_avg_pars$theta - rql_plus
  model_avg_pars$expected_mean = model_avg_pars$expected_mean - rql_plus
  rql_avg_pars <- model_avg_pars
  rql_exp <- setNames(model_avg_pars$expected_mean, rownames(model_avg_pars))
  print(head(model_avg_pars))
  
  plot_data <- melt(model_avg_pars)
  RQL_plot_data <- plot_data
  plot_data$tip_state <- factor(plot_data$tip_state)
  plot_data$tip_state <- fct_relevel(plot_data$tip_state, c("other verts","lizards","snakes","fish","eels","mollusks","worms","bugs","crusts"))
  prql = ggplot(plot_data, aes(x = tip_state, y = value, color = tip_state)) +
    geom_point(size = 3, shape = 21) +
    stat_summary(fun=mean,geom="point",aes(group=1, size = 2)) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", aes(group=1), width = 0.15, color = "black") +
    theme_classic() +
    facet_wrap(~variable, scales = "free") + scale_x_discrete(labels=c("Nonelongate
tetrapods", "Lizards", "Elongate
squamates","Nonelongate
fishes", 
                                                                       "Anguilliform
vertebrates","Gastropods","Earthworms","Arhtropods","Crustaceans"), name = NULL) +
    scale_color_manual(breaks=c("other verts","lizards","snakes","fish","eels","mollusks","worms","bugs","crusts"), guide = "none",
                       values = c("gray50","forestgreen","mediumpurple2","royalblue2","turquoise2","saddlebrown","lightpink2", "goldenrod2","tomato2"))+
    ggtitle("Relative quadrate length")
  print(prql)
  
  ggsave(plot = prql,
         filename = paste("figs/",ror,"RQL_averaged_results.png",sep=""),
         height = 6,
         width = 9)
  
  write_rds(rql_avg_pars, paste("hOUwie/",ror,"_RQL_results.RDS",sep=""))
  
}

## phylogenetic ANOVA tests on expected values
{
  rql.panova <- phylANOVA(phy, setNames(rql_avg_pars$tip_state, rownames(rql_avg_pars)), setNames(rql_avg_pars$expected_mean, rownames(rql_avg_pars)))
  ma.panova <- phylANOVA(phy, setNames(ma_avg_pars$tip_state, rownames(ma_avg_pars)), setNames(ma_avg_pars$expected_mean, rownames(ma_avg_pars)))
  
  (rql.panova$Pt<0.05)["other verts",c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs")]
  (ma.panova$Pt<0.05)["other verts",c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs")]
  
  rql.theta.sig = names((rql.panova$Pt<0.05)["other verts",c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs")])[(rql.panova$Pt<0.05)["other verts",c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs")]]
  ma.theta.sig = names((ma.panova$Pt<0.05)["other verts",c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs")])[(ma.panova$Pt<0.05)["other verts",c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs")]]
}

## plotting the results of the optima and expected values together for both traits
{
  states <- rql_avg_pars$tip_state
  colors <- ifelse(states=="other verts", "gray50",
                   ifelse(states=="snakes","mediumpurple2",
                          ifelse(states=="worms","lightpink2",
                                 ifelse(states=="eels","turquoise2",
                                        ifelse(states=="mollusks","saddlebrown",
                                               ifelse(states=="crusts","tomato2",
                                                      ifelse(states=="bugs","goldenrod2",
                                                             ifelse(states=="fish","royalblue2",
                                                                    ifelse(states=="lizards", "forestgreen", NA)))))))))
  ## RQL plot
  {    
    plot_data <- melt(rql_avg_pars)
    plot_data$value[plot_data$variable=="expected_mean"|plot_data$variable=="theta"] <- exp(plot_data$value[plot_data$variable=="expected_mean"|plot_data$variable=="theta"])
    plot_data$tip_state <- factor(plot_data$tip_state)
    plot_data$tip_state <- fct_relevel(plot_data$tip_state, c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"))
    theta_data <- unique(plot_data[plot_data$variable=="theta",])
    exp_mean <- plot_data[plot_data$variable=="expected_mean",]
    exp_mean <- rbind(exp_mean, theta_data)
    sig.dat = data.frame(diet = rql.theta.sig, yval = NA, label = "*")
    for(i in 1:nrow(sig.dat)){
      sig.dat$yval[i] = ifelse(mean(exp_mean$value[exp_mean$tip_state==sig.dat$diet[i]])<mean(exp_mean$value[exp_mean$tip_state=="other verts"]),
                               min(exp_mean$value[exp_mean$tip_state==sig.dat$diet[i]]) - (max(exp_mean$value)-min(exp_mean$value))/10,
                               max(exp_mean$value[exp_mean$tip_state==sig.dat$diet[i]]) + (max(exp_mean$value)-min(exp_mean$value))/10
      )        
    }
    
    rql.p = ggplot(exp_mean[exp_mean$variable=="expected_mean",], aes(x=tip_state, fill = tip_state, y = value)) + 
      geom_hline(yintercept = mean(exp_mean$value[exp_mean$tip_state=="other verts" & exp_mean$variable=="expected_mean"]),
                 col = "gray50", lty = 2) + 
      geom_quasirandom(pch = 21, size = 4, alpha = 0.7) +
      scale_color_manual(breaks=c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"), guide = "none",
                         values = c("gray50","forestgreen","mediumpurple2","royalblue2","turquoise2","tomato2","saddlebrown","lightpink2", "goldenrod2"),
                         name=NA,labels=rep(NA,9)) +
      theme_light() + theme(axis.text.x = element_text(angle = 30, hjust = 1), 
                            plot.title = element_text(hjust = 0.5), 
                            axis.title = element_text(size=15), 
                            title=element_text(size=15),
                            panel.border = element_blank(),
                            axis.line = element_line(color = "gray60"),
                            axis.line.y.right = element_blank(),
                            axis.line.x.top = element_blank(),
                            legend.spacing = unit(20, "lines")
      ) +  
      ggtitle("Relative quadrate length") +
      geom_point(data=exp_mean[exp_mean$variable=="theta",], 
                 aes(x=tip_state, fill = tip_state, y = value), 
                 alpha = 0.7, pch = 23, size = 8, color = "black", stroke=1) +
      scale_fill_manual(breaks=c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"), guide = "none",
                        values = c("gray50","forestgreen","mediumpurple2","royalblue2","turquoise2","tomato2","saddlebrown","lightpink2", "goldenrod2"),
                        name=NA,labels=rep(NA,9)) +
      scale_x_discrete(name=NULL, labels = c("Generalist", "Saurophage", "Ophiophage", "Piscivore", "Anguillivore", 
                                             "Crustacivore", "Molluscivore", "Vermivore", "Insectivore")) +
      scale_y_continuous(name="") + 
      geom_text(
        data = sig.dat,
        aes(x = diet, y = yval, label = label),
        inherit.aes = FALSE, size = 8, col = "gray30"
      )
  }
  
  ## MA plot
  {
    plot_data <- melt(ma_avg_pars)
    plot_data$value[plot_data$variable=="expected_mean"|plot_data$variable=="theta"] <- exp(plot_data$value[plot_data$variable=="expected_mean"|plot_data$variable=="theta"])
    plot_data$tip_state <- fct_relevel(plot_data$tip_state, c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"))
    theta_data <- unique(plot_data[plot_data$variable=="theta",])
    exp_mean <- plot_data[plot_data$variable=="expected_mean",]
    exp_mean <- rbind(exp_mean, theta_data)
    sig.dat = data.frame(diet = ma.theta.sig, yval = NA, label = "*")
    for(i in 1:nrow(sig.dat)){
      sig.dat$yval[i] = ifelse(mean(exp_mean$value[exp_mean$tip_state==sig.dat$diet[i]])<mean(exp_mean$value[exp_mean$tip_state=="other verts"]),
                               min(exp_mean$value[exp_mean$tip_state==sig.dat$diet[i]]) - (max(exp_mean$value)-min(exp_mean$value))/10,
                               max(exp_mean$value[exp_mean$tip_state==sig.dat$diet[i]]) + (max(exp_mean$value)-min(exp_mean$value))/10
      )        
    }
    
    ma.p = ggplot(exp_mean[exp_mean$variable=="expected_mean",], aes(x=tip_state, fill = tip_state, y = value)) + 
      geom_hline(yintercept = mean(exp_mean$value[exp_mean$tip_state=="other verts" & exp_mean$variable=="expected_mean"]),
                 col = "gray50", lty = 2) + 
      geom_quasirandom(pch = 21, size = 4, alpha = 0.7) +
      scale_color_manual(breaks=c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"), guide = "none",
                         values = c("gray50","forestgreen","mediumpurple2","royalblue2","turquoise2","tomato2","saddlebrown","lightpink2", "goldenrod2"),
                         name=NA,labels=rep(NA,9)) +
      theme_light() + theme(axis.text.x = element_text(angle = 30, hjust = 1), 
                            plot.title = element_text(hjust = 0.5), 
                            axis.title = element_text(size=15), 
                            title=element_text(size=15),
                            panel.border = element_blank(), # Remove the default panel border
                            axis.line = element_line(color = "gray60"), # Add axis lines for bottom and left
                            axis.line.y.right = element_blank(),       # Remove right axis line
                            axis.line.x.top = element_blank()          # Remove top axis line
      ) +  
      ggtitle("Mechanical advantage") +
      geom_point(data=exp_mean[exp_mean$variable=="theta",], 
                 aes(x=tip_state, fill = tip_state, y = value), 
                 alpha = 0.7, pch = 23, size = 8, color = "black", stroke=1) +
      scale_fill_manual(breaks=c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"), guide = "none",
                        values = c("gray50","forestgreen","mediumpurple2","royalblue2","turquoise2","tomato2","saddlebrown","lightpink2", "goldenrod2"),
                        name=NA,labels=rep(NA,9)) +
      scale_x_discrete(name=NULL, labels = c("Generalist", "Saurophage", "Ophiophage", "Piscivore", "Anguillivore", 
                                             "Crustacivore", "Molluscivore", "Vermivore", "Insectivore")) +
      scale_y_continuous(name=" ")  + 
      geom_text(
        data = sig.dat,
        aes(x = diet, y = yval, label = label),
        inherit.aes = FALSE, size = 8, col = "gray30"
      )
  }
  
  ## merge together and add legend
  {
    legend_plot <- ggplot() +
      geom_point(aes(x = 1, y = 1, shape = "Expected tip value"), size = 4, fill = "black", color = "black") +
      geom_point(aes(x = 1, y = 0.8, shape = "evzBLANK"), size = 8, fill = "white", color = "white") +
      geom_point(aes(x = 1, y = 0.8, shape = "Evolutionary optima"), size = 8, fill = "black", color = "black") +
      scale_shape_manual(
        values = c("Expected tip value" = 21, "Evolutionary optima" = 23, "evzBLANK" = NA),
        labels = c("Evolutionary
optima","","Expected
tip values"),
        guide = guide_legend(override.aes = list(size = 4, fill = "white"))
      ) +
      theme_void() + # Remove axes and gridlines
      theme(legend.position = "right",
            legend.spacing = unit(10, "lines")) +
      guides(shape = guide_legend(title = ""))
    
    # Extract the legend from the legend plot
    suppressWarnings(custom_legend <- get_legend(legend_plot))
    
    ma.p = plot_grid(ma.p, custom_legend, rel_widths = c(5, 1))
    
    ggp = plot_grid(rql.p, ma.p, rel_widths = c(5,6.075709))
    ggp = grid.arrange(ggp,
                       bottom = textGrob("Diet          .", gp = gpar(fontsize = 15)))
  }
  
  print(ggp)
  
  ggsave(plot = ggp,
         filename = paste("figs/hOUwie_results_optima_",ror,".png",sep=""),
         bg="white",
         units="in",
         height=6,
         width=12,
         dpi = 1200
  )
}

## plotting optima and expected values for poster
{
  states <- rql_avg_pars$tip_state
  colors <- ifelse(states=="other verts", "gray75",
                   ifelse(states=="snakes","mediumpurple2",
                          ifelse(states=="worms","lightpink2",
                                 ifelse(states=="eels","turquoise2",
                                        ifelse(states=="mollusks","saddlebrown",
                                               ifelse(states=="crusts","tomato2",
                                                      ifelse(states=="bugs","goldenrod2",
                                                             ifelse(states=="fish","royalblue2",
                                                                    ifelse(states=="lizards", "forestgreen", NA)))))))))
  ## RQL plot
  {    
    plot_data <- melt(rql_avg_pars)
    plot_data$value[plot_data$variable=="expected_mean"|plot_data$variable=="theta"] <- exp(plot_data$value[plot_data$variable=="expected_mean"|plot_data$variable=="theta"])
    plot_data$tip_state <- factor(plot_data$tip_state)
    plot_data$tip_state <- fct_relevel(plot_data$tip_state, c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"))
    theta_data <- unique(plot_data[plot_data$variable=="theta",])
    exp_mean <- plot_data[plot_data$variable=="expected_mean",]
    exp_mean <- rbind(exp_mean, theta_data)
    sig.dat = data.frame(diet = rql.theta.sig, yval = NA, label = "*")
    for(i in 1:nrow(sig.dat)){
      sig.dat$yval[i] = ifelse(mean(exp_mean$value[exp_mean$tip_state==sig.dat$diet[i]])<mean(exp_mean$value[exp_mean$tip_state=="other verts"]),
                               min(exp_mean$value[exp_mean$tip_state==sig.dat$diet[i]]) - (max(exp_mean$value)-min(exp_mean$value))/10,
                               max(exp_mean$value[exp_mean$tip_state==sig.dat$diet[i]]) + (max(exp_mean$value)-min(exp_mean$value))/10
      )        
    }
    
    rql.p = ggplot(exp_mean[exp_mean$variable=="expected_mean",], aes(x=tip_state, fill = tip_state, y = value)) + 
      geom_hline(yintercept = mean(exp_mean$value[exp_mean$tip_state=="other verts" & exp_mean$variable=="expected_mean"]),
                 col = "gray75", lty = 2) + 
      geom_quasirandom(pch = 21, size = 4, alpha = 0.7) +
      scale_color_manual(breaks=c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"), guide = "none",
                         values = c("gray75","forestgreen","mediumpurple2","royalblue2","turquoise2","tomato2","saddlebrown","lightpink2", "goldenrod2"),
                         name=NA,labels=rep(NA,9)) +
      theme(
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent') #transparent legend panel
      ) + theme(axis.text.x = element_text(angle = 30, hjust = 1,color = "white"), 
                            plot.title = element_text(hjust = 0.5,color = "white"), 
                            axis.title = element_text(size=15,color = "white"), 
                            title=element_text(size=15,color = "white"),
                            panel.border = element_blank(),
                            axis.line = element_line(color = "white"),
                            axis.line.y.right = element_blank(),
                            axis.line.x.top = element_blank(),
                            legend.spacing = unit(20, "lines")
      ) +  
      ggtitle("Relative quadrate length") +
      geom_point(data=exp_mean[exp_mean$variable=="theta",], 
                 aes(x=tip_state, fill = tip_state, y = value), 
                 alpha = 0.7, pch = 23, size = 8, color = "black", stroke=1) +
      scale_fill_manual(breaks=c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"), guide = "none",
                        values = c("gray75","forestgreen","mediumpurple2","royalblue2","turquoise2","tomato2","saddlebrown","lightpink2", "goldenrod2"),
                        name=NA,labels=rep(NA,9)) +
      scale_x_discrete(name=NULL, labels = c("Generalist", "Saurophage", "Ophiophage", "Piscivore", "Anguillivore", 
                                             "Crustacivore", "Molluscivore", "Vermivore", "Insectivore")) +
      scale_y_continuous(name="") + 
      geom_text(
        data = sig.dat,
        aes(x = diet, y = yval, label = label),
        inherit.aes = FALSE, size = 8, col = "white"
      )
  }
  
  ## MA plot
  {
    plot_data <- melt(ma_avg_pars)
    plot_data$value[plot_data$variable=="expected_mean"|plot_data$variable=="theta"] <- exp(plot_data$value[plot_data$variable=="expected_mean"|plot_data$variable=="theta"])
    plot_data$tip_state <- fct_relevel(plot_data$tip_state, c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"))
    theta_data <- unique(plot_data[plot_data$variable=="theta",])
    exp_mean <- plot_data[plot_data$variable=="expected_mean",]
    exp_mean <- rbind(exp_mean, theta_data)
    sig.dat = data.frame(diet = ma.theta.sig, yval = NA, label = "*")
    for(i in 1:nrow(sig.dat)){
      sig.dat$yval[i] = ifelse(mean(exp_mean$value[exp_mean$tip_state==sig.dat$diet[i]])<mean(exp_mean$value[exp_mean$tip_state=="other verts"]),
                               min(exp_mean$value[exp_mean$tip_state==sig.dat$diet[i]]) - (max(exp_mean$value)-min(exp_mean$value))/10,
                               max(exp_mean$value[exp_mean$tip_state==sig.dat$diet[i]]) + (max(exp_mean$value)-min(exp_mean$value))/10
      )        
    }
    
    ma.p = ggplot(exp_mean[exp_mean$variable=="expected_mean",], aes(x=tip_state, fill = tip_state, y = value)) + 
      geom_hline(yintercept = mean(exp_mean$value[exp_mean$tip_state=="other verts" & exp_mean$variable=="expected_mean"]),
                 col = "gray75", lty = 2) + 
      geom_quasirandom(pch = 21, size = 4, alpha = 0.7) +
      scale_color_manual(breaks=c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"), guide = "none",
                         values = c("gray75","forestgreen","mediumpurple2","royalblue2","turquoise2","tomato2","saddlebrown","lightpink2", "goldenrod2"),
                         name=NA,labels=rep(NA,9)) +
      theme(
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent') #transparent legend panel
      ) + theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "white"), 
                            plot.title = element_text(hjust = 0.5, color = "white"), 
                            axis.title = element_text(size=15, color = "white"), 
                            title=element_text(size=15, color = "white"),
                            panel.border = element_blank(), # Remove the default panel border
                            axis.line = element_line(color = "gray95"), # Add axis lines for bottom and left
                            axis.line.y.right = element_blank(),       # Remove right axis line
                            axis.line.x.top = element_blank()          # Remove top axis line
      ) +  
      ggtitle("Mechanical advantage") +
      geom_point(data=exp_mean[exp_mean$variable=="theta",], 
                 aes(x=tip_state, fill = tip_state, y = value), 
                 alpha = 0.7, pch = 23, size = 8, color = "black", stroke=1) +
      scale_fill_manual(breaks=c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"), guide = "none",
                        values = c("gray75","forestgreen","mediumpurple2","royalblue2","turquoise2","tomato2","saddlebrown","lightpink2", "goldenrod2"),
                        name=NA,labels=rep(NA,9)) +
      scale_x_discrete(name=NULL, labels = c("Generalist", "Saurophage", "Ophiophage", "Piscivore", "Anguillivore", 
                                             "Crustacivore", "Molluscivore", "Vermivore", "Insectivore")) +
      scale_y_continuous(name=" ")  + 
      geom_text(
        data = sig.dat,
        aes(x = diet, y = yval, label = label),
        inherit.aes = FALSE, size = 8, col = "white"
      )
    
    
  }
  
  ## merge together and add legend
  {
    legend_plot <- ggplot() +
      geom_point(aes(x = 1, y = 1, shape = "Expected tip value"), size = 4, fill = "white", color = "white") +
      geom_point(aes(x = 1, y = 0.8, shape = "evzBLANK"), size = 8, fill = "white", color = "white") +
      geom_point(aes(x = 1, y = 0.8, shape = "Evolutionary optima"), size = 8, fill = "white", color = "white") +
      scale_shape_manual(
        values = c("Expected tip value" = 21, "Evolutionary optima" = 23, "evzBLANK" = NA),
        labels = c("Evolutionary
optima","","Expected
tip values"),
        guide = guide_legend(override.aes = list(size = 4, fill = "white"))
      ) +
      theme_void() + # Remove axes and gridlines
      theme(legend.position = "right",
            legend.spacing = unit(10, "lines")) +
      guides(shape = guide_legend(title = ""))
    
    # Extract the legend from the legend plot
    suppressWarnings(custom_legend <- get_legend(legend_plot))
    
    ma.p = plot_grid(ma.p, custom_legend, rel_widths = c(5, 1))
    
    ggp = plot_grid(rql.p, ma.p, rel_widths = c(5,6.075709))
    ggp = grid.arrange(ggp,
                       bottom = textGrob("Diet          .", gp = gpar(fontsize = 15))
                       )
  }
  
  print(ggp)
  
  ggsave(plot = ggp,
         filename = paste("figs/poster_hOUwie_results_optima_",ror,".png",sep=""),
         bg="transparent",
         units="in",
         height=6,
         width=12,
         dpi = 600
  )
}

## phylogenetic ANOVA tests on expected variances
{
  rql.var.panova <- phylANOVA(phy, setNames(rql_avg_pars$tip_state, rownames(rql_avg_pars)), setNames(rql_avg_pars$expected_var, rownames(rql_avg_pars)))
  ma.var.panova <- phylANOVA(phy, setNames(ma_avg_pars$tip_state, rownames(ma_avg_pars)), setNames(ma_avg_pars$expected_var, rownames(ma_avg_pars)))
  
  (rql.var.panova$Pt<0.05)["other verts",c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs")]
  (ma.var.panova$Pt<0.05)["other verts",c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs")]
  
  rql.sigma.sig = names((rql.var.panova$Pt<0.05)["other verts",c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs")])[(rql.var.panova$Pt<0.05)["other verts",c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs")]]
  ma.sigma.sig = names((ma.var.panova$Pt<0.05)["other verts",c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs")])[(ma.var.panova$Pt<0.05)["other verts",c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs")]]
}

## plotting the results of the expected variances
{
  ## sigma squared plot for RQL
  {
    plot_data <- melt(rql_avg_pars)
    plot_data$value[plot_data$variable=="expected_mean"|plot_data$variable=="theta"] <- exp(plot_data$value[plot_data$variable=="expected_mean"|plot_data$variable=="theta"])
    plot_data$tip_state <- factor(plot_data$tip_state)
    plot_data$tip_state <- fct_relevel(plot_data$tip_state, c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"))
    theta_data <- unique(plot_data[plot_data$variable=="sigma.sq",])
    exp_mean <- plot_data[plot_data$variable=="expected_var",]
    sig.dat = data.frame(diet = rql.sigma.sig, yval = NA, label = "*")
    for(i in 1:nrow(sig.dat)){
      sig.dat$yval[i] = ifelse(mean(exp_mean$value[exp_mean$tip_state==sig.dat$diet[i]])<mean(exp_mean$value[exp_mean$tip_state=="other verts"]),
                               min(exp_mean$value[exp_mean$tip_state==sig.dat$diet[i]]) - (max(exp_mean$value)-min(exp_mean$value))/10,
                               max(exp_mean$value[exp_mean$tip_state==sig.dat$diet[i]]) + (max(exp_mean$value)-min(exp_mean$value))/10
      )
      
    }
    
    rql.var.p = ggplot(exp_mean[exp_mean$variable=="expected_var",], aes(x=tip_state, fill = tip_state, y = value)) + 
      geom_hline(yintercept = mean(exp_mean$value[exp_mean$tip_state=="other verts" & exp_mean$variable=="expected_var"]),
                 col = "gray50", lty = 2) + 
      geom_quasirandom(pch = 24, size = 4, alpha = 0.7) +
      scale_color_manual(breaks=c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"), guide = "none",
                         values = c("gray50","forestgreen","mediumpurple2","royalblue2","turquoise2","tomato2","saddlebrown","lightpink2", "goldenrod2"),
                         name=NA,labels=rep(NA,9)) +
      theme_light() + theme(axis.text.x = element_text(angle = 30, hjust = 1), 
                            plot.title = element_text(hjust = 0.5), 
                            axis.title = element_text(size=15), 
                            title=element_text(size=15),
                            panel.border = element_blank(),
                            axis.line = element_line(color = "gray60"),
                            axis.line.y.right = element_blank(),
                            axis.line.x.top = element_blank(),
                            legend.spacing = unit(20, "lines")
      ) +  
      ggtitle("Relative quadrate length") +
      scale_fill_manual(breaks=c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"), guide = "none",
                        values = c("gray50","forestgreen","mediumpurple2","royalblue2","turquoise2","tomato2","saddlebrown","lightpink2", "goldenrod2"),
                        name=NA,labels=rep(NA,9)) +
      scale_x_discrete(name=NULL, labels = c("Generalist", "Saurophage", "Ophiophage", "Piscivore", "Anguillivore", 
                                             "Crustacivore", "Molluscivore", "Vermivore", "Insectivore")) +
      scale_y_continuous(name="") + 
      geom_text(
        data = sig.dat,
        aes(x = diet, y = yval, label = label),
        inherit.aes = FALSE, size = 8, col = "gray30"
      )
  }
  print(rql.var.p)
  
  ## sigma squared plot for MA
  {
    plot_data <- melt(ma_avg_pars)
    plot_data$value[plot_data$variable=="expected_mean"|plot_data$variable=="theta"] <- exp(plot_data$value[plot_data$variable=="expected_mean"|plot_data$variable=="theta"])
    plot_data$tip_state <- factor(plot_data$tip_state)
    plot_data$tip_state <- fct_relevel(plot_data$tip_state, c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"))
    theta_data <- unique(plot_data[plot_data$variable=="sigma.sq",])
    exp_mean <- plot_data[plot_data$variable=="expected_var",]
    exp_mean <- rbind(exp_mean)
    sig.dat = data.frame(diet = ma.sigma.sig, yval = NA, label = "*")
    for(i in 1:nrow(sig.dat)){
      sig.dat$yval[i] = ifelse(mean(exp_mean$value[exp_mean$tip_state==sig.dat$diet[i]])<mean(exp_mean$value[exp_mean$tip_state=="other verts"]),
                               min(exp_mean$value[exp_mean$tip_state==sig.dat$diet[i]]) - (max(exp_mean$value)-min(exp_mean$value))/10,
                               max(exp_mean$value[exp_mean$tip_state==sig.dat$diet[i]]) + (max(exp_mean$value)-min(exp_mean$value))/10
      )
    }
    
    ma.var.p = ggplot(exp_mean[exp_mean$variable=="expected_var",], aes(x=tip_state, fill = tip_state, y = value)) + 
      geom_hline(yintercept = mean(exp_mean$value[exp_mean$tip_state=="other verts" & exp_mean$variable=="expected_var"]),
                 col = "gray50", lty = 2) + 
      geom_quasirandom(pch = 24, size = 4, alpha = 0.7) +
      scale_color_manual(breaks=c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"), guide = "none",
                         values = c("gray50","forestgreen","mediumpurple2","royalblue2","turquoise2","tomato2","saddlebrown","lightpink2", "goldenrod2"),
                         name=NA,labels=rep(NA,9)) +
      theme_light() + theme(axis.text.x = element_text(angle = 30, hjust = 1), 
                            plot.title = element_text(hjust = 0.5), 
                            axis.title = element_text(size=15), 
                            title=element_text(size=15),
                            panel.border = element_blank(),
                            axis.line = element_line(color = "gray60"),
                            axis.line.y.right = element_blank(),
                            axis.line.x.top = element_blank(),
                            legend.spacing = unit(20, "lines")
      ) +  
      ggtitle("Mechanical advantage") +
      scale_fill_manual(breaks=c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"), guide = "none",
                        values = c("gray50","forestgreen","mediumpurple2","royalblue2","turquoise2","tomato2","saddlebrown","lightpink2", "goldenrod2"),
                        name=NA,labels=rep(NA,9)) +
      scale_x_discrete(name=NULL, labels = c("Generalist", "Saurophage", "Ophiophage", "Piscivore", "Anguillivore", 
                                             "Crustacivore", "Molluscivore", "Vermivore", "Insectivore")) +
      scale_y_continuous(name="") + 
      geom_text(
        data = sig.dat,
        aes(x = diet, y = yval, label = label),
        inherit.aes = FALSE, size = 8, col = "gray30"
      )
  }
  print(ma.var.p)
  
  ## merge them together and add legend
  {
    legend_plot <- ggplot() +
      geom_point(aes(x = 1, y = 0.8, shape = "EVT"), size = 4, fill = "black", color = "black") +
      scale_shape_manual(
        values = c("EVT" = 24),
        labels = c("Expected 
tip variances"),
        guide = guide_legend(override.aes = list(size = 4, fill = "white"))
      ) +
      theme_void() + # Remove axes and gridlines
      theme(legend.position = "right",
            legend.spacing = unit(10, "lines")) +
      guides(shape = guide_legend(title = ""))
    
    # Extract the legend from the legend plot
    suppressWarnings(custom_legend <- get_legend(legend_plot))
    
    ma.var.p = plot_grid(ma.var.p, custom_legend, rel_widths = c(5, 1))
    
    ggpvar = plot_grid(rql.var.p, ma.var.p, rel_widths = c(5,6.075709))
    ggpvar = grid.arrange(ggpvar,
                       bottom = textGrob("Diet          .", gp = gpar(fontsize = 15)))
  }
  
  print(ggpvar)
  
  ggsave(plot = ggpvar,
                  filename = paste("figs/hOUwie_results_variance_",ror,".png",sep=""),
                  bg="white",
                  units="in",
                  height=6,
                  width=12,
                  dpi = 1200
  )
}

## heat map of model supports
ror2 = "ratio" # ""ratio" or "resid"
{
  mtyp = c(null="null",
           free="9cats",
           g8 = "fish-crusts",
           g11 = "snakes-eels",
           g9 = "FC-SE",
           g10 = "FC-SW",
           g1 = "bigbulky",
           g2 = "ingestionratio",
           g3 = "nolimbs",
           g6 = "snakes-eels",
           g7 = "snakes-worms",
           g5 = "shapesize"
  )
  matab = read.csv(paste("hOUwie/",ror2,"/model_table_MA.csv",sep=""),row.names=1)
  rqltab = read.csv(paste("hOUwie/",ror2,"/model_table_RQL.csv",sep=""),row.names=1)
  mmod = setNames(matab$AICcwt,rownames(matab))
  rmod = setNames(rqltab$AICcwt,rownames(rqltab))
  mmod = mmod[mmod>0.001]
  rmod = rmod[rmod>0.001]

  goodgroups = unique(c(word(names(mmod),sep="_",2,2),word(names(rmod),sep="_",2,2)))
  gg = unique(c("null","9cats",mtyp[mtyp %in% goodgroups]))
  
  mmod = setNames(matab$AICcwt,rownames(matab))
  rmod = setNames(rqltab$AICcwt,rownames(rqltab))
  
  fill_mat = function(gg, models=c("bm1","ou1","bmv","ouv","oum","oumv"), vec){
    mat = matrix(nrow=length(gg), ncol = 6)
    rownames(mat)=gg
    colnames(mat)=models
    
    for (row in rownames(mat)) {
      for (col in colnames(mat)) {
        # Use grepl to find the matching name in the vector
        matching_name <- names(vec)[grepl(row,names(vec))&grepl(paste(col,"\\b",sep=""),names(vec))]
        
        # If a match is found, fill the matrix cell
        if (length(matching_name) > 0) {
          mat[row, col] <- vec[matching_name]
        }
      }
    }
  return(mat)
  }
  mmat = fill_mat(gg, vec = mmod)
  rmat = fill_mat(gg, vec = rmod)

  ylabs = 
                 rev(c("Null","No constraints","(Piscivores +
Crustacivores)",
                       "(Ophiophages +
Anguillivores)",
                       "(Piscivores +
Crustacivores) & 
(Ophiophages +
Anguillivores)",
                       "(Piscivores +
Crustacivores) & 
(Ophiophages +
Vermivores)","(Generalists +
Piscivores + 
Crustacivores) &
(Ophiophages +
Vermivores +
Anguillivores)"))
  
  mat_df <- as.data.frame(as.table(mmat))
  colnames(mat_df) <- c("Row", "Column", "Value")
  mat_df$Row <- factor(mat_df$Row, levels = rev(levels(factor(mat_df$Row))))
  maheat <- ggplot(mat_df, aes(x = Column, y = Row, fill = Value)) +
    geom_tile() +                
    scale_fill_gradient(low = "white", high = "#F2685D", limits = c(0, 1)) +
    theme_minimal() +
    scale_y_discrete(labels=ylabs) +
    scale_x_discrete(position = "top",labels=c("BM1", "OU1", "BMV", "OUV", "OUM", "OUMV")) + 
    labs(title = "", x = "", y = "", fill = "AICc weight
(MA models)  ") + 
    theme(
      axis.text.x = element_text(size = 14), 
      axis.text.y = element_text(size = 14),
      legend.position = "bottom"
    )

  mat_df <- as.data.frame(as.table(rmat))
  colnames(mat_df) <- c("Row", "Column", "Value")
  mat_df$Row <- factor(mat_df$Row, levels = rev(levels(factor(mat_df$Row))))
  rqlheat <- ggplot(mat_df, aes(x = Column, y = Row, fill = Value)) +
    geom_tile() +                
    scale_fill_gradient(low = "white", high = "#29AA81", limits = c(0, 1)) +
    theme_minimal() +
    scale_y_discrete(labels=ylabs) + 
    scale_x_discrete(position = "top",labels=c("BM1", "OU1", "BMV", "OUV", "OUM", "OUMV")) + 
    labs(title = "", x = "", y = "", fill = "AICc weight
(RQL models)  ") + 
    theme(
      axis.text.x = element_text(size = 14), 
      axis.text.y = element_text(size = 14),
      legend.position = "bottom"
    )

heatp = maheat + rqlheat  
print(heatp)

ggsave(filename = "figs/heat_plot_raw.png",
       plot = heatp,
       bg="white",
       units="in",
       height=6,
       width=12,
       dpi = 1200)

}

## supplementary plot ratio vs resid
{
  ram = readRDS("hOUwie/ratio_MA_results.RDS")
  rar = readRDS("hOUwie/ratio_RQL_results.RDS")
  rem = readRDS("hOUwie/residual_MA_results.RDS")
  rer = readRDS("hOUwie/residual_RQL_results.RDS")
  
  compdat = data.frame(row.names = rownames(ram),
                       diet = ram$tip_state,
                       ram = dat$MA,
                       rar = dat$RQL,
                       rem = dat$inlever,
                       rer = dat$quad,
                       ramexp = ram$expected_mean,
                       rarexp = rar$expected_mean,
                       remexp = rem$expected_mean,
                       rerexp = rer$expected_mean,
                       ramvar = ram$expected_var,
                       rarvar = rar$expected_var,
                       remvar = rem$expected_var,
                       rervar = rer$expected_var
  )
  
  {
    x = compdat$rar; y = compdat$rer
    x_lim <- c(min(x), max(x))
    y_lim <- c(min(y), max(y))
    rawr = ggplot(compdat,aes(x = rar, y = rer)) + geom_abline(slope = (y_lim[2] - y_lim[1]) / (x_lim[2] - x_lim[1]), intercept = y_lim[1] - (y_lim[2] - y_lim[1]) / (x_lim[2] - x_lim[1]) * x_lim[1], color = "black", linetype = "dashed", size = 1) +
      geom_point(aes(fill = diet), alpha = 0.7, size = 3, pch = 21) + 
      scale_fill_manual(breaks=c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"), guide = "none",
                        values = c("gray50","forestgreen","mediumpurple2","royalblue2","turquoise2","tomato2","saddlebrown","lightpink2", "goldenrod2")) + 
      scale_x_continuous(name="RQL (ratio)") +
      scale_y_continuous(name="RQL (residual)") + 
      theme_light() + theme(plot.margin = margin(b = 0), plot.title = element_text(hjust = 0.5),axis.title = element_text(size=15),title=element_text(size=15)) +  
      ggtitle("RQL data")
    
    x = compdat$ram; y = compdat$rem
    x_lim <- c(min(x), max(x))
    y_lim <- c(min(y), max(y))
    rawm = ggplot(compdat, aes(x = ram, y = rem)) + geom_abline(slope = (y_lim[2] - y_lim[1]) / (x_lim[2] - x_lim[1]), intercept = y_lim[1] - (y_lim[2] - y_lim[1]) / (x_lim[2] - x_lim[1]) * x_lim[1], color = "black", linetype = "dashed", size = 1) +
      geom_point(aes(fill = diet), alpha = 0.7, size = 3, pch = 21) + 
      scale_fill_manual(breaks=c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"), guide = "none",
                        values = c("gray50","forestgreen","mediumpurple2","royalblue2","turquoise2","tomato2","saddlebrown","lightpink2", "goldenrod2")) + 
      scale_x_continuous(name="MA (ratio)") +
      scale_y_continuous(name="MA (residual)") +
      theme_light() + theme(plot.margin = margin(b = 0), plot.title = element_text(hjust = 0.5),axis.title = element_text(size=15),title=element_text(size=15)) +  
      ggtitle("MA data")
    
    x = compdat$rarexp; y = compdat$rerexp
    x_lim <- c(min(x), max(x))
    y_lim <- c(min(y), max(y))
    expr = ggplot(compdat, aes(x = rarexp, y = rerexp)) + geom_abline(slope = (y_lim[2] - y_lim[1]) / (x_lim[2] - x_lim[1]), intercept = y_lim[1] - (y_lim[2] - y_lim[1]) / (x_lim[2] - x_lim[1]) * x_lim[1], color = "black", linetype = "dashed", size = 1) +
      geom_point(aes(fill = diet), alpha = 0.7, size = 3, pch = 23) + 
      scale_fill_manual(breaks=c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"), guide = "none",
                        values = c("gray50","forestgreen","mediumpurple2","royalblue2","turquoise2","tomato2","saddlebrown","lightpink2", "goldenrod2")) + 
      scale_x_continuous(name="Expected RQL (ratio)") +
      scale_y_continuous(name="Expected RQL (residual)") + 
      theme_light() + theme(plot.margin = margin(b = 0), plot.title = element_text(hjust = 0.5),axis.title = element_text(size=15),title=element_text(size=15)) +  
      ggtitle("Expected RQL tip values")
    
    x = compdat$ramexp; y = compdat$remexp
    x_lim <- c(min(x), max(x))
    y_lim <- c(min(y), max(y))
    expm = ggplot(compdat, aes(x = ramexp, y = remexp)) + geom_abline(slope = (y_lim[2] - y_lim[1]) / (x_lim[2] - x_lim[1]), intercept = y_lim[1] - (y_lim[2] - y_lim[1]) / (x_lim[2] - x_lim[1]) * x_lim[1], color = "black", linetype = "dashed", size = 1) +
      geom_point(aes(fill = diet), alpha = 0.7, size = 3, pch = 23) + 
      scale_fill_manual(breaks=c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"), guide = "none",
                        values = c("gray50","forestgreen","mediumpurple2","royalblue2","turquoise2","tomato2","saddlebrown","lightpink2", "goldenrod2")) + 
      scale_x_continuous(name="Expected MA (ratio)") +
      scale_y_continuous(name="Expected MA (residual)") +
      theme_light() + theme(plot.margin = margin(b = 0), plot.title = element_text(hjust = 0.5),axis.title = element_text(size=15),title=element_text(size=15)) +  
      ggtitle("Expected MA tip values")
    
    x = compdat$rarvar; y = compdat$rervar
    x_lim <- c(min(x), max(x))
    y_lim <- c(min(y), max(y))
    varr = ggplot(compdat, aes(x = rarvar, y = rervar)) + geom_abline(slope = (y_lim[2] - y_lim[1]) / (x_lim[2] - x_lim[1]), intercept = y_lim[1] - (y_lim[2] - y_lim[1]) / (x_lim[2] - x_lim[1]) * x_lim[1], color = "black", linetype = "dashed", size = 1) +
      geom_point(aes(fill = diet), alpha = 0.7, size = 3, pch = 24) + 
      scale_fill_manual(breaks=c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"), guide = "none",
                        values = c("gray50","forestgreen","mediumpurple2","royalblue2","turquoise2","tomato2","saddlebrown","lightpink2", "goldenrod2")) + 
      scale_x_continuous(name="Expected RQL variance (ratio)") +
      scale_y_continuous(name="Expected RQL variance (residual)") + 
      theme_light() + theme(plot.margin = margin(b = 0), plot.title = element_text(hjust = 0.5),axis.title = element_text(size=15),title=element_text(size=15)) +  
      ggtitle("Expected RQL variances")
    
    
    x = compdat$ramvar; y = compdat$remvar
    x_lim <- c(min(x), max(x))
    y_lim <- c(min(y), max(y))
    varm = ggplot(compdat, aes(x = ramvar, y = remvar)) + geom_abline(slope = (y_lim[2] - y_lim[1]) / (x_lim[2] - x_lim[1]), intercept = y_lim[1] - (y_lim[2] - y_lim[1]) / (x_lim[2] - x_lim[1]) * x_lim[1], color = "black", linetype = "dashed", size = 1) +
      geom_point(aes(fill = diet), alpha = 0.7, size = 3, pch = 24) + 
      scale_fill_manual(breaks=c("other verts","lizards","snakes","fish","eels","crusts","mollusks","worms","bugs"), guide = "none",
                        values = c("gray50","forestgreen","mediumpurple2","royalblue2","turquoise2","tomato2","saddlebrown","lightpink2", "goldenrod2")) + 
      scale_x_continuous(name="Expected MA variance (ratio)") +
      scale_y_continuous(name="Expected MA variance (residual)") +
      theme_light() + theme(plot.margin = margin(b = 0), plot.title = element_text(hjust = 0.5),axis.title = element_text(size=15),title=element_text(size=15)) +  
      ggtitle("Expected MA variances")
    
    suppgrid <- (rawr | rawm) / (expr | expm) / (varr | varm) +
      plot_annotation(tag_levels = 'a', tag_prefix = " ", tag_suffix = ')')
  }
  
  print(suppgrid)
  
  ggsave(plot = suppgrid,
         filename = "figs/supplementary_ratio_vs_residual.png",
         bg="white",
         units="in",
         height=12,
         width=8,
         dpi = 1200
  )
}


mean(ma_avg_pars$expected_mean[ma_avg_pars$tip_state=="fish"])
mean(ma_avg_pars$expected_mean[ma_avg_pars$tip_state=="other verts"])

elaps = extract.clade(phy, node = getMRCA(phy, tip = c("Aipysurus_laevis","Hydrophis_platurus")))
plot(elaps,show.tip.label = T, cex = 0.5)
axisPhylo()
