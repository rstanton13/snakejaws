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

## load data
make.data.new <- TRUE
if(make.data.new){
  # load data
  view.skull.data <- F
  {
    ## make the get.skull.data function
    get.skull.data <- function(input,type,abnormal,scaler){
      if(is.na(abnormal)){
        if(type == "points"|type == "json"){
          path.to.file <- paste("Landmarks/round2/",input,".json",sep="")
          json <- read_json(path.to.file)
          pts <- json$markups[[1]]$controlPoints
          pts.v <- list()
          for(v in 1:length(pts)){
            coord <- unlist(pts[[v]]$position)
            orie <- unlist(pts[[v]]$orientation)[c(1,5,9)]
            coord <- coord*orie
            pts.v[[v]] <- coord
          }
          pts.v <- matrix(unlist(pts.v),nrow=length(pts.v),ncol=3,byrow=T)
          pts.v <- data.frame(pts.v)
        }
        else if(type == "text"){
          fill.vec <- rep(NA,21*3)
          spl <- str_split(input, "title: ",simplify = TRUE)
          for(j in 1:21){
            spl.j <- spl[j]
            ext <- sub(".*position: ", "", spl.j) 
            ext <- sub("scale.*", "", ext)
            ext <- gsub("[^0-9-] ", "", ext)
            ext <- gsub(",","",ext)
            ext <- as.numeric(str_split(ext," ",simplify=TRUE))
            ext <- na.omit(ext)
            fill.vec[seq(j*3-2,j*3)] <- ext
          }
          m <- matrix(fill.vec, nrow=21, byrow=TRUE)
          pts.v <- data.frame(m)
        }
        
        ## making the data from the points
        {
          colnames(pts.v)=c("x","y","z")
          rownames(pts.v) <- c("kill1","kill2","lhole","lcrest","ldent","lquad.pmp","lquad.amp","lsupra.pmp","lsupra.amp",
                               "kill3","kill4","rhole","rcrest","rdent","rquad.pmp","rquad.amp","rsupra.pmp","rsupra.amp",
                               "parietal.pmp","parietal.amp","premax")
          left.quad.base <- data.frame(x=mean(pts.v[1:2,1]),
                                       y=mean(pts.v[1:2,2]),
                                       z=mean(pts.v[1:2,3]),
                                       row.names = "lquad.base")
          right.quad.base <- data.frame(x=mean(pts.v[10:11,1]),
                                        y=mean(pts.v[10:11,2]),
                                        z=mean(pts.v[10:11,3]),
                                        row.names = "rquad.base")
          left.quad.tip <- data.frame(x=mean(pts.v[6:7,1]),
                                      y=mean(pts.v[6:7,2]),
                                      z=mean(pts.v[6:7,3]),
                                      row.names = "lquad.tip")
          right.quad.tip <- data.frame(x=mean(pts.v[15:16,1]),
                                       y=mean(pts.v[15:16,2]),
                                       z=mean(pts.v[15:16,3]),
                                       row.names = "rquad.tip")
          pts.v2 <- rbind(pts.v[which(!grepl("kill",rownames(pts.v))),],left.quad.base,right.quad.base,left.quad.tip, right.quad.tip)
          
          d <- as.matrix(dist(pts.v2))
          
          skull.length <- d["parietal.pmp","premax"]*sin(pi/2 - acos( 
            ( (d["parietal.pmp","parietal.amp"])^2 +
                (d["parietal.pmp","premax"])^2 - 
                (d["parietal.amp","premax"])^2)/(2*d["parietal.pmp","parietal.amp"]*d["parietal.pmp","premax"]))
          )
          
          left.fossa <- d["lhole","lquad.base"]
          left.in.lever <- d["lcrest","lquad.base"]
          left.quad <- d["lquad.base","lquad.tip"]
          left.mand <- d["lquad.base","ldent"]
          left.supra <- d["lsupra.pmp","lsupra.amp"]
          right.fossa <- d["rhole","rquad.base"]
          right.in.lever <- d["rcrest","rquad.base"]
          right.quad <- d["rquad.base","rquad.tip"]
          right.mand <- d["rquad.base","rdent"]
          right.supra <- d["rsupra.pmp","rsupra.amp"]
          
          dat <- data.frame(left=c(left.mand, left.fossa, left.in.lever, left.supra, left.quad),
                            right=c(right.mand, right.fossa, right.in.lever, right.supra, right.quad),
                            row.names = c("mandible","old.inlever","inlever","supratemporal", "quadrate"))
          dat$mean <- rowMeans(dat)
          dat <- rbind(dat,data.frame(left=skull.length,
                                      right=skull.length,
                                      mean=skull.length,
                                      row.names = "skull.length")) 
        }
        
        ## scaling and other changes depending on abnormality note and scaling factor presence
        {
          if(!is.na(scaler)){
            dat <- (dat/5)*scaler
          }
          if(type == "text"){
            dat <- dat*1000
          }
        }
        return(dat)
      }
      else if(!is.na(abnormal)){
        if(abnormal == "quadrate articulation messed up" & type == "text"){
          fill.vec <- rep(NA,23*3)
          spl <- str_split(input, "title: ",simplify = TRUE)
          for(j in 1:23){
            spl.j <- spl[j]
            ext <- sub(".*position: ", "", spl.j) 
            ext <- sub("scale.*", "", ext)
            ext <- gsub("[^0-9-] ", "", ext)
            ext <- gsub(",","",ext)
            ext <- as.numeric(str_split(ext," ",simplify=TRUE))
            ext <- na.omit(ext)
            fill.vec[seq(j*3-2,j*3)] <- ext
          }
          m <- matrix(fill.vec, nrow=23, byrow=TRUE)
          pts.v <- data.frame(m)
          ## making the data from the points
          {
            colnames(pts.v)=c("x","y","z")
            rownames(pts.v) <- c("kill1","kill2","lhole","lcrest","ldent","lquad.pmp","lquad.amp","lsupra.pmp","lsupra.amp",
                                 "kill3","kill4","rhole","rcrest","rdent","rquad.pmp","rquad.amp","rsupra.pmp","rsupra.amp",
                                 "parietal.pmp","parietal.amp","premax",
                                 "left_articulation", "right_articulation"
            )
            left.quad.base <- data.frame(x=mean(pts.v[1:2,1]),
                                         y=mean(pts.v[1:2,2]),
                                         z=mean(pts.v[1:2,3]),
                                         row.names = "lquad.base")
            right.quad.base <- data.frame(x=mean(pts.v[10:11,1]),
                                          y=mean(pts.v[10:11,2]),
                                          z=mean(pts.v[10:11,3]),
                                          row.names = "rquad.base")
            left.quad.tip <- data.frame(x=mean(pts.v[6:7,1]),
                                        y=mean(pts.v[6:7,2]),
                                        z=mean(pts.v[6:7,3]),
                                        row.names = "lquad.tip")
            right.quad.tip <- data.frame(x=mean(pts.v[15:16,1]),
                                         y=mean(pts.v[15:16,2]),
                                         z=mean(pts.v[15:16,3]),
                                         row.names = "rquad.tip")
            pts.v2 <- rbind(pts.v[which(!grepl("kill",rownames(pts.v))),],left.quad.base,right.quad.base,left.quad.tip, right.quad.tip)
            
            d <- as.matrix(dist(pts.v2))
            
            skull.length <- d["parietal.pmp","premax"]*sin(pi/2 - acos( 
              ( (d["parietal.pmp","parietal.amp"])^2 +
                  (d["parietal.pmp","premax"])^2 - 
                  (d["parietal.amp","premax"])^2)/(2*d["parietal.pmp","parietal.amp"]*d["parietal.pmp","premax"]))
            )
            
            left.fossa <- d["lhole","left_articulation"]
            left.in.lever <- d["lcrest","left_articulation"]
            left.quad <- d["lquad.base","lquad.tip"]
            left.mand <- d["left_articulation","ldent"]
            left.supra <- d["lsupra.pmp","lsupra.amp"]
            right.fossa <- d["rhole","right_articulation"]
            right.in.lever <- d["rcrest","right_articulation"]
            right.quad <- d["rquad.base","rquad.tip"]
            right.mand <- d["right_articulation","rdent"]
            right.supra <- d["rsupra.pmp","rsupra.amp"]
            
            dat <- data.frame(left=c(left.mand, left.fossa, left.in.lever, left.supra, left.quad),
                              right=c(right.mand, right.fossa, right.in.lever, right.supra, right.quad),
                              row.names = c("mandible","old.inlever","inlever","supratemporal", "quadrate"))
            dat$mean <- rowMeans(dat)
            dat <- rbind(dat,data.frame(left=skull.length,
                                        right=skull.length,
                                        mean=skull.length,
                                        row.names = "skull.length")) 
          }
        }
        else if(type == "points"|type == "json"){
          path.to.file <- paste("Landmarks/round2/",input,".json",sep="")
          json <- read_json(path.to.file)
          pts <- json$markups[[1]]$controlPoints
          pts.v <- list()
          for(v in 1:length(pts)){
            coord <- unlist(pts[[v]]$position)
            orie <- unlist(pts[[v]]$orientation)[c(1,5,9)]
            coord <- coord*orie
            pts.v[[v]] <- coord
          }
          pts.v <- matrix(unlist(pts.v),nrow=length(pts.v),ncol=3,byrow=T)
          pts.v <- data.frame(pts.v)
          ## making the data from the points
          {
            colnames(pts.v)=c("x","y","z")
            rownames(pts.v) <- c("kill1","kill2","lhole","lcrest","ldent","lquad.pmp","lquad.amp","lsupra.pmp","lsupra.amp",
                                 "kill3","kill4","rhole","rcrest","rdent","rquad.pmp","rquad.amp","rsupra.pmp","rsupra.amp",
                                 "parietal.pmp","parietal.amp","premax")
            left.quad.base <- data.frame(x=mean(pts.v[1:2,1]),
                                         y=mean(pts.v[1:2,2]),
                                         z=mean(pts.v[1:2,3]),
                                         row.names = "lquad.base")
            right.quad.base <- data.frame(x=mean(pts.v[10:11,1]),
                                          y=mean(pts.v[10:11,2]),
                                          z=mean(pts.v[10:11,3]),
                                          row.names = "rquad.base")
            left.quad.tip <- data.frame(x=mean(pts.v[6:7,1]),
                                        y=mean(pts.v[6:7,2]),
                                        z=mean(pts.v[6:7,3]),
                                        row.names = "lquad.tip")
            right.quad.tip <- data.frame(x=mean(pts.v[15:16,1]),
                                         y=mean(pts.v[15:16,2]),
                                         z=mean(pts.v[15:16,3]),
                                         row.names = "rquad.tip")
            pts.v2 <- rbind(pts.v[which(!grepl("kill",rownames(pts.v))),],left.quad.base,right.quad.base,left.quad.tip, right.quad.tip)
            
            d <- as.matrix(dist(pts.v2))
            
            skull.length <- d["parietal.pmp","premax"]*sin(pi/2 - acos( 
              ( (d["parietal.pmp","parietal.amp"])^2 +
                  (d["parietal.pmp","premax"])^2 - 
                  (d["parietal.amp","premax"])^2)/(2*d["parietal.pmp","parietal.amp"]*d["parietal.pmp","premax"]))
            )
            
            left.fossa <- d["lhole","lquad.base"]
            left.in.lever <- d["lcrest","lquad.base"]
            left.quad <- d["lquad.base","lquad.tip"]
            left.mand <- d["lquad.base","ldent"]
            left.supra <- d["lsupra.pmp","lsupra.amp"]
            right.fossa <- d["rhole","rquad.base"]
            right.in.lever <- d["rcrest","rquad.base"]
            right.quad <- d["rquad.base","rquad.tip"]
            right.mand <- d["rquad.base","rdent"]
            right.supra <- d["rsupra.pmp","rsupra.amp"]
            
            dat <- data.frame(left=c(left.mand, left.fossa, left.in.lever, left.supra, left.quad),
                              right=c(right.mand, right.fossa, right.in.lever, right.supra, right.quad),
                              row.names = c("mandible","old.inlever","inlever","supratemporal", "quadrate"))
            dat$mean <- rowMeans(dat)
            dat <- rbind(dat,data.frame(left=skull.length,
                                        right=skull.length,
                                        mean=skull.length,
                                        row.names = "skull.length")) 
          }
        }
        else if(type == "text"){
          fill.vec <- rep(NA,21*3)
          spl <- str_split(input, "title: ",simplify = TRUE)
          for(j in 1:21){
            spl.j <- spl[j]
            ext <- sub(".*position: ", "", spl.j) 
            ext <- sub("scale.*", "", ext)
            ext <- gsub("[^0-9-] ", "", ext)
            ext <- gsub(",","",ext)
            ext <- as.numeric(str_split(ext," ",simplify=TRUE))
            ext <- na.omit(ext)
            fill.vec[seq(j*3-2,j*3)] <- ext
          }
          m <- matrix(fill.vec, nrow=21, byrow=TRUE)
          pts.v <- data.frame(m)
          ## making the data from the points
          {
            colnames(pts.v)=c("x","y","z")
            rownames(pts.v) <- c("kill1","kill2","lhole","lcrest","ldent","lquad.pmp","lquad.amp","lsupra.pmp","lsupra.amp",
                                 "kill3","kill4","rhole","rcrest","rdent","rquad.pmp","rquad.amp","rsupra.pmp","rsupra.amp",
                                 "parietal.pmp","parietal.amp","premax")
            left.quad.base <- data.frame(x=mean(pts.v[1:2,1]),
                                         y=mean(pts.v[1:2,2]),
                                         z=mean(pts.v[1:2,3]),
                                         row.names = "lquad.base")
            right.quad.base <- data.frame(x=mean(pts.v[10:11,1]),
                                          y=mean(pts.v[10:11,2]),
                                          z=mean(pts.v[10:11,3]),
                                          row.names = "rquad.base")
            left.quad.tip <- data.frame(x=mean(pts.v[6:7,1]),
                                        y=mean(pts.v[6:7,2]),
                                        z=mean(pts.v[6:7,3]),
                                        row.names = "lquad.tip")
            right.quad.tip <- data.frame(x=mean(pts.v[15:16,1]),
                                         y=mean(pts.v[15:16,2]),
                                         z=mean(pts.v[15:16,3]),
                                         row.names = "rquad.tip")
            pts.v2 <- rbind(pts.v[which(!grepl("kill",rownames(pts.v))),],left.quad.base,right.quad.base,left.quad.tip, right.quad.tip)
            
            d <- as.matrix(dist(pts.v2))
            
            skull.length <- d["parietal.pmp","premax"]*sin(pi/2 - acos( 
              ( (d["parietal.pmp","parietal.amp"])^2 +
                  (d["parietal.pmp","premax"])^2 - 
                  (d["parietal.amp","premax"])^2)/(2*d["parietal.pmp","parietal.amp"]*d["parietal.pmp","premax"]))
            )
            
            left.fossa <- d["lhole","lquad.base"]
            left.in.lever <- d["lcrest","lquad.base"]
            left.quad <- d["lquad.base","lquad.tip"]
            left.mand <- d["lquad.base","ldent"]
            left.supra <- d["lsupra.pmp","lsupra.amp"]
            right.fossa <- d["rhole","rquad.base"]
            right.in.lever <- d["rcrest","rquad.base"]
            right.quad <- d["rquad.base","rquad.tip"]
            right.mand <- d["rquad.base","rdent"]
            right.supra <- d["rsupra.pmp","rsupra.amp"]
            
            dat <- data.frame(left=c(left.mand, left.fossa, left.in.lever, left.supra, left.quad),
                              right=c(right.mand, right.fossa, right.in.lever, right.supra, right.quad),
                              row.names = c("mandible","old.inlever","inlever","supratemporal", "quadrate"))
            dat$mean <- rowMeans(dat)
            dat <- rbind(dat,data.frame(left=skull.length,
                                        right=skull.length,
                                        mean=skull.length,
                                        row.names = "skull.length")) 
          }
        }
        
        ## scaling
        {
          if(!is.na(scaler)){
            dat <- (dat/5)*scaler
          }
          if(type == "text"){
            dat <- dat*1000
          }
        }
        
        ## changes because of abnormality
        {
          if(grepl("3=4", abnormal)){
            dat["inlever",] <- dat["old.inlever",]
          }
          if(grepl("no supra", abnormal)){
            dat["supratemporal",] <- NA
          }
          if(grepl("skip", abnormal)){
            dat[1:nrow(dat),] <- NA
          }
          if(grepl("left mandible broken", abnormal)){
            dat["mandible",c("left","mean")] <- dat["mandible","right"]
          }
          if(grepl("right mandible broken", abnormal)){
            dat["mandible",c("right","mean")] <- dat["mandible","left"]
          }
          if(grepl("left side only", abnormal)){
            dat[,c("right","mean")] <- dat[,"left"]
          }
          if(grepl("right side only", abnormal)){
            dat[,c("left","mean")] <- dat[,"right"]
          }
          if(grepl("no quads", abnormal)){
            dat["quadrate",] <- NA
          }
          if(grepl("no left quad", abnormal)){
            dat["quadrate",c("left","mean")] <- dat["quadrate","right"]
          }
          if(grepl("no right quad", abnormal)){
            dat["quadrate",c("right","mean")] <- dat["quadrate","left"]
          }
          if(grepl("no left supra", abnormal)){
            dat["supratemporal",c("left","mean")] <- dat["supratemporal","right"]
          }
          if(grepl("no right supra", abnormal)){
            dat["supratemporal",c("right","mean")] <- dat["supratemporal","left"]
          }
        }
        return(dat)
      }
    }
    ## setting stuff up (loading in the main sheet, making blank data frames, etc) and filling in the data
    suppressWarnings({
      sheet <- readxl::read_xlsx("input_sheets/speicmen_data.xlsx")
      endpt <- length(sheet$type)
      sheet <- sheet[1:endpt,]
      sheet <- sheet[!is.na(sheet$round2),]
      rownames(sheet) <- sheet$Specimen
      
      if(exists("n.runs")){
        if(all(sheet.copy$Specimen %in% sheet$Specimen) & all(sheet$Specimen %in% sheet.copy$Specimen)){
          sheet.copy <- sheet.copy[rownames(sheet),]
          rownames(sheet.copy) <- sheet.copy$Specimen
        }
        n.runs <- n.runs + 1
        if(identical(sheet.copy$round2,sheet$round2)&
           identical(sheet.copy$Notes,sheet$Notes)&
           identical(sheet.copy$scale_factor,sheet$scale_factor)){
          sheet.changed <- FALSE
          cat("The data sheet has not changed since the last time the code ran.")
        }else{
          sheet.changed <- TRUE
          a1 <- sheet
          a2 <- sheet.copy
          new.rows <- sqldf('SELECT * FROM a1 EXCEPT SELECT * FROM a2')
          old.rows <- sqldf('SELECT * FROM a2 EXCEPT SELECT * FROM a1')
          rm("a1","a2")
          
          if(nrow(new.rows)>0){
            sheet <- new.rows
            additional.dat <- data.frame(species=sheet$Species,
                                         specimen = sheet$Specimen,
                                         mand=rep(NA,length(sheet$Species)),
                                         old.inlever=rep(NA,length(sheet$Species)),
                                         inlever=rep(NA,length(sheet$Species)),
                                         supra=rep(NA,length(sheet$Species)),
                                         quad=rep(NA,length(sheet$Species)),
                                         skull=rep(NA,length(sheet$Species))
            )
            
            for(i in 1:nrow(sheet)){
              rep.dat <- get.skull.data(input = sheet$round2[i],
                                        type = sheet$type[i],
                                        abnormal = sheet$Notes[i],
                                        scaler = sheet$scale_factor[i]
              )
              additional.dat[i,"mand"] <- rep.dat$mean[1]
              additional.dat[i,"old.inlever"] <- rep.dat$mean[2]
              additional.dat[i,"inlever"] <- rep.dat$mean[3]
              additional.dat[i,"supra"] <- rep.dat$mean[4]
              additional.dat[i,"quad"] <- rep.dat$mean[5]
              additional.dat[i,"skull"] <- rep.dat$mean[6]
              print(noquote(paste(i, gsub("_", " ",additional.dat$species[i]),additional.dat$specimen[i], sep = ", ")))
            }
            additional.dat <- additional.dat[rowSums(is.na(additional.dat[,3:8]))!=6,]
            rm(rep.dat)
          }
          if(nrow(old.rows)>0){
            skull.dat <- skull.dat[-which(skull.dat$specimen %in% old.rows$Specimen),]
            sheet.copy <- sheet.copy[-which(sheet.copy$Specimen%in%old.rows$Specimen),]
          }
          if(nrow(new.rows)>0){
            additional.dat$MA = additional.dat$inlever/additional.dat$mand
            additional.dat$rql = additional.dat$quad/additional.dat$skull
            skull.dat <- rbind(skull.dat, additional.dat)
            sheet.copy <- rbind(sheet, sheet.copy)
          }
          rownames(sheet.copy) <- sheet.copy$Specimen
        }
      }else{
        n.runs <- 1
        skull.dat <- data.frame(species=sheet$Species,
                                specimen = sheet$Specimen,
                                mand=rep(NA,length(sheet$Species)),
                                old.inlever=rep(NA,length(sheet$Species)),
                                inlever=rep(NA,length(sheet$Species)),
                                supra=rep(NA,length(sheet$Species)),
                                quad=rep(NA,length(sheet$Species)),
                                skull=rep(NA,length(sheet$Species))
        )
        rm(endpt)
        for(i in 1:nrow(sheet)){ 
          rep.dat <- get.skull.data(input = sheet$round2[i],
                                    type = sheet$type[i],
                                    abnormal = sheet$Notes[i],
                                    scaler = sheet$scale_factor[i]
          )
          skull.dat[i,"mand"] <- rep.dat$mean[1]
          skull.dat[i,"old.inlever"] <- rep.dat$mean[2]
          skull.dat[i,"inlever"] <- rep.dat$mean[3]
          skull.dat[i,"supra"] <- rep.dat$mean[4]
          skull.dat[i,"quad"] <- rep.dat$mean[5]
          skull.dat[i,"skull"] <- rep.dat$mean[6]
          print(noquote(paste(i, gsub("_", " ",skull.dat$species[i]),skull.dat$specimen[i], sep = ", ")))
        }
        skull.dat <- skull.dat[rowSums(is.na(skull.dat[,3:8]))!=6,]
        rm(rep.dat)
      }
    })
    
    cat((paste("There are currently",length(unique(skull.dat$species)), "species in the dataset")))
    
    if(view.skull.data) View(skull.dat)
    
    ## plotting variables regressed agaisnt scalers
    {
      p1 <- ggplot(skull.dat,
                   aes(x = skull, y = mand) ) +
        geom_point(color="black") + 
        scale_x_log10(name = "Skull length") + 
        scale_y_log10(name = "Mandible length") + 
        geom_smooth(method = "lm", color = "yellow") +
        ggtitle("Mandible scaling")
      
      p2 <- ggplot(skull.dat,
                   aes(x = mand, y = inlever) ) +
        geom_point(color="black") + 
        scale_x_log10(name = "Mandible length") + 
        scale_y_log10(name = "In-lever length") + 
        geom_smooth(method = "lm", color = "green") +
        ggtitle("In-lever scaling")
      
      p3 <- ggplot(skull.dat,
                   aes(x = skull, y = supra) ) +
        geom_point(color="black") + 
        scale_x_log10(name = "Skull length") + 
        scale_y_log10(name = "Supratemporal length") + 
        geom_smooth(method = "lm", color = "blue") +
        ggtitle("Supratemporal scaling")
      
      p4 <- ggplot(skull.dat,
                   aes(x = skull, y = quad) ) +
        geom_point(color="black") + 
        scale_x_log10(name = "Skull length") + 
        scale_y_log10(name = "Quadrate length") + 
        geom_smooth(method = "lm", color = "purple") +
        ggtitle("Quadrate scaling")
      
      p5 <- ggplot(skull.dat,
                   aes(x = mand, y = inlever) ) +
        geom_point(color="black") + 
        scale_x_log10() + scale_y_log10() + 
        geom_smooth(method = "lm", color = "red") +
        ggtitle("Mechancial advantage scaling")
      
      if(n.runs==1){
        sheet.copy <- sheet
        suppressWarnings(suppressMessages(print(p1+p2+p3+p4)))
      }else if(sheet.changed) {
        suppressWarnings(suppressMessages(print(p1+p2+p3+p4)))
      }
    }
  }
  
  ## setting up data
  {
    skull.dat$MA <- skull.dat$inlever/skull.dat$mand
    skull.dat$rql <- skull.dat$quad/skull.dat$skull
    phy <- read.tree("trees/Title2024_macrostomata_names_replaced.tre")
    skl.dat <- skull.dat[,c("species","inlever","mand","quad","skull")]
    skl.dat <- na.omit(skl.dat[skl.dat$species %in% phy$tip.label,])
    skl.dat$MA <- skl.dat$inlever/skl.dat$mand
    skl.dat$RQL <- skl.dat$quad/skl.dat$skull
    skl.dat[,2:7] <- log(skl.dat[,2:7])
    phy <- keep.tip(phy, unique(skl.dat$species))
    dat <- aggregate(skl.dat[,2:7], by = list(skl.dat$species), FUN = "mean")
    colnames(dat)[1] <- "species"
    dat <- dat[match(phy$tip.label, dat$species),]

    eco <- read_xlsx("input_sheets/snake_species_ecology_discrete.xlsx")
    eco <- eco[eco$species %in% dat$species,]
    eco <- eco[match(dat$species, eco$species),]
    eco$diet[eco$diet=="bird eggs"|eco$diet=="reptile eggs"] <- "other verts"
    eco$diet[eco$diet=="fish eggs"] <- NA
    eco$diet[eco$diet=="unknown"] <- NA
    dat$diet <- eco$diet
    dat <- na.omit(dat)
    tree <- keep.tip(phy, dat$species)
    dat <- dat[match(tree$tip.label, dat$species),]
    diet <- setNames(dat$diet,dat$species)
    par(bg="white")
    dat$color <- ifelse(dat$diet=="other verts", "gray50",
                        ifelse(dat$diet=="snakes","mediumpurple2",
                               ifelse(dat$diet=="worms","lightpink2",
                                      ifelse(dat$diet=="eels","royalblue2",
                                             ifelse(dat$diet=="mollusks","saddlebrown",
                                                    ifelse(dat$diet=="crusts","tomato2",
                                                           ifelse(dat$diet=="bugs","goldenrod2",
                                                                  ifelse(dat$diet=="fish","turquoise2",
                                                                         ifelse(dat$diet=="lizards", "forestgreen", NA)))))))))
    skl.dat <- skl.dat[skl.dat$species %in% dat$species,]
    all(skl.dat$species %in% phy$tip.label)
    
    inlever_resid = phyl.resid(keep.tip(phy, unique(skl.dat$species)), setNames(skl.dat$mand,skl.dat$species), setNames(skl.dat$inlever,skl.dat$species))
    quadrate_resid = phyl.resid(keep.tip(phy, unique(skl.dat$species)), setNames(skl.dat$skull,skl.dat$species), setNames(skl.dat$quad,skl.dat$species))
    resid_data = data.frame(species=rownames(inlever_resid$resid), inlever = inlever_resid$resid[,1], quad = quadrate_resid$resid[,1])
    
    sp <- dat$species
    dat$ma_var <- NA
    dat$rql_var <- NA
    
    n_i = setNames(rep(NA,length(sp)),sp)
    ma_vars = setNames(rep(NA,length(sp)),sp)
    rql_vars = setNames(rep(NA,length(sp)),sp)

    for(i in 1:length(sp)){
        s.i <- skl.dat[skl.dat$species==sp[i],]
        n_i[sp[i]] = nrow(s.i)
        ma_vars[sp[i]] = var(s.i$MA)
        rql_vars[sp[i]] = var(s.i$RQL)
    }

    MAv = na.omit(ma_vars)
    RQLv = na.omit(rql_vars) 
    ni = n_i[names(MAv)]
    
    MA_sigma_w = ( sum(MAv*(ni-1)) )/( sum(ni-1) )
    RQL_sigma_w = ( sum(RQLv*(ni-1)) )/( sum(ni-1) )
    
    all(dat$species==names(n_i))
    dat$ma_var = MA_sigma_w/n_i
    dat$rql_var = RQL_sigma_w/n_i
    
    ## see https://doi.org/10.1111/jbi.14292 for formula used to derive this error measure
    
    resid_data$inlever_error = lm(resid_data$inlever~dat$MA)$coeff[2]*dat$ma_var
    resid_data$quad_error = lm(resid_data$quad~dat$RQL)$coeff[2]*dat$rql_var
    
    write.csv(dat, "hOUwie/skull_dat.csv")
    
    resid_data$diet = dat$diet
    resid_data$color = dat$color
    write.csv(resid_data, "hOUwie/resid_dat.csv")
    
  }

}else{
  dat <- read.csv("hOUwie/skull_dat.csv")
}


write.csv(skull.dat, "c:/Users/riley/Downloads/skull_data.csv")
