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
  library(patchwork)
  library(png)
  `%nin%` <- Negate(`%in%`)
  options(scipen=999)
}

## load data
dat <- read.csv("hOUwie/skull_dat.csv",row.names=1)

## load tree and hOUwie ancestral state reconstruction
{
  reconmod <- readRDS("hOUwie/ratio/MA/MA_9cats_cd_oum.RDS")
  tree <- reconmod$phy
  diet <- setNames(reconmod$data$diet,reconmod$data$species)
  
  hrecon = readRDS("hOUwie/hrecon.RDS")
  colnames(hrecon)=sort(unique(reconmod$data$diet))
  colnames(hrecon)[7]="other"
}

## set up
{
  state.colors <- c(other = "black",
                    lizards="forestgreen",
                    snakes = "mediumpurple2",
                    fish="royalblue2",
                    eels = "turquoise2",
                    mollusks = "saddlebrown",
                    worms = "lightpink2",
                    bugs = "goldenrod2",
                    crusts="tomato2"
  )
  
  diet.data <- data.frame(worms=ifelse(diet=="worms",1,0),
                          bugs=ifelse(diet=="bugs",1,0),
                          other=ifelse(diet=="other verts",1,0),
                          snakes=ifelse(diet=="snakes",1,0),
                          mollusks=ifelse(diet=="mollusks",1,0),
                          crusts=ifelse(diet=="crusts",1,0),
                          eels=ifelse(diet=="eels",1,0),
                          fish=ifelse(diet=="fish",1,0),
                          lizards=ifelse(diet=="lizards",1,0)
  )
  
  allStates<-rbind(diet.data[tree$tip.label,colnames(hrecon)],
                   hrecon)
  
  map = reorder(tree, by= "clade")
  map<-paintSubTree(map,node=Ntip(map)+1,
                    state=names(which(hrecon[1,]==max(hrecon[1,]))))
  for(i in 1:nrow(map$edge)){
    states<-sapply(map$edge[i,],function(x,aa) 
      names(aa)[which(aa[x,]==max(aa[x,]))],
      aa=allStates)
    if(length(unique(states))!=1)
      map<-paintSubTree(map,node=map$edge[i,2],state=states[2],
                        stem=0.5)
  }
  
  eco <- readxl::read_xlsx("input_sheets/snake_species_ecology_discrete.xlsx")
  vipers <- eco$species[eco$family == "Viperidae" & eco$species %in% tree$tip.label]
  elapids <- eco$species[eco$family == "Elapidae" & eco$species %in% tree$tip.label]
  dipsadids <- eco$species[eco$family == "Dipsadidae" & eco$species %in% tree$tip.label]
  natricids <- eco$species[eco$family == "Natricidae" & eco$species %in% tree$tip.label]
  colubrids <- eco$species[eco$family == "Colubridae" & eco$species %in% tree$tip.label]
  
  viper.mrca <- getMRCA(tree, vipers)
  elapid.mrca <- getMRCA(tree, elapids)
  dipsadid.mrca <- getMRCA(tree, dipsadids)
  natricid.mrca <- getMRCA(tree, natricids)
  colubrid.mrca <- getMRCA(tree, c("Coluber_constrictor","Ahaetulla_prasina"))
  
  legend.text <- c("Nonelongate tetrapods","Lizards","Elongate squamates","Nonelognate fishes","Anguillform vertebrates","Gastropods","Earthworms","Arthropods","Crustaceans")
}

## make phylogeny trait plot and diet tree with families
{
  ## dietary reconstruction with families
  pdf("figs/tree_dietary_states_with_families.pdf")
  {
    par(bg="white")
    plot(map,state.colors,ftype="off",lwd=1.5,type="fan",part=0.5)
    nodelabels(pie=hrecon,cex=0.15,piecol = state.colors[order(names(state.colors))])
    add.simmap.legend(leg=legend.text, colors=state.colors,prompt=FALSE, x=0.9*par()$usr[1],
                      y=1.1*max(nodeHeights(tree)),fsize=0.7)
    text("Diet",x=0.86*par()$usr[1], y=1.15*max(nodeHeights(tree)))
    arc.cladelabels(tree, "Viperidae", node = viper.mrca, cex = 0.8, mark.node=FALSE)
    arc.cladelabels(tree, "Elapidae", node = elapid.mrca, cex = 0.8, mark.node=FALSE)
    arc.cladelabels(tree, "Dipsadidae", node = dipsadid.mrca, cex = 0.8, mark.node=FALSE)
    arc.cladelabels(tree, "Natricidae", node = natricid.mrca, cex = 0.8, mark.node=FALSE)
    arc.cladelabels(tree, "Colubridae", node = colubrid.mrca, cex = 0.8, mark.node=FALSE)
    
    par(bg="white")
    plot(map,state.colors,ftype="i",fsize=0.175,offset=2,lwd=1.5,type="fan",part=1)
    nodelabels(pie=hrecon,cex=0.15,piecol = state.colors[order(names(state.colors))])
    add.simmap.legend(leg=legend.text, colors=state.colors,prompt=FALSE, x=0.9*par()$usr[1],
                      y=1.1*max(nodeHeights(tree)),fsize=0.7)
    text("Diet",x=0.86*par()$usr[1], y=1.17*max(nodeHeights(tree)))
    
    par(bg="white")
    plot(map,state.colors,ftype="i",fsize=0.175,offset=2,lwd=1.5,type="fan",part=1)
    nodelabels(pie=hrecon,cex=0.15,piecol = state.colors[order(names(state.colors))])
    nodelabels(bg="black",col="white",frame="circle",cex=0.1)
  }
  dev.off()
  
  dat <- read.csv("hOUwie/skull_dat.csv",row.names=1)
  dat = dat[!grepl("Acrochordus",dat$species),]
  dat$color[dat$diet=="other verts"] = "gray75"
  dat$color[dat$diet=="lizards"] = "forestgreen"
  dat$color[dat$diet=="mollusks"] = "saddlebrown"
  dat$color[dat$diet=="fish"] = "royalblue2"
  dat$color[dat$diet=="eels"] = "turquoise2"
  dat$color[dat$diet=="worms"] = "lightpink2"
  dat$color[dat$diet=="snakes"] = "mediumpurple2"
  dat$color[dat$diet=="crusts"] = "tomato2"
  dat$color[dat$diet=="bugs"] = "goldenrod2"
  
  ## dietary reconstruction with MA and RQL plotted
  {
    par(bg="white")
    trait1 = setNames((dat$RQL),dat$species)[map$tip.label]
    trait1[trait1>quantile(trait1,0.99)] <- quantile(trait1, 0.99)
    trait2 = setNames((dat$MA), dat$species)[map$tip.label]
    trait2[trait2<quantile(trait2,0.01)] <- quantile(trait2, 0.01)
    
    trait2 <- (trait2-min(trait2))/(max(trait2)-min(trait2))
    trait1 <- (trait1-min(trait1))/(max(trait1)-min(trait1))
    
    all(names(trait2)==names(trait1))
    map2 <- map
    all(map2$tip.label==map$tip.label)
    
    trait1 = 2+7*trait1
    trait2 = 2+7*trait2
    
    desat <- function(cols, sat=0.5) {
      X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
      hsv(X[1,], X[2,], X[3,])
    }
    hot_color_func <- viridis::magma
    hot_color_func <- colorRampPalette(c("#000004FF","#150E37FF", "#3B0F70FF", "#641A80FF", "#8C2981FF", "#B63679FF", "#DE4968FF", "#F76F5CFF","#FE9F6D","#FECF6E","#FFFF6F"))
    cold_color_func <- colorRampPalette(c("black","#000040", "navyblue", "#183487", "#31688EFF", "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF","#B4DE2CFF","#D8E228"))
    MA.cols <- data.frame(MA=seq(quantile(trait2,0.01),quantile(trait2,0.99),length.out=length(trait2)),col=hot_color_func(length(trait2)))
    quad.cols <- data.frame(rql=seq(quantile(trait1,0.01),quantile(trait1,0.99),length.out=length(trait1)),col=cold_color_func(length(trait1)))
    wider=1
    bg.thinner=0.4
    dx.mult <- 1.005
    dy.mult <- 1.005
    diet.legend = FALSE
    pdf("figs/tree_reconstruction_with_cont_data.pdf")
    {
      h <- max(nodeHeights(map2))
      state.colors["other"] <- "black"
      capture.output(plot(map2, ftype="off",colors=state.colors, 
                          xlim = c(-1.25*h,1.25*h),ylim = c(0,0.75*h), 
                          lwd = 1, type = "fan", part=0.5))
      obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
      sw <- strwidth("l")
      w <- (par()$usr[4] - par()$usr[3])/(max(c(max(trait1)/max(nodeHeights(tree)), 1)) * length(tree$tip.label))
      
      for(i in 1:length(trait1)){
        {
          theta <- atan(obj$yy[i]/obj$xx[i])
          
          s <- ifelse(obj$xx[i] > 0, 1, -1)
          dx <- dx.mult *s * h * cos(theta) + s * cos(theta) * sw
          dy <- dy.mult * s * h * sin(theta) + s * sin(theta) * sw
          
          x1 <- s * max(trait1) * cos(theta) + dx - (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
          x2 <- s * max(trait1) * cos(theta) + dx + (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
          x3 <- x2 - s * trait1[i] * cos(theta)
          x4 <- x1 - s * trait1[i] * cos(theta)
          x5 <- x2 + s * trait2[i] * cos(theta)
          x6 <- x1 + s * trait2[i] * cos(theta)
          
          y1 <- s * max(trait1) * sin(theta) + dy + (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
          y2 <- s * max(trait1) * sin(theta) + dy - (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
          y3 <- y2 - s * trait1[i] * sin(theta)
          y4 <- y1 - s * trait1[i] * sin(theta)
          y5 <- y2 + s * trait2[i] * sin(theta)
          y6 <- y1 + s * trait2[i] * sin(theta)
          
          x1w <- s * max(trait1) * cos(theta) + dx - wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
          x2w <- s * max(trait1) * cos(theta) + dx + wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
          y1w <- s * max(trait1) * sin(theta) + dy + wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
          y2w <- s * max(trait1) * sin(theta) + dy - wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
          
          x01 <- x2w - s * max(trait1) * cos(theta)
          x02 <- x1w - s * max(trait1) * cos(theta)
          x03 <- x2w + s * max(trait1) * cos(theta)
          x04 <- x1w + s * max(trait1) * cos(theta)
          y01 <- y2w - s * max(trait1) * sin(theta)
          y02 <- y1w - s * max(trait1) * sin(theta)
          y03 <- y2w + s * max(trait1) * sin(theta)
          y04 <- y1w + s * max(trait1) * sin(theta)
        }
        
        ## background state colors
        polygon(c(x02, x01, x03, x04), c(y02, y01, y03, y04), 
                col = desat(dat$color[i],bg.thinner), border = desat(dat$color[i],bg.thinner), lwd = 0.01)
        ## trait 1: Quadrate
        #polygon(c(x1, x2, x3, x4), c(y1, y2, y3, y4), 
        #        col = quad.cols$col[which(abs(trait1[i]-quad.cols$rql)==min(abs(trait1[i]-quad.cols$rql)))], border = "black",lwd=0.005)
        ## trait 2: MA
        #polygon(c(x1, x2, x5, x6), c(y1, y2, y5, y6), 
        #        col = MA.cols$col[which(abs(trait2[i]-MA.cols$MA)==min(abs(trait2[i]-MA.cols$MA)))], border = "black",lwd=0.005)
        
        if(i+1<=Ntip(map2)){
          j <- i+1
          
          if(theta >= 0){
            x1.j <- x02
            y1.j <- y02
            x2.j <- x04
            y2.j <- y04
          }else if(theta < 0){
            x1.j <- x01
            y1.j <- y01
            x2.j <- x03
            y2.j <- y03
          }
          
          theta <- atan(obj$yy[j]/obj$xx[j])
          
          if(theta >= 0){
            {
              s <- ifelse(obj$xx[j] > 0, 1, -1)
              dx <- dx.mult *s * h * cos(theta) + s * cos(theta) * sw
              dy <- dy.mult * s * h * sin(theta) + s * sin(theta) * sw
              
              x1 <- s * max(trait1) * cos(theta) + dx - (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
              x2 <- s * max(trait1) * cos(theta) + dx + (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
              x3 <- x2 - s * trait1[j] * cos(theta)
              x4 <- x1 - s * trait1[j] * cos(theta)
              x5 <- x2 + s * trait2[j] * cos(theta)
              x6 <- x1 + s * trait2[j] * cos(theta)
              
              y1 <- s * max(trait1) * sin(theta) + dy + (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
              y2 <- s * max(trait1) * sin(theta) + dy - (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
              y3 <- y2 - s * trait1[j] * sin(theta)
              y4 <- y1 - s * trait1[j] * sin(theta)
              y5 <- y2 + s * trait2[j] * sin(theta)
              y6 <- y1 + s * trait2[j] * sin(theta)
              
              x1w <- s * max(trait1) * cos(theta) + dx - wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
              x2w <- s * max(trait1) * cos(theta) + dx + wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
              y1w <- s * max(trait1) * sin(theta) + dy + wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
              y2w <- s * max(trait1) * sin(theta) + dy - wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
              
              x01 <- x2w - s * max(trait1) * cos(theta)
              x02 <- x1w - s * max(trait1) * cos(theta)
              x03 <- x2w + s * max(trait1) * cos(theta)
              x04 <- x1w + s * max(trait1) * cos(theta)
              y01 <- y2w - s * max(trait1) * sin(theta)
              y02 <- y1w - s * max(trait1) * sin(theta)
              y03 <- y2w + s * max(trait1) * sin(theta)
              y04 <- y1w + s * max(trait1) * sin(theta)
            }
            
            x3.j <- x01
            y3.j <- y01
            x4.j <- x03
            y4.j <- y03
            
            x5.j <- mean(c(x1.j, x3.j))
            x6.j <- mean(c(x2.j, x4.j))
            y5.j <- mean(c(y1.j, y3.j))
            y6.j <- mean(c(y2.j, y4.j))
          }else if(theta<0){
            {
              s <- ifelse(obj$xx[j] > 0, 1, -1)
              dx <- dx.mult *s * h * cos(theta) + s * cos(theta) * sw
              dy <- dy.mult * s * h * sin(theta) + s * sin(theta) * sw
              
              x1 <- s * max(trait1) * cos(theta) + dx - (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
              x2 <- s * max(trait1) * cos(theta) + dx + (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
              x3 <- x2 - s * trait1[j] * cos(theta)
              x4 <- x1 - s * trait1[j] * cos(theta)
              x5 <- x2 + s * trait2[j] * cos(theta)
              x6 <- x1 + s * trait2[j] * cos(theta)
              
              y1 <- s * max(trait1) * sin(theta) + dy + (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
              y2 <- s * max(trait1) * sin(theta) + dy - (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
              y3 <- y2 - s * trait1[j] * sin(theta)
              y4 <- y1 - s * trait1[j] * sin(theta)
              y5 <- y2 + s * trait2[j] * sin(theta)
              y6 <- y1 + s * trait2[j] * sin(theta)
              
              x1w <- s * max(trait1) * cos(theta) + dx - wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
              x2w <- s * max(trait1) * cos(theta) + dx + wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
              y1w <- s * max(trait1) * sin(theta) + dy + wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
              y2w <- s * max(trait1) * sin(theta) + dy - wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
              
              x01 <- x2w - s * max(trait1) * cos(theta)
              x02 <- x1w - s * max(trait1) * cos(theta)
              x03 <- x2w + s * max(trait1) * cos(theta)
              x04 <- x1w + s * max(trait1) * cos(theta)
              y01 <- y2w - s * max(trait1) * sin(theta)
              y02 <- y1w - s * max(trait1) * sin(theta)
              y03 <- y2w + s * max(trait1) * sin(theta)
              y04 <- y1w + s * max(trait1) * sin(theta)
            }
            
            x3.j <- x02
            y3.j <- y02
            x4.j <- x04
            y4.j <- y04
            
            x5.j <- mean(c(x1.j, x3.j))
            x6.j <- mean(c(x2.j, x4.j))
            y5.j <- mean(c(y1.j, y3.j))
            y6.j <- mean(c(y2.j, y4.j))
            
          }
          polygon(c(x2.j,x1.j,x5.j,x6.j), c(y2.j,y1.j,y5.j,y6.j), col=desat(dat$color[i],bg.thinner), border=desat(dat$color[i],bg.thinner), lwd = 0.01)
          polygon(c(x4.j,x3.j,x5.j,x6.j), c(y4.j,y3.j,y5.j,y6.j), col=desat(dat$color[j],bg.thinner), border=desat(dat$color[i],bg.thinner), lwd = 0.01)
        }
        
        {
          theta <- atan(obj$yy[i]/obj$xx[i])
          
          s <- ifelse(obj$xx[i] > 0, 1, -1)
          dx <- dx.mult *s * h * cos(theta) + s * cos(theta) * sw
          dy <- dy.mult * s * h * sin(theta) + s * sin(theta) * sw
          
          x1 <- s * max(trait1) * cos(theta) + dx - (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
          x2 <- s * max(trait1) * cos(theta) + dx + (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
          x3 <- x2 - s * trait1[i] * cos(theta)
          x4 <- x1 - s * trait1[i] * cos(theta)
          x5 <- x2 + s * trait2[i] * cos(theta)
          x6 <- x1 + s * trait2[i] * cos(theta)
          
          y1 <- s * max(trait1) * sin(theta) + dy + (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
          y2 <- s * max(trait1) * sin(theta) + dy - (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
          y3 <- y2 - s * trait1[i] * sin(theta)
          y4 <- y1 - s * trait1[i] * sin(theta)
          y5 <- y2 + s * trait2[i] * sin(theta)
          y6 <- y1 + s * trait2[i] * sin(theta)
          
          x1w <- s * max(trait1) * cos(theta) + dx - wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
          x2w <- s * max(trait1) * cos(theta) + dx + wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
          y1w <- s * max(trait1) * sin(theta) + dy + wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
          y2w <- s * max(trait1) * sin(theta) + dy - wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
          
          x01 <- x2w - s * max(trait1) * cos(theta)
          x02 <- x1w - s * max(trait1) * cos(theta)
          x03 <- x2w + s * max(trait1) * cos(theta)
          x04 <- x1w + s * max(trait1) * cos(theta)
          y01 <- y2w - s * max(trait1) * sin(theta)
          y02 <- y1w - s * max(trait1) * sin(theta)
          y03 <- y2w + s * max(trait1) * sin(theta)
          y04 <- y1w + s * max(trait1) * sin(theta)
        }
        ## trait 1: Quadrate
        polygon(c(x1, x2, x3, x4), c(y1, y2, y3, y4), 
                col = quad.cols$col[which(abs(trait1[i]-quad.cols$rql)==min(abs(trait1[i]-quad.cols$rql)))], border = "black",lwd=0.01)
        ## trait 2: MA
        polygon(c(x1, x2, x5, x6), c(y1, y2, y5, y6), 
                col = MA.cols$col[which(abs(trait2[i]-MA.cols$MA)==min(abs(trait2[i]-MA.cols$MA)))], border = "black",lwd=0.01)
      }
      hot_legend_image <- as.raster(matrix(rev(hot_color_func(100)), ncol=100))
      rasterImage(hot_legend_image, -h ,-7.5, -5,-10)
      polygon(col = NA, border = "black", lwd=0.25, x=c(-h,-h,-5,-5), y = c(-10,-7.5,-7.5,-10))
      text(x=mean(c(-h,-5)),y=-5.5,"Mechanical advantage (MA)",cex=0.75)
      text(x=-h+5, y = -11.5, cex = 0.5, "More force-modified")
      text(x=-10, y = -11.5, cex = 0.5, "Less force-modified")
      arrows(-22.5,-11.5,-h+17.5,-11.5, lwd=0.75, angle = 30, length=0.05, code = 3, col = "black")
      
      cold_legend_image <- as.raster(matrix(rev(cold_color_func(100)), ncol=100))
      rasterImage(cold_legend_image, h,-7.5,5,-10)
      polygon(col = NA, border = "black", lwd=0.25, x=c(h,h,5,5), y = c(-10,-7.5,-7.5,-10))
      text(x=mean(c(h,5)),y=-5.5,"Relative quadrate length (RQL)",cex=0.75)
      text(x=h-5, y = -11.5, cex = 0.5, "Larger gape")
      text(x=10, y = -11.5, cex = 0.5, "Smaller gape")
      arrows(20,-11.5,h-15,-11.5, lwd=0.75, angle = 30, length=0.05, code = 3, col = "black")
      
      if(diet.legend){
        legend(leg=c("Non-elongate
tetrapods","Lizards","Elongate
squamates"), 
               fill =state.colors[c("other","lizards","snakes")], ncol = 1, x=-19, y = 25, cex=0.6, text.width = 16, box.lty = 0, bg = NA, horiz = TRUE
        )
        legend(leg=c("Non-elongate
fishes","Anguilliform
vertebrates","Crustaceans"), 
               fill =state.colors[c("fish","eels","crusts")], ncol = 1, x=-19, y = 16.5, cex=0.6, text.width = 16, box.lty = 0, bg = NA, horiz = TRUE
        )
        legend(leg=c("Gastropods","Earthworms","Arthropods"), 
               fill =c(state.colors[c("mollusks","worms","bugs")]), ncol = 1, x=-19, y = 7, cex=0.6, text.width = 16, box.lty = 0, bg = NA, horiz = TRUE,
               border = c("black")
        )
        legend(leg=c(NA, NA, NA), 
               fill = NA, border = NA, ncol = 1, x=-19, y = 29, cex=0.65, text.width = 16, box.lty = 0, bg = NA, horiz = TRUE, title = "Diet", title.cex = 0.8
        )
      }
      
      bgcol = rgb(255,255,255, max = 255, alpha = 175)
      pie.adj = NULL
      nodelabels("1",node=viper.mrca,frame="circle",bg=bgcol,cex=0.5,adj=pie.adj)  ## 1 = viperids
      nodelabels("2",node=getMRCA(map, c("Brachyorrhos_albus","Homalopsis_buccata")),frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 2 = homalopsids
      nodelabels("3",node=getMRCA(map, c("Atractaspis_bibronii","Micrurus_fulvius")),frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 3 = elapoids
      nodelabels("4",node=getMRCA(map, c("Langaha_madagascariensis","Duberria_lutrix")),frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 4 = psuedoxyrhophiids
      nodelabels("5",node=getMRCA(map, c("Atractaspis_bibronii","Aparallactus_werneri")),frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 5 = atractaspididsa
      nodelabels("6",node=elapid.mrca,frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ##  6 = elapids
      nodelabels("7",node=getMRCA(map, c("Coluber_constrictor","Carphophis_amoenus")),frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 7 = colubroids
      nodelabels("8",node=dipsadid.mrca,frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 8 = dipsadids
      nodelabels("9",node=natricid.mrca,frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 9 = natricids
      nodelabels(":0",node=colubrid.mrca,frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 10 = colubrids
      
      points(x=rep(0,5), y = c(seq(17,4,by=-3.25)), pch = 21, cex = 1.5)
      text(x=rep(0,5), y = c(seq(17,4,by=-3.25)), cex = 0.5, c(1:5))
      text(x=rep(0,5), y = c(seq(17,4,by=-3.25)), cex = 0.5, c("Viperidae","Homalopsidae","Elapoidea","Pseudoxyrhophiidae","Atractaspididiae"), pos=4)
      points(x=rep(23,5), y = c(seq(17,4,by=-3.25)), pch = 21, cex = 1.5)
      text(x=rep(23,5), y = c(seq(17,4,by=-3.25)), cex = 0.5, c(6:10))
      text(x=rep(23,5), y = c(seq(17,4,by=-3.25)), cex = 0.5, c("Elapidae","Colubroidea","Dipsadidae","Natricidae","Colubridae"), pos=4)
    }
    dev.off()
  }
}

## make phylogeny trait plot and diet tree with families, merging colors states weighted by their probabilities at a specific node
{
  ## set up
  {
    dat <- read.csv("hOUwie/skull_dat.csv",row.names=1)
    dat = dat[!grepl("Acrochordus",dat$species),]
    dat$color[dat$diet=="other verts"] = "gray50"
    dat$color[dat$diet=="lizards"] = "forestgreen"
    dat$color[dat$diet=="mollusks"] = "saddlebrown"
    dat$color[dat$diet=="fish"] = "royalblue2"
    dat$color[dat$diet=="eels"] = "turquoise2"
    dat$color[dat$diet=="worms"] = "lightpink2"
    dat$color[dat$diet=="snakes"] = "mediumpurple2"
    dat$color[dat$diet=="crusts"] = "tomato2"
    dat$color[dat$diet=="bugs"] = "goldenrod2"
    
    
    state_colors <- c(other = "gray50",
                      lizards="forestgreen",
                      snakes = "mediumpurple2",
                      fish="royalblue2",
                      eels = "turquoise2",
                      mollusks = "saddlebrown",
                      worms = "lightpink2",
                      bugs = "goldenrod2",
                      crusts="tomato2"
    )
    
    diet.data <- data.frame(worms=ifelse(diet=="worms",1,0),
                            bugs=ifelse(diet=="bugs",1,0),
                            other=ifelse(diet=="other verts",1,0),
                            snakes=ifelse(diet=="snakes",1,0),
                            mollusks=ifelse(diet=="mollusks",1,0),
                            crusts=ifelse(diet=="crusts",1,0),
                            eels=ifelse(diet=="eels",1,0),
                            fish=ifelse(diet=="fish",1,0),
                            lizards=ifelse(diet=="lizards",1,0)
    )
    
    allStates<-rbind(diet.data[tree$tip.label,colnames(hrecon)],
                     hrecon)
    
    map = reorder(tree, by= "clade")
    map<-paintSubTree(map,node=Ntip(map)+1,
                      state=names(which(hrecon[1,]==max(hrecon[1,]))))
    for(i in 1:nrow(map$edge)){
      states<-sapply(map$edge[i,],function(x,aa) 
        names(aa)[which(aa[x,]==max(aa[x,]))],
        aa=allStates)
      if(length(unique(states))!=1)
        map<-paintSubTree(map,node=map$edge[i,2],state=states[2],
                          stem=0.5)
    }
    
    eco <- readxl::read_xlsx("input_sheets/snake_species_ecology_discrete.xlsx")
    vipers <- eco$species[eco$family == "Viperidae" & eco$species %in% tree$tip.label]
    elapids <- eco$species[eco$family == "Elapidae" & eco$species %in% tree$tip.label]
    dipsadids <- eco$species[eco$family == "Dipsadidae" & eco$species %in% tree$tip.label]
    natricids <- eco$species[eco$family == "Natricidae" & eco$species %in% tree$tip.label]
    colubrids <- eco$species[eco$family == "Colubridae" & eco$species %in% tree$tip.label]
    
    viper.mrca <- getMRCA(tree, vipers)
    elapid.mrca <- getMRCA(tree, elapids)
    dipsadid.mrca <- getMRCA(tree, dipsadids)
    natricid.mrca <- getMRCA(tree, natricids)
    colubrid.mrca <- getMRCA(tree, c("Coluber_constrictor","Ahaetulla_prasina"))
    
    legend.text <- c("Nonelongate tetrapods","Lizards","Elongate squamates","Nonelognate fishes","Anguillform vertebrates","Gastropods","Earthworms","Arthropods","Crustaceans")
    
    phy <- ladderize(reconmod$phy)
    num_nodes <- Nnode(phy) + Ntip(phy)  # Total number of nodes (internal + tips)
    node_probs <- hrecon
    
    #node_probs <- matrix(runif(num_nodes * num_states), nrow = num_nodes, ncol = num_states)
    node_probs <- t(apply(node_probs, 1, function(x) x / sum(x)))  # Normalize to sum to 1 per node
    weighted_avg_color <- function(colors, weights) {
      # Convert each color to RGB and multiply by the corresponding weight
      rgb_vals <- rowSums(sapply(seq_along(colors), function(i) {
        col2rgb(colors[i]) * weights[i]
      })) / 255  
      # Create final color in hex format
      rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3])
    }
    
    node_colors=vector()
    colors = state_colors[colnames(node_probs)]
    for(i in 1:nrow(node_probs)){
      node_colors[i] = weighted_avg_color(colors, as.vector(node_probs[i,]))
    }
    node_colors=c(dat$color, node_colors)
    
    
    map2 <- map
    map3 = map2
    me = map3$mapped.edge
    nc_df = data.frame(names=node_colors,colors=node_colors)
    nc_df = nc_df[!duplicated(nc_df),]
    rownames(nc_df)=nc_df$names
    nc_df$names=NULL
    
    me = as.data.frame(matrix(data=0,nrow=nrow(me),ncol=(9+nrow(nc_df))))
    colnames(me)=c(colnames(map3$mapped.edge),nc_df$colors)
    rownames(me)=rownames(map3$mapped.edge)
    
    for(i in 1:nrow(me)){
      start_node=as.numeric(str_split(rownames(me)[i],",")[[1]][1])
      if(map3$edge[i,2]<=Ntip(phy)){
        this_state = dat$diet[dat$species==phy$tip.label[map3$edge[i,2]]]
        if(this_state=="other verts") this_state = "other"
        me[i,this_state] = map3$edge.length[i]/2
        me[i,node_colors[start_node]] = map3$edge.length[i]/2
      }else{
        end_node=as.numeric(str_split(rownames(me)[i],",")[[1]][2])
        if(node_colors[start_node]!=node_colors[end_node]){
          me[i,node_colors[start_node]] = map3$edge.length[i]/2
          me[i,node_colors[end_node]] = map3$edge.length[i]/2
        }else me[i,node_colors[end_node]] = map3$edge.length[i]
      }
    }
    
    mps = map3$maps
    
    for(i in 1:nrow(me)){
      this_row = me[i,]
      start_node=as.numeric(str_split(rownames(this_row),",")[[1]][1])
      cn=colnames(this_row)[this_row!=0]
      this_row = this_row[,this_row!=0]
      the_vector = setNames(c(t(this_row)), cn)
      if(names(the_vector)[1]!=node_colors[start_node]) the_vector = rev(the_vector)
      mps[[i]] = the_vector
    }
    
    map3$mapped.edge=me
    map3$maps=mps
    state_colors2=c(state_colors, nc_df$colors)
    names(state_colors2)=c(names(state_colors),nc_df$colors)
  }
  
  ## dietary reconstruction with families
  pdf("figs/tree_dietary_states_with_families_weighted_state_colors.pdf")
  {
    par(bg="white")
    plot(map3,state_colors2,ftype="off",lwd=1.5,type="fan",part=0.5)
    nodelabels(pie=hrecon,cex=0.15,piecol = state.colors[order(names(state.colors))])
    add.simmap.legend(leg=legend.text, colors=state.colors,prompt=FALSE, x=0.9*par()$usr[1],
                      y=1.1*max(nodeHeights(tree)),fsize=0.7)
    text("Diet",x=0.86*par()$usr[1], y=1.15*max(nodeHeights(tree)))
    arc.cladelabels(tree, "Viperidae", node = viper.mrca, cex = 0.8, mark.node=FALSE)
    arc.cladelabels(tree, "Elapidae", node = elapid.mrca, cex = 0.8, mark.node=FALSE)
    arc.cladelabels(tree, "Dipsadidae", node = dipsadid.mrca, cex = 0.8, mark.node=FALSE)
    arc.cladelabels(tree, "Natricidae", node = natricid.mrca, cex = 0.8, mark.node=FALSE)
    arc.cladelabels(tree, "Colubridae", node = colubrid.mrca, cex = 0.8, mark.node=FALSE)
    
    par(bg="white")
    plot(map3,state_colors2,ftype="i",fsize=0.175,offset=2,lwd=1.5,type="fan",part=1)
    nodelabels(pie=hrecon,cex=0.15,piecol = state.colors[order(names(state.colors))])
    add.simmap.legend(leg=legend.text, colors=state.colors,prompt=FALSE, x=0.9*par()$usr[1],
                      y=1.1*max(nodeHeights(tree)),fsize=0.7)
    text("Diet",x=0.86*par()$usr[1], y=1.17*max(nodeHeights(tree)))
    
    par(bg="white")
    plot(map3,state_colors2,ftype="i",fsize=0.175,offset=2,lwd=1.5,type="fan",part=1)
    nodelabels(pie=hrecon,cex=0.15,piecol = state.colors[order(names(state.colors))])
    nodelabels(bg="black",col="white",frame="circle",cex=0.1)
  }
  dev.off()
  
  
  
  dat <- read.csv("hOUwie/skull_dat.csv",row.names=1)
  dat = dat[!grepl("Acrochordus",dat$species),]
  dat$color[dat$diet=="other verts"] = "gray90"
  dat$color[dat$diet=="lizards"] = "forestgreen"
  dat$color[dat$diet=="mollusks"] = "saddlebrown"
  dat$color[dat$diet=="fish"] = "royalblue2"
  dat$color[dat$diet=="eels"] = "turquoise2"
  dat$color[dat$diet=="worms"] = "lightpink2"
  dat$color[dat$diet=="snakes"] = "mediumpurple2"
  dat$color[dat$diet=="crusts"] = "tomato2"
  dat$color[dat$diet=="bugs"] = "goldenrod2"
  
  ## dietary reconstruction with MA and RQL plotted
  ## set up
  {
    par(bg="white")
    trait1 = setNames((dat$RQL),dat$species)[map$tip.label]
    trait1[trait1>quantile(trait1,0.99)] <- quantile(trait1, 0.99)
    trait2 = setNames((dat$MA), dat$species)[map$tip.label]
    trait2[trait2<quantile(trait2,0.01)] <- quantile(trait2, 0.01)
    
    trait2 <- (trait2-min(trait2))/(max(trait2)-min(trait2))
    trait1 <- (trait1-min(trait1))/(max(trait1)-min(trait1))
    
    all(names(trait2)==names(trait1))
    map2 <- map
    all(map2$tip.label==map$tip.label)
    
    trait1 = 2+7*trait1
    trait2 = 2+7*trait2
    
    desat <- function(cols, sat=0.5) {
      X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
      hsv(X[1,], X[2,], X[3,])
    }
    hot_color_func <- viridis::magma
    hot_color_func <- colorRampPalette(c("#000004FF","#150E37FF", "#3B0F70FF", "#641A80FF", "#8C2981FF", "#B63679FF", "#DE4968FF", "#F76F5CFF","#FE9F6D","#FECF6E","#FFFF6F"))
    cold_color_func <- colorRampPalette(c("black","#000040", "navyblue", "#183487", "#31688EFF", "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF","#B4DE2CFF","#D8E228"))
    MA.cols <- data.frame(MA=seq(quantile(trait2,0.01),quantile(trait2,0.99),length.out=length(trait2)),col=hot_color_func(length(trait2)))
    quad.cols <- data.frame(rql=seq(quantile(trait1,0.01),quantile(trait1,0.99),length.out=length(trait1)),col=cold_color_func(length(trait1)))
    wider=1
    bg.thinner=0.5
    dx.mult <- 1.005
    dy.mult <- 1.005
    trait.border = 0.25
  }
  
  pdf("figs/tree_reconstruction_with_cont_data_weighted_state_colors.pdf",
      width = 13,
      height = 13
      )
  {
    h <- max(nodeHeights(map3))
    state.colors["other"] <- "gray50"
    capture.output(plot(map3, ftype="off",colors=state_colors2, 
                        xlim = c(-2.25*h,2.25*h),ylim = c(0,h), 
                        lwd = 1, type = "fan", part=0.5))
    obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    sw <- strwidth("l")
    w <- (par()$usr[4] - par()$usr[3])/(max(c(max(trait1)/max(nodeHeights(tree)), 1)) * length(tree$tip.label))*0.55
    
    skulls = paste("figs/skulls/",list.files("figs/skulls"),sep="")
    sp = gsub(".png","",word(skulls,3,3,sep="/"))
    sp = dat$species[dat$species %in% sp]
    skulls = paste("figs/skulls/",sp,".png",sep="")
    skdat = read.csv("figs/skulls_for_figure.csv")
    
    dat2 = dat
    dat2$color[dat2$diet=="other verts"] = "gray50"

    r = 1.7*h
    angles = pi/16 + seq(from=0, to = 7*pi/8, length.out = 9)
    xc = cos(angles)*r
    yc = sin(angles)*r
    yc[c(2,8)]=0.46*r
    yc[c(3,7)]=0.695*r
    yc[c(4,6)]=0.907*r
    xc = sqrt(r^2-yc^2)
    xc = xc^1.01
    xc[6:9]=-1*xc[6:9]
    xc
    
    for(i in 1:9){
      sk = readPNG(skulls[i])
      s = sp[i]
      catnum = skdat$specimen[skdat$species==s]
      dietstate = skdat$diet[skdat$diet==s]
      sx = ncol(sk)/2/sqrt(3000)
      sy = nrow(sk)/2/sqrt(3000)
      cx = xc[i]
      cy = yc[i]
      rasterImage(sk, 
                  xleft = cx - sx, xright = cx + sx,
                  ybottom = cy - sy, ytop = cy + sy
      )
      
      stext = gsub("_"," ",s)
      spdat = skdat[skdat$species==s,]
      
      text(x=cx, y = cy - 1.3*sy, labels = bquote(italic(.(stext))),col=dat2$color[dat$species==s],
           cex=0.6)
      text(x=cx, y = cy - 1.3*sy - 2.5, labels = spdat$specimen,col=dat2$color[dat2$species==s],
           cex=0.45)
      text(x=cx, y = cy - 1.3*sy - 5, labels = spdat$diet,col=dat2$color[dat2$species==s],
           cex=0.6)
      
      px0 = cos((which(map3$tip.label==s)-ifelse((which(map3$tip.label==s)/Ntip(map3))<0.5, +0.625, -0))/Ntip(map3)*pi)*1.315*h
      py0 = sin((which(map3$tip.label==s)-ifelse((which(map3$tip.label==s)/Ntip(map3))<0.5, +0.625, -0))/Ntip(map3)*pi)*1.315*h
      if(i >= 7){
        px1 = cx + strwidth(stext, cex = 0.6)/2.2
        py1 = cy - 1.3*sy - 2.5
      }else if(i >= 4){
        px1 = cx
        py1 = cy - 1.3*sy - 6.5
      }else if(i <= 3){
        px1 = cx - strwidth(stext, cex = 0.6)/2.2
        py1 = cy - 1.3*sy - 2.5
      }
      segments(x0 = px0, y0 = py0,
               y1 = py1, x1 = px1, 
               col = dat2$color[dat2$species==s])
    }
    
    for(i in 1:length(trait1)){
      {
        theta <- atan(obj$yy[i]/obj$xx[i])
        
        s <- ifelse(obj$xx[i] > 0, 1, -1)
        dx <- dx.mult *s * h * cos(theta) + s * cos(theta) * sw
        dy <- dy.mult * s * h * sin(theta) + s * sin(theta) * sw
        
        x1 <- s * max(trait1) * cos(theta) + dx - (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
        x2 <- s * max(trait1) * cos(theta) + dx + (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
        x3 <- x2 - s * trait1[i] * cos(theta)
        x4 <- x1 - s * trait1[i] * cos(theta)
        x5 <- x2 + s * trait2[i] * cos(theta)
        x6 <- x1 + s * trait2[i] * cos(theta)
        
        y1 <- s * max(trait1) * sin(theta) + dy + (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
        y2 <- s * max(trait1) * sin(theta) + dy - (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
        y3 <- y2 - s * trait1[i] * sin(theta)
        y4 <- y1 - s * trait1[i] * sin(theta)
        y5 <- y2 + s * trait2[i] * sin(theta)
        y6 <- y1 + s * trait2[i] * sin(theta)
        
        x1w <- s * max(trait1) * cos(theta) + dx - wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
        x2w <- s * max(trait1) * cos(theta) + dx + wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
        y1w <- s * max(trait1) * sin(theta) + dy + wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
        y2w <- s * max(trait1) * sin(theta) + dy - wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
        
        x01 <- x2w - s * max(trait1) * cos(theta)
        x02 <- x1w - s * max(trait1) * cos(theta)
        x03 <- x2w + s * max(trait1) * cos(theta)
        x04 <- x1w + s * max(trait1) * cos(theta)
        y01 <- y2w - s * max(trait1) * sin(theta)
        y02 <- y1w - s * max(trait1) * sin(theta)
        y03 <- y2w + s * max(trait1) * sin(theta)
        y04 <- y1w + s * max(trait1) * sin(theta)
      }
      
      ## background state colors
      polygon(c(x02, x01, x03, x04), c(y02, y01, y03, y04), 
              col = desat(dat$color[i],bg.thinner), border = desat(dat$color[i],bg.thinner), lwd = 0.01)
      ## trait 1: Quadrate
      #polygon(c(x1, x2, x3, x4), c(y1, y2, y3, y4), 
      #        col = quad.cols$col[which(abs(trait1[i]-quad.cols$rql)==min(abs(trait1[i]-quad.cols$rql)))], border = "black",lwd=0.005)
      ## trait 2: MA
      #polygon(c(x1, x2, x5, x6), c(y1, y2, y5, y6), 
      #        col = MA.cols$col[which(abs(trait2[i]-MA.cols$MA)==min(abs(trait2[i]-MA.cols$MA)))], border = "black",lwd=0.005)
      
      if(i+1<=Ntip(map2)){
        j <- i+1
        
        if(theta >= 0){
          x1.j <- x02
          y1.j <- y02
          x2.j <- x04
          y2.j <- y04
        }else if(theta < 0){
          x1.j <- x01
          y1.j <- y01
          x2.j <- x03
          y2.j <- y03
        }
        
        theta <- atan(obj$yy[j]/obj$xx[j])
        
        if(theta >= 0){
          {
            s <- ifelse(obj$xx[j] > 0, 1, -1)
            dx <- dx.mult *s * h * cos(theta) + s * cos(theta) * sw
            dy <- dy.mult * s * h * sin(theta) + s * sin(theta) * sw
            
            x1 <- s * max(trait1) * cos(theta) + dx - (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
            x2 <- s * max(trait1) * cos(theta) + dx + (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
            x3 <- x2 - s * trait1[j] * cos(theta)
            x4 <- x1 - s * trait1[j] * cos(theta)
            x5 <- x2 + s * trait2[j] * cos(theta)
            x6 <- x1 + s * trait2[j] * cos(theta)
            
            y1 <- s * max(trait1) * sin(theta) + dy + (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
            y2 <- s * max(trait1) * sin(theta) + dy - (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
            y3 <- y2 - s * trait1[j] * sin(theta)
            y4 <- y1 - s * trait1[j] * sin(theta)
            y5 <- y2 + s * trait2[j] * sin(theta)
            y6 <- y1 + s * trait2[j] * sin(theta)
            
            x1w <- s * max(trait1) * cos(theta) + dx - wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
            x2w <- s * max(trait1) * cos(theta) + dx + wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
            y1w <- s * max(trait1) * sin(theta) + dy + wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
            y2w <- s * max(trait1) * sin(theta) + dy - wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
            
            x01 <- x2w - s * max(trait1) * cos(theta)
            x02 <- x1w - s * max(trait1) * cos(theta)
            x03 <- x2w + s * max(trait1) * cos(theta)
            x04 <- x1w + s * max(trait1) * cos(theta)
            y01 <- y2w - s * max(trait1) * sin(theta)
            y02 <- y1w - s * max(trait1) * sin(theta)
            y03 <- y2w + s * max(trait1) * sin(theta)
            y04 <- y1w + s * max(trait1) * sin(theta)
          }
          
          x3.j <- x01
          y3.j <- y01
          x4.j <- x03
          y4.j <- y03
          
          x5.j <- mean(c(x1.j, x3.j))
          x6.j <- mean(c(x2.j, x4.j))
          y5.j <- mean(c(y1.j, y3.j))
          y6.j <- mean(c(y2.j, y4.j))
        }else if(theta<0){
          {
            s <- ifelse(obj$xx[j] > 0, 1, -1)
            dx <- dx.mult *s * h * cos(theta) + s * cos(theta) * sw
            dy <- dy.mult * s * h * sin(theta) + s * sin(theta) * sw
            
            x1 <- s * max(trait1) * cos(theta) + dx - (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
            x2 <- s * max(trait1) * cos(theta) + dx + (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
            x3 <- x2 - s * trait1[j] * cos(theta)
            x4 <- x1 - s * trait1[j] * cos(theta)
            x5 <- x2 + s * trait2[j] * cos(theta)
            x6 <- x1 + s * trait2[j] * cos(theta)
            
            y1 <- s * max(trait1) * sin(theta) + dy + (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
            y2 <- s * max(trait1) * sin(theta) + dy - (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
            y3 <- y2 - s * trait1[j] * sin(theta)
            y4 <- y1 - s * trait1[j] * sin(theta)
            y5 <- y2 + s * trait2[j] * sin(theta)
            y6 <- y1 + s * trait2[j] * sin(theta)
            
            x1w <- s * max(trait1) * cos(theta) + dx - wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
            x2w <- s * max(trait1) * cos(theta) + dx + wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
            y1w <- s * max(trait1) * sin(theta) + dy + wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
            y2w <- s * max(trait1) * sin(theta) + dy - wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
            
            x01 <- x2w - s * max(trait1) * cos(theta)
            x02 <- x1w - s * max(trait1) * cos(theta)
            x03 <- x2w + s * max(trait1) * cos(theta)
            x04 <- x1w + s * max(trait1) * cos(theta)
            y01 <- y2w - s * max(trait1) * sin(theta)
            y02 <- y1w - s * max(trait1) * sin(theta)
            y03 <- y2w + s * max(trait1) * sin(theta)
            y04 <- y1w + s * max(trait1) * sin(theta)
          }
          
          x3.j <- x02
          y3.j <- y02
          x4.j <- x04
          y4.j <- y04
          
          x5.j <- mean(c(x1.j, x3.j))
          x6.j <- mean(c(x2.j, x4.j))
          y5.j <- mean(c(y1.j, y3.j))
          y6.j <- mean(c(y2.j, y4.j))
          
        }
        polygon(c(x2.j,x1.j,x5.j,x6.j), c(y2.j,y1.j,y5.j,y6.j), col=desat(dat$color[i],bg.thinner), border=desat(dat$color[i],bg.thinner), lwd = 0.01)
        polygon(c(x4.j,x3.j,x5.j,x6.j), c(y4.j,y3.j,y5.j,y6.j), col=desat(dat$color[j],bg.thinner), border=desat(dat$color[i],bg.thinner), lwd = 0.01)
      }
      
      {
        theta <- atan(obj$yy[i]/obj$xx[i])
        
        s <- ifelse(obj$xx[i] > 0, 1, -1)
        dx <- dx.mult *s * h * cos(theta) + s * cos(theta) * sw
        dy <- dy.mult * s * h * sin(theta) + s * sin(theta) * sw
        
        x1 <- s * max(trait1) * cos(theta) + dx - (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
        x2 <- s * max(trait1) * cos(theta) + dx + (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
        x3 <- x2 - s * trait1[i] * cos(theta)
        x4 <- x1 - s * trait1[i] * cos(theta)
        x5 <- x2 + s * trait2[i] * cos(theta)
        x6 <- x1 + s * trait2[i] * cos(theta)
        
        y1 <- s * max(trait1) * sin(theta) + dy + (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
        y2 <- s * max(trait1) * sin(theta) + dy - (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
        y3 <- y2 - s * trait1[i] * sin(theta)
        y4 <- y1 - s * trait1[i] * sin(theta)
        y5 <- y2 + s * trait2[i] * sin(theta)
        y6 <- y1 + s * trait2[i] * sin(theta)
        
        x1w <- s * max(trait1) * cos(theta) + dx - wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
        x2w <- s * max(trait1) * cos(theta) + dx + wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
        y1w <- s * max(trait1) * sin(theta) + dy + wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
        y2w <- s * max(trait1) * sin(theta) + dy - wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
        
        x01 <- x2w - s * max(trait1) * cos(theta)
        x02 <- x1w - s * max(trait1) * cos(theta)
        x03 <- x2w + s * max(trait1) * cos(theta)
        x04 <- x1w + s * max(trait1) * cos(theta)
        y01 <- y2w - s * max(trait1) * sin(theta)
        y02 <- y1w - s * max(trait1) * sin(theta)
        y03 <- y2w + s * max(trait1) * sin(theta)
        y04 <- y1w + s * max(trait1) * sin(theta)
      }
      ## trait 1: Quadrate
      polygon(c(x1, x2, x3, x4), c(y1, y2, y3, y4), 
              col = quad.cols$col[which(abs(trait1[i]-quad.cols$rql)==min(abs(trait1[i]-quad.cols$rql)))], border = "black",lwd=trait.border)
      ## trait 2: MA
      polygon(c(x1, x2, x5, x6), c(y1, y2, y5, y6), 
              col = MA.cols$col[which(abs(trait2[i]-MA.cols$MA)==min(abs(trait2[i]-MA.cols$MA)))], border = "black",lwd=trait.border)
    }
    hot_legend_image <- as.raster(matrix(rev(hot_color_func(100)), ncol=100))
    rasterImage(hot_legend_image, -h ,-7.5, -5,-10)
    polygon(col = NA, border = "black", lwd=0.5, x=c(-h,-h,-5,-5), y = c(-10,-7.5,-7.5,-10))
    text(x=mean(c(-h,-5)),y=-5.5,"Mechanical advantage (MA)",cex=0.8)
    text(x=-h+5, y = -11.5, cex = 0.55, "More force-modified")
    text(x=-10, y = -11.5, cex = 0.55, "Less force-modified")
    arrows(-22.5,-11.5,-h+17.5,-11.5, lwd=0.75, angle = 30, length=0.05, code = 3, col = "black")
    
    cold_legend_image <- as.raster(matrix(rev(cold_color_func(100)), ncol=100))
    rasterImage(cold_legend_image, h,-7.5,5,-10)
    polygon(col = NA, border = "black", lwd=0.5, x=c(h,h,5,5), y = c(-10,-7.5,-7.5,-10))
    text(x=mean(c(h,5)),y=-5.5,"Relative quadrate length (RQL)",cex=0.8)
    text(x=h-5, y = -11.5, cex = 0.55, "Larger gape")
    text(x=10, y = -11.5, cex = 0.55, "Smaller gape")
    arrows(20,-11.5,h-15,-11.5, lwd=0.75, angle = 30, length=0.05, code = 3, col = "black")
    
    diet.legend = readPNG("figs/phylogeny_diet_legend.png")
    dil.x = ncol(diet.legend)/2/sqrt(3000)
    dil.y = nrow(diet.legend)/2/sqrt(3000)
    rasterImage(diet.legend,
                xleft = -1*dil.x,
                xright = dil.x,
                ybottom = -25-dil.y,
                ytop=-25+dil.y
                )
    rect(xleft = -1*dil.x,
         xright = dil.x,
         ybottom = -25-dil.y,
         ytop=-25++1.225*dil.y+2.5,
         border = "white",
         lwd = 0.5
         )
    
    text("Main prey item in diet",
         x=0, y = -25+dil.y,
         cex = 0.8, col = "black"
    )
    
    bgcol = rgb(255,255,255, max = 255, alpha = 175)
    pie.adj = NULL
    nodelabels("1",node=viper.mrca,frame="circle",bg=bgcol,cex=0.5,adj=pie.adj)  ## 1 = viperids
    nodelabels("2",node=getMRCA(map3, c("Brachyorrhos_albus","Homalopsis_buccata")),frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 2 = homalopsids
    nodelabels("3",node=getMRCA(map3, c("Atractaspis_bibronii","Micrurus_fulvius")),frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 3 = elapoids
    nodelabels("4",node=getMRCA(map3, c("Langaha_madagascariensis","Duberria_lutrix")),frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 4 = psuedoxyrhophiids
    nodelabels("5",node=getMRCA(map3, c("Atractaspis_bibronii","Aparallactus_werneri")),frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 5 = atractaspididsa
    nodelabels("6",node=elapid.mrca,frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ##  6 = elapids
    nodelabels("7",node=getMRCA(map3, c("Coluber_constrictor","Carphophis_amoenus")),frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 7 = colubroids
    nodelabels("8",node=dipsadid.mrca,frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 8 = dipsadids
    nodelabels("9",node=natricid.mrca,frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 9 = natricids
    nodelabels(" ",node=colubrid.mrca,frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 10 = colubrids
    
    points(x=rep(0,5), y = c(seq(18,5,by=-3.25)), pch = 21, cex = 1.5)
    text(x=rep(0,5), y = c(seq(18,5,by=-3.25)), cex = 0.5, c(1:5))
    text(x=rep(0,5), y = c(seq(18,5,by=-3.25)), cex = 0.5, c("Viperidae","Homalopsidae","Elapoidea","Pseudoxyrhophiidae","Atractaspididiae"), pos=4)
    points(x=rep(23,5), y = c(seq(18,5,by=-3.25)), pch = 21, cex = 1.5)
    text(x=rep(23,5), y = c(seq(18,5,by=-3.25)), cex = 0.5, c(6:10))
    text(x=rep(23,5), y = c(seq(18,5,by=-3.25)), cex = 0.5, c("Elapidae","Colubroidea","Dipsadidae","Natricidae","Colubridae"), pos=4)
    
  }
  dev.off()
    
  png("figs/tree_reconstruction_with_cont_data_weighted_state_colors.png",
      width = 13,
      height = 13,
      units = "in",
      res = 1200)
  {
    h <- max(nodeHeights(map3))
    state.colors["other"] <- "gray50"
    capture.output(plot(map3, ftype="off",colors=state_colors2, 
                        xlim = c(-2.25*h,2.25*h),ylim = c(0,h), 
                        lwd = 1, type = "fan", part=0.5))
    obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    sw <- strwidth("l")
    w <- (par()$usr[4] - par()$usr[3])/(max(c(max(trait1)/max(nodeHeights(tree)), 1)) * length(tree$tip.label))*0.55
    
    skulls = paste("figs/skulls/",list.files("figs/skulls"),sep="")
    sp = gsub(".png","",word(skulls,3,3,sep="/"))
    sp = dat$species[dat$species %in% sp]
    skulls = paste("figs/skulls/",sp,".png",sep="")
    skdat = read.csv("figs/skulls_for_figure.csv")
    
    dat2 = dat
    dat2$color[dat2$diet=="other verts"] = "gray50"
    
    r = 1.7*h
    angles = pi/16 + seq(from=0, to = 7*pi/8, length.out = 9)
    xc = cos(angles)*r
    yc = sin(angles)*r
    yc[c(2,8)]=0.46*r
    yc[c(3,7)]=0.695*r
    yc[c(4,6)]=0.907*r
    xc = sqrt(r^2-yc^2)
    xc = xc^1.01
    xc[6:9]=-1*xc[6:9]
    xc
    
    for(i in 1:9){
      sk = readPNG(skulls[i])
      s = sp[i]
      catnum = skdat$specimen[skdat$species==s]
      dietstate = skdat$diet[skdat$diet==s]
      sx = ncol(sk)/2/sqrt(3000)
      sy = nrow(sk)/2/sqrt(3000)
      cx = xc[i]
      cy = yc[i]
      rasterImage(sk, 
                  xleft = cx - sx, xright = cx + sx,
                  ybottom = cy - sy, ytop = cy + sy
      )
      
      stext = gsub("_"," ",s)
      spdat = skdat[skdat$species==s,]
      
      text(x=cx, y = cy - 1.3*sy, labels = bquote(italic(.(stext))),col=dat2$color[dat$species==s],
           cex=0.6)
      text(x=cx, y = cy - 1.3*sy - 2.5, labels = spdat$specimen,col=dat2$color[dat2$species==s],
           cex=0.45)
      text(x=cx, y = cy - 1.3*sy - 5, labels = spdat$diet,col=dat2$color[dat2$species==s],
           cex=0.6)
      
      px0 = cos((which(map3$tip.label==s)-ifelse((which(map3$tip.label==s)/Ntip(map3))<0.5, +0.625, -0))/Ntip(map3)*pi)*1.315*h
      py0 = sin((which(map3$tip.label==s)-ifelse((which(map3$tip.label==s)/Ntip(map3))<0.5, +0.625, -0))/Ntip(map3)*pi)*1.315*h
      if(i >= 7){
        px1 = cx + strwidth(stext, cex = 0.6)/2.2
        py1 = cy - 1.3*sy - 2.5
      }else if(i >= 4){
        px1 = cx
        py1 = cy - 1.3*sy - 6.5
      }else if(i <= 3){
        px1 = cx - strwidth(stext, cex = 0.6)/2.2
        py1 = cy - 1.3*sy - 2.5
      }
      segments(x0 = px0, y0 = py0,
               y1 = py1, x1 = px1, 
               col = dat2$color[dat2$species==s])
    }
    
    for(i in 1:length(trait1)){
      {
        theta <- atan(obj$yy[i]/obj$xx[i])
        
        s <- ifelse(obj$xx[i] > 0, 1, -1)
        dx <- dx.mult *s * h * cos(theta) + s * cos(theta) * sw
        dy <- dy.mult * s * h * sin(theta) + s * sin(theta) * sw
        
        x1 <- s * max(trait1) * cos(theta) + dx - (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
        x2 <- s * max(trait1) * cos(theta) + dx + (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
        x3 <- x2 - s * trait1[i] * cos(theta)
        x4 <- x1 - s * trait1[i] * cos(theta)
        x5 <- x2 + s * trait2[i] * cos(theta)
        x6 <- x1 + s * trait2[i] * cos(theta)
        
        y1 <- s * max(trait1) * sin(theta) + dy + (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
        y2 <- s * max(trait1) * sin(theta) + dy - (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
        y3 <- y2 - s * trait1[i] * sin(theta)
        y4 <- y1 - s * trait1[i] * sin(theta)
        y5 <- y2 + s * trait2[i] * sin(theta)
        y6 <- y1 + s * trait2[i] * sin(theta)
        
        x1w <- s * max(trait1) * cos(theta) + dx - wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
        x2w <- s * max(trait1) * cos(theta) + dx + wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
        y1w <- s * max(trait1) * sin(theta) + dy + wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
        y2w <- s * max(trait1) * sin(theta) + dy - wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
        
        x01 <- x2w - s * max(trait1) * cos(theta)
        x02 <- x1w - s * max(trait1) * cos(theta)
        x03 <- x2w + s * max(trait1) * cos(theta)
        x04 <- x1w + s * max(trait1) * cos(theta)
        y01 <- y2w - s * max(trait1) * sin(theta)
        y02 <- y1w - s * max(trait1) * sin(theta)
        y03 <- y2w + s * max(trait1) * sin(theta)
        y04 <- y1w + s * max(trait1) * sin(theta)
      }
      
      ## background state colors
      polygon(c(x02, x01, x03, x04), c(y02, y01, y03, y04), 
              col = desat(dat$color[i],bg.thinner), border = desat(dat$color[i],bg.thinner), lwd = 0.01)
      ## trait 1: Quadrate
      #polygon(c(x1, x2, x3, x4), c(y1, y2, y3, y4), 
      #        col = quad.cols$col[which(abs(trait1[i]-quad.cols$rql)==min(abs(trait1[i]-quad.cols$rql)))], border = "black",lwd=0.005)
      ## trait 2: MA
      #polygon(c(x1, x2, x5, x6), c(y1, y2, y5, y6), 
      #        col = MA.cols$col[which(abs(trait2[i]-MA.cols$MA)==min(abs(trait2[i]-MA.cols$MA)))], border = "black",lwd=0.005)
      
      if(i+1<=Ntip(map2)){
        j <- i+1
        
        if(theta >= 0){
          x1.j <- x02
          y1.j <- y02
          x2.j <- x04
          y2.j <- y04
        }else if(theta < 0){
          x1.j <- x01
          y1.j <- y01
          x2.j <- x03
          y2.j <- y03
        }
        
        theta <- atan(obj$yy[j]/obj$xx[j])
        
        if(theta >= 0){
          {
            s <- ifelse(obj$xx[j] > 0, 1, -1)
            dx <- dx.mult *s * h * cos(theta) + s * cos(theta) * sw
            dy <- dy.mult * s * h * sin(theta) + s * sin(theta) * sw
            
            x1 <- s * max(trait1) * cos(theta) + dx - (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
            x2 <- s * max(trait1) * cos(theta) + dx + (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
            x3 <- x2 - s * trait1[j] * cos(theta)
            x4 <- x1 - s * trait1[j] * cos(theta)
            x5 <- x2 + s * trait2[j] * cos(theta)
            x6 <- x1 + s * trait2[j] * cos(theta)
            
            y1 <- s * max(trait1) * sin(theta) + dy + (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
            y2 <- s * max(trait1) * sin(theta) + dy - (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
            y3 <- y2 - s * trait1[j] * sin(theta)
            y4 <- y1 - s * trait1[j] * sin(theta)
            y5 <- y2 + s * trait2[j] * sin(theta)
            y6 <- y1 + s * trait2[j] * sin(theta)
            
            x1w <- s * max(trait1) * cos(theta) + dx - wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
            x2w <- s * max(trait1) * cos(theta) + dx + wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
            y1w <- s * max(trait1) * sin(theta) + dy + wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
            y2w <- s * max(trait1) * sin(theta) + dy - wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
            
            x01 <- x2w - s * max(trait1) * cos(theta)
            x02 <- x1w - s * max(trait1) * cos(theta)
            x03 <- x2w + s * max(trait1) * cos(theta)
            x04 <- x1w + s * max(trait1) * cos(theta)
            y01 <- y2w - s * max(trait1) * sin(theta)
            y02 <- y1w - s * max(trait1) * sin(theta)
            y03 <- y2w + s * max(trait1) * sin(theta)
            y04 <- y1w + s * max(trait1) * sin(theta)
          }
          
          x3.j <- x01
          y3.j <- y01
          x4.j <- x03
          y4.j <- y03
          
          x5.j <- mean(c(x1.j, x3.j))
          x6.j <- mean(c(x2.j, x4.j))
          y5.j <- mean(c(y1.j, y3.j))
          y6.j <- mean(c(y2.j, y4.j))
        }else if(theta<0){
          {
            s <- ifelse(obj$xx[j] > 0, 1, -1)
            dx <- dx.mult *s * h * cos(theta) + s * cos(theta) * sw
            dy <- dy.mult * s * h * sin(theta) + s * sin(theta) * sw
            
            x1 <- s * max(trait1) * cos(theta) + dx - (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
            x2 <- s * max(trait1) * cos(theta) + dx + (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
            x3 <- x2 - s * trait1[j] * cos(theta)
            x4 <- x1 - s * trait1[j] * cos(theta)
            x5 <- x2 + s * trait2[j] * cos(theta)
            x6 <- x1 + s * trait2[j] * cos(theta)
            
            y1 <- s * max(trait1) * sin(theta) + dy + (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
            y2 <- s * max(trait1) * sin(theta) + dy - (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
            y3 <- y2 - s * trait1[j] * sin(theta)
            y4 <- y1 - s * trait1[j] * sin(theta)
            y5 <- y2 + s * trait2[j] * sin(theta)
            y6 <- y1 + s * trait2[j] * sin(theta)
            
            x1w <- s * max(trait1) * cos(theta) + dx - wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
            x2w <- s * max(trait1) * cos(theta) + dx + wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
            y1w <- s * max(trait1) * sin(theta) + dy + wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
            y2w <- s * max(trait1) * sin(theta) + dy - wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
            
            x01 <- x2w - s * max(trait1) * cos(theta)
            x02 <- x1w - s * max(trait1) * cos(theta)
            x03 <- x2w + s * max(trait1) * cos(theta)
            x04 <- x1w + s * max(trait1) * cos(theta)
            y01 <- y2w - s * max(trait1) * sin(theta)
            y02 <- y1w - s * max(trait1) * sin(theta)
            y03 <- y2w + s * max(trait1) * sin(theta)
            y04 <- y1w + s * max(trait1) * sin(theta)
          }
          
          x3.j <- x02
          y3.j <- y02
          x4.j <- x04
          y4.j <- y04
          
          x5.j <- mean(c(x1.j, x3.j))
          x6.j <- mean(c(x2.j, x4.j))
          y5.j <- mean(c(y1.j, y3.j))
          y6.j <- mean(c(y2.j, y4.j))
          
        }
        polygon(c(x2.j,x1.j,x5.j,x6.j), c(y2.j,y1.j,y5.j,y6.j), col=desat(dat$color[i],bg.thinner), border=desat(dat$color[i],bg.thinner), lwd = 0.01)
        polygon(c(x4.j,x3.j,x5.j,x6.j), c(y4.j,y3.j,y5.j,y6.j), col=desat(dat$color[j],bg.thinner), border=desat(dat$color[i],bg.thinner), lwd = 0.01)
      }
      
      {
        theta <- atan(obj$yy[i]/obj$xx[i])
        
        s <- ifelse(obj$xx[i] > 0, 1, -1)
        dx <- dx.mult *s * h * cos(theta) + s * cos(theta) * sw
        dy <- dy.mult * s * h * sin(theta) + s * sin(theta) * sw
        
        x1 <- s * max(trait1) * cos(theta) + dx - (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
        x2 <- s * max(trait1) * cos(theta) + dx + (w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
        x3 <- x2 - s * trait1[i] * cos(theta)
        x4 <- x1 - s * trait1[i] * cos(theta)
        x5 <- x2 + s * trait2[i] * cos(theta)
        x6 <- x1 + s * trait2[i] * cos(theta)
        
        y1 <- s * max(trait1) * sin(theta) + dy + (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
        y2 <- s * max(trait1) * sin(theta) + dy - (w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
        y3 <- y2 - s * trait1[i] * sin(theta)
        y4 <- y1 - s * trait1[i] * sin(theta)
        y5 <- y2 + s * trait2[i] * sin(theta)
        y6 <- y1 + s * trait2[i] * sin(theta)
        
        x1w <- s * max(trait1) * cos(theta) + dx - wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
        x2w <- s * max(trait1) * cos(theta) + dx + wider*(w/2) * cos(pi/2 - theta) - s * min(0,min(trait1)) * cos(theta)
        y1w <- s * max(trait1) * sin(theta) + dy + wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
        y2w <- s * max(trait1) * sin(theta) + dy - wider*(w/2) * sin(pi/2 - theta) - s * min(0,min(trait1)) * sin(theta)
        
        x01 <- x2w - s * max(trait1) * cos(theta)
        x02 <- x1w - s * max(trait1) * cos(theta)
        x03 <- x2w + s * max(trait1) * cos(theta)
        x04 <- x1w + s * max(trait1) * cos(theta)
        y01 <- y2w - s * max(trait1) * sin(theta)
        y02 <- y1w - s * max(trait1) * sin(theta)
        y03 <- y2w + s * max(trait1) * sin(theta)
        y04 <- y1w + s * max(trait1) * sin(theta)
      }
      ## trait 1: Quadrate
      polygon(c(x1, x2, x3, x4), c(y1, y2, y3, y4), 
              col = quad.cols$col[which(abs(trait1[i]-quad.cols$rql)==min(abs(trait1[i]-quad.cols$rql)))], border = "black",lwd=trait.border)
      ## trait 2: MA
      polygon(c(x1, x2, x5, x6), c(y1, y2, y5, y6), 
              col = MA.cols$col[which(abs(trait2[i]-MA.cols$MA)==min(abs(trait2[i]-MA.cols$MA)))], border = "black",lwd=trait.border)
    }
    hot_legend_image <- as.raster(matrix(rev(hot_color_func(100)), ncol=100))
    rasterImage(hot_legend_image, -h ,-7.5, -5,-10)
    polygon(col = NA, border = "black", lwd=0.5, x=c(-h,-h,-5,-5), y = c(-10,-7.5,-7.5,-10))
    text(x=mean(c(-h,-5)),y=-5.5,"Mechanical advantage (MA)",cex=0.8)
    text(x=-h+5, y = -11.5, cex = 0.55, "More force-modified")
    text(x=-10, y = -11.5, cex = 0.55, "Less force-modified")
    arrows(-22.5,-11.5,-h+17.5,-11.5, lwd=0.75, angle = 30, length=0.05, code = 3, col = "black")
    
    cold_legend_image <- as.raster(matrix(rev(cold_color_func(100)), ncol=100))
    rasterImage(cold_legend_image, h,-7.5,5,-10)
    polygon(col = NA, border = "black", lwd=0.5, x=c(h,h,5,5), y = c(-10,-7.5,-7.5,-10))
    text(x=mean(c(h,5)),y=-5.5,"Relative quadrate length (RQL)",cex=0.8)
    text(x=h-5, y = -11.5, cex = 0.55, "Larger gape")
    text(x=10, y = -11.5, cex = 0.55, "Smaller gape")
    arrows(20,-11.5,h-15,-11.5, lwd=0.75, angle = 30, length=0.05, code = 3, col = "black")
    
    diet.legend = readPNG("figs/phylogeny_diet_legend.png")
    dil.x = ncol(diet.legend)/2/sqrt(3000)
    dil.y = nrow(diet.legend)/2/sqrt(3000)
    rasterImage(diet.legend,
                xleft = -1*dil.x,
                xright = dil.x,
                ybottom = -25-dil.y,
                ytop=-25+dil.y
    )
    rect(xleft = -1*dil.x,
         xright = dil.x,
         ybottom = -25-dil.y,
         ytop=-25++1.225*dil.y+2.5,
         border = "white",
         lwd = 0.5
    )
    
    text("Main prey item in diet",
         x=0, y = -25+dil.y,
         cex = 0.8, col = "black"
    )
    
    bgcol = rgb(255,255,255, max = 255, alpha = 175)
    pie.adj = NULL
    nodelabels("1",node=viper.mrca,frame="circle",bg=bgcol,cex=0.5,adj=pie.adj)  ## 1 = viperids
    nodelabels("2",node=getMRCA(map3, c("Brachyorrhos_albus","Homalopsis_buccata")),frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 2 = homalopsids
    nodelabels("3",node=getMRCA(map3, c("Atractaspis_bibronii","Micrurus_fulvius")),frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 3 = elapoids
    nodelabels("4",node=getMRCA(map3, c("Langaha_madagascariensis","Duberria_lutrix")),frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 4 = psuedoxyrhophiids
    nodelabels("5",node=getMRCA(map3, c("Atractaspis_bibronii","Aparallactus_werneri")),frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 5 = atractaspididsa
    nodelabels("6",node=elapid.mrca,frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ##  6 = elapids
    nodelabels("7",node=getMRCA(map3, c("Coluber_constrictor","Carphophis_amoenus")),frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 7 = colubroids
    nodelabels("8",node=dipsadid.mrca,frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 8 = dipsadids
    nodelabels("9",node=natricid.mrca,frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 9 = natricids
    nodelabels(" ",node=colubrid.mrca,frame="circle",bg=bgcol,cex=0.5,adj=pie.adj) ## 10 = colubrids
    
    points(x=rep(0,5), y = c(seq(18,5,by=-3.25)), pch = 21, cex = 1.5)
    text(x=rep(0,5), y = c(seq(18,5,by=-3.25)), cex = 0.5, c(1:5))
    text(x=rep(0,5), y = c(seq(18,5,by=-3.25)), cex = 0.5, c("Viperidae","Homalopsidae","Elapoidea","Pseudoxyrhophiidae","Atractaspididiae"), pos=4)
    points(x=rep(23,5), y = c(seq(18,5,by=-3.25)), pch = 21, cex = 1.5)
    text(x=rep(23,5), y = c(seq(18,5,by=-3.25)), cex = 0.5, c(6:10))
    text(x=rep(23,5), y = c(seq(18,5,by=-3.25)), cex = 0.5, c("Elapidae","Colubroidea","Dipsadidae","Natricidae","Colubridae"), pos=4)
    
  }
  dev.off()
    
}