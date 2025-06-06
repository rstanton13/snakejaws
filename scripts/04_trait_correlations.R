## set wd and load packages
{
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
  `%nin%` <- Negate(`%in%`)
  options(scipen=999)
  library(slouch)
  library(nlme)
  library(extrafont)
}

## load data
{
  dat <- read.csv("hOUwie/skull_dat.csv",row.names=1)
  phy <- read.tree("trees/Title2024_macrostomata_names_replaced.tre")
  phy <- keep.tip(phy, dat$species)
}

## slouch analysis
{
  par(mfrow=c(1,1))
  
  null.fit <- slouch.fit(phy,
                         species=dat$species,
                         response=dat$MA
  )
  
  test.fit <- slouch.fit(phy,
                         species=dat$species,
                         response=dat$MA,
                         random.cov = dat$RQL
  )
  
  print(c(null=null.fit$modfit$AICc, test = test.fit$modfit$AICc))
  ## test model is better by a lot by AICc
}

## PGLS analysis
{
  spp <- dat$species
  corBM <- corBrownian(phy=phy, form = ~spp)
  pgls.test <- gls(MA ~ RQL, correlation = corBM, data = dat, method = "ML")
  pgls.null <- gls(MA ~ 1, correlation = corBM, data = dat, method = "ML")
  
  print(c(null=AIC(pgls.null), test = AIC(pgls.test)))
  ## test model is better by a lot by AIC
}

## plot for optima
ror = "ratio" ## "ratio" or "residual"
{
  rql_avg_pars = readRDS(paste("hOUwie/",ror,"_RQL_results.RDS",sep=""))
  ma_avg_pars  = readRDS(paste("hOUwie/",ror,"_MA_results.RDS",sep=""))
  
  optd = data.frame(r = unique(rql_avg_pars$theta), m = unique(ma_avg_pars$theta))
  optd$diet = unique(rql_avg_pars$tip_state)
  lmf = summary(lm(optd$m ~ optd$r))
  
  p=ggplot(optd,aes(x=m,y=r)) + geom_smooth(method="lm",color="black") + 
    geom_point(data = optd, aes(x=m, y = r, fill=diet),pch=23,size=4, alpha = 0.7) +
    scale_fill_manual(
      breaks = c("other verts", "lizards", "snakes", "fish", "eels", "crusts", "mollusks", "worms", "bugs"),
      guide = "none", 
      values = c("gray50", "forestgreen", "mediumpurple2", "royalblue2", "turquoise2", "tomato2", "saddlebrown", "lightpink2", "goldenrod2"),
      name = NA, 
      labels = rep(NA, 9)
    ) + theme_light() + 
    scale_x_continuous(name="MA optima") +
    scale_y_continuous(name="RQL optima") + 
    ggtitle(paste("p =",lmf$coefficients[2,4]))
  print(p)
}

## plot for PGLS and SLOUCH
{
  # Create a new data frame for the custom legend lines
  legend_data <- data.frame(
    x = c(min(dat$RQL), max(dat$RQL)),  # Range for RQL
    y = c(pgls.test$coefficients[1] + pgls.test$coefficients[2] * min(dat$RQL), 
          pgls.test$coefficients[1] + pgls.test$coefficients[2] * max(dat$RQL)),  # Calculate corresponding y values
    line_type = "PGLS"
  )
  
  # Add another data frame for the 'SLOUCH' line
  legend_data_slo <- rbind(legend_data, data.frame(
    x = c(min(dat$RQL), max(dat$RQL)),  # Range for RQL
    y = c(test.fit$beta_primary$coefficients_bias_corr[1,1] + test.fit$beta_primary$coefficients_bias_corr[2,1] * min(dat$RQL), 
          test.fit$beta_primary$coefficients_bias_corr[1,1] + test.fit$beta_primary$coefficients_bias_corr[2,1] * max(dat$RQL)),  # Corresponding y values for SLOUCH
    line_type = "SLOUCH"
  ))
  
  #legend_data_slo <- rbind(legend_data_slo, data.frame(
  #  x = c(min(dat$RQL), max(dat$RQL)),  # Range for RQL
  #  y = c(lmf$coefficients[1,1] + lmf$coefficients[2,1] * min(dat$RQL), 
  #        lmf$coefficients[1,1] + lmf$coefficients[2,1] * max(dat$RQL)),  # Corresponding y values for optima
  #  line_type = "Optima"
  #))
  
  # Plot the data with custom lines and the legend
  p2 <- ggplot(dat, aes(x = RQL, y = MA, color = diet)) + 
    geom_point(alpha = 0.7) +
    
    # Customizing color scale for 'diet' without showing it in the legend
    scale_color_manual(
      breaks = c("other verts", "lizards", "snakes", "fish", "eels", "crusts", "mollusks", "worms", "bugs"),
      guide = "none", 
      values = c("gray50", "forestgreen", "mediumpurple2", "royalblue2", "turquoise2", "tomato2", "saddlebrown", "lightpink2", "goldenrod2"),
      name = NA, 
      labels = rep(NA, 9)
    ) +
    
    geom_point(data = optd, aes(x = r, y = m, fill=diet),pch=23, color = "black", size=4, alpha = 0.7) +
    
    scale_fill_manual(
      breaks = c("other verts", "lizards", "snakes", "fish", "eels", "crusts", "mollusks", "worms", "bugs"),
      guide = "none", 
      values = c("gray50", "forestgreen", "mediumpurple2", "royalblue2", "turquoise2", "tomato2", "saddlebrown", "lightpink2", "goldenrod2"),
      name = NA, 
      labels = rep(NA, 9)
    ) +
    # Add the 'PGLS' and 'SLOUCH' lines
    geom_line(data = legend_data_slo, aes(x = x, y = y, linetype = line_type), size = 1, color = "black") +
    
    # Labels for the x and y axis
    scale_x_continuous(name = "log(RQL)") + 
    scale_y_continuous(name = "log(MA)") + 
    
    # Title for the plot (empty here)
    ggtitle("") + 
    
    # Using minimal theme
    theme_minimal() + 
    
    # Center the title
    theme(plot.title = element_text(hjust = 0.5)) + 
    
    # Adding a custom legend for the lines (SLOUCH and PGLS)
    guides(
      linetype = guide_legend(title = NULL, 
                              labels = c("PGLS", "SLOUCH", "Optima"), 
                              override.aes = list(color = "black", size = 1))
    ) +
    theme(
      axis.line = element_blank(),  # Remove all axis lines
      axis.line.x = element_line(color = "gray80", size = 0.250),  # Dark line at the bottom
      axis.line.y = element_line(color = "gray80", size = 0.250),  # Dark line on the left
      axis.ticks = element_line(color = "gray80"),  # Dark ticks
      axis.ticks.length = unit(0.2, "cm"),  # Set tick length
      axis.text = element_text(size = 12),  # Adjust axis text size
      axis.title = element_text(size = 14)  # Adjust axis title size
    )
  
  # Display the plot
  print(p2)
  
  ggsave(plot=p2, 
         filename = "figs/inverse_correlation.png",
         bg="white",
         units="in",
         height=6,
         width=7,
         dpi=1200
  )
  
}

## plot for PGLS and SLOUCH for the poster
{
  # Create a new data frame for the custom legend lines
  legend_data <- data.frame(
    x = c(min(dat$RQL), max(dat$RQL)),  # Range for RQL
    y = c(pgls.test$coefficients[1] + pgls.test$coefficients[2] * min(dat$RQL), 
          pgls.test$coefficients[1] + pgls.test$coefficients[2] * max(dat$RQL)),  # Calculate corresponding y values
    line_type = "PGLS"
  )
  
  # Add another data frame for the 'SLOUCH' line
  legend_data_slo <- rbind(legend_data, data.frame(
    x = c(min(dat$RQL), max(dat$RQL)),  # Range for RQL
    y = c(test.fit$beta_primary$coefficients_bias_corr[1,1] + test.fit$beta_primary$coefficients_bias_corr[2,1] * min(dat$RQL), 
          test.fit$beta_primary$coefficients_bias_corr[1,1] + test.fit$beta_primary$coefficients_bias_corr[2,1] * max(dat$RQL)),  # Corresponding y values for SLOUCH
    line_type = "SLOUCH"
  ))
  
  #legend_data_slo <- rbind(legend_data_slo, data.frame(
  #  x = c(min(dat$RQL), max(dat$RQL)),  # Range for RQL
  #  y = c(lmf$coefficients[1,1] + lmf$coefficients[2,1] * min(dat$RQL), 
  #        lmf$coefficients[1,1] + lmf$coefficients[2,1] * max(dat$RQL)),  # Corresponding y values for optima
  #  line_type = "Optima"
  #))
  
  # Plot the data with custom lines and the legend
  p2 <- ggplot(dat, aes(x = RQL, y = MA, color = diet)) + 
    geom_point(alpha = 0.7) +
    
    # Customizing color scale for 'diet' without showing it in the legend
    scale_color_manual(
      breaks = c("other verts", "lizards", "snakes", "fish", "eels", "crusts", "mollusks", "worms", "bugs"),
      guide = "none", 
      values = c("gray75", "forestgreen", "mediumpurple2", "royalblue2", "turquoise2", "tomato2", "saddlebrown", "lightpink2", "goldenrod2"),
      name = NA, 
      labels = rep(NA, 9)
    ) +
    
    geom_point(data = optd, aes(x = r, y = m, fill=diet),pch=23, color = "black", size=4, alpha = 0.7) +
    
    scale_fill_manual(
      breaks = c("other verts", "lizards", "snakes", "fish", "eels", "crusts", "mollusks", "worms", "bugs"),
      guide = "none", 
      values = c("gray75", "forestgreen", "mediumpurple2", "royalblue2", "turquoise2", "tomato2", "saddlebrown", "lightpink2", "goldenrod2"),
      name = NA, 
      labels = rep(NA, 9)
    ) +
    # Add the 'PGLS' and 'SLOUCH' lines
    geom_line(data = legend_data_slo, aes(x = x, y = y, linetype = line_type), size = 1, color = "white") +
    
    # Labels for the x and y axis
    scale_x_continuous(name = "log(RQL)") + 
    scale_y_continuous(name = "log(MA)") + 
    
    # Title for the plot (empty here)
    ggtitle("") + 
    
    # Using minimal theme
    theme(
      panel.background = element_rect(fill='transparent'), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      legend.background = element_rect(fill='transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent'), #transparent legend panel
      plot.title = element_text(hjust = 0.5)
    ) + 
    
    # Adding a custom legend for the lines (SLOUCH and PGLS)
    guides(
      linetype = guide_legend(title = NULL, 
                              labels = c("PGLS", "SLOUCH", "Optima"), 
                              override.aes = list(color = "white", size = 1))
    ) +
    theme(
      axis.line = element_blank(),  # Remove all axis lines
      axis.line.x = element_line(color = "white", size = 0.250),  # Dark line at the bottom
      axis.line.y = element_line(color = "white", size = 0.250),  # Dark line on the left
      axis.ticks = element_line(color = "white"),  # Dark ticks
      axis.ticks.length = unit(0.2, "cm"),  # Set tick length
      axis.text = element_text(size = 12),  # Adjust axis text size
      axis.title = element_text(size = 14)  # Adjust axis title size
    )
  
  # Display the plot
  print(p2)
  
  ggsave(plot=p2, 
         filename = "figs/poster_inverse_correlation.png",
         bg="transparent",
         units="in",
         height=6,
         width=9,
         dpi=600
  )
  
}

## pgls
t_stat = pgls.test$coefficients[2]/(pgls.test$varBeta[2,2]/593)
df <- 592
p_value <- 2 * (1 - pt(abs(t_stat), df))
p_value

## slouch
s2 = (sqrt(593)*(test.fit$beta_primary$coefficients[2,2])^2)
t_stat = test.fit$beta_primary$coefficients[2,1]/(s2/593)
df <- 592
p_value <- 2 * (1 - pt(abs(t_stat), df))
p_value
