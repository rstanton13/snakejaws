library(tidyverse)
library(squamatabase)
library(ggbeeswarm)

data(diet)
d=diet[!is.na(diet$predator_mass) & !is.na(diet$prey_mass),]

d = subset(d, select=c(predator, predator_taxon, prey, 
                       prey_taxon, prey_mass, prey_age, predator_mass, 
                       prey_count, event_outcome))
d$prey_age=tolower(d$prey_age)
d = d[!grepl("Boidae|Cylindrophiidae",d$predator_taxon),]
d=d[d$event_outcome=="predation_interrupted_by_observer"|d$event_outcome=="prey_eaten",]

d$prey_taxon[grepl("Reptilia",d$prey_taxon) & grepl("egg",d$prey_age)]="rept_egg"
d$prey_taxon[grepl("Serpentes|Amphisbaenidae|Ophiodes",d$prey_taxon)]="long_squam"
d$prey_taxon[grepl("Squamata",d$prey_taxon)]="lizard"
d$prey_taxon[grepl("Sirenidae|Amphiumidae|Anguilliformes|Synrbanchiformes|Gymnotiformes|Petromyzontiformes",d$prey_taxon)]="eel"
d$prey_taxon[grepl("Amphibia",d$prey_taxon) & grepl("neonate|larvae|tadpole|juvenile|subadult",d$prey_age)]="amphib_larvae"
d$prey_taxon[grepl("Amphibia",d$prey_taxon) & grepl("egg",d$prey_age)]="amphib_egg"
d$prey_taxon[grepl("Anura",d$prey_taxon)]="frog"
d$prey_taxon[grepl("Caudata",d$prey_taxon)]="salamander"
d$prey_taxon[grepl("Gymnophiona",d$prey_taxon)]="caecilian"
d$prey_taxon[grepl("Mammalia",d$prey_taxon)]="mammal"
d$prey_taxon[grepl("Actinopterygii",d$prey_taxon) & grepl("egg",d$prey_age)]="fish_egg"
d$prey_taxon[grepl("Actinopterygii",d$prey_taxon)]="fish"
d$prey_taxon[grepl("Aves",d$prey_taxon) & grepl("egg",d$prey_age)]="bird_egg"
d$prey_taxon[grepl("Aves",d$prey_taxon)]="bird"
d$prey_taxon[grepl("Gastropoda",d$prey_taxon)]="gastropod"
d$prey_taxon[grepl("Haplotaxida",d$prey_taxon)]="worm"
d$prey_taxon[grepl("Testudines",d$prey_taxon)]="turtle"
d$prey_taxon[grepl("Decapoda",d$prey_taxon)]="crust"
d$prey_taxon[grepl("Arthropoda",d$prey_taxon)]="bug"



d$rpm = (d$prey_mass/d$predator_mass)/d$prey_count
d$log_rpm = log(d$rpm)

d = d[grepl("long_squam|lizard|bird|mammal|fish|frog|worm|bug|caecilian|gastropod",d$prey_taxon),]


d$prey_taxon = factor(d$prey_taxon, levels = c("bird","mammal","frog","salamander","lizard","caecilian","long_squam","fish","gastropod","worm","bug"),
                      labels = c("Birds","Mammals","Frogs","Salamanders","Lizards","Caecilians","Elongate
squamates", "Non-elongate
fishes", "Gastropods","Earthworms","Arthropods"))

rpmplot = ggplot(data=d, aes(x=prey_taxon, y = rpm, color = prey_taxon)) +
  geom_boxplot(outlier.shape = NA) +
  labs(
    x = "Type of prey",
    y = "RPM"
  ) + 
  geom_quasirandom(alpha=0.7) +
  scale_color_manual(values=c("gray50","gray50","gray50","forestgreen","mediumpurple2","mediumpurple2",
                              "royalblue2","saddlebrown","lightpink2","goldenrod2"))+
  theme_light() + theme(axis.text.x = element_text(angle = 30, hjust = 1), 
                        plot.title = element_text(hjust = 0.5), 
                        axis.title = element_text(size=15), 
                        title=element_text(size=15),
                        panel.border = element_blank(),
                        axis.line = element_line(color = "gray60"),
                        axis.line.y.right = element_blank(),
                        axis.line.x.top = element_blank(),
                        legend.spacing = unit(20, "lines"),
                        legend.position = "none"
  )

ggsave(plot = rpmplot,
       filename = "figs/RPM_plot.png",
       bg="white",
       units="in",
       height=6,
       width=9,
       dpi = 1200
)

