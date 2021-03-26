library(tidyverse)
library(extrafont)
library(ggpubr)
#if(!exists('flathead_mort')) source('data/data_prep.R')  
#devtools::install_github("kassambara/ggpubr")
if(!exists('flathead_mort')) source('data/data_prep.R')
# clean data and calculate frequency
blue_combined$lcat20 = lencat(blue_combined$Length, w=20)
blue_length =  count(blue_combined,lcat20, Gear)

total = blue_length %>% 
  group_by(Gear) %>% 
  summarise(total = sum(n))
blue_length = merge(blue_length, total)
blue_length$freq = (blue_length$n/blue_length$total)*100
blue_length$species = 'Blue catfish'


channel_combined$lcat20 = lencat(channel_combined$Length, w=20)
Channel_length =  count(channel_combined,lcat20, Gear)

total = Channel_length %>% 
  group_by(Gear) %>% 
  summarise(total = sum(n))
Channel_length = merge(Channel_length, total)
Channel_length$freq = (Channel_length$n/Channel_length$total)*100
Channel_length$species = 'Channel catfish'


flathead_combined$lcat20 = lencat(flathead_combined$Length, w=20)
Flathead_length =  count(flathead_combined,lcat20, Gear)

total = Flathead_length %>% 
  group_by(Gear) %>% 
  summarise(total = sum(n))
Flathead_length = merge(Flathead_length, total)
Flathead_length$freq = (Flathead_length$n/Flathead_length$total)*100
Flathead_length$species = 'Flathead catfish'

rm(list = ls()[!ls() %in% c('Flathead_length', 'Channel_length', 'blue_length' )])

all_length = rbind(blue_length, Channel_length, Flathead_length)
blue_length = droplevels(blue_length[-which(blue_length$Gear %in% c('Unknown')),])
Flathead_length = droplevels(Flathead_length[-which(Flathead_length$Gear %in% c('Unknown')),])


#Blue plot
blue_plot = ggplot(blue_length, aes(lcat20, freq, fill = Gear)) +
  geom_bar(stat="identity", position = position_dodge(preserve = "single"),colour="black",size=0.25) +
  scale_y_continuous(name = "",limits=c(0,25),expand = c(0, 0)) +
  scale_x_continuous(name ="", breaks = c(0,100,200,
                                          300,400,500,
                                          600,700,800,
                                          900,1000,1100,1200),limits = c(0,1250), expand = c(0, 0)) +
  theme(axis.title.y = element_text(family="Arial")) + 
  scale_fill_grey(start = 1, end = 0) + theme_bw() +
  theme(panel.grid = element_line(colour = "transparent")) + 
  theme(legend.position = "none") +
  annotate("text", x=1000, y=20, label= "Blue Catfish")

#Channel plot
channel_plot = ggplot(Channel_length, aes(lcat20, freq, fill = Gear)) +
  geom_bar(stat="identity", position = position_dodge(preserve = "single"),colour="black",size=0.25) +
  scale_y_continuous(name = "",limits=c(0,40),expand = c(0, 0)) +
  scale_x_continuous(name ="", breaks = c(0,100,200,
                                                           300,400,500,
                                                           600,700,800,
                                                           900,1000, 1100, 1200),limits = c(0,1250), expand = c(0, 0)) +
  theme(axis.title.y = element_text(family="Arial")) + 
  scale_fill_grey(start = 1, end = 0) + theme_bw() +
  theme(panel.grid = element_line(colour = "transparent")) + 
  theme(legend.position = "none") +
  annotate("text", x=1000, y=33, label= "Channel Catfish")

#Flathead plot
flathead_plot = ggplot(Flathead_length, aes(lcat20, freq, fill = Gear)) +
  geom_bar(stat="identity", position = position_dodge(preserve = "single"),colour="black",size=0.25) +
  scale_y_continuous(name = "",limits=c(0,15),expand = c(0, 0)) +
  scale_x_continuous(name ="Total Length (mm)", breaks = c(0,100,200,
                                                           300,400,500,
                                                           600,700,800,
                                                           900,1000,1100,1200),limits = c(0,1250), expand = c(0, 0)) +
  theme(axis.title.y = element_text(family="Arial")) + 
  scale_fill_grey(start = 1, end = 0) + theme_bw() +
  theme(panel.grid = element_line(colour = "transparent")) + 
  theme(legend.position = "none") +
  annotate("text", x=1000, y=12, label= "Flathead Catfish")

#arrange and export plots

fig = ggarrange(blue_plot, channel_plot, flathead_plot,
          ncol = 1, nrow = 3)
png("results\\Figure 1.png",width = 6.95, height = 4.89,units = 'in', res = 1080,bg = "white")
annotate_figure(fig,left = text_grob("Relative frequency (percent of gear total catch)",family = "Arial", color = "Black", rot = 90))
dev.off()
