library(FSA)
library(tidyr)
library(dplyr)
library(nnet)


#### clean and split data ####
aged = data.frame(read.csv(file="data/Catfish_age.csv"))
noaged = data.frame(read.csv(file="data/Catfish_noage.csv"))

aged$Weight = as.numeric(levels(aged$Weight))[aged$Weight]
aged$Length = as.numeric(levels(aged$Length))[aged$Length]
noaged$Age = as.numeric(noaged$Age)
aged = aged %>% drop_na(Length,Age,Weight)
aged= droplevels(aged[-which(aged$Fish %in% c('C3107','F222','F35','C3122','F3332','F3117','F3121', 'F37','F2621')),])

aged = aged %>%mutate(lcat10 = lencat(Length, w=10))
blue_age = aged %>% subset(Species == "Blue") %>% droplevels()
channel_age = aged %>% subset(Species == "Channel") %>% droplevels()
flathead_age = aged %>% subset(Species == "Flathead") %>% droplevels()


noaged$Weight = as.numeric(levels(noaged$Weight))[noaged$Weight]
noaged$Length = as.numeric(levels(noaged$Length))[noaged$Length]
noaged = noaged %>% drop_na(Length, Weight)
noaged = noaged %>%mutate(lcat10 = lencat(Length, w=10))
blue_noage = noaged %>% subset(Species == "Blue") %>% droplevels()
channel_noage = noaged %>% subset(Species == "Channel") %>% droplevels()
flathead_noage = noaged %>% subset(Species == "Flathead") %>% droplevels()

#### make traditional (raw prop. table) age-length keys  ####

alk_freq_blue =  xtabs(~lcat10+Age, data = blue_age)
alk_freq_chan =  xtabs(~lcat10+Age, data = channel_age)
alk_freq_flat =  xtabs(~lcat10+Age, data = flathead_age)

alk_T_blue = prop.table(alk_freq_blue, margin = 1)
alk_T_chan = prop.table(alk_freq_chan, margin = 1)
alk_T_flat = prop.table(alk_freq_flat, margin = 1)

#### make smoothed modeled (multinomial estimated) age-length keys  ####
blue_mlr = multinom(Age~lcat10, data = blue_age, maxit = 500)
channel_mlr = multinom(Age~lcat10, data = channel_age, maxit = 500)
flathead_mlr = multinom(Age~lcat10, data = flathead_age, maxit = 500)

blue_lens = seq(10, 1020, 10)# set based on min and max length for both age and noage
channel_lens = seq(40, 800, 10)
flathead_lens = seq(40, 4190, 10)

alk_sm_blue = predict(blue_mlr,data.frame(lcat10=blue_lens),type = 'probs')
row.names(alk_sm_blue) = blue_lens
alk_sm_chan = predict(channel_mlr,data.frame(lcat10=channel_lens),type = 'probs')
row.names(alk_sm_chan) = channel_lens
alk_sm_flat = predict(flathead_mlr,data.frame(lcat10=flathead_lens),type = 'probs')
row.names(alk_sm_flat) = flathead_lens

#### apply key of choice and remove junk ####

blue_sm_aged = alkIndivAge(alk_sm_blue, Age~Length, data=blue_noage)
blue_noage$Age = blue_sm_aged$Age

channel_sm_aged = alkIndivAge(alk_sm_chan, Age~Length, data=channel_noage)
channel_noage$Age = channel_sm_aged$Age

flathead_sm_aged = alkIndivAge(alk_sm_flat, Age~Length, data=flathead_noage)
flathead_noage$Age = flathead_sm_aged$Age

remove(list = c("alk_sm_blue", "alk_sm_chan", "alk_sm_flat",
                "blue_mlr", "blue_sm_aged", "channel_mlr",
                "channel_sm_aged", "flathead_mlr", "flathead_sm_aged",
                "alk_freq_blue","alk_T_blue", "blue_lens",
                "alk_freq_chan","alk_T_chan", "channel_lens",
                "alk_freq_flat","alk_T_flat", "flathead_lens"))

### Combine, clean  and make jags datasets ####

blue_combined = rbind(blue_age,blue_noage)
channel_combined = rbind(channel_age, channel_noage)
flathead_combined = rbind(flathead_age, flathead_noage)

remove(list = c("blue_age","blue_noage","channel_age", "channel_noage", "flathead_age", "flathead_noage"))

blue_growth =  blue_combined %>% select('Age', 'Weight','Length')
blue_growth$lWeight = log10(blue_growth$Weight)
blue_growth$lLength = log10(blue_growth$Length)
channel_growth = channel_combined %>% select('Age', 'Weight','Length')
channel_growth$lWeight = log10(channel_growth$Weight)
channel_growth$lLength = log10(channel_growth$Length)
flathead_growth = flathead_combined %>% select('Age', 'Weight','Length')
flathead_growth$lWeight = log10(flathead_growth$Weight)
flathead_growth$lLength = log10(flathead_growth$Length)
# remove unknown gears and set first age to age with highest catch
blue_mort = blue_combined %>% group_by(Gear,Age) %>% count()
blue_mort = droplevels(blue_mort[-which(blue_mort$Gear == "Unknown"),])
blue_mort = droplevels(blue_mort[-which(blue_mort$Age < 3),])
blue_mort = droplevels(blue_mort[-which(blue_mort$Age < 5 & blue_mort$Gear == "Trotline"),])
blue_mort = blue_mort %>% spread(Gear, n)

channel_mort = channel_combined %>% group_by(Gear,Age) %>% count()
channel_mort = droplevels(channel_mort[-which(channel_mort$Age < 6 & channel_mort$Gear == "e60"),])
channel_mort = droplevels(channel_mort[-which(channel_mort$Age < 6 & channel_mort$Gear == "Hoopnet"),])
channel_mort = droplevels(channel_mort[-which(channel_mort$Age < 6 & channel_mort$Gear == "Trotline"),])
channel_mort = channel_mort %>% spread(Gear, n)

flathead_mort = flathead_combined %>% group_by(Gear,Age) %>% count()
flathead_mort = droplevels(flathead_mort[-which(flathead_mort$Gear == "Unknown"),])
flathead_mort = droplevels(flathead_mort[-which(flathead_mort$Gear == "Trotline"),])
flathead_mort = droplevels(flathead_mort[-which(flathead_mort$Age < 4 & flathead_mort$Gear == "e15"),])
flathead_mort = droplevels(flathead_mort[-which(flathead_mort$Age < 2 & flathead_mort$Gear == "e60"),])
flathead_mort = droplevels(flathead_mort[-which(flathead_mort$Age < 5 & flathead_mort$Gear == "Hoopnet"),])
flathead_mort = flathead_mort %>% spread(Gear, n)

# Detach packages to avoid masking in other script(s)
detach("package:tidyr", unload = TRUE)

