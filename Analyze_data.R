##Load libraries & data ----
library(ggplot2)
library(dplyr)

all_data <- read.csv("./Clean_data/cleaned_merged_data.csv",
                     stringsAsFactors = F)

#Okabe and Ito 2008 colorblind-safe qualitative color scale ----
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

#Drop NA Density values
all_data <- all_data[!is.na(all_data$Density),]

#Recode phage presence
all_data$Phage_presence[all_data$Phage_presence == 1] <- "Phage present"
all_data$Phage_presence[all_data$Phage_presence == 0] <- "Phage absent"

#Reorder Pop levels
all_data$Pop <- factor(all_data$Pop,
                       levels = c("P_aeruginosa", "PT7",
                                  "PA14", "AB", "BC", "SA", "DMS3vir",
                                  "Klebsiella", "Pseudomonas", "Staphylococcus",
                                  "Kleb_phage"))

#Reorder competitor community levels
all_data$Bact_community <- 
  factor(all_data$Bact_community,
         levels = c("K", "P", "S", "K+P+S", 
                    "Pa", "Pa+Sa", "Pa+Sm", "Pa+Sa+Sm",
                    "PA", "PA+AB", "PA+BC", "PA+SA", "PA+AB+BC+SA", 
                    "PA+surfacemutant"))
all_data$Focal_strain <- 
  factor(all_data$Focal_strain,
         levels = c("Klebsiella", "PAO1", "lasR", "PA14"))


#Get summarized values
all_data <- group_by(all_data, Study, Focal_strain, Bact_community,
                     Phage_presence, Time_day, Rep_pop, Pop)
all_data_repsum <- summarize(all_data,
                          dens_rep_avg = mean(Density, na.rm = TRUE),
                          dens_rep_sd = sd(Density, na.rm = TRUE),
                          dens_rep_n = sum(sapply(Density, is.numeric)))

##TODO check that repsum is correct!

all_data_repsum <- group_by(all_data_repsum,
                            Study, Focal_strain, Bact_community,
                            Phage_presence, Time_day, Pop)
all_data_popsum <- summarize(all_data_repsum,
                             dens_pop_avg = mean(dens_rep_avg),
                             dens_pop_sd = sd(dens_rep_avg, na.rm = TRUE),
                             dens_pop_n = sum(sapply(dens_rep_avg, is.numeric)))

#Plots ----
dir.create("./Plots", showWarnings = FALSE)

#Mumford
temp_repsum <- all_data_repsum[all_data_repsum$Study == "Mumford_Friman", ]
temp_popsum <- all_data_popsum[all_data_popsum$Study == "Mumford_Friman", ]
tiff("./Plots/Mumford_Friman.tiff", width = 10, height = 6,
     units = "in", res = 300)
ggplot(data = temp_repsum, 
       aes(x = Time_day, y = dens_rep_avg+1, color = Pop)) +
  geom_point(alpha = 0.5, 
             #position = position_dodge(width = 1.5)
             ) +
  # geom_line(aes(group = interaction(Rep_pop, Competitor_community)),
  #           alpha = 0.2, position = position_dodge(width = 1.5)) +
  geom_point(data = temp_popsum, aes(y = dens_pop_avg+1),
             #position = position_dodge(width = 1.5), 
             size = 3) +
  geom_line(data = temp_popsum, aes(y = dens_pop_avg+1)) +
  geom_hline(yintercept = 1, lty = 2) +
  facet_grid(Phage_presence ~ Focal_strain + Bact_community) +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = my_cols) +
  labs(x = "Time (days)", y = "Density + 1", color = "Population", 
       title = "Focal Host Genotype", subtitle = "Community Treatment") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text.x = element_text(size = 11.5),
        strip.text.y = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 13),
        plot.subtitle = element_text(size = 13))
dev.off()

#Alseth
temp_repsum <- all_data_repsum[all_data_repsum$Study == "Alseth_etal", ]
temp_popsum <- all_data_popsum[all_data_popsum$Study == "Alseth_etal", ]
tiff("./Plots/Alseth_etal.tiff", width = 10, height = 5,
     units = "in", res = 300)
ggplot(data = temp_repsum,
       aes(x = Time_day, y = dens_rep_avg+1,
           color = Pop)) +
  geom_point(alpha = 0.5) +
  geom_point(data = temp_popsum, aes(y = dens_pop_avg+1), size = 3) +
  geom_line(data = temp_popsum, aes(y = dens_pop_avg+1), lwd = 1) +
  geom_hline(yintercept = 1, lty = 2) +
  facet_grid(Phage_presence~Bact_community) +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = my_cols) +
  labs(x = "Time (days)", y = "Density + 1",
       color = "Population", subtitle = "Community Treatment") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        plot.subtitle = element_text(size = 13))
dev.off()

#Johnke
temp_repsum <- all_data_repsum[all_data_repsum$Study == "Johnke_etal", ]
temp_popsum <- all_data_popsum[all_data_popsum$Study == "Johnke_etal", ]
tiff("./Plots/Johnke_etal.tiff", width = 10, height = 5,
     units = "in", res = 300)
ggplot(data = temp_repsum,
       aes(x = Time_day, y = dens_rep_avg+1,
           color = Pop)) +
  geom_point(alpha = 0.5) +
  geom_point(data = temp_popsum, aes(y = dens_pop_avg+1), size = 3) +
  geom_line(data = temp_popsum, aes(y = dens_pop_avg+1), lwd = 1) +
  geom_hline(yintercept = 1, lty = 2) +
  facet_grid(Phage_presence~Bact_community) +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = my_cols) +
  labs(x = "Time (days)", y = "Density + 1",
       color = "Population", subtitle = "Community Treatment") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        plot.subtitle = element_text(size = 13))
dev.off()

##Stats ----

model_mum_bact <- lm(log10(dens_rep_avg+1) ~ Bact_community*Focal_strain*Phage_presence,
                 all_data_repsum[all_data_repsum$Time_day > 0 &
                                   all_data_repsum$Study == "Mumford_Friman" &
                                   all_data_repsum$Pop == "P_aeruginosa", ])
summary(model_mum_bact)

model_mum_phg <- lm(log10(dens_rep_avg) ~ Bact_community*Time_day*Focal_strain,
                all_data_repsum[all_data_repsum$Time_day > 0 &
                                all_data_repsum$Study == "Mumford_Friman" &
                                all_data_repsum$Pop == "PT7", ])
summary(model_mum_phg)
##Only sig effect is that phage declines over time

model_alseth_bact <- lm(log10(dens_rep_avg) ~ Bact_community*Time_day,
                    all_data_repsum[all_data_repsum$Time_day > 0 &
                                      all_data_repsum$Study == "Alseth_etal" &
                                      all_data_repsum$Pop == "PA14", ])
summary(model_alseth_bact)
##PA dens is suppressed by           


model_alseth_phg <- lm(log10(dens_rep_avg) ~ Bact_community*Time_day,
                       all_data_repsum[all_data_repsum$Time_day > 0 &
                                         all_data_repsum$Study == "Alseth_etal" &
                                         all_data_repsum$Pop == "DMS3vir", ])
summary(model_alseth_phg)
##Phages decline faster over time with BC


model_jon_bact <- lm(log10(dens_rep_avg+1) ~ Bact_community*Phage_presence*Time_day,
             all_data_repsum[all_data_repsum$Time_day > 0 &
                              all_data_repsum$Study == "Johnke_etal" &
                              all_data_repsum$Pop == "Klebsiella", ])
summary(model_jon_bact)

model_jon_phg <- lm(log10(dens_rep_avg) ~ Bact_community*Time_day,
                 all_data_repsum[all_data_repsum$Time_day > 0 &
                                   all_data_repsum$Study == "Johnke_etal" &
                                   all_data_repsum$Pop == "Kleb_phage", ])
summary(model_jon_phg)
             
             
             
             
             
             