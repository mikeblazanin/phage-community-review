##Load libraries & data ----
library(ggplot2)

all_data <- read.csv("./Clean_data/cleaned_merged_data.csv",
                     stringsAsFactors = F)

#Okabe and Ito 2008 colorblind-safe qualitative color scale ----
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

#Drop NA Density values
all_data <- all_data[!is.na(all_data$Density),]

#Recode phage presence
all_data$PhagePresence[all_data$PhagePresence == 1] <- "Phage present"
all_data$PhagePresence[all_data$PhagePresence == 0] <- "Phage absent"

#Reorder Pop levels
all_data$Pop <- factor(all_data$Pop,
                       levels = c("P_aeruginosa", "PT7",
                                  "PA14", "AB", "BC", "SA", "DMS3vir"))

#Reorder competitor community levels
all_data$Competitor_community <- 
  factor(all_data$Competitor_community,
         levels = c("None", "Sa", "Sm", "Sa+Sm",
                    "AB", "BC", "SA", "AB+BC+SA", "surfacemutant"))
                                        

#Get summarized values
all_data <- group_by(all_data, Study, Focal_strain, Competitor_community,
                     PhagePresence, Time_day, Rep_pop, Pop)
all_data_repsum <- summarize(all_data,
                          dens_rep_avg = mean(Density, na.rm = TRUE),
                          dens_rep_sd = sd(Density, na.rm = TRUE),
                          dens_rep_n = sum(sapply(Density, is.numeric)))

##TODO check that repsum is correct!

all_data_repsum <- group_by(all_data_repsum,
                            Study, Focal_strain, Competitor_community,
                            PhagePresence, Time_day, Pop)
all_data_popsum <- summarize(all_data_repsum,
                             dens_pop_avg = mean(dens_rep_avg),
                             dens_pop_sd = sd(dens_rep_avg, na.rm = TRUE),
                             dens_pop_n = sum(sapply(dens_rep_avg, is.numeric)))


#Plots ----
temp_repsum <- all_data_repsum[all_data_repsum$Study == "Mumford_Friman", ]
temp_popsum <- all_data_popsum[all_data_popsum$Study == "Mumford_Friman", ]
ggplot(data = temp_repsum, 
       aes(x = Time_day, y = dens_rep_avg+1, color = Competitor_community)) +
  geom_point(alpha = 0.2, position = position_dodge(width = 1.5)) +
  # geom_line(aes(group = interaction(Rep_pop, Competitor_community)),
  #           alpha = 0.2, position = position_dodge(width = 1.5)) +
  geom_point(data = temp_popsum, aes(y = dens_pop_avg+1)) +
  geom_line(data = temp_popsum, aes(y = dens_pop_avg+1)) +
  facet_grid(Focal_strain + Pop ~ PhagePresence, scales = "free_y") +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = my_cols)

temp_repsum <- all_data_repsum[all_data_repsum$Study == "Alseth_etal", ]
temp_popsum <- all_data_popsum[all_data_popsum$Study == "Alseth_etal", ]
ggplot(data = temp_repsum,
       aes(x = Time_day, y = dens_rep_avg+1,
           color = Competitor_community)) +
  geom_point(alpha = 0.2) +
  geom_point(data = temp_popsum, aes(y = dens_pop_avg+1)) +
  geom_line(data = temp_popsum, aes(y = dens_pop_avg+1), lwd = 1) +
  facet_grid(Pop~PhagePresence, scales = "free") +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = my_cols)
