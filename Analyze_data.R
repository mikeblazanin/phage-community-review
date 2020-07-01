##Load libraries & data ----
library(ggplot2)

all_data <- read.csv("./Clean_data/cleaned_merged_data.csv",
                     stringsAsFactors = F)

#Okabe and Ito 2008 colorblind-safe qualitative color scale ----
my_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000")
scales::show_col(my_cols)

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
                             dens_pop_n = sum(is.numeric(dens_rep_avg)))


#Plots ----
ggplot(data = all_data_repsum[all_data_repsum$Study == "Mumford_Friman", ], 
       aes(x = Time_day, y = dens_rep_avg+1,
                            color = Competitor_community)) +
  geom_point(alpha = 0.5, 
             position = position_dodge(width = 1.5)) +
  facet_grid(Pop ~ PhagePresence, scales = "free_y") +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = my_cols)

ggplot(data = all_data[all_data$Study == "Alseth_etal", ],
       aes(x = Time_day, y = Density,
           color = Competitor_community)) +
  geom_point() +
  facet_grid(Pop~PhagePresence, scales = "free") +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = my_cols)
