##Load libraries & data ----
library(ggplot2)
library(dplyr)
library(rjags)

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
                          dens_rep_avg = 10**mean(log10(Density+1), na.rm = TRUE),
                          dens_rep_sd = sd(Density, na.rm = TRUE),
                          dens_rep_n = sum(sapply(Density, is.numeric)))

all_data_repsum <- group_by(all_data_repsum,
                            Study, Focal_strain, Bact_community,
                            Phage_presence, Time_day, Pop)
all_data_popsum <- summarize(all_data_repsum,
                             dens_pop_avg = 10**mean(log10(dens_rep_avg+1)),
                             dens_pop_sd = sd(dens_rep_avg, na.rm = TRUE),
                             dens_pop_n = sum(sapply(dens_rep_avg, is.numeric)))

#Add log10-transformed col
all_data_repsum$dens_rep_avg_log10 <- log10(all_data_repsum$dens_rep_avg+1)

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
  scale_color_manual(values = my_cols[c(8, 5)]) +
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

#Simplified plot of bact only for talks
temp_repsum <- all_data_repsum[all_data_repsum$Study == "Mumford_Friman" &
                                 all_data_repsum$Focal_strain == "PAO1" &
                                 all_data_repsum$Pop == "P_aeruginosa" &
                                 all_data_repsum$Phage_presence == "Phage absent", ]
temp_popsum <- all_data_popsum[all_data_popsum$Study == "Mumford_Friman"  &
                                 all_data_popsum$Focal_strain == "PAO1" &
                                 all_data_popsum$Pop == "P_aeruginosa" &
                                 all_data_popsum$Phage_presence == "Phage absent", ]
tiff("./Plots/Mumford_Friman_simplified.tiff", width = 5, height = 4,
     units = "in", res = 300)
ggplot(data = temp_repsum, 
       aes(x = Bact_community, y = dens_rep_avg+1)) +
  geom_point(alpha = 0.5) +
  geom_point(data = temp_popsum, aes(y = dens_pop_avg+1),
             size = 4) +
  geom_hline(yintercept = 1, lty = 2) +
  facet_grid(Phage_presence ~ .) +
  scale_y_continuous(trans = "log10") +
  labs(x = "Community Treatment", y = "Bacterial Density + 1") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 16),
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
  scale_color_manual(values = my_cols[c(8, 1:3, 5)]) +
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
  scale_color_manual(values = my_cols[c(8, 1, 3, 5)]) +
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

#Simplified plot for talks
temp_repsum <- all_data_repsum[all_data_repsum$Study == "Johnke_etal" &
                                 all_data_repsum$Pop %in% c("Klebsiella", "Kleb_phage"), ]
temp_popsum <- all_data_popsum[all_data_popsum$Study == "Johnke_etal" &
                                 all_data_popsum$Pop %in% c("Klebsiella", "Kleb_phage"), ]
tiff("./Plots/Johnke_etal_simplified.tiff", width = 6, height = 4,
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
  scale_color_manual(values = my_cols[c(8, 1, 3, 5)]) +
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

#Mumford ----
mum_data_bact <- all_data_repsum[all_data_repsum$Time_day > 0 &
                              all_data_repsum$Study == "Mumford_Friman" &
                              all_data_repsum$Pop == "P_aeruginosa", ]
mum_data_bact$geno_com_phg <- paste(mum_data_bact$Focal_strain,
                                    mum_data_bact$Bact_community,
                                    mum_data_bact$Phage_presence,
                                    sep = "_")

#Note: Mumford did all their stats so I don't have to re-do them
model_mum_bact <- lm(log10(dens_rep_avg+1) ~ 
                       Bact_community*Focal_strain*Phage_presence,
                 data = mum_data_bact)
summary(model_mum_bact)

#Bayesian test of alternate ecological hypotheses
mum_bact_model <- "
  model{
    for(i in 1:length(dens_rep_avg_log10)){
      dens_rep_avg_log10[i] ~ dnorm(mu[i], tau)
      mu[i] <- geno_com_phg_int[geno_com_phg[i]]
    }
    tau <- 1 / (sig * sig)
    
    sig ~ dunif(0, 100)
    for(i in 1:n_geno_com_phg){
      geno_com_phg_int[i] ~ dnorm(0, 1.0E-3)
    }
  }"

#build the model (this includes, by default, 1000 adaptation steps)
mum_bact_jgsmdl <- 
  jags.model(file = textConnection(mum_bact_model),
             data = list(dens_rep_avg_log10 = mum_data_bact$dens_rep_avg_log10,
                         geno_com_phg = as.factor(mum_data_bact$geno_com_phg),
                         n_geno_com_phg = length(unique(mum_data_bact$geno_com_phg))),
             inits = list(sig = 3, 
                          geno_com_phg_int = rep(8, length(unique(mum_data_bact$geno_com_phg)))),
             n.adapt = 1000)
#burn-in 1000 steps
update(mum_bact_jgsmdl, n.iter = 1000)
#collect samples
mum_bact_jgs_samples <- 
  coda.samples(mum_bact_jgsmdl, 
               variable.names = c("sig", "geno_com_phg_int"), 
               50000)
mum_bact_jgs_samples_df <- as.data.frame(mum_bact_jgs_samples[[1]])

#Plot outcomes (check that has reached convergence)
traceplot(mum_bact_jgs_samples)
densplot(mum_bact_jgs_samples)

summary(mum_bact_jgs_samples)

colnames(mum_bact_jgs_samples_df)[1:16] <-
  paste(levels(as.factor(mum_data_bact$geno_com_phg)),
        "_int", sep = "")

mum_bact_contrasts <- 
  data.frame(matrix(data = c(9, 11, 1, 1,
                             9, 11, 2, 1,
                             9, 15, 1, 1,
                             9, 15, 2, 1,
                             9, 13, 1, 1,
                             9, 13, 3, 1,
                             
                             1, 3, 1, 1,
                             1, 3, 2, 1,
                             1, 7, 1, 1,
                             1, 7, 2, 1,
                             1, 5, 1, 1,
                             1, 5, 3, 1,
                             
                             10, 12, 1, 1,
                             10, 12, 2, 1,
                             10, 16, 1, 1,
                             10, 16, 2, 1,
                             10, 14, 1, 1,
                             10, 14, 3, 1,
                             
                             2, 4, 1, 1,
                             2, 4, 2, 1,
                             2, 8, 1, 1,
                             2, 8, 2, 1,
                             2, 6, 1, 1,
                             2, 6, 3, 1),
                    byrow = TRUE, ncol = 4))
colnames(mum_bact_contrasts) <- c("first_int", "second_int",
                             "divisor_first", "divisor_second")

for (i in 1:nrow(mum_bact_contrasts)) {
  mum_bact_jgs_samples_df[, paste("contrast", i, sep = "_")] <-
    (mum_bact_jgs_samples_df[, mum_bact_contrasts[i, 2]] - log10(mum_bact_contrasts[i, 4])) -
    (mum_bact_jgs_samples_df[, mum_bact_contrasts[i, 1]] - log10(mum_bact_contrasts[i, 3]))
}
                             
for (i in grep("contrast", colnames(mum_bact_jgs_samples_df))) {
  hist(mum_bact_jgs_samples_df[, i],
       main = paste(colnames(mum_bact_jgs_samples_df)[i],
                    mean(mum_bact_jgs_samples_df[, i] < 0), "< 0"),
       plot = TRUE)
  
}

mum_bact_contrasts_res <- mum_bact_contrasts
mum_bact_contrasts_res$first_int <- colnames(mum_bact_jgs_samples_df)[mum_bact_contrasts_res$first_int]
mum_bact_contrasts_res$second_int <- colnames(mum_bact_jgs_samples_df)[mum_bact_contrasts_res$second_int]
colnames(mum_bact_contrasts_res)[3:4] <- c("Adjusted", "Diff < 0")
mum_bact_contrasts_res$Adjusted <- mum_bact_contrasts_res$Adjusted != 1
mum_bact_contrasts_res$`Diff < 0` <- 
  apply(mum_bact_jgs_samples_df[, grep("contrast", colnames(mum_bact_jgs_samples_df))],
        MARGIN = 2,
        FUN = function(x) {mean(x < 0)})
mum_bact_contrasts_res <- tidyr::pivot_wider(mum_bact_contrasts_res,
                                        names_from = Adjusted,
                                        values_from = `Diff < 0`,
                                        names_prefix = "Adjusted=")

#Plotting simulated means against ea other to see if there's a pattern
for (i in 1:nrow(mum_bact_contrasts)) {
  if(mum_bact_contrasts[i, 3] == 1 &
     mum_bact_contrasts[i, 4] == 1) {
    print(plot(mum_bact_jgs_samples_df[, mum_bact_contrasts[i, 1]],
               mum_bact_jgs_samples_df[, mum_bact_contrasts[i, 2]],
               xlab = colnames(mum_bact_jgs_samples_df)[mum_bact_contrasts[i, 1]],
               ylab = colnames(mum_bact_jgs_samples_df)[mum_bact_contrasts[i, 2]],
               type = "p", pch = 20, col = rgb(0, 0, 0, 0.1), cex = 0.5,
               main = cor(mum_bact_jgs_samples_df[, mum_bact_contrasts[i, 1]],
                          mum_bact_jgs_samples_df[, mum_bact_contrasts[i, 2]])))
  }
}


##Mumford phage
mum_data_phg <- all_data_repsum[all_data_repsum$Time_day > 0 &
                                  all_data_repsum$Study == "Mumford_Friman" &
                                  all_data_repsum$Pop == "PT7", ]
model_mum_phg <- lm(log10(dens_rep_avg) ~ Bact_community*Time_day*Focal_strain,
                data = mum_data_phg)
summary(model_mum_phg)

#Bayesian
mum_data_phg$geno_com <- paste(mum_data_phg$Focal_strain,
                               mum_data_phg$Bact_community,
                               sep = "_")
mum_phg_model <- "
  model{
    for(i in 1:length(dens_rep_avg_log10)){
      dens_rep_avg_log10[i] ~ dnorm(mu[i], tau)
      mu[i] <- geno_com_int[geno_com[i]] + geno_com_slp[geno_com[i]]*Time_day[i]
    }
    tau <- 1 / (sig * sig)
    
    sig ~ dunif(0, 100)
    for(i in 1:n_geno_com_trt){
      geno_com_int[i] ~ dnorm(0, 1.0E-3)
      geno_com_slp[i] ~ dnorm(0, 1.0E-3)
    }
  }"

#build the model (this includes, by default, 1000 adaptation steps)
mum_phg_jgsmdl <- 
  jags.model(file = textConnection(mum_phg_model),
             data = list(dens_rep_avg_log10 = mum_data_phg$dens_rep_avg_log10,
                         Time_day = mum_data_phg$Time_day,
                         geno_com = as.factor(mum_data_phg$geno_com),
                         n_geno_com_trt = length(unique(mum_data_phg$geno_com))),
             inits = list(sig = 3, geno_com_int = rep(8, 8), 
                          geno_com_slp = rep(0, 8)),
             n.adapt = 1000)
#burn-in 1000 steps
update(mum_phg_jgsmdl, n.iter = 1000)
#collect samples
mum_phg_jgs_samples <- 
  coda.samples(mum_phg_jgsmdl, 
               variable.names = c("sig", "geno_com_int", "geno_com_slp"), 
               50000)
mum_phg_jgs_samples_df <- as.data.frame(mum_phg_jgs_samples[[1]])

colnames(mum_phg_jgs_samples_df)[1:8] <-
  paste(levels(as.factor(mum_data_phg$geno_com)),
        "_int", sep = "")

summary(mum_phg_jgs_samples)

mum_phg_contrasts <- 
  data.frame(matrix(data = c(5, 6, 1, 1,
                             5, 6, 2, 1,
                             5, 8, 1, 1,
                             5, 8, 2, 1,
                             5, 7, 1, 1,
                             5, 7, 3, 1,
                             
                             1, 2, 1, 1,
                             1, 2, 2, 1,
                             1, 4, 1, 1,
                             1, 4, 2, 1,
                             1, 3, 1, 1,
                             1, 3, 3, 1),
                    byrow = TRUE, ncol = 4))
colnames(mum_phg_contrasts) <- c("first_int", "second_int",
                                  "divisor_first", "divisor_second")

for (i in 1:nrow(mum_phg_contrasts)) {
  mum_phg_jgs_samples_df[, paste("contrast", i, sep = "_")] <-
    (mum_phg_jgs_samples_df[, mum_phg_contrasts[i, 2]] - log10(mum_phg_contrasts[i, 4])) -
    (mum_phg_jgs_samples_df[, mum_phg_contrasts[i, 1]] - log10(mum_phg_contrasts[i, 3]))
}

mum_phg_contrasts_res <- mum_phg_contrasts
mum_phg_contrasts_res$first_int <- colnames(mum_phg_jgs_samples_df)[mum_phg_contrasts_res$first_int]
mum_phg_contrasts_res$second_int <- colnames(mum_phg_jgs_samples_df)[mum_phg_contrasts_res$second_int]
colnames(mum_phg_contrasts_res)[3:4] <- c("Adjusted", "Diff < 0")
mum_phg_contrasts_res$Adjusted <- mum_phg_contrasts_res$Adjusted != 1
mum_phg_contrasts_res$`Diff < 0` <- 
  apply(mum_phg_jgs_samples_df[, grep("contrast", colnames(mum_phg_jgs_samples_df))],
        MARGIN = 2,
        FUN = function(x) {mean(x < 0)})
mum_phg_contrasts_res <- tidyr::pivot_wider(mum_phg_contrasts_res,
                                             names_from = Adjusted,
                                             values_from = `Diff < 0`,
                                             names_prefix = "Adjusted=")

#Alseth ----
alseth_data_bact <- all_data_repsum[all_data_repsum$Time_day > 0 &
                                      all_data_repsum$Study == "Alseth_etal" &
                                      all_data_repsum$Pop == "PA14", ]
model_alseth_bact <- lm(log10(dens_rep_avg) ~ Bact_community*Time_day,
                        alseth_data_bact)
summary(model_alseth_bact)

#What is effect of comm on phage density          
alseth_data_phg <- all_data_repsum[all_data_repsum$Time_day > 0 &
                                     all_data_repsum$Study == "Alseth_etal" &
                                     all_data_repsum$Pop == "DMS3vir", ]
model_alseth_phg <- lm(log10(dens_rep_avg) ~ Bact_community*Time_day,
                       alseth_data_phg)
summary(model_alseth_phg)

#Bayesian
alseth_phg_model <- "
  model{
    for(i in 1:length(dens_rep_avg_log10)){
      dens_rep_avg_log10[i] ~ dnorm(mu[i], tau)
      mu[i] <- com_int[Bact_community[i]] + com_slp[Bact_community[i]]*Time_day[i]
    }
    tau <- 1 / (sig * sig)
    
    sig ~ dunif(0, 100)
    for(i in 1:n_com_trt){
      com_int[i] ~ dnorm(0, 1.0E-3)
      com_slp[i] ~ dnorm(0, 1.0E-3)
    }
  }"

#build the model (this includes, by default, 1000 adaptation steps)
alseth_phg_jgsmdl <- 
  jags.model(file = textConnection(alseth_phg_model),
             data = list(dens_rep_avg_log10 = alseth_data_phg$dens_rep_avg_log10,
                         Time_day = alseth_data_phg$Time_day,
                         Bact_community = droplevels(as.factor(alseth_data_phg$Bact_community)),
                         n_com_trt = length(unique(alseth_data_phg$Bact_community))),
             inits = list(sig = 3, com_int = rep(8, 6), 
                          com_slp = rep(0, 6)),
             n.adapt = 1000)
#burn-in 1000 steps
update(alseth_phg_jgsmdl, n.iter = 1000)
#collect samples
alseth_phg_jgs_samples <- 
  coda.samples(alseth_phg_jgsmdl, 
               variable.names = c("sig", "com_int", "com_slp"), 
               50000)
alseth_phg_jgs_samples_df <- as.data.frame(alseth_phg_jgs_samples[[1]])

colnames(alseth_phg_jgs_samples_df)[1:6] <-
  paste(levels(droplevels(as.factor(alseth_data_phg$Bact_community))),
        "_int", sep = "")

summary(alseth_phg_jgs_samples)

alseth_phg_contrasts <- 
  data.frame(matrix(data = c(1, 2, 1, 1,
                             1, 2, 2, 1,
                             1, 3, 1, 1,
                             1, 3, 2, 1,
                             1, 4, 1, 1,
                             1, 4, 2, 1,
                             1, 5, 1, 1,
                             1, 5, 4, 1),
                    byrow = TRUE, ncol = 4))
colnames(alseth_phg_contrasts) <- c("first_int", "second_int",
                                 "divisor_first", "divisor_second")

for (i in 1:nrow(alseth_phg_contrasts)) {
  alseth_phg_jgs_samples_df[, paste("contrast", i, sep = "_")] <-
    (alseth_phg_jgs_samples_df[, alseth_phg_contrasts[i, 2]] - log10(alseth_phg_contrasts[i, 4])) -
    (alseth_phg_jgs_samples_df[, alseth_phg_contrasts[i, 1]] - log10(alseth_phg_contrasts[i, 3]))
}

alseth_phg_contrasts_res <- alseth_phg_contrasts
alseth_phg_contrasts_res$first_int <- colnames(alseth_phg_jgs_samples_df)[alseth_phg_contrasts_res$first_int]
alseth_phg_contrasts_res$second_int <- colnames(alseth_phg_jgs_samples_df)[alseth_phg_contrasts_res$second_int]
colnames(alseth_phg_contrasts_res)[3:4] <- c("Adjusted", "Diff < 0")
alseth_phg_contrasts_res$Adjusted <- alseth_phg_contrasts_res$Adjusted != 1
alseth_phg_contrasts_res$`Diff < 0` <- 
  apply(alseth_phg_jgs_samples_df[, grep("contrast", colnames(alseth_phg_jgs_samples_df))],
        MARGIN = 2,
        FUN = function(x) {mean(x < 0)})
alseth_phg_contrasts_res <- tidyr::pivot_wider(alseth_phg_contrasts_res,
                                            names_from = Adjusted,
                                            values_from = `Diff < 0`,
                                            names_prefix = "Adjusted=")

#Johnke ----
johnke_data_bact <- all_data_repsum[all_data_repsum$Time_day > 0 &
                                   all_data_repsum$Study == "Johnke_etal" &
                                   all_data_repsum$Pop == "Klebsiella", ]
johnke_data_bact$com_phg <- paste(johnke_data_bact$Bact_community,
                                   johnke_data_bact$Phage_presence,
                                   sep = "_")

#what is effect of comm on Kleb dens relative to monoculture
model_johnke_bact <- lm(log10(dens_rep_avg+1) ~ com_phg*Time_day,
             johnke_data_bact)
summary(model_johnke_bact)

johnke_data_bact2 <- johnke_data_bact
johnke_data_bact2$com_phg <- relevel(as.factor(johnke_data_bact2$com_phg), 
                                     ref = "K_Phage present")
model_johnke_bact2 <- 
  lm(log10(dens_rep_avg+1) ~ com_phg*Time_day,
     johnke_data_bact2)
summary(model_johnke_bact2)

#Bayesian
johnke_bact_model <- "
  model{
    for(i in 1:length(dens_rep_avg_log10)){
      dens_rep_avg_log10[i] ~ dnorm(mu[i], tau)
      mu[i] <- com_phg_int[com_phg[i]] + com_phg_slp[com_phg[i]]*Time_day[i]
    }
    tau <- 1 / (sig * sig)
    
    sig ~ dunif(0, 100)
    for(i in 1:n_com_phg_trt){
      com_phg_int[i] ~ dnorm(0, 1.0E-3)
      com_phg_slp[i] ~ dnorm(0, 1.0E-3)
    }
  }"

#build the model (this includes, by default, 1000 adaptation steps)
johnke_bact_jgsmdl <- jags.model(file = textConnection(johnke_bact_model),
                 data = list(dens_rep_avg_log10 = johnke_data_bact$dens_rep_avg_log10,
                             Time_day = johnke_data_bact$Time_day,
                             com_phg = as.factor(johnke_data_bact$com_phg),
                             n_com_phg_trt = length(unique(johnke_data_bact$com_phg))),
                 inits = list(sig = 3, com_phg_int = c(8, 8, 8, 8), 
                              com_phg_slp = c(0, 0, 0, 0)),
                 n.adapt = 1000)
#burn-in 1000 steps
update(johnke_bact_jgsmdl, n.iter = 1000)
#collect samples
johnke_bact_jgs_samples <- 
    coda.samples(johnke_bact_jgsmdl, 
                 variable.names = c("sig", "com_phg_int", "com_phg_slp"), 
                 50000)
johnke_bact_jgs_samples_df <- as.data.frame(johnke_bact_jgs_samples[[1]])

#Plot outcomes (check that has reached convergence)
traceplot(johnke_bact_jgs_samples)
densplot(johnke_bact_jgs_samples)

summary(johnke_bact_jgs_samples)

levels(as.factor(johnke_data_bact$com_phg))
#1 - K solo, no phg
#2 - K solo, phg
#3 - KPS, no phg
#4 - KPS, phg

#positive = coculture higher than monoculture (phage absent)
johnke_bact_jgs_samples_df$diff_nophg <- 
  johnke_bact_jgs_samples_df$`com_phg_int[3]` -
  johnke_bact_jgs_samples_df$`com_phg_int[1]`
#positive = coculture higher than expected from 1/3 mono (phage -)
johnke_bact_jgs_samples_df$diff_prop_nophg <- 
  johnke_bact_jgs_samples_df$`com_phg_int[3]` -
  (johnke_bact_jgs_samples_df$`com_phg_int[1]` - log10(3)) 
#positive = coculture higher than mono (phage present)
johnke_bact_jgs_samples_df$diff_phg <- 
  johnke_bact_jgs_samples_df$`com_phg_int[4]` -
  johnke_bact_jgs_samples_df$`com_phg_int[2]`
#positive = coculture higher than expected from 1/3 mono (phage +)
johnke_bact_jgs_samples_df$diff_prop_phg <- 
  johnke_bact_jgs_samples_df$`com_phg_int[4]` -
  (johnke_bact_jgs_samples_df$`com_phg_int[2]`- log10(3))

for (i in 10:13) {
  hist(johnke_bact_jgs_samples_df[, i],
       main = colnames(johnke_bact_jgs_samples_df)[i])
}

apply(johnke_bact_jgs_samples_df[, 10:13],
      MARGIN = 2, FUN = function(x) {mean(x < 0)})

temp <- tidyr::pivot_longer(johnke_bact_jgs_samples_df,
                            cols = 10:13,
                            names_to = "contrast",
                            values_to = "obs")
ggplot(data = temp,
       mapping = aes(x = -obs)) +
  geom_histogram() +
  facet_grid(contrast ~ .) +
  theme_bw() +
  geom_vline(xintercept = 0, lty = 2) +
  labs(y = "MCMC Difference in Means Count",
       x = "<-- Community suppresses more than null  |  Community suppresses less than null -->")

#Phage
johnke_data_phg <- all_data_repsum[all_data_repsum$Time_day > 0 &
                                     all_data_repsum$Study == "Johnke_etal" &
                                     all_data_repsum$Pop == "Kleb_phage", ]

#what is effect of comm on phg dens relative to monoculture
model_jon_phg <- lm(log10(dens_rep_avg) ~ Bact_community*Time_day,
                    johnke_data_phg)
summary(model_jon_phg)

johnke_phg_model <- "
  model{
    for(i in 1:length(dens_rep_avg_log10)){
      dens_rep_avg_log10[i] ~ dnorm(mu[i], tau)
      mu[i] <- com_int[Bact_community[i]] + com_slp[Bact_community[i]]*Time_day[i]
    }
    tau <- 1 / (sig * sig)
    
    sig ~ dunif(0, 100)
    for(i in 1:n_com_trt){
      com_int[i] ~ dnorm(0, 1.0E-3)
      com_slp[i] ~ dnorm(0, 1.0E-3)
    }
  }"

#build the model (this includes, by default, 1000 adaptation steps)
johnke_phg_jgsmdl <- 
  jags.model(file = textConnection(johnke_phg_model),
             data = list(dens_rep_avg_log10 = johnke_data_phg$dens_rep_avg_log10,
                         Time_day = johnke_data_phg$Time_day,
                         Bact_community = droplevels(johnke_data_phg$Bact_community),
                         n_com_trt = length(unique(johnke_data_phg$Bact_community))),
             inits = list(sig = 3, com_int = c(8, 8), com_slp = c(0, 0)),
             n.adapt = 1000)
#burn-in 1000 steps
update(johnke_phg_jgsmdl, n.iter = 1000)
#collect samples
johnke_phg_jgs_samples <- 
  coda.samples(johnke_phg_jgsmdl, 
               variable.names = c("sig", "com_int", "com_slp"), 
               30000)
johnke_phg_jgs_samples_df <- as.data.frame(johnke_phg_jgs_samples[[1]])

#Plot outcomes (check that has reached convergence)
traceplot(johnke_phg_jgs_samples)
densplot(johnke_phg_jgs_samples)

summary(johnke_phg_jgs_samples)

levels(droplevels(johnke_data_phg$Bact_community))
#1 - K solo, phg
#2 - KPS, phg

#positive = higher density in comm than monoculture (+ phage)
johnke_phg_jgs_samples_df$diff_phg <- 
  johnke_phg_jgs_samples_df$`com_int[2]` -
  johnke_phg_jgs_samples_df$`com_int[1]`
#positive = higher density in comm than expected from monoculture (+ phage)
johnke_phg_jgs_samples_df$diff_prop_phg <- 
  johnke_phg_jgs_samples_df$`com_int[2]` -
  (johnke_phg_jgs_samples_df$`com_int[1]` - log10(3))

for (i in 6:7) {
  hist(johnke_phg_jgs_samples_df[, i],
       main = colnames(johnke_phg_jgs_samples_df)[i])
}

apply(johnke_phg_jgs_samples_df[, 6:7],
      MARGIN = 2, FUN = function(x) {mean(x < 0)})


##Make citation network graph
cits <- read.csv("./Clean_data/citation_network.csv", 
                 stringsAsFactors = FALSE, row.names = 1)
#Note that format is the rows cite the columns

row.names(cits) <- gsub(x = row.names(cits), 
                            pattern = "(\\s|\\.)+", 
                            replacement = "_", fixed = FALSE)
colnames(cits) <- gsub(x = colnames(cits),
                       pattern = "(\\.)+",
                       replacement = "_", fixed = FALSE)

# net <- graph_from_incidence_matrix(cits, directed = TRUE,
#                                    mode = "out", add.names = TRUE)
net <- graph_from_adjacency_matrix(as.matrix(cits),
                                   mode = "directed",
                                   add.colnames = TRUE,
                                   add.rownames = TRUE)

set.seed(12)
png("./Plots/citation_network.png",
    width = 4, height = 4, units = "in", res = 300)
print(plot(net, 
           vertex.label = row.names(cits), vertex.label.color = "black",
           vertex.color = "black", vertex.size = 5,
           edge.color = "gray50",
           vertex.label.dist = 2, vertex.label.cex = 0.7,
           edge.arrow.size = .5, edge.arrow.width = 1.5,
           edge.width = 2, margin = c(0, 0, 0, 0)))
dev.off()

##Simulated density data ----
set.seed(1)
sim_data <- data.frame(trt = rep(c("no_comm", "comm"), each = 10),
                       dens_log10 = c(rnorm(10, mean = 9, sd = 0.6),
                                      rnorm(10, mean = 9 - log10(8), sd = 0.6)))
ggplot(data = sim_data, aes(x = trt, y = dens_log10)) +
  geom_point()

#Bayesian
sim_bact_model <- "
  model{
    for(i in 1:length(dens_log10)){
      dens_log10[i] ~ dnorm(mu[i], tau)
      mu[i] <- trt_int[trt[i]]
    }
    tau <- 1 / (sig * sig)
    
    sig ~ dunif(0, 100)
    for(i in 1:n_trt){
      trt_int[i] ~ dnorm(0, 1.0E-3)
    }
  }"

#build the model (this includes, by default, 1000 adaptation steps)
sim_bact_jgsmdl <- jags.model(file = textConnection(sim_bact_model),
                                 data = list(dens_log10 = sim_data$dens_log10,
                                             trt = as.factor(sim_data$trt),
                                             n_trt = length(unique(sim_data$trt))),
                                 inits = list(sig = 3, trt_int = c(8, 8)),
                                 n.adapt = 1000)
#burn-in 1000 steps
update(sim_bact_jgsmdl, n.iter = 1000)
#collect samples
sim_bact_jgs_samples <- 
  coda.samples(sim_bact_jgsmdl, 
               variable.names = c("sig", "trt_int"), 
               50000)
sim_bact_jgs_samples_df <- as.data.frame(sim_bact_jgs_samples[[1]])
colnames(sim_bact_jgs_samples_df)[2:3] <- c("comm", "no_comm")
sim_bact_jgs_samples_lng <- 
  tidyr::pivot_longer(sim_bact_jgs_samples_df,
                    cols = c(2, 3),
                    names_to = "trt",
                    values_to = "Bayesian_mean")
sim_bact_jgs_samples_lng$trt <- factor(sim_bact_jgs_samples_lng$trt,
                                       levels = c("no_comm", "comm"))

tiff("./Plots/Bayesian_example.tiff", width = 5, height = 3, units = "in", res = 300)
ggplot(data = sim_bact_jgs_samples_lng, aes(y = Bayesian_mean)) +
  geom_histogram() +
  geom_point(data = sim_data, mapping = aes(x = 0, y = dens_log10),
             color = "red", size = 2) +
  facet_grid(~as.factor(trt),
             labeller = as_labeller(c("no_comm" = "Community Absent",
                                         "comm" = "Community Present"))) +
  geom_hline(yintercept = mean(sim_bact_jgs_samples_df$no_comm), lty = 2) +
  geom_hline(yintercept = mean(sim_bact_jgs_samples_df$no_comm)-log10(3), lty = 3) +
  #geom_hline(yintercept = mean(sim_bact_jgs_samples_df$no_comm)-log10(6), lty = 3) +
  theme_bw() +
  labs(y = "log10 Bacterial Density") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
