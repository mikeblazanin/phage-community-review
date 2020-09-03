##Load libraries ----
library(tidyr)
library(dplyr)

##Clean up Mumford data ----

#Read in
mum_phg <- read.csv("./Raw_data/Mumford & Friman/Phagedensitydata.csv",
                    stringsAsFactors = F)
mum_bact <- read.csv("./Raw_data/Mumford & Friman/Pseudomonasdensitydata.csv",
                     stringsAsFactors = F)

##Make a community column in bact that is same as phage version
##  (no Pseudomonas included)
mum_bact$Competitor.Community <- gsub("PAO\\+?", "",
                                      gsub("lasR\\+?", "", 
                                           gsub("^PAO$", 0,
                                                gsub("^lasR$", 0, 
                                                     mum_bact$Community))))

#Make a phage presence column in phage
mum_phg$PhagePresence <- 1

#Rename Pseudomonas column in mum_bact
mum_bact$Pseudomonas <- gsub("^PAO$", "PAO1", mum_bact$Pseudomonas)

#Convert timepoints in mum_phg to actual days
mum_phg$Time_day <- mum_phg$Time*4

#Add Time_day to mum_bact (only measured bact on last day)
mum_bact$Time_day <- 16

#Add in Rep_pop column (assuming pops are listed in consistent order)
mum_bact$Rep_pop <- as.character(rep_len(1:5, nrow(mum_bact)))
mum_phg$Rep_pop <- as.character(rep_len(1:5, nrow(mum_phg)))

#Join dataframes together
mum_data <- full_join(mum_phg, mum_bact,
                      by = c("Pseudomonas", "Competitor.Community",
                             "PhagePresence", "Time_day", "Rep_pop"))

#Pivot densities into one column
mum_data <- pivot_longer(mum_data,
                         cols = c("PhagesPerML", "PseudomonasCFUperml"),
                         names_to = "Pop",
                         values_to = "Density")

##Tidy up

#Shared columns: Time_day, Pseudomonas, Competitor.Community,
#                 PhagePresence, Pop, Density, Rep_pop
#Uniq cols: Competitor community (for phg includes only non-PA, for bact includes all)
#           PhagePresence (absent from phage dataframe)
#           QS (absent from phage dataframe) - but just a shorthand for PAO vs lasR
#           Treatmentcode - absent from bact, 
#                           just shorthand for Pseudomonas - Competitor combinations
mum_data <- select(mum_data, Pseudomonas, Competitor.Community,
                   PhagePresence, Time_day, Pop, Density, Rep_pop)
mum_data$Pop[mum_data$Pop == "PhagesPerML"] <- "PT7"
mum_data$Pop[mum_data$Pop == "PseudomonasCFUperml"] <- "P_aeruginosa"
mum_data$Competitor.Community[mum_data$Competitor.Community == 0] <- "None"
colnames(mum_data)[1:2] <- c("Focal_strain", "Competitor_community")
mum_data <- mum_data[complete.cases(mum_data),]

#Add study info
mum_data$Study <- "Mumford_Friman"

##Clean up Alseth data ----

#Read in files
alseth_phg <- read.csv("./Raw_data/Alseth et al/Phage titres.csv",
                       stringsAsFactors = FALSE)
alseth_bact_list <- list(NA, NA, NA, NA)
i <- 1
for (bact_fil in c("qPCR results Evolution Experiment no. 4 T0.csv",
                   "qPCR results Evolution Experiment no. 4 T1.csv",
                   "qPCR results Evolution Experiment no. 4 T2.csv",
                   "qPCR results Evolution Experiment no. 4 T3.csv")) {
  alseth_bact_list[[i]] <- 
    read.csv(paste("./Raw_data/Alseth et al/", bact_fil, sep = ""),
             stringsAsFactors = FALSE)
  i <- i+1
}

#Add Repeat to T0
alseth_bact_list[[1]]$Rep_pop <- "Inoc"
#Convert non T0 repeat to character
for (i in 2:length(alseth_bact_list)) {
  alseth_bact_list[[i]]$Rep_pop <- as.character(alseth_bact_list[[i]]$Repeat)
}

#Add timepoint info to bact dataframes & merge to one dataframe
for (i in 1:length(alseth_bact_list)) {
  alseth_bact_list[[i]] <- cbind(alseth_bact_list[[i]],
                                 data.frame(Time_day = i-1))
}
alseth_bact <- Reduce(
  function(x, y, ...) {full_join(x, y, ...)},
  alseth_bact_list
)
alseth_bact <- alseth_bact[alseth_bact$Treatment != "", ]

#All bacterial-density measures are from with-phage treatments
alseth_bact$PhagePresence <- 1

#Clean up Treatment
alseth_bact$Treatment <- gsub(" and ", "+", alseth_bact$Treatment)
alseth_bact$Treatment <- gsub("MC", "AB+BC+SA", alseth_bact$Treatment)

#Create Pop column from Primer.set (which bact is being measured)
alseth_bact$Pop <- sapply(alseth_bact$Primer.set,
                          FUN = function(x) {strsplit(x, split = " ")[[1]][1]})

#Create Density column (qPCR was measured on 500 uL samples)
alseth_bact$Density <- alseth_bact$Quantity*2

#Add repeat information to phage
alseth_phg$Rep_pop <- as.character(rep_len(1:6, len = nrow(alseth_phg)))

#Pivot phage into tidy format
alseth_phg <- pivot_longer(alseth_phg,
                           cols = c("T0", "T1", "T2", "T3"),
                           names_to = "Time_day",
                           values_to = "Density")

#Cleanup phage for merge
alseth_phg$Time_day <- as.numeric(gsub("T", "", alseth_phg$Time_day))
alseth_phg$Treatment <- gsub(" ", "", alseth_phg$Treatment)
alseth_phg$Treatment <- gsub("MC", "AB+BC+SA", alseth_phg$Treatment)
alseth_phg$Pop <- "DMS3vir"
alseth_phg$PhagePresence <- 1

#Merge
alseth_data <- full_join(alseth_bact, alseth_phg)

#Split Treatment into Focal strain & Competitor Community
alseth_data$Focal_strain <- sapply(alseth_data$Treatment,
                                   function(x) {strsplit(x, "\\+")[[1]][1]})
alseth_data$Competitor_community <- 
  sapply(alseth_data$Treatment, function(x) {
    temp <- strsplit(x, "\\+")[[1]]
    if (length(temp) <= 1) {return("None")
    } else {return(paste(temp[2:length(temp)], collapse = "+"))
    }
  })
                                           

#Drop columns we don't need
alseth_data <- select(alseth_data, Focal_strain, Competitor_community,
                      PhagePresence, Time_day, Rep_pop, Pop, Density)

#Add study name
alseth_data$Study <- "Alseth_etal"

##Clean up Johnke data ----

#Read in file
johnke_data <- read.csv("./Raw_data/Johnke et al/Density_data.csv",
                        stringsAsFactors = F)

#Drop treatments we don't want
johnke_data <- johnke_data[johnke_data$info %in% 
                             c("no predator, CI ", "no predator, CII",
                               "Phage, I", "Phage, II"), ]

#Rename & relevel variables
johnke_data$PhagePresence <- 
  sapply(johnke_data$info, FUN = function(x) {strsplit(x, split = ", ")[[1]][1]})
johnke_data$PhagePresence <- recode(johnke_data$PhagePresence,
                                    Phage = 1, `no predator` = 0)

johnke_data$Competitor_community <- 
  sapply(johnke_data$info, FUN = function(x) {strsplit(x, split = ", ")[[1]][2]})
johnke_data$Competitor_community <- recode(johnke_data$Competitor_community,
                                           `CI ` = "None", I = "None",
                                           CII = "P+S", II = "P+S")

johnke_data$Time_day <- johnke_data$time..h/24
johnke_data$Focal_strain <- "Klebsiella"
johnke_data$Pop <- recode(johnke_data$prey.predator,
                          Phage = "Kleb_phage")
johnke_data$Density <- johnke_data$cell.count.ml
johnke_data$Rep_pop <- as.character(1:nrow(johnke_data))
johnke_data$Study <- "Johnke_etal"

#Drop unneeded columns
johnke_data <- johnke_data[ , c("Focal_strain", "Competitor_community",
                                "PhagePresence", "Time_day", "Pop", "Density",
                                "Rep_pop", "Study")]

##Merge studies together ----
all_data <- full_join(full_join(mum_data, alseth_data), johnke_data)

#Output clean data
write.csv(all_data,
          "./Clean_data/cleaned_merged_data.csv",
          row.names = FALSE)
