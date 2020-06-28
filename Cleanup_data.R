library(tidyr)
library(dplyr)

##Read in Mumford data ----

mum_phg <- read.csv("./Mumford & Friman data/Raw/Phagedensitydata.csv",
                    stringsAsFactors = F)
mum_bact <- read.csv("./Mumford & Friman data/Raw/Pseudomonasdensitydata.csv",
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

#Join dataframes together
mum_data <- full_join(mum_phg, mum_bact,
                      by = c("Pseudomonas", "Competitor.Community",
                             "PhagePresence", "Time_day"))

#Pivot densities into one column
mum_data <- pivot_longer(mum_data,
                         cols = c("PhagesPerML", "PseudomonasCFUperml"),
                         names_to = "Pop",
                         values_to = "Density")

##Tidy up

#Shared columns: Time_day, Pseudomonas, Competitor.Community,
#                 PhagePresence, Pop, Density
#Uniq cols: Competitor community (for phg includes only non-PA, for bact includes all)
#           PhagePresence (absent from phage dataframe)
#           QS (absent from phage dataframe) - but just a shorthand for PAO vs lasR
#           Treatmentcode - absent from bact, 
#                           just shorthand for Pseudomonas - Competitor combinations
mum_data <- select(mum_data, Pseudomonas, Competitor.Community,
                   PhagePresence, Time_day, Pop, Density)
mum_data$Pop[mum_data$Pop == "PhagesPerML"] <- "PT7"
mum_data$Pop[mum_data$Pop == "PseudomonasCFUperml"] <- "P_aeruginosa"
mum_data$Competitor.Community[mum_data$Competitor.Community == 0] <- "None"
mum_data <- mum_data[complete.cases(mum_data),]

#Output clean data
write.csv(mum_data,
          "./Clean_data/mumford_friman_data.csv",
          row.names = FALSE)
