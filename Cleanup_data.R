library(tidyr)

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


#Shared columns: Time, Pseudomonas, Competitor.presence
#Uniq cols: Competitor community (for phg includes only non-PA, for bact includes all)
#           PhagePresence (absent from phage dataframe)
#           QS (absent from phage dataframe) - but just a shorthand for PAO vs lasR
#           Treatmentcode - absent from bact, 
#                           just shorthand for Pseudomonas - Competitor combinations
#           