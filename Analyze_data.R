library(ggplot2)

mum_data <- read.csv("./Clean_data/mumford_friman_data.csv",
                     stringsAsFactors = F)

ggplot(data = mum_data, aes(x = Time_day, y = Density+1,
                            color = Competitor.Community)) +
  geom_point(alpha = 0.5, 
             position = position_dodge(width = 1.5)) +
  facet_grid(Pop ~ PhagePresence, scales = "free_y") +
  scale_y_continuous(trans = "log10")
