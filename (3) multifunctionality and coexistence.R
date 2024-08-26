# R Code attached to "Stably coexisting communities deliver higher ecosystem multifunctionality" ####
## code author: Caroline Daniel
## contact: caroline.daniel@unibe.ch

### Part 3 - Multifunctionality and coexistence mechanisms analysis in control

# source previous file (cascading sourcing with (1))
source("(2) individual functions net effects.R")

#calculate the percentage of each function compared to max
data_multi_plot <- data_multi

# biomass
max_biomass <- head(data_multi_plot[order(-data_multi_plot$biomass), ], 5)
max_biomass <- mean(max_biomass$biomass)
data_multi_plot$relative.biomass <- data_multi_plot$biomass / max_biomass

# beta_glucosidase
max_betaglucosidase <- head(data_multi_plot[order(-data_multi_plot$betaglucosidase), ], 5)
max_betaglucosidase <- mean(max_betaglucosidase$betaglucosidase)
data_multi_plot$relative.betaglucosidase <- data_multi_plot$betaglucosidase / max_betaglucosidase

# phosphatase
max_phosphatase <- head(data_multi_plot[order(-data_multi_plot$phosphatase), ], 5)
max_phosphatase <- mean(max_phosphatase$phosphatase)
data_multi_plot$relative.phosphatase <- data_multi_plot$phosphatase / max_phosphatase

# root_biomass
max_root_biomass <- head(data_multi_plot[order(-data_multi_plot$root_biomass), ], 5)
max_root_biomass <- mean(max_root_biomass$root_biomass)
data_multi_plot$relative.root.biomass <- data_multi_plot$root_biomass / max_root_biomass

# decomposition
max_decomposition <- head(data_multi_plot[order(-data_multi_plot$decomposition), ], 5)
max_decomposition <- mean(max_decomposition$decomposition)
data_multi_plot$relative.decomposition <- data_multi_plot$decomposition / max_decomposition

# fungi_damage
data_multi_plot$fungi_damage <- data_multi_plot$fungi_damage + abs(min(data_multi_plot$fungi_damage))
max_fungi_damage <- head(data_multi_plot[order(-data_multi_plot$fungi_damage), ], 5)
max_fungi_damage <- mean(max_fungi_damage$fungi_damage)
data_multi_plot$relative.fungi.damage <- data_multi_plot$fungi_damage / max_fungi_damage

# herbivory_damage
data_multi_plot$herbivory_damage <- data_multi_plot$herbivory_damage + abs(min(data_multi_plot$herbivory_damage))
max_herbivory_damage <- head(data_multi_plot[order(-data_multi_plot$herbivory_damage), ], 5)
max_herbivory_damage <- mean(max_herbivory_damage$herbivory_damage)
data_multi_plot$relative.herbivory.damage <- data_multi_plot$herbivory_damage / max_herbivory_damage

#compute percentage threshold for multifunctionality
thresholds <- c(0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8)
data_multi_threshold <- matrix(ncol = length(thresholds), nrow = nrow(data_multi_plot))
data_multi_threshold <- as.data.frame(data_multi_threshold)
data_multi_threshold$plot <- data_multi_plot$plot
data_multi_threshold$nitrogen <- data_multi_plot$nitrogen
data_multi_threshold$species <- data_multi$species
data_multi_threshold$structural.niche <- data_multi$omega
data_multi_threshold$structural.fitness <- data_multi$non.logged.theta
data_multi_threshold$indirect.interactions <- data_multi$differential
data_multi_threshold$min.distance <- data_multi$non.logged.min.distance

for (j in 1:length(thresholds)) {
  step <- thresholds[j]
  step_vector <- NULL
  for (i in 1:nrow(data_multi_plot)) {
    # how many function surpass the threshold?
    percentage_multi <- sum(data_multi_plot[i, 35:41] >= step, na.rm = T) / 7
    step_vector[i] <- percentage_multi
  }
  data_multi_threshold[, j] <- step_vector
}

# paste the name of the threshold in the columns of the data frame
colnames(data_multi_threshold)[1:9] <- paste("threshold", as.character(thresholds), sep = "_")
save_multi<-data_multi_threshold
data_multi_threshold[,c(1:9,13:16)]<-scale(data_multi_threshold[,c(1:9,13:16)], center=T)

corrplot(cor(cbind(all_functions[,3:9], data_multi_threshold[,1:9])), method="number")
# little plot
par(mfrow = c(3, 3), mar=c(0,0,0,0))
hist(data_multi_threshold$threshold_0.4)
hist(data_multi_threshold$threshold_0.45)
hist(data_multi_threshold$threshold_0.5)
hist(data_multi_threshold$threshold_0.55)
hist(data_multi_threshold$threshold_0.6)
hist(data_multi_threshold$threshold_0.65)
hist(data_multi_threshold$threshold_0.7)
hist(data_multi_threshold$threshold_0.75)
hist(data_multi_threshold$threshold_0.8)
