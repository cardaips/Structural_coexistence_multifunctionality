# R Code attached to "Stably coexisting communities deliver higher ecosystem multifunctionality" ####
## code author: Caroline Daniel
## contact: caroline.daniel@unibe.ch

### Part 3 - Multifunctionality and coexistence mechanisms analysis in control

# source previous files (cascading sourcing with (1))
source("(2) individual functions net effects.R")

# calculate the percentage of each function compared to max ####
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

# compute percentage threshold for multifunctionality ####
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
save_multi <- data_multi_threshold

corrplot(cor(cbind(all_functions[, 3:9], data_multi_threshold[, 1:9])),
  method = "number"
) # looks ok, but beta-glucosidase <-> high levels of multifunctionality = 0.6

# little plot checks
par(mfrow = c(3, 3), mar = c(0, 0, 0, 0))
hist(data_multi_threshold$threshold_0.4)
hist(data_multi_threshold$threshold_0.45)
hist(data_multi_threshold$threshold_0.5)
hist(data_multi_threshold$threshold_0.55)
hist(data_multi_threshold$threshold_0.6)
hist(data_multi_threshold$threshold_0.65)
hist(data_multi_threshold$threshold_0.7)
hist(data_multi_threshold$threshold_0.75)
hist(data_multi_threshold$threshold_0.8)

# creating the reliability data to test how a different number of functions changes multifunctionality values ####
# relative name of each function is registered
func_vec <- colnames(data_multi_plot[, 35:41])

# loop init
reliability_data <- NULL

for (l in 2:7) {
  # combination of i number of functions
  combinations <- combn(func_vec, l)

  # creation of an empty data frame that will be repeated for each number of column of combination
  data_multi_threshold_empty <- matrix(ncol = length(thresholds), nrow = nrow(data_multi_plot))
  data_multi_threshold_empty <- as.data.frame(data_multi_threshold_empty)
  data_multi_threshold_empty <- cbind(data_multi_threshold_empty, data_multi_threshold[, 10:16])
  colnames(data_multi_threshold_empty)[1:9] <- paste("threshold", as.character(thresholds), sep = "_")
  # save the combination level
  data_multi_threshold_empty$combination <- l
  data_multi_threshold_empty <- cbind(data_multi_threshold_empty, pres_matrix)
  data_multi_threshold_empty$iteration <- NA
  data_multi_threshold_full <- data_multi_threshold_empty


  for (k in 1:ncol(combinations)) {
    for (j in 1:length(thresholds)) {
      step <- thresholds[j]
      step_vector <- NULL
      for (i in 1:nrow(data_multi_plot)) {
        # I need to create a subset of data_multi_plot data frame with the combination in column k
        func <- data_multi_plot[i, 35:41]
        func <- func[, combinations[, k]]
        # how many function surpass the threshold?
        percentage_multi <- sum(func >= step, na.rm = T) / l
        step_vector[i] <- percentage_multi
      }
      data_multi_threshold_full[, j] <- step_vector
    }
    data_multi_threshold_full$iteration <- k
    reliability_data <- rbind(reliability_data, data_multi_threshold_full)
  }
} # done !

# now we need to transform the data in the long format
data_multi_threshold <- cbind(data_multi_threshold, pres_matrix)
long_data_multi_threshold <- melt(data_multi_threshold, measure.vars = c(1:9))

# threshold must be continuous
long_data_multi_threshold[c("threshold", "threshold.continuous")] <- str_split_fixed(long_data_multi_threshold$variable, "_", 2)
long_data_multi_threshold$threshold.continuous <- as.numeric(long_data_multi_threshold$threshold.continuous)

# multifunctionality only in control
long_data_multi_threshold_control <- subset(long_data_multi_threshold, long_data_multi_threshold$nitrogen == 0)
rownames(long_data_multi_threshold_control) <- NULL
long_data_multi_threshold_control$nitrogen <- NULL
pres_matrix_long_control <- long_data_multi_threshold_control[, 7:18]

long_data_multi_threshold_control_scaled <- long_data_multi_threshold_control
long_data_multi_threshold_control_scaled[, c(3:6, 22)] <- scale(long_data_multi_threshold_control[, c(3:6, 22)], center = T)

# now we can run the multimembership model in control
# min distance
formula <- "value ~ threshold.continuous*min.distance + (1 | species)"

long_multi_model <- multimembership_model(formula, pres_matrix_long_control, long_data_multi_threshold_control_scaled)
dist_multi_plot <- plot_multi(long_multi_model)

dist_multi_plot

# then all coexistence mechanisms
formula <- "value ~ threshold.continuous * (structural.niche * indirect.interactions + structural.fitness + structural.fitness:structural.niche) + (1 | species)"

long_multi_model_coex <- multimembership_model(formula, pres_matrix_long_control, long_data_multi_threshold_control_scaled)
coex_multi_plot <- plot_multi(long_multi_model_coex)

coex_multi_plot

# predictions and figure plotting ####
#minimum distance to exclusion first
new.data<-expand.grid(min.distance=seq(min(long_data_multi_threshold_control_scaled$min.distance), 
                                       max(long_data_multi_threshold_control_scaled$min.distance),
                                       length.out=50),
                      threshold.continuous=mean(long_data_multi_threshold_control_scaled$threshold.continuous),
                      species="A")

formula <- "value ~ threshold.continuous*min.distance + (1 | species)"
long_multi_model_pred <- multimembership_model(formula, pres_matrix_long_control, long_data_multi_threshold_control)

pred_dist<-predict_multifunctionality_dist(model=long_multi_model, new.data=new.data)


plotdist<-ggplot(data=pred_dist, aes(x=min.distance,y=fit))+
  geom_line(size=0.8)+
  #geom_point(data=long_data_multi_threshold,aes(x=structural.niche,y=indirect.interactions), size = 1.5, alpha = 0.035)+
  #scale_fill_distiller(palette= "YlGnBu", direction = -1, name = "multifunctionality")+
  ylab("predicted.multifunctionality")+
  xlab("minimum.distance.to.exclusion")+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.2)+
  theme_classic()
plotdist


# new.data is the dataset I use for predictions
new.data <- expand.grid(
  structural.niche = seq(min(long_data_multi_threshold_control$structural.niche),
    max(long_data_multi_threshold_control$structural.niche),
    length.out = 50
  ),
  structural.fitness = seq(min(long_data_multi_threshold_control$structural.fitness),
    max(long_data_multi_threshold_control$structural.fitness),
    length.out = 50
  ),
  threshold.continuous = mean(long_data_multi_threshold_control$threshold.continuous),
  indirect.interactions = mean(long_data_multi_threshold_control$indirect.interactions),
  species = "A"
)

# prediction model - unscaled values
formula <- "value ~ threshold.continuous * (structural.niche * indirect.interactions + structural.fitness + structural.fitness:structural.niche) + (1 | species)"
long_multi_model_coex_pred <- multimembership_model(formula, pres_matrix_long_control, long_data_multi_threshold_control)

pred <- predict_multifunctionality_coex(model = long_multi_model_coex_pred, new.data = new.data)

plot1 <- ggplot(data = pred, aes(x = structural.niche, y = structural.fitness)) +
  geom_tile(aes(fill = fit)) +
  geom_point(data = long_data_multi_threshold_control, aes(x = structural.niche, y = structural.fitness), size = 1.5, alpha = 0.035) +
  scale_fill_distiller(palette = "YlGnBu", direction = -1, name = "multifunctionality") +
  theme_classic()
plot1
