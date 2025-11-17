# R Code attached to "Stably coexisting communities deliver higher ecosystem multifunctionality" ####
## code author: Caroline Daniel
## contact: caroline.daniel@unibe.ch

### Part 2 - compute net effects for individual functions and analyse the
### relationship with structural coexistence measure

# source the data from the 3-species communities ####
source("(1) individual functions analysis.R")

# load data on monocultures in order to calculate the net effect for each function
biomass_mono <- read.table("data/2022_biomass_may_august_mono.txt", sep = "\t", header = T)
phosphatase_mono <- read.table("data/2022_phosphatase_may_august_mono.txt", sep = "\t", header = T)
beta_glucosidase_mono <- read.table("data/2022_betaglucosidase_may_august_mono.txt", sep = "\t", header = T)
damage_mono <- read.table("data/202205_herbi_patho_ostermundigen_mono.txt", sep = "\t", header = T, fill = T)
root_biomass_mono <- read.table("data/202208_root_biomass_mono.txt", sep = "\t", header = T)
decomposition_mono <- read.table("data/202304_decomposition_mono.txt", sep = "\t", header = T)
cover_june_mono <- subset(cover_june, cover_june$plot >= 49)

# keep one value per monoculture for root biomass
root_biomass_mono <- root_biomass_mono %>%
  group_by(plot, nitrogen, species) %>%
  summarise_at(vars(mg.g), list(mean))

# do the same for damage
damage_mono$perc.cover <- cover_june_mono$percentage.cover
damage_mono$fungi.cwm <- damage_mono$perc_dmg_fung * damage_mono$perc.cover
damage_mono$herbivory.cwm <- damage_mono$perc_dmg * damage_mono$perc.cover

# calculate the decomposition per monoculture plot
decomposition_mono$difference <- decomposition_mono$mono.weight.before - decomposition_mono$mono.weight.after

# average per monoculture the enzymatic activity
beta_glucosidase_mono$mono.mean.beta.glucosidase <- (beta_glucosidase_mono$mono.betaglucosidase_spring + beta_glucosidase_mono$mono.betaglucosidase_summer) / 2
phosphatase_mono$mono.mean.phosphatase <- (phosphatase_mono$mono.phosphatase_spring + phosphatase_mono$mono.phosphatase_summer) / 2

# build a big all function predicted dataset to compare it with all function dataset observed
all_functions_predicted <- data.frame(plot = biomass$plot, nitrogen = biomass$nitrogen, biomass = NA, betaglucosidase = NA, phosphatase = NA, root_biomass = NA, decomposition = NA, fungi_damage = NA, herbivory_damage = NA)

for (i in 1:nrow(all_functions_predicted)) {
  # subset of species for each function
  species_focus <- na.omit(subset(plot_information, plot_information$plot == all_functions_predicted$plot[i]))
  species_focus$plot <- NULL
  species_predicted_plot <- NULL

  for (j in 1:length(species_focus)) {
    cover_focus <- subset(cover_mean, cover_mean$plot == all_functions_predicted$plot[i])
    cover_focus <- subset(cover_focus, cover_focus$nitrogen == all_functions_predicted$nitrogen[i])
    cover_focus_mean <- subset(cover_focus$mean.percentage.cover, cover_focus$species == as.character(species_focus[j]))
    cover_focus_june <- subset(cover_focus$percentage.cover.x, cover_focus$species == as.character(species_focus[j]))
    cover_focus_august <- subset(cover_focus$percentage.cover.y, cover_focus$species == as.character(species_focus[j]))


    # plot number corresponding to each species
    plotnr <- unique(biomass_mono$plot[biomass_mono$species == as.character(species_focus[j])])

    # functions sampled in both periods
    biomass_focus <- subset(biomass_mono, biomass_mono$plot == plotnr)
    biomass_focus <- subset(biomass_focus, biomass_focus$nitrogen == all_functions_predicted$nitrogen[i])
    biomass_predicted <- biomass_focus$mean.mono.biomass * 0.33

    betagluco_focus <- subset(beta_glucosidase_mono, beta_glucosidase_mono$plot == plotnr)
    betagluco_focus <- subset(betagluco_focus, betagluco_focus$nitrogen == all_functions_predicted$nitrogen[i])
    betagluco_predicted <- betagluco_focus$mono.betaglucosidase_spring * 0.33

    phosphatase_focus <- subset(phosphatase_mono, phosphatase_mono$plot == plotnr)
    phosphatase_focus <- subset(phosphatase_focus, phosphatase_focus$nitrogen == all_functions_predicted$nitrogen[i])
    phosphatase_predicted <- phosphatase_focus$mono.phosphatase_spring * 0.33

    decomposition_focus <- subset(decomposition_mono, decomposition_mono$plot == plotnr)
    decomposition_focus <- subset(decomposition_focus, decomposition_focus$nitrogen == all_functions_predicted$nitrogen[i])
    decomposition_predicted <- mean(decomposition_focus$difference) * 0.33

    # functions sampled in june
    damage_focus <- subset(damage_mono, damage_mono$plot == plotnr)
    damage_focus <- subset(damage_focus, damage_focus$nitrogen == all_functions_predicted$nitrogen[i])
    herbivory_predicted <- damage_focus$herbivory.cwm * 0.33
    fungi_predicted <- damage_focus$fungi.cwm * 0.33

    # function sampled in august
    root_biomass_focus <- subset(root_biomass_mono, root_biomass_mono$plot == plotnr)
    root_biomass_focus <- subset(root_biomass_focus, root_biomass_focus$nitrogen == all_functions_predicted$nitrogen[i])
    root_biomass_predicted <- root_biomass_focus$mg.g * 0.33

    species_predicted <- data.frame(biomass = biomass_predicted, betaglucosidase = betagluco_predicted, phosphatase = phosphatase_predicted, root_biomass = root_biomass_predicted, decomposition = decomposition_predicted, fungi = fungi_predicted, herbivory = herbivory_predicted, species = as.character(species_focus[j]))
    species_predicted_plot <- rbind(species_predicted_plot, species_predicted)
  }

  plot_predicted <- colSums(species_predicted_plot[1:7])
  plot_predicted <- c(all_functions_predicted$plot[i], all_functions_predicted$nitrogen[i], plot_predicted)
  all_functions_predicted[i, ] <- plot_predicted
}

# I want to log it so I need to replace 0 with very small measure
all_functions_predicted$herbivory_damage[all_functions_predicted$herbivory_damage == 0] <- 0.01

# log what is not normally distributed
all_functions_predicted$fungi_damage <- log(all_functions_predicted$fungi_damage)
all_functions_predicted$herbivory_damage <- log(all_functions_predicted$herbivory_damage)

# build the net effect dataset for each 48 plot by comparing with predicted values from monocultures ####
net_effect <- all_functions[3:9] - all_functions_predicted[3:9]
# quick check of correlations
corrplot(cor(net_effect),
  method = "number", type = "upper", title = "b.",
  mar = c(0, 0, 1, 0)
) # looks ok
net_effect$plot <- all_functions$plot
net_effect$nitrogen <- all_functions$nitrogen
net_effect <- merge(net_effect, structural_coexistence_all)
net_effect <- merge(net_effect, plot_information, all = T)

net_effect <- cbind(net_effect, pres_matrix)

# model with multimembership random factor for individual functions net effects ####
## individual functions in control, them multimembership models with minimum distance to exclusion
### functions in control
net_effect_control <- subset(net_effect, net_effect$nitrogen == 0)
rownames(net_effect_control) <- NULL
net_effect_control$species <- rep(LETTERS[1:12], length.out = nrow(net_effect_control))

# I need to scale the data just before running the models
net_effect_control_scaled <- net_effect_control
net_effect_control_scaled[, 3:19] <- scale(net_effect_control[, 3:19], center = T)

# run the model for all functions with multimembership random factor
function_names <- colnames(net_effect_control_scaled[, 3:9])
all_models <- NULL

for (i in 1:length(function_names)) {
  formula <- paste(function_names[i], "~ non.logged.min.distance + (1 | species)", sep = " ")
  multi_model <- multimembership_model(formula, pres_matrix_control, net_effect_control_scaled)
  all_models <- c(all_models, multi_model)
}

names(all_models) <- function_names
net_effect_dist_plot <- plot_multi(all_models)
net_effect_dist_plot <- net_effect_dist_plot +
  scale_x_continuous(breaks = c(-0.3, 0, 0.3))
net_effect_dist_plot

# now not only minimum distance to exclusion but models with all coexistence mechanisms
all_models <- NULL

for (i in 1:length(function_names)) {
  formula <- paste(function_names[i], "~ (omega * differential) + non.logged.theta + non.logged.theta:omega + (1 | species)", sep = " ")
  multi_model <- multimembership_model(formula, pres_matrix_control, net_effect_control_scaled)
  all_models <- c(all_models, multi_model)
}

names(all_models) <- function_names
net_effect_coex_plot <- plot_multi(all_models)
net_effect_coex_plot

# now not only main effects to check which coex mechanism is more driving (even if not significant)
all_models <- NULL

for (i in 1:length(function_names)) {
  formula <- paste(function_names[i], "~ omega + differential + non.logged.theta + (1 | species)", sep = " ")
  multi_model <- multimembership_model(formula, pres_matrix_control, net_effect_control_scaled)
  all_models <- c(all_models, multi_model)
}

names(all_models) <- function_names
net_effect_coex_plot_main <- plot_multi(all_models)
net_effect_coex_plot_main

model_net_effect_estimates <- lapply(all_models, function(m) {
  coef <- fixef(m)
  return(coef)
})
net_effect_estimates <- do.call(rbind, model_net_effect_estimates)
net_effect_estimates <- as.data.frame(net_effect_estimates)
net_effect_estimates$functions <- row.names(net_effect_estimates)

net_effect_estimates_plot <- ggplot(net_effect_estimates, aes(x = omega, y = differential, label = functions)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "blue") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "blue") +
  geom_text(vjust = -0.8, size = 4) + # Labels slightly above the points
  labs(
    x = "Estimate of ND",
    y = "Estimate of ID",
    title = "Estimates of ND vs. ID, individual net effect"
  ) +
  xlim(-0.3, 0.5) +
  ylim(-0.3, 0.25) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )
