# R Code attached to "Stably coexisting communities deliver higher ecosystem multifunctionality" ####
## code author: Caroline Daniel
## contact: caroline.daniel@unibe.ch

### Part 1 - compute individual functions and analyse the relationship with
### structural coexistence measure

# renv loading and library ####
# renv::restore()
library(dplyr) # for efficient data wrangling
library(tidyr) # for efficient data wrangling
library(corrplot) # to check correlations between functions
library(modelsummary) # to easily present model outputs
library(ggplot2) # beautiful figures, of course
library(lme4) # linear model packages
library(devtools) # to manipulate lme4 in order to hack the random factors
library(reshape) # just another efficient data wrangling package
library(ghibli) # nice colors package
library(broom.mixed) # mixed model package
library(broom) # normal linear model package
library(stringr) # easy strings manipulation
library(effects) # aesthetic package
library(ggpubr) # complements ggplot2
library(scales) # complements ggplot2
library(styler) # standardized code syntax, make it pretty!
library(vegan) # to get the Shannon diversity easily
library(flextable) # make beautiful table output easily
library(officer) # save and write office documents from R
library(purrr) # extract some model estimates for plotting

# loading function from anisoFun package manually
## These functions come from Allen-Perkin et al. 2023, Ecology Letters,  doi: https://doi.org/10.1111/ele.14291
files.sources <- list.files(paste(getwd(), "/anisoFun functions", sep = ""))
files.sources <- paste("anisoFun functions/", files.sources, sep = "")
sapply(files.sources, source)
# source the multimembership random factor model function
source("multimembership_function.R")

# we need to compile the following functions in a common dataset:

## biomass in June and August 2022
## beta-glucosidase and phosphatase in June and August 2022
## herbivory and pathogen damage in June 2022
## root biomass in August 2022
## decomposition in April 2023

# loading plot information and experimental design ####
structural_coexistence <- read.table("data/structural_coexistence_experimental_design.txt", sep = "\t", header = T)
plot_information <- read.table("data/plot_info.txt", sep = "\t", header = T)

# loading cover data and putting it in long format
# loading cover of June 2022
cover_june_22 <- read.table("data/202206_%cover.txt", sep = "\t", header = T)
cover_june_22[is.na(cover_june_22)] <- 0
cover_june_22[cover_june_22 == "-"] <- NA
# calculate relative cover
cover_june_22[, 6:20] <- cover_june_22[, 6:20] / (cover_june_22$total.cover + cover_june_22$bare.ground)
# loading cover of August 2022
cover_august_22 <- read.table("data/202207_%cover.txt", sep = "\t", header = T)
cover_june_22[cover_june_22 == "x"] <- NA
cover_august_22[, 6:20] <- cover_august_22[, 6:20] / (cover_august_22$total.cover + cover_august_22$bare.ground)
cover_august_22[cover_august_22 == "x"] <- NA

# making them both in long format
cover_june <- data.frame(plot = NA, nitrogen = NA, species = NA, percentage.cover = NA)
sp <- na.omit(unique(c(cover_june_22$species.1, cover_june_22$species.2, cover_june_22$species.3)))

for (i in 1:nrow(cover_june_22)) {
  for (j in 1:length(sp)) {
    sp_focus <- sp[j]
    cover_focus <- cover_june_22[i, sp_focus]
    if (cover_focus != 0) {
      vec <- c(cover_june_22$plot_ID[i], cover_june_22$Treatment[i], sp_focus, cover_focus)
      cover_june <- rbind(cover_june, vec)
    }
  }
}

cover_june <- cover_june[-1, ]
cover_june$nitrogen[cover_june$nitrogen == "C"] <- 0
cover_june$nitrogen[cover_june$nitrogen == "N"] <- 1
cover_june$percentage.cover <- as.numeric(cover_june$percentage.cover)
cover_june$plot <- as.numeric(cover_june$plot)
cover_june$nitrogen <- as.numeric(cover_june$nitrogen)

# manual correction of small mistakes
cover_june$species[cover_june$plot == 10 & cover_june$species == "Cb"] <- "To"
cover_june$species[cover_june$plot == 20 & cover_june$species == "Cb"] <- "Dc"
cover_june$species[cover_june$plot == 42 & cover_june$species == "Fr"] <- "Hl"
cover_august <- data.frame(plot = NA, nitrogen = NA, species = NA, percentage.cover = NA)

for (i in 1:nrow(cover_august_22)) {
  for (j in 1:length(sp)) {
    sp_focus <- sp[j]
    cover_focus <- cover_august_22[i, sp_focus]
    if (cover_focus != 0) {
      vec <- c(cover_august_22$plot_ID[i], cover_august_22$treatment[i], sp_focus, cover_focus)
      cover_august <- rbind(cover_august, vec)
    }
  }
}

cover_august <- cover_august[-1, ]
cover_august$nitrogen[cover_august$nitrogen == "C"] <- 0
cover_august$nitrogen[cover_august$nitrogen == "N"] <- 1
cover_august$percentage.cover <- as.numeric(cover_august$percentage.cover)
cover_august$plot <- as.numeric(cover_august$plot)
cover_august$nitrogen <- as.numeric(cover_august$nitrogen)

# manual correction of small mistakes
cover_august$species[cover_august$plot == 10 & cover_august$species == "Cb"] <- "To"
cover_august$species[cover_august$plot == 20 & cover_august$species == "Cb"] <- "Dc"
cover_august$species[cover_august$plot == 42 & cover_august$species == "Fr"] <- "Hl"
cover_mean <- merge(cover_june, cover_august, by = c("plot", "nitrogen", "species"), all = T)
cover_mean[is.na(cover_mean) == T] <- 0
cover_mean$mean.percentage.cover <- (cover_mean$percentage.cover.x + cover_mean$percentage.cover.y) / 2

# counting how many extinctions
extinctions <- 0
plot_extinctions <- NULL
species_extinct <- NULL

for (i in 1:nrow(cover_june_22)) {
  foc_sp <- cover_june_22[i, 3:5]
  foc_sp <- as.vector(unlist(foc_sp))

  foc_june <- cover_june_22[i, which(colnames(cover_june_22) %in% foc_sp)]
  foc_aug <- cover_august_22[i, which(colnames(cover_august_22) %in% foc_sp)]

  for (j in 1:length(foc_june)) {
    if (foc_june[j] & foc_aug[j] == 0) {
      extinctions <- extinctions + 1
      plot_extinctions <- rbind(plot_extinctions, cover_june_22$plot_ID[i])
      species_extinct <- rbind(species_extinct, names(foc_june[j]))
    } else {
      extinctions <- extinctions
    }
  }
}

species_extinct <- as.data.frame(species_extinct)
colnames(species_extinct) <- "species"
species_extinct$plot <- plot_extinctions[, 1]

species_extinct <- species_extinct %>%
  regulartable() %>%
  autofit()
word_table <- read_docx() %>%
  body_add_flextable(species_extinct) %>%
  print(target = "extinction table.docx")

# loading data for triplet plots per function ####
# biomass
biomass <- read.table("data/2022_biomass_may_august.txt", sep = "\t", header = T)

# enzymatic activity
phosphatase <- read.table("data/2022_phosphatase_may_august.txt", sep = "\t", header = T)
beta_glucosidase <- read.table("data/2022_betaglucosidase_may_august.txt", sep = "\t", header = T)

# herbivory and pathogen damage
damage <- read.table("data/202205_herbivory_pathogen.txt", sep = "\t", header = T, fill = T)
damage <- subset(damage, damage$plot <= 48)
damage$date <- NULL
damage$comments <- NULL
damage$atribute <- NULL
damage$season <- "spring"

# root biomass
root_biomass <- read.table("data/202208_root_biomass.txt", sep = "\t", header = T)

# decomposition
decomposition <- read.table("data/202304_decomposition.txt", sep = "\t", header = T)

## calculate one value for each plot for each function

### biomass is done

### enzymatic activity
beta_glucosidase <- beta_glucosidase %>%
  mutate(mean_beta_glucosidase = (betaglucosidase_spring + betaglucosidase_summer) / 2)

phosphatase <- phosphatase %>%
  mutate(mean_phosphatase = (phosphatase_spring + phosphatase_summer) / 2)

### herbivory and pathogen damage
# make a new table to summarise plot and treatment
plot_damage <- NULL
plot_damage$plot <- biomass$plot
plot_damage$nitrogen <- biomass$nitrogen
plot_damage <- as.data.frame(plot_damage)
plot_damage$fungi_damage <- NA
plot_damage$herbivory_damage <- NA

for (i in 1:nrow(plot_damage)) {
  plot_selection <- subset(damage, damage$plot == plot_damage$plot[i])
  plot_selection <- subset(plot_selection, plot_selection$nitrogen == plot_damage$nitrogen[i])
  plot_selection <- plot_selection[order(plot_selection$species), ]
  plot_selection <- na.omit(plot_selection)
  plot_cover_selection <- subset(cover_june, cover_june$plot == plot_damage$plot[i])
  plot_cover_selection <- subset(plot_cover_selection, plot_cover_selection$nitrogen == plot_damage$nitrogen[i])
  plot_cover_selection <- plot_cover_selection[order(plot_cover_selection$species), ]
  plot_cover_selection <- subset(plot_cover_selection, (plot_cover_selection$species %in% plot_selection$species))
  plot_cover_selection <- plot_cover_selection[order(plot_cover_selection$species), ]
  plot_selection <- subset(plot_selection, (plot_selection$species %in% plot_cover_selection$species))
  plot_selection <- plot_selection[order(plot_selection$species), ]
  fungi_damage <- sum(plot_selection$perc_dmg_fung * plot_cover_selection$percentage.cover)
  herbi_damage <- sum(plot_selection$perc_dmg * plot_cover_selection$percentage.cover)
  plot_damage$fungi_damage[i] <- fungi_damage
  plot_damage$herbivory_damage[i] <- herbi_damage
}

# we need to log so I need to replace the 0 by a small percentage
plot_damage$herbivory_damage[plot_damage$herbivory_damage == 0] <- 0.01

### root biomass is already done

### decomposition
decomposition$difference <- decomposition$weight_before - decomposition$weight_after
plot_decomposition <- decomposition %>%
  group_by(plot.nr, nitrogen) %>%
  summarize_at(vars(difference), list(mean))
colnames(plot_decomposition) <- c("plot", "nitrogen", "delta_weight")

## check correlations between functions
all_functions <- data.frame(
  plot = biomass$plot, nitrogen = biomass$nitrogen, biomass = biomass$mean_biomass,
  betaglucosidase = beta_glucosidase$betaglucosidase_spring, phosphatase = phosphatase$phosphatase_spring,
  root_biomass = root_biomass$root_dry_weight, decomposition = plot_decomposition$delta_weight,
  fungi_damage = plot_damage$fungi_damage, herbivory_damage = plot_damage$herbivory_damage
)
corrplot(cor(all_functions[, 3:9]),
  method = "number", type = "upper", title = "a.",
  mar = c(0, 0, 1, 0)
) # looks ok

## transform function to have a normal distribution
# log what is not normally distributed
all_functions$fungi_damage <- log(all_functions$fungi_damage)
all_functions$herbivory_damage <- log(all_functions$herbivory_damage)
all_functions$root_biomass <- log(all_functions$root_biomass)
all_functions$nitrogen <- as.factor(all_functions$nitrogen)

# final functions data.frame creation ####
# load image from structural_coexistence_control_nitrogen.qmd
load("data/structural_coexistence_image.Rdata")
data_multi <- merge(all_functions, structural_coexistence_all)
data_multi <- merge(data_multi, plot_information, all = T)

data_all_triplets <- read.table("data/all_triplets_experiment.txt", header = T, sep = "\t")

# here I plot the structural coexistence
plota <- ggplot(data_all_triplets) +
  geom_point(aes(x = Omega, y = differential), alpha = 0.15) +
  geom_point(data = structural_coexistence_all, aes(x = omega, y = differential, color = non.logged.min.distance), size = 2) +
  xlab("structural niche differences") +
  ylab("indirect interactions") +
  scale_color_gradient2(low = "darkblue", mid = "beige", high = "darkred") +
  theme_classic()

plotb <- ggplot(data_all_triplets) +
  geom_point(aes(x = Omega, y = theta), alpha = 0.15) +
  geom_point(data = structural_coexistence_all, aes(x = omega, y = non.logged.theta, color = non.logged.min.distance), size = 2) +
  xlab("structural niche differences") +
  ylab("structural fitness differences") +
  scale_color_gradient2(low = "darkblue", mid = "beige", high = "darkred", name = "min. distance to exclusion") +
  theme_classic()

ggarrange(plotb, plota, common.legend = T)

# plot for concept figure 1

ggplot(data_all_triplets) +
  geom_point(
    data = subset(structural_coexistence_all, nitrogen == 0),
    aes(x = omega, y = differential, color = non.logged.min.distance),
    size = 3
  ) +
  geom_smooth(
    data = subset(structural_coexistence_all, nitrogen == 0),
    aes(x = omega, y = differential),
    method = "lm",
    se = TRUE,
    linetype = "dotted",
    color = "black",
    fill = "grey50",
    alpha = 0.1
  ) +
  xlab("structural niche differences") +
  ylab("indirect interactions") +
  scale_color_gradient2(
    name = "min. distance to exclusion",
    low = "darkblue", mid = "grey", high = "darkred"
  ) +
  theme_classic()


# model with multimembership random factor for individual functions ####
### Preparation for fancy models
# create a fake species variable, we have 12 different species
Species_fake <- rep(LETTERS[1:12], length.out = nrow(data_multi), each = 2)
data_multi$species <- Species_fake

# create a presence/absence matrix
sp <- unique(c(data_multi$sp1, data_multi$sp2, data_multi$sp3))
pres_matrix <- matrix(ncol = length(sp))
pres_matrix <- as.data.frame(pres_matrix)
colnames(pres_matrix) <- sp

for (i in 1:nrow(data_multi)) {
  for (j in 1:length(pres_matrix)) {
    if (data_multi$sp1[i] == colnames(pres_matrix)[j] | data_multi$sp2[i] == colnames(pres_matrix)[j] | data_multi$sp3[i] == colnames(pres_matrix)[j]) {
      pres_matrix[i, j] <- 1
    }
  }
}

pres_matrix[is.na(pres_matrix)] <- 0
data_multi <- cbind(data_multi, pres_matrix)

## individual functions in control, them multimembership models with minimum distance to exclusion
### functions in control
data_multi_control <- subset(data_multi, data_multi$nitrogen == 0)
rownames(data_multi) <- NULL
data_multi_control$species <- rep(LETTERS[1:12], length.out = nrow(data_multi_control))
pres_matrix_control <- data_multi_control[, 24:35]

# I need to scale the data just before running the models
data_multi_control_scaled <- data_multi_control
data_multi_control_scaled[, 3:19] <- scale(data_multi_control[, 3:19], center = T)

# run the model for all functions with multimembership random factor
function_names <- colnames(data_multi_control_scaled[, 3:9])
all_models <- NULL

for (i in 1:length(function_names)) {
  formula <- paste(function_names[i], "~ non.logged.min.distance + (1 | species)", sep = " ")
  multi_model <- multimembership_model(formula, pres_matrix_control, data_multi_control_scaled)
  all_models <- c(all_models, multi_model)
}

names(all_models) <- function_names
dist_plot <- plot_multi(all_models)
dist_plot <- dist_plot +
  scale_x_continuous(breaks = c(-0.3, 0, 0.3))
dist_plot

# now not only minimum distance to exclusion but models with all coexistence mechanisms
all_models <- NULL

for (i in 1:length(function_names)) {
  formula <- paste(function_names[i], "~ (omega * differential) + non.logged.theta + non.logged.theta:omega + (1 | species)", sep = " ")
  multi_model <- multimembership_model(formula, pres_matrix_control, data_multi_control_scaled)
  all_models <- c(all_models, multi_model)
}

names(all_models) <- function_names
coex_plot <- plot_multi(all_models)
coex_plot

# now not only main effects to check which coex mechanism is more driving (even if not significant)
all_models <- NULL

for (i in 1:length(function_names)) {
  formula <- paste(function_names[i], "~ omega + differential + non.logged.theta + (1 | species)", sep = " ")
  multi_model <- multimembership_model(formula, pres_matrix_control, data_multi_control_scaled)
  all_models <- c(all_models, multi_model)
}
names(all_models) <- function_names
coex_plot_main <- plot_multi(all_models)
coex_plot_main

model_multi_estimates <- lapply(all_models, function(m) {
  coef <- fixef(m)
  return(coef)
})
multi_estimates <- do.call(rbind, model_multi_estimates)
multi_estimates <- as.data.frame(multi_estimates)
multi_estimates$functions <- row.names(multi_estimates)

multi_estimates_plot <- ggplot(multi_estimates, aes(x = omega, y = differential, label = functions)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "blue") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "blue") +
  geom_text(vjust = -0.8, size = 4) + # Labels slightly above the points
  labs(
    x = "Estimate of ND",
    y = "Estimate of ID",
    title = "Estimates of ND vs. ID, individual functions"
  ) +
  xlim(-0.3, 0.5) +
  ylim(-0.3, 0.25) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )
