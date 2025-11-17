# Time to extinction simulations ####
## code author:Oscar Godoy, Caroline Daniel
## contact: caroline.daniel@unibe.ch
# Measuring time to extinction, eveness and predicted abundances for communities of 3 species

## we extract to vectors here a list of species and its instrinsic growth rates (lambda)
library(dplyr)
library(ggplot2)
alpha_control <- read.table("data/species alphas and lambdas/biommatrix_june_2020_control.txt", header = T, sep = "\t", row.names = 1)
SE_alpha_control <- read.table("data/species alphas and lambdas/SE_biommatrix_june_2020_control.txt", header = T, sep = "\t", row.names = 1)
alpha_control <- alpha_control / SE_alpha_control
alpha_control <- as.data.frame(alpha_control)
lambda_control <- read.table("data/species alphas and lambdas/biomintrinsic_june_2020_control.txt", header = T, sep = "\t")
lambda_control <- as.data.frame(lambda_control)

# removing facilitation
alpha_control[alpha_control > 0] <- (-0.00001)

species <- rownames(alpha_control)
dat_triplets <- read.table("data/plot_info.txt", header = T, sep = "\t", row.names = 1)
triplets <- data.frame(dat_triplets)
colnames(triplets) <- c("sp1", "sp2", "sp3")
rownames(triplets) <- NULL

# this loop simulate the triplets
yearext <- data.frame(triplets = "", splim = "", year_to_extinction = "")
evenness_data <- data.frame(triplet = "", evenness.22 = "")
all_plots <- NULL

for (i in 1:nrow(triplets)) {
  extinction <- 0
  print(i)

  subsp <- select(alpha_control, triplets[i, 1], triplets[i, 2], triplets[i, 3])
  subsp <- subset(subsp, rownames(subsp) %in% triplets[i, 1] | rownames(subsp) %in% triplets[i, 2] | rownames(subsp) %in% triplets[i, 3])

  sublambda <- lambda_control %>% filter(row.names(lambda_control) %in% c(triplets[i, 1], triplets[i, 2], triplets[i, 3]))
  sublambda <- sublambda[, 1]

  species <- rownames(subsp)
  nyears <- 100
  N <- matrix(NA, nrow = nyears, ncol = length(species))
  colnames(N) <- species
  N[1, ] <- 1


  limy <- data.frame(year = NA, sp = NA)

  for (t in 2:nyears) {
    for (j in 1:(length(species))) {
      if (N[t - 1, j] != Inf & N[t - 1, j] != -Inf & is.na(N[t - 1, j]) == FALSE) {
        # Beverton-Holt equation:
        N[t, j] <- N[t - 1, j] * sublambda[j] / (1 - (sum(N[t - 1, ] * subsp[j, ])))

        if (N[t - 1, j] >= 0.001 & N[t, j] < 0.001 & is.na(N[t, j]) == FALSE) {
          limy <- rbind(limy, c(colnames(N)[j], t))
          extinction <- 1
        }
      }
    }
  }
  limy <- limy[-1, ]
  yearext <- rbind(yearext, c(paste(triplets[i, 1], triplets[i, 2], triplets[i, 3], sep = "_"), limy[1, 1], limy[1, 2]))

  # calculation of relative abundances and evenness
  N_df <- as.data.frame(N)
  year15 <- N_df[15, ]
  year15[year15 < 0] <- 0.001
  p <- year15 / sum(year15)
  # Calculate Shannon diversity
  H <- -sum(p * log(p))
  # Calculate Pielou's Evenness
  S <- length(year22)
  evenness <- H / log(S)
  evenness[is.na(evenness) == T] <- 0

  year15 <- p
  # pivot to long format
  year15_long <- year15 %>%
    pivot_longer(
      cols = everything(),
      names_to = "species",
      values_to = "value"
    ) %>%
    mutate(plot = i)

  # bind to the final dataframe
  all_plots <- bind_rows(all_plots, year15_long)

  # save evenness
  evenness_data[i, ] <- c(yearext[i + 1, 1], evenness)
}

yearext <- yearext[-1, ]
yearext$year_to_extinction <- as.numeric(as.character(yearext$year_to_extinction))

# yearext has NAN when there has been no extinction
write.table(yearext, "data/years_to_extinction_triplets.txt", sep = "\t", row.names = FALSE)
write.table(evenness_data, "data/bh_predicted_eveness_15y_2022.txt", sep = "\t", row.names = FALSE)
write.table(all_plots, "data/abundances_predictedBH_all_plots.txt", sep = "\t", row.names = FALSE)
