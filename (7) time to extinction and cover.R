# R Code attached to "Stably coexisting communities deliver higher ecosystem multifunctionality" ####
## code author: Caroline Daniel
## contact: caroline.daniel@unibe.ch

### Part 7 - Compare predictions of time to extinction with cover data
### in June and August 2022

# source the cover data from the 3-species communities ####
source("(2) individual functions net effects.R")

# create time point 0 with sown diversity
cover_mean_c <- subset(cover_mean, cover_mean$nitrogen == 0)
cover_mean_c$percentage.cover.o <- 0.333
cover_mean_c <- subset(cover_mean_c, cover_mean_c$plot <= 48)

# load data from code from Oscar
year_ext <- read.table("data/years_to_extinction_triplets.txt", sep = "\t", header = T)
year_ext_se <- read.table("data/years_to_extinction_triplets_SE.txt", sep = "\t", header = T)
pred_evenness_22 <- read.table("data/predicted_eveness_2022.txt", sep = "\t", header = T)
bhpred_evenness_22 <- read.table("data/bh_predicted_eveness_2022.txt", sep = "\t", header = T)

# look at the evolution of cover per species and see if it relates to the count of extinctions
count_lim_sp <- year_ext %>%
  count(splim, name = "count") %>%
  arrange(desc(count))
count_lim_sp$species <- count_lim_sp$splim
count_lim_sp$splim <- NULL

cover_mean_c$delta.june <- cover_mean_c$percentage.cover.x - cover_mean_c$percentage.cover.o
cover_mean_c$delta.august <- cover_mean_c$percentage.cover.y - cover_mean_c$percentage.cover.o
cover_mean_c$delta.mean <- (cover_mean_c$delta.june + cover_mean_c$delta.august) / 2

cover_mean_c$delta.june.aug <- cover_mean_c$percentage.cover.x - cover_mean_c$percentage.cover.y

delta_cover_sp <- cover_mean_c %>%
  group_by(species) %>%
  summarize(delta.mean = mean(delta.mean, na.rm = TRUE))

comp_splim_delta <- merge(count_lim_sp, delta_cover_sp, by.y = "species")
comp_splim_delta <- comp_splim_delta %>%
  arrange(desc(count))

ggplot(comp_splim_delta, aes(x = count, y = delta.mean)) +
  geom_point() +
  geom_smooth(method = "lm")
model <- lm(delta.mean ~ count, data = comp_splim_delta)
summary(model)

# let's look at predicted vs. observed evenness (a little bit of aggregation!)

obs_evenness_22 <- data.frame(triplets = pred_evenness_22$triplet, evenness = "")
for (i in 1:nrow(obs_evenness_22)) {
  plot_focus <- subset(cover_mean_c, cover_mean_c$plot == i)
  vec <- plot_focus$percentage.cover.x
  p <- vec / sum(vec)
  # Calculate Shannon diversity
  H <- -sum(p * log(p))
  # Calculate Pielou's Evenness
  S <- length(vec)
  evenness <- H / log(S)
  evenness[is.na(evenness) == T] <- 0
  obs_evenness_22[i, 2] <- evenness
}

obs_evenness_22$pred.evenness <- pred_evenness_22$evenness.22
obs_evenness_22$bhpred.evenness <- bhpred_evenness_22$evenness.22

obs_evenness_22$evenness <- as.numeric(obs_evenness_22$evenness)
obs_evenness_22$pred.evenness <- as.numeric(obs_evenness_22$pred.evenness)
obs_evenness_22$bhpred.evenness <- as.numeric(obs_evenness_22$bhpred.evenness)

ggplot(obs_evenness_22, aes(x = evenness, y = pred.evenness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", alpha = .5) + # 1:1 line
  theme_minimal()

ggplot(obs_evenness_22, aes(x = bhpred.evenness, y = evenness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", alpha = .5) + # 1:1 line
  theme_minimal()

model <- lm(evenness ~ bhpred.evenness, data = obs_evenness_22)
summary(model)

# some exploration

comp_all_cover <- NULL
for (i in 1:length(sp)) {
  spfocus <- sp[i]
  sp_cover <- subset(cover_mean_c, cover_mean_c$species == spfocus)
  sp_year_ext <- subset(year_ext, grepl(spfocus, year_ext$triplets) == T)
  sp_year_ext$plot <- rownames(sp_year_ext)

  comp_sp <- merge(sp_year_ext, sp_cover, by.y = "plot") %>%
    select(plot, splim, year_to_extinction, delta.mean) %>%
    arrange(desc(year_to_extinction))

  comp_sp[is.na(comp_sp)] <- 100
  comp_sp$same <- comp_sp$splim == spfocus

  comp_all_cover <- rbind(comp_all_cover, comp_sp)
}

ggplot(comp_all_cover, aes(x = year_to_extinction, y = delta.mean, color = same)) +
  geom_point()
