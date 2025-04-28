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
bhpred_evenness_22 <- read.table("data/bh_predicted_eveness_15y_2022.txt", sep = "\t", header = T)

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
  geom_point(size = 2, alpha = 0.7) +  # Slightly larger points with transparency
  geom_smooth(method = "lm", color = "black", fill = "gray80", alpha = 0.3) +  # Subtle shaded regression line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +  # 1:1 reference line
  labs(x = "Predicted Evenness (BH model)", 
       y = "Observed Evenness",
       title = "Observed vs. Predicted Evenness") +
  theme_minimal(base_size = 14) +  # Minimal theme with slightly larger text
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Centered, bold title
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),  # Softer gridlines
    panel.grid.minor = element_blank()  # Remove minor gridlines
  )

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

all_plots<-read.table("data/abundances_predictedBH_all_plots.txt", sep = "\t", header = T)
cover_june_c<-subset(cover_june, cover_june$nitrogen==0)
cover_june_c$nitrogen<-NULL
all_plots$percentage.cover<-all_plots$value
all_plots$value<-NULL

all_plots <- all_plots %>%
  mutate(predicted.percentage.cover = percentage.cover) %>%
  select(plot, species, predicted.percentage.cover)

cover_june_c <- cover_june_c %>%
  mutate(observed.percentage.cover = percentage.cover) %>%
  select(plot, species, observed.percentage.cover)

merged_cover_june <- left_join(all_plots, cover_june_c, by = c("plot", "species"))
merged_cover_june[is.na(merged_cover_june)] <- 0
merged_cover_june$delta.observed.percentage.cover<-merged_cover_june$observed.percentage.cover-0.33

merged_cover_june_no_extinction<-subset(merged_cover_june, merged_cover_june$predicted.percentage.cover>=0.001)

ggplot(merged_cover_june_no_extinction, aes(x = predicted.percentage.cover, y = observed.percentage.cover)) +
  geom_point(size = 2, alpha = 0.7) +  # Slightly larger, semi-transparent points
  geom_smooth(method = "lm", color = "black", fill = "gray80", alpha = 0.3) +  # Subtle regression line
  labs(
    x = "Predicted Percentage Cover",
    y = "Observed Percentage Cover",
    title = "Observed vs Predicted Cover (June)"
  ) +
  theme_minimal(base_size = 14) +  # Minimal clean theme, bigger text
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Centered, bold title
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),  # Soft major gridlines
    panel.grid.minor = element_blank()  # No minor gridlines
  )

model <- lm(observed.percentage.cover ~ predicted.percentage.cover, data = merged_cover_june_no_extinction)
summary(model)

ggplot(merged_cover_june_no_extinction, aes(x = predicted.percentage.cover, y = observed.percentage.cover, color = species)) +
  geom_point() +                              # scatter points
  geom_smooth(method = "lm", se = TRUE, aes(group = species)) +  # separate regression line for each species
  labs(
    x = "Predicted Percentage Cover",
    y = "Observed Percentage Cover",
    title = "Observed vs Predicted Cover (June) by Species"
  ) +
  theme_minimal()

model_per_species_tidy <- merged_cover_june_no_extinction %>%
  group_by(species) %>%
  do(tidy(lm(observed.percentage.cover ~ predicted.percentage.cover, data = .)))

View(model_per_species_tidy)
