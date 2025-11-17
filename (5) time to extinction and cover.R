# R Code attached to "Stably coexisting communities deliver higher ecosystem multifunctionality" ####
## code author: Caroline Daniel
## contact: caroline.daniel@unibe.ch

### Part 5 - Compare predictions of time to extinction with cover data
### in June 2022

# source the cover data from the 3-species communities ####
source("(2) individual functions net effects.R")

# create time point 0 with sown diversity
cover_mean_c <- subset(cover_mean, cover_mean$nitrogen == 0)
cover_mean_c$percentage.cover.o <- 0.333
cover_mean_c <- subset(cover_mean_c, cover_mean_c$plot <= 48)

# load data from code modified from Oscar's code: "time_to_extinction_tool_no_SE.R"
year_ext <- read.table("data/years_to_extinction_triplets.txt", sep = "\t", header = T)
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

# let's look at predicted vs. observed evenness (a little bit of aggregation!)

obs_evenness_22 <- data.frame(triplets = bhpred_evenness_22$triplet, evenness = "")
for (i in 1:nrow(obs_evenness_22)) {
  plot_focus <- subset(cover_mean_c, cover_mean_c$plot == i)
  vec <- plot_focus$percentage.cover.x
  p <- vec / sum(vec)
  # calculate Shannon diversity
  H <- -sum(p * log(p))
  # calculate Pielou's Evenness
  S <- length(vec)
  evenness <- H / log(S)
  evenness[is.na(evenness) == T] <- 0
  obs_evenness_22[i, 2] <- evenness
}

obs_evenness_22$bhpred.evenness <- bhpred_evenness_22$evenness.22

obs_evenness_22$evenness <- as.numeric(obs_evenness_22$evenness)
obs_evenness_22$bhpred.evenness <- as.numeric(obs_evenness_22$bhpred.evenness)


ggplot(obs_evenness_22, aes(x = bhpred.evenness, y = evenness)) +
  geom_point(size = 2, alpha = 0.2) + # Slightly larger points with transparency
  geom_smooth(method = "lm", color = "black", fill = "gray80", alpha = 0.3) + # Subtle shaded regression line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", alpha = 0.5) + # 1:1 reference line
  labs(
    x = "Predicted Evenness (BH model)",
    y = "Observed Evenness",
    title = "Observed vs. Predicted Evenness"
  ) +
  theme_minimal(base_size = 14) + # Minimal theme with slightly larger text
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Centered, bold title
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )

model <- lm(evenness ~ bhpred.evenness, data = obs_evenness_22)
summary(model)
