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

corrplot(cor(cbind(all_functions[, 3:9], data_multi_threshold[, 1:9]))[1:7,],
  method = "number", type = "upper") # looks ok, but beta-glucosidase <-> high levels of multifunctionality = 0.6

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

# let's have a look at how the individual thresholds keep the effect consistent or not
models_list <- list()
for (i in 1:length(unique(long_data_multi_threshold_control$threshold.continuous))){
  
  threshold <- unique(long_data_multi_threshold_control$threshold.continuous)[i]
  long_data_multi <- long_data_multi_threshold_control[long_data_multi_threshold_control$threshold.continuous == threshold,]
  
  pres_temp<- long_data_multi[,7:18]
  long_data_multi$species <- rep(LETTERS[1:12], length.out = nrow(long_data_multi))
  formula <- "value ~ min.distance + (1 | species)"
  
  temp_multi_model <- multimembership_model(formula, pres_temp, long_data_multi)
  #print(summary(temp_multi_model))
  dist_multi_temp <- plot_multi(temp_multi_model)
  #print(dist_multi_temp + ggtitle(i))
  
  models_list[[i]] <- temp_multi_model
}

fixed_df <- map2_dfr(
  models_list,  seq_along(models_list),
  ~ {
    est <- fixef(.x)["min.distance"]
    ci <- confint(.x, parm = "min.distance", method = "profile")
    data.frame(
      model = .y,
      estimate = est,
      lci = ci[1],
      uci = ci[2]
    )
  }
)

fixed_df$model<-as.character(unique(long_data_multi_threshold_control$threshold.continuous))

catplot_multi <- ggplot(fixed_df, aes(x = factor(model), y = estimate, color = factor(model))) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "blue")+
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.2) +
  labs(
    title = "a.",
    x = "Model (multifunctionality threshold)",
    y = "Estimate with 95% CI"
  ) +
  theme_minimal(base_size = 14) +
  labs(color = "Model threshold") +
  coord_flip() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) #seems the trend is the same but the power to test each threshold is lower

# then all coexistence mechanisms
formula <- "value ~ threshold.continuous * (structural.niche * indirect.interactions + structural.fitness + structural.fitness:structural.niche) + (1 | species)"

long_multi_model_coex <- multimembership_model(formula, pres_matrix_long_control, long_data_multi_threshold_control_scaled)
coex_multi_plot <- plot_multi(long_multi_model_coex)

coex_multi_plot

# predictions and figure plotting ####
#minimum distance to exclusion first
new.data<-expand.grid(min.distance=seq(min(long_data_multi_threshold_control$min.distance), 
                                       max(long_data_multi_threshold_control$min.distance),
                                       length.out=50),
                      threshold.continuous=mean(long_data_multi_threshold_control$threshold.continuous),
                      species="A")

formula <- "value ~ threshold.continuous*min.distance + (1 | species)"
long_multi_model_pred <- multimembership_model(formula, pres_matrix_long_control, long_data_multi_threshold_control)

# extracting residuals to plot the data points in the prediction
formula <- "value ~  (1 | species)"
formula <- "value ~ threshold.continuous + (1 | species)"
long_multi_model_resid <- multimembership_model(formula, pres_matrix_long_control, long_data_multi_threshold_control)

long_data_multi_threshold_control$resid <- resid(long_multi_model_resid)
long_data_multi_threshold_control$corrected_y <-long_data_multi_threshold_control$value + long_data_multi_threshold_control$resid

# here it's scaled values, I use the model that ran for the effect sizes
pred_dist<-predict_multifunctionality_dist(model=long_multi_model_pred, new.data=new.data)

plotdist<-ggplot(data=pred_dist, aes(x=min.distance,y=fit))+
  geom_line(size=0.8)+
  #geom_point(data= long_data_multi_threshold_control, aes(x=min.distance, y=value), alpha = 0.05)+
  #geom_point(data=long_data_multi_threshold,aes(x=structural.niche,y=indirect.interactions), size = 1.5, alpha = 0.035)+
  #scale_fill_distiller(palette= "YlGnBu", direction = -1, name = "multifunctionality")+
  ylab("predicted multifunctionality")+
  xlab("minimum distance to exclusion")+
  ylim(-0.5,1.5)+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.2)+
  theme_classic()

plotdist_points <- plotdist +
  geom_point(
    data = long_data_multi_threshold_control,
    aes(x = min.distance, y = corrected_y),  
    inherit.aes = FALSE,
    size = 1.5,
    alpha = 0.1,
    color = "black"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
  ggtitle("a.") +
  theme(plot.title = element_text(size = 12))
plotdist_points

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

pred_nd_fd <- predict_multifunctionality_coex(model = long_multi_model_coex_pred, new.data = new.data)

plot_nd_fd <- ggplot(data = pred_nd_fd, aes(x = structural.niche, y = structural.fitness)) +
  geom_tile(aes(fill = fit)) +
  #geom_point(data = long_data_multi_threshold_control, aes(x = structural.niche, y = structural.fitness), size = 1.5, alpha = 0.035) +
  scale_fill_distiller(palette = "YlGnBu", direction = -1, name = "multifunctionality", limits = c(0.25,0.63)) +
  labs(x = "niche differences", 
       y = "fitness differences") +
  theme_classic()
plot_nd_fd

# new.data with ND and ID to predict
new.data <- expand.grid(
  structural.niche = seq(min(long_data_multi_threshold_control$structural.niche),
                         max(long_data_multi_threshold_control$structural.niche),
                         length.out = 50
  ),
  indirect.interactions = seq(min(long_data_multi_threshold_control$indirect.interactions),
                              max(long_data_multi_threshold_control$indirect.interactions),
                              length.out = 50
  ),
  threshold.continuous = mean(long_data_multi_threshold_control$threshold.continuous),
  structural.fitness = mean(long_data_multi_threshold_control$structural.fitness),
  species = "A"
)

# prediction model - unscaled values
formula <- "value ~ threshold.continuous * (structural.niche * indirect.interactions + structural.fitness + structural.fitness:structural.niche) + (1 | species)"
long_multi_model_coex_pred <- multimembership_model(formula, pres_matrix_long_control, long_data_multi_threshold_control)

pred_nd_id_multi <- predict_multifunctionality_coex(model = long_multi_model_coex_pred, new.data = new.data)

plot_nd_id_multi <- ggplot(data = pred_nd_id_multi, aes(x = structural.niche, y = indirect.interactions)) +
  geom_tile(aes(fill = fit)) +
  #geom_point(data = long_data_multi_threshold_control, aes(x = structural.niche, y = indirect.interactions), size = 1.5, alpha = 0.035) +
  scale_fill_distiller(palette = "YlGnBu", direction = -1, name = "multifunctionality", limits = c(0.25,0.63)) +
  labs(x = "niche differences", 
       y = "indirect interactions") +
  theme_classic()
plot_nd_id_multi

fig.3a <- ggarrange(plot_nd_fd, plot_nd_id_multi, common.legend = T, legend = "right")

#reliability of multifunctionality ####

#first we need to transform the data in the long format
long_reliability_data<-melt(reliability_data, measure.vars = c(1:9))

#threshold must be continuous
long_reliability_data[c("threshold" ,"threshold.continuous")]<-str_split_fixed(long_reliability_data$variable, "_", 2)
long_reliability_data$threshold.continuous<-as.numeric(long_reliability_data$threshold.continuous)
pres_matrix_long_reliability<-long_reliability_data[,9:20]
long_reliability_data_control<-subset(long_reliability_data, long_reliability_data$nitrogen==0)
rownames(long_reliability_data_control)<-NULL
long_reliability_data_control[,c(4:7,23,25)]<-scale(long_reliability_data_control[,c(4:7,23,25)], center = T)
pres_matrix_long_reliability_control<-long_reliability_data_control[,9:20]

#check overall first, it corresponds more or less to the overall estimate with small CI
#min distance only in control
#now we can run the multimembership model
# Omega and differential only
lmod <- lFormula(value ~ threshold.continuous * min.distance  + (1 | species), data = long_reliability_data_control)
# I modified the random factor by introducing "speciesmat" which is a presence/absence matrix of PaNDiv species in my communities of 5 plant species.
lmod$reTrms$Zt <- lmod$reTrms$Ztlist[[1]] <- Matrix(t(pres_matrix_long_reliability_control))
# preparing to run the optimisation of the lmer
devfun <- do.call(mkLmerDevfun, lmod)
# optimising per se
opt <- optimizeLmer(devfun)
# run the optimisation and save the model into a variable called m1, m2,... to m5 -> important !
long_reliability_threshold_control_dist <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)

modelplot(long_reliability_threshold_control_dist, coef_omit = "SD|Intercept",
          coef_rename = c("min.distance" = "min. distance to exclusion",
                          "threshold.continuous" = "multifunctionality threshold")) +
  aes(shape = ifelse(p.value > 0.05,
                     "Not significant",
                     "Significant"
  )) +
  scale_shape_manual(values = c(4, 19)) +
  theme_bw() +
  theme(legend.title = element_blank())

#loop init
coef_data<-NULL
vec_coef<-NULL
#loop saving all combination level coef separately
for (i in 2:7){
  combinations<-combn(func_vec,i)
  rint<-ncol(combinations)
  sub_reliability_data<-subset(long_reliability_data_control, long_reliability_data_control$combination==i)
  
  for (j in 1:rint){
    sub_reliability_data_iteration<-subset(sub_reliability_data,sub_reliability_data$iteration==j)
    pres_matrix_iteration<-sub_reliability_data_iteration[9:20]
    sub_reliability_data_iteration$species<-rep(LETTERS[1:12], length.out = nrow(sub_reliability_data_iteration))
    
    sub_reliability_data_iteration[,c(4:7,23,25)]<-scale(sub_reliability_data_iteration[,c(4:7,23,25)], center = T)
    lmod <- lFormula(value ~ threshold.continuous * min.distance  + (1 | species), data =   sub_reliability_data_iteration)
    # I modified the random factor by introducing "speciesmat" which is a presence/absence matrix of PaNDiv   species in my communities of 5 plant species.
    lmod$reTrms$Zt <- lmod$reTrms$Ztlist[[1]] <- Matrix(t(pres_matrix_iteration))
    # preparing to run the optimisation of the lmer
    devfun <- do.call(mkLmerDevfun, lmod)
    # optimising per se
    opt <- optimizeLmer(devfun)
    # run the optimisation and save the model into a variable called m1, m2,... to m5 -> important !
    model_iteration_control_dist <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
    coef_iteration<-as.data.frame(coef(model_iteration_control_dist)[1])[1,3]
    vec_coef<-rbind(vec_coef,coef_iteration)
    rownames(vec_coef)<-NULL
    vec_coef<-as.data.frame(vec_coef)
    vec_coef$V2<-i
  }
  coef_data<-rbind(coef_data,vec_coef)
  vec_coef<-NULL
}

colnames(coef_data)<-c("min.distance.coef","number.of.functions")

reliability_plot<-ggplot(data=coef_data, aes(number.of.functions,min.distance.coef))+
  geom_point(alpha=0.1)+
  geom_smooth(se=F,color="#A2CD5A")+
  stat_summary(geom = "point",
               fun.y = "mean",
               col = "black",
               size = 3)+
  ylim(c(-0.1,0.3))+
  ylab("Coef. min. distance to exclusion")+
  xlab("Number of functions included")+
  ggtitle("Overall multifunctionality")+
  theme_classic()
reliability_plot
