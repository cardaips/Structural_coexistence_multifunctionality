# R Code attached to "Stably coexisting communities deliver higher ecosystem multifunctionality" ####
## code author: Caroline Daniel
## contact: caroline.daniel@unibe.ch

### Part 4 - Net effect of multifunctionality and coexistence mechanisms analysis in control

# source previous files (cascading sourcing with (1))
source("(3) multifunctionality and coexistence.R")

# calculate the percentage of each function compared to max ####
data_net_effect_plot <- net_effect

# biomass
max_biomass <- head(data_net_effect_plot[order(-data_net_effect_plot$biomass), ], 5)
max_biomass <- mean(max_biomass$biomass)
data_net_effect_plot$relative.biomass <- data_net_effect_plot$biomass / max_biomass

# beta_glucosidase
max_betaglucosidase <- head(data_net_effect_plot[order(-data_net_effect_plot$betaglucosidase), ], 5)
max_betaglucosidase <- mean(max_betaglucosidase$betaglucosidase)
data_net_effect_plot$relative.betaglucosidase <- data_net_effect_plot$betaglucosidase / max_betaglucosidase

# phosphatase
max_phosphatase <- head(data_net_effect_plot[order(-data_net_effect_plot$phosphatase), ], 5)
max_phosphatase <- mean(max_phosphatase$phosphatase)
data_net_effect_plot$relative.phosphatase <- data_net_effect_plot$phosphatase / max_phosphatase

# root_biomass
max_root_biomass <- head(data_net_effect_plot[order(-data_net_effect_plot$root_biomass), ], 5)
max_root_biomass <- mean(max_root_biomass$root_biomass)
data_net_effect_plot$relative.root.biomass <- data_net_effect_plot$root_biomass / max_root_biomass

# decomposition
max_decomposition <- head(data_net_effect_plot[order(-data_net_effect_plot$decomposition), ], 5)
max_decomposition <- mean(max_decomposition$decomposition)
data_net_effect_plot$relative.decomposition <- data_net_effect_plot$decomposition / max_decomposition

# fungi_damage
data_net_effect_plot$fungi_damage <- data_net_effect_plot$fungi_damage + abs(min(data_net_effect_plot$fungi_damage))
max_fungi_damage <- head(data_net_effect_plot[order(-data_net_effect_plot$fungi_damage), ], 5)
max_fungi_damage <- mean(max_fungi_damage$fungi_damage)
data_net_effect_plot$relative.fungi.damage <- data_net_effect_plot$fungi_damage / max_fungi_damage

# herbivory_damage
data_net_effect_plot$herbivory_damage <- data_net_effect_plot$herbivory_damage + abs(min(data_net_effect_plot$herbivory_damage))
max_herbivory_damage <- head(data_net_effect_plot[order(-data_net_effect_plot$herbivory_damage), ], 5)
max_herbivory_damage <- mean(max_herbivory_damage$herbivory_damage)
data_net_effect_plot$relative.herbivory.damage <- data_net_effect_plot$herbivory_damage / max_herbivory_damage

# compute percentage threshold for multifunctionality ####
thresholds <- c(0, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8)
data_net_effect_threshold <- matrix(ncol = length(thresholds), nrow = nrow(data_net_effect_plot))
data_net_effect_threshold <- as.data.frame(data_net_effect_threshold)
data_net_effect_threshold$plot <- data_net_effect_plot$plot
data_net_effect_threshold$nitrogen <- data_net_effect_plot$nitrogen
data_net_effect_threshold$species <- data_multi$species
data_net_effect_threshold$structural.niche <- data_multi$omega
data_net_effect_threshold$structural.fitness <- data_multi$non.logged.theta
data_net_effect_threshold$indirect.interactions <- data_multi$differential
data_net_effect_threshold$min.distance <- data_multi$non.logged.min.distance

for (j in 1:length(thresholds)) {
  step <- thresholds[j]
  step_vector <- NULL
  for (i in 1:nrow(data_net_effect_plot)) {
    # how many function surpass the threshold?
    percentage_multi <- sum(data_net_effect_plot[i, 35:41] >= step, na.rm = T) / 7
    step_vector[i] <- percentage_multi
  }
  data_net_effect_threshold[, j] <- step_vector
}

# paste the name of the threshold in the columns of the data frame
colnames(data_net_effect_threshold)[1:10] <- paste("threshold", as.character(thresholds), sep = "_")
save_net_effect<-data_net_effect_threshold
plot_net<-ggplot(data_net_effect_threshold, aes(x=threshold_0))+
  geom_histogram(binwidth=0.1, color="grey18", fill= "grey18")+
  xlab("percentage of functions with a positive net effect")+
  ylab("number of plots")+
  ggtitle("Net multifunctionality, threshold = 0")+
  theme_classic()
plot_net

corrplot(cor(cbind(all_functions[,3:9], data_net_effect_threshold[,1:10])), method="number")

# little plots
par(mfrow = c(3, 3), mar=c(0.3,0.3,.3,0.3))
hist(data_net_effect_threshold$threshold_0.4)
hist(data_net_effect_threshold$threshold_0.45)
hist(data_net_effect_threshold$threshold_0.5)
hist(data_net_effect_threshold$threshold_0.55)
hist(data_net_effect_threshold$threshold_0.6)
hist(data_net_effect_threshold$threshold_0.65)
hist(data_net_effect_threshold$threshold_0.7)
hist(data_net_effect_threshold$threshold_0.75)
hist(data_net_effect_threshold$threshold_0.8)

# creating the reliability data to test how a different number of functions changes multifunctionality values ####
#relative name of each function is registered
func_vec<-colnames(data_net_effect_plot[,35:41])

#loop init
reliability_data_net_effect<-NULL

thresholds_rel<-thresholds[thresholds!=0]
for (l in 2:7){
  #combination of i number of functions
  combinations<-combn(func_vec,l)
  
  #creation of an empty data frame that will be repeated for each number of column of combination
  data_net_effect_threshold_empty <- matrix(ncol = length(thresholds_rel), nrow = nrow(data_net_effect_plot))
  data_net_effect_threshold_empty <- as.data.frame(data_net_effect_threshold_empty)
  data_net_effect_threshold_empty<-cbind(data_net_effect_threshold_empty, data_net_effect_threshold[,11:17])
  colnames(data_net_effect_threshold_empty)[1:9] <- paste("threshold", as.character(thresholds_rel), sep = "_")
  #save the combination level
  data_net_effect_threshold_empty$combination<-l
  data_net_effect_threshold_empty<-cbind(data_net_effect_threshold_empty, pres_matrix)
  data_net_effect_threshold_empty$iteration<-NA
  data_net_effect_threshold_full<-data_net_effect_threshold_empty
  
  
  for (k in 1:ncol(combinations)){
    
    for (j in 1:length(thresholds_rel)) {
      step <- thresholds[j]
      step_vector <- NULL
      for (i in 1:nrow(data_net_effect_plot)) {
        # I need to create a subset of data_net_effect_plot data frame with the combination in column k
        func<-data_net_effect_plot[i,35:41]
        func<-func[,combinations[,k]]
        # how many function surpass the threshold?
        percentage_multi <- sum(func >= step, na.rm = T) / l
        step_vector[i] <- percentage_multi
      }
      data_net_effect_threshold_full[, j] <- step_vector
    }
    data_net_effect_threshold_full$iteration<-k
    reliability_data_net_effect<-rbind(reliability_data_net_effect, data_net_effect_threshold_full)
  }
} #done!

# now we need to transform the data in the long format

data_net_effect_threshold <- cbind(data_net_effect_threshold, pres_matrix)
data_net_effect_threshold$threshold_0<-NULL
long_data_net_effect_threshold <- melt(data_net_effect_threshold, measure.vars = c(1:9))

# threshold must be continuous
long_data_net_effect_threshold[c("threshold", "threshold.continuous")] <- str_split_fixed(long_data_net_effect_threshold$variable, "_", 2)
long_data_net_effect_threshold$threshold.continuous <- as.numeric(long_data_net_effect_threshold$threshold.continuous)

# multifunctionality only in control
long_data_net_effect_threshold_control <- subset(long_data_net_effect_threshold, long_data_net_effect_threshold$nitrogen == 0)
rownames(long_data_net_effect_threshold_control) <- NULL
long_data_net_effect_threshold_control$nitrogen <- NULL
pres_matrix_long_control <- long_data_net_effect_threshold_control[, 7:18]

long_data_net_effect_threshold_control_scaled <- long_data_net_effect_threshold_control
long_data_net_effect_threshold_control_scaled[, c(3:6, 22)] <- scale(long_data_net_effect_threshold_control[, c(3:6, 22)], center = T)

# now we can run the multimembership model in control
# min distance
formula <- "value ~ threshold.continuous*min.distance + (1 | species)"

long_net_effect_model <- multimembership_model(formula, pres_matrix_long_control, long_data_net_effect_threshold_control_scaled)
dist_net_effect_plot <- plot_multi(long_net_effect_model)

dist_net_effect_plot <- dist_net_effect_plot +
  ggtitle("b.")+
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    strip.text = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12)
  )

dist_multi_plot <- dist_multi_plot +
  ggtitle("a.")+
  theme(
    axis.title.x = element_blank(),
    strip.text = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12)
  )

figs1 <- ggarrange(dist_multi_plot,dist_net_effect_plot, common.legend = T,
          legend = "right", ncol = 2, widths = c(2, 1))

annotate_figure(figs1,
                bottom = text_grob("Coefficient estimates and 95% confidence intervals", size = 14)
)

# then all coexistence mechanisms
formula <- "value ~ threshold.continuous * (structural.niche * indirect.interactions + structural.fitness + structural.fitness:structural.niche) + (1 | species)"

long_net_effect_model_coex <- multimembership_model(formula, pres_matrix_long_control, long_data_net_effect_threshold_control_scaled)
coex_net_effect_plot <- plot_multi(long_net_effect_model_coex)

coex_net_effect_plot

# predictions and figure plotting for net effects####
#minimum distance to exclusion first
new.data<-expand.grid(min.distance=seq(min(long_data_net_effect_threshold_control$min.distance), 
                                       max(long_data_net_effect_threshold_control$min.distance),
                                       length.out=50),
                      threshold.continuous=mean(long_data_net_effect_threshold_control$threshold.continuous),
                      species="A")

formula <- "value ~ threshold.continuous*min.distance + (1 | species)"
long_net_effect_model_pred <- multimembership_model(formula, pres_matrix_long_control, long_data_net_effect_threshold_control)

# extracting residuals to plot the data points in the prediction

formula <- "value ~ threshold.continuous + (1 | species)"
long_net_effect_model_resid <- multimembership_model(formula, pres_matrix_long_control, long_data_net_effect_threshold_control)

long_data_net_effect_threshold_control$resid <- resid(long_net_effect_model_resid)
long_data_net_effect_threshold_control$corrected_y <-long_data_net_effect_threshold_control$value + long_data_net_effect_threshold_control$resid

pred_net_effect_dist<-predict_multifunctionality_dist(model=long_net_effect_model_pred, new.data=new.data)

plotdist_net_effect<-ggplot(data=pred_net_effect_dist, aes(x=min.distance,y=fit))+
  geom_line(size=0.8)+
  #geom_point(data=long_data_multi_threshold,aes(x=structural.niche,y=indirect.interactions), size = 1.5, alpha = 0.035)+
  #scale_fill_distiller(palette= "YlGnBu", direction = -1, name = "multifunctionality")+
  ylab("predicted net effect")+
  xlab("minimum distance to exclusion")+
  ylim(-0.5,1.5)+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.2)+
  theme_classic()

plotdist_net_effect_points <- plotdist_net_effect +
  geom_point(
    data = long_data_net_effect_threshold_control,
    aes(x = min.distance, y = corrected_y),  
    inherit.aes = FALSE,
    size = 1.5,
    alpha = 0.1,
    color = "black"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue") +
  ggtitle("b.") +
  theme(plot.title = element_text(size = 12))
plotdist_net_effect_points

fig.2 <- ggarrange(plotdist_points,plotdist_net_effect_points)

# new.data with ND and FD to predict
new.data <- expand.grid(
  structural.niche = seq(min(long_data_net_effect_threshold_control$structural.niche),
                         max(long_data_net_effect_threshold_control$structural.niche),
                         length.out = 50
  ),
  structural.fitness = seq(min(long_data_net_effect_threshold_control$structural.fitness),
                           max(long_data_net_effect_threshold_control$structural.fitness),
                           length.out = 50
  ),
  threshold.continuous = mean(long_data_net_effect_threshold_control$threshold.continuous),
  indirect.interactions = mean(long_data_net_effect_threshold_control$indirect.interactions),
  species = "A"
)

# prediction model - unscaled values
formula <- "value ~ threshold.continuous * (structural.niche * indirect.interactions + structural.fitness + structural.fitness:structural.niche) + (1 | species)"
long_net_effect_model_coex_pred <- multimembership_model(formula, pres_matrix_long_control, long_data_net_effect_threshold_control)

pred_nd_fd_net_effect <- predict_multifunctionality_coex(model = long_net_effect_model_coex_pred, new.data = new.data)

plot_nd_fd_net_effect <- ggplot(data = pred_nd_fd_net_effect, aes(x = structural.niche, y = structural.fitness)) +
  geom_tile(aes(fill = fit)) +
  geom_point(data = long_data_multi_threshold_control, aes(x = structural.niche, y = structural.fitness), size = 1.5, alpha = 0.035) +
  scale_fill_distiller(palette = "YlGnBu", direction = -1, name = "net effect", limits = c(-0.1,0.5)) +
  theme_classic()
plot_nd_fd_net_effect

# new.data with ND and ID to predict
new.data <- expand.grid(
  structural.niche = seq(min(long_data_net_effect_threshold_control$structural.niche),
                         max(long_data_net_effect_threshold_control$structural.niche),
                         length.out = 50
  ),
  indirect.interactions = seq(min(long_data_net_effect_threshold_control$indirect.interactions),
                           max(long_data_net_effect_threshold_control$indirect.interactions),
                           length.out = 50
  ),
  threshold.continuous = mean(long_data_net_effect_threshold_control$threshold.continuous),
  structural.fitness = mean(long_data_net_effect_threshold_control$structural.fitness),
  species = "A"
)

# prediction model - unscaled values
formula <- "value ~ threshold.continuous * (structural.niche * indirect.interactions + structural.fitness + structural.fitness:structural.niche) + (1 | species)"
long_net_effect_model_coex_pred <- multimembership_model(formula, pres_matrix_long_control, long_data_net_effect_threshold_control)

pred_nd_id_net_effect <- predict_multifunctionality_coex(model = long_net_effect_model_coex_pred, new.data = new.data)

plot_nd_id_net_effect <- ggplot(data = pred_nd_id_net_effect, aes(x = structural.niche, y = indirect.interactions)) +
  geom_tile(aes(fill = fit)) +
  geom_point(data = long_data_multi_threshold_control, aes(x = structural.niche, y = indirect.interactions), size = 1.5, alpha = 0.035) +
  scale_fill_distiller(palette = "YlGnBu", direction = -1, name = "net effect", limits = c(-0.1,0.5)) +
  theme_classic()
plot_nd_id_net_effect

fig.3b <- ggarrange(plot_nd_fd_net_effect, plot_nd_id_net_effect, common.legend = T, legend = "right")

# reliability plot ####
#first we need to transform the data in the long format
long_reliability_data_net_effect<-melt(reliability_data_net_effect, measure.vars = c(1:9))

#threshold must be continuous
long_reliability_data_net_effect[c("threshold" ,"threshold.continuous")]<-str_split_fixed(long_reliability_data_net_effect$variable, "_", 2)
long_reliability_data_net_effect$threshold.continuous<-as.numeric(long_reliability_data_net_effect$threshold.continuous)
pres_matrix_long_reliability<-long_reliability_data_net_effect[,9:20]
long_reliability_data_net_effect_control<-subset(long_reliability_data_net_effect, long_reliability_data_net_effect$nitrogen==0)
rownames(long_reliability_data_net_effect_control)<-NULL
long_reliability_data_net_effect_control[,c(4:7,23,25)]<-scale(long_reliability_data_net_effect_control[,c(4:7,23,25)], center = T)
pres_matrix_long_reliability_control<-long_reliability_data_net_effect_control[,9:20]

#check overall first, it corresponds more or less to the overall estimate with small CI
#min distance only in control
#now we can run the multimembership model
# Omega and differential only
lmod <- lFormula(value ~ threshold.continuous * min.distance  + (1 | species), data = long_reliability_data_net_effect_control)
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
  sub_reliability_data<-subset(long_reliability_data_net_effect_control, long_reliability_data_net_effect_control$combination==i)
  
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

reliability_plot_net_effect<-ggplot(data=coef_data, aes(number.of.functions,min.distance.coef))+
  geom_point(alpha=0.1)+
  geom_smooth(se=F,color="#A2CD5A")+
  stat_summary(geom = "point",
               fun.y = "mean",
               col = "black",
               size = 3)+
  ylim(c(-0.1,0.3))+
  ylab("Coef. min. distance to exclusion")+
  xlab("Number of functions included")+
  ggtitle("Net effect")+
  theme_classic()
reliability_plot_net_effect
