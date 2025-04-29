# R Code attached to "Stably coexisting communities deliver higher ecosystem multifunctionality" ####
## code author: Caroline Daniel
## contact: caroline.daniel@unibe.ch

### Part 6 - Look at invasion data depending on minimum distance to exclusion, ND, FD and ID
### in June 2023

# source the cover data from the 3-species communities ####
source("(2) individual functions net effects.R")

#reading invasion data from June 2023 to calculate sp richness and evenness of invasion depending on coexistence values
invasion_data<- read.table("data/202306_%cover_invaded_plots.txt", sep = "\t", header = T)
invasion_data_c <- subset(invasion_data, invasion_data$Treatment=="C")
invasion_data_c$Treatment<-NULL
rownames(invasion_data_c)<-NULL
invasion_data_c<-cbind(plot_information[,2:4], invasion_data_c)

#getting only the necessary coexistence measure but for control
structural_coexistence_c <- subset(structural_coexistence, structural_coexistence$nitrogen==0)
structural_coexistence_c$nitrogen<-NULL
rownames(structural_coexistence_c)<-NULL

#merge my table and prepare for the loop where I will add 4 measures of richness and evenness for weeds
invasion_data_coex_c<-cbind(invasion_data_c, structural_coexistence_c)
invasion_data_coex_c$cover.target<-NA
invasion_data_coex_c$cover.weeds<-NA
invasion_data_coex_c$richness.weeds<-NA
invasion_data_coex_c$shannon.weeds<-NA

v<-NULL
v2<-NULL
c<-0
invasion_data_coex_full<-NULL

for (i in 1:nrow(invasion_data_coex_c)){
print(i)
target.plot<-invasion_data_coex_c[i,]
target.species<-target.plot[,1:3]

for (j in 11:51){
  if (colnames(target.plot)[j] %in% target.species){
    cover<-target.plot[1,j]
    v<-c(v, cover)
  }
  else{
    cover2<-target.plot[1,j]
    v2<-c(v2, cover2)
    if(is.na(cover2)==F){
    c<-(c+1)
    }
  }
}

target.plot$cover.weeds<-sum(na.omit(v2))
target.plot$cover.target<-sum(na.omit(v))
target.plot$richness.weeds<-c
target.plot$shannon.weeds<-diversity(na.omit(v2), index = "shannon")/log(c)
v<-NULL 
v2<-NULL
c<-0

invasion_data_coex_full<-rbind(invasion_data_coex_full,target.plot)
}

invasion_data_coex_full$shannon.weeds[is.na(invasion_data_coex_full$shannon.weeds)] <- 0

#initialising the multimembership model
invasion_data_coex_full <- cbind(invasion_data_coex_full,pres_matrix_control)
#log for normality
invasion_data_coex_full$cover.weeds<-log(invasion_data_coex_full$cover.weeds)
#scale for comparison
invasion_data_coex_full_scaled<-invasion_data_coex_full
invasion_data_coex_full_scaled[,c(53,54,56,57,58:61)]<-scale(invasion_data_coex_full_scaled[,c(53,54,56,57,58:61)],center=T)
Species_fake_c <- rep(LETTERS[1:12], length.out = nrow(invasion_data_coex_full_scaled), each = 2)
invasion_data_coex_full_scaled$species <- Species_fake_c

invasion_names<-colnames(invasion_data_coex_full_scaled[,58:61])
all_models <- NULL

# minimum distance to exclusion only
for (i in 1:length(invasion_names)) {
  formula <- paste(invasion_names[i], "~ min.distance + (1 | species)", sep = " ")
  multi_model <- multimembership_model(formula, pres_matrix_control, invasion_data_coex_full_scaled)
  all_models <- c(all_models, multi_model)
}

names(all_models) <- invasion_names
inv_dist_plot <- plot_multi(all_models)
inv_dist_plot

# all coexistence variables
all_models <- NULL

for (i in 1:length(invasion_names)) {
  formula <- paste(invasion_names[i], "~ (omega * differential) + theta + theta:omega + (1 | species)", sep = " ")
  multi_model <- multimembership_model(formula, pres_matrix_control, invasion_data_coex_full_scaled)
  all_models <- c(all_models, multi_model)
}

names(all_models) <- invasion_names
inv_coex_plot <- plot_multi(all_models)
inv_coex_plot

