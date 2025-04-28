#%%
### Measuring time to extinction
##we extract to vectors here a list of species and its instrinsic growth rates (lambda)
library(dplyr)
library(ggplot2)
alpha_control<-read.table("data/species alphas and lambdas/biommatrix_june_2020_control.txt",header = T, sep = "\t", row.names = 1)
SE_alpha_control<-read.table("data/species alphas and lambdas/SE_biommatrix_june_2020_control.txt",header = T, sep = "\t", row.names = 1)
alpha_control<-alpha_control/SE_alpha_control
alpha_control<-as.data.frame(alpha_control)
lambda_control<-read.table("data/species alphas and lambdas/biomintrinsic_june_2020_control.txt", header = T, sep = "\t")
lambda_control<-as.data.frame(lambda_control)

alpha_control[alpha_control>0]<-(-0.00001)

species<-rownames(alpha_control)
#%%
dat_triplets<-read.table("data/plot_info.txt",header = T, sep = "\t", row.names = 1)
#triplets<-as.data.frame(strsplit(dat_triplets, "_"))##run coex_plot_... cluster analysis for dat_triplets
#triplets<-data.frame(t(triplets))
triplets<-data.frame(dat_triplets)
colnames(triplets)<-c("sp1", "sp2", "sp3")
rownames(triplets)<-NULL

#%%
##this loop simulates the whole system for nyears
# nyears = 1000
# N<-matrix(NA, nrow=nyears, ncol = length(species))
# N[1,]<-rep(1, length(species))
# lambda<-lambda_control[,1]
# #this is a LV function
# for (t in 2:nyears)  {
#  for (j in 1:(length(species))){
   
#    N[t,j]<-N[t-1,j]*lambda[j] + sum(N[t-1,]*alpha_control[j,])
   
#  }
# }

#%% this loop simulate triplets

yearext<-data.frame(triplets = "", splim = "", year_to_extinction = "")
evenness_data<-data.frame(triplet="", evenness.22="")
all_plots<-NULL

for (i in 1:nrow(triplets)){
    extinction <- 0
    print(i)
    
    subsp<-select(alpha_control, triplets[i,1], triplets[i,2], triplets[i,3])
    subsp<-subset(subsp, rownames(subsp) %in% triplets[i,1]|rownames(subsp) %in% triplets[i,2]|rownames(subsp) %in% triplets[i,3])

    sublambda<-lambda_control %>% filter(row.names(lambda_control) %in% c(triplets[i,1] , triplets[i,2] , triplets[i,3]))
    sublambda<-sublambda[,1]

    species<-rownames(subsp)
    nyears <- 100
    N<-matrix(NA, nrow=nyears, ncol = length(species))
    colnames(N)<-species
    N[1,]<-1


    limy<-data.frame(year = NA, sp = NA)

    for (t in 2:nyears)  {
        #if(extinction == 1){
        #    break
        #}else{
        for (j in 1:(length(species))){
            if (N[t-1,j]!= Inf & N[t-1,j]!= -Inf & is.na(N[t-1,j])== FALSE){
                #Lotka-Volterra
                #N[t,j]<-N[t-1,j]*sublambda[j] + sum(N[t-1,]*subsp[j,])
                #BH
                N[t,j]<-N[t-1,j]*sublambda[j]/(1-(sum(N[t-1,]*subsp[j,])))

                if(N[t-1,j]>=0.001 & N[t,j]<0.001 & is.na(N[t,j])==FALSE) {
                    limy<-rbind(limy, c(colnames(N)[j], t))
                    extinction <- 1
                }
            #}  
        }}
    }
    limy<-limy[-1,]
    yearext<-rbind(yearext, c(paste(triplets[i,1], triplets[i,2], triplets[i,3], sep="_"), limy[1,1], limy[1,2] ))
    
    #calculation of relative abundances and evenness
    N_df <- as.data.frame(N)
    year15<-N_df[15,]
    year15[year15<0]<-0.001
    p <- year15 / sum(year15)
    # Calculate Shannon diversity 
    H <- -sum(p * log(p))
    # Calculate Pielou's Evenness 
    S <- length(year22)
    evenness <- H / log(S)
    evenness[is.na(evenness)==T]<-0
    
    year15<-p
    # Pivot to long format
    year15_long <- year15 %>%
      pivot_longer(
        cols = everything(),    
        names_to = "species",
        values_to = "value"
      ) %>%
      mutate(plot = i)            
    
    # Bind to the final dataframe
    all_plots <- bind_rows(all_plots, year15_long)
    
    #save evenness 
   evenness_data[i,]<-c(yearext[i+1,1],evenness)
    
    #if(extinction == 1){
       # N <- N[1:(as.numeric(limy[1,2])+1),]
       # N_df <- as.data.frame(N)
       # N_df$Year <- 1:(as.numeric(limy[1,2])+1)
       # N_melt <- reshape2::melt(N_df, id.vars = "Year", variable.name = "Species", value.name = "Population")
##
       # myplot <- ggplot(N_melt, aes(x = Year, y = Population, color = Species)) + 
       #     geom_line(linewidth = 1.5) +
       #     ggtitle(paste(triplets[i,1], triplets[i,2], triplets[i,3], sep=" ")) +
       #     theme_bw() +
       #     theme(
       #         plot.title = element_text(size = 20),
       #         axis.title.x = element_text(size = 16),
       #         axis.title.y = element_text(size = 16),
       #         axis.text.x = element_text(size = 14),
       #         axis.text.y = element_text(size = 14),
       #         legend.title = element_text(size = 16),
       #         legend.text = element_text(size = 14)
       #    )
        #ggsave(paste(triplets[i,1], triplets[i,2], triplets[i,3], ".pdf", sep="_") , plot = myplot, path = "figures/withoutSE", dpi = 300, width = 10, height = 7)
    #}
}
#%%
yearext<-yearext[-1,]
yearext$year_to_extinction<-as.numeric(as.character(yearext$year_to_extinction))

#yearext has NAN when there has been no extinction
write.table(yearext, "data/years_to_extinction_triplets.txt", sep = "\t", row.names = FALSE)
write.table(evenness_data, "data/bh_predicted_eveness_15y_2022.txt", sep = "\t", row.names = FALSE)
write.table(all_plots, "data/abundances_predictedBH_all_plots.txt", sep = "\t", row.names = FALSE)


#%% plots
# number of NaNs
print(paste(sum(is.na(yearext$year_to_extinction)),"triplets did not go extinct",sep = " "))
# histogram on years of extinction
p1 <- ggplot(yearext, aes(x = year_to_extinction)) + geom_histogram(bins = 30) + theme_bw() + theme(
    plot.title = element_text(size = 20),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
) + ggtitle("Years to extinction of triplets") + xlab("Years") + ylab("Counts") 
print(p1)
#ggsave("figures/years_to_extinction_triplets.pdf", plot = p1, dpi = 300, width = 10, height = 7)

# histogram of the times the species go extinct
yearext_non_na <- yearext[!is.na(yearext$year_to_extinction), ]
p2 <- ggplot(yearext_non_na, aes(x = splim)) + 
    geom_bar() + 
    theme_bw() + theme(
    plot.title = element_text(size = 20),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16), 
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
) + ggtitle("Number of extinctions") + xlab("Species") + ylab("Counts")
print(p2)

#ggsave("figures/number_of_extinctions_triplets.pdf", plot = p2, dpi = 300, width = 10, height = 7)
