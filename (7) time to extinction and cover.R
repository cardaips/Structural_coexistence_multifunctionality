# R Code attached to "Stably coexisting communities deliver higher ecosystem multifunctionality" ####
## code author: Caroline Daniel
## contact: caroline.daniel@unibe.ch

### Part 7 - Compare predictions of time to extinction with cover data
### in June and August 2022

# source the cover data from the 3-species communities ####
source("(2) individual functions net effects.R")

# create time point 0 with sown diversity
cover_mean_c<-subset(cover_mean, cover_mean$nitrogen==0)
cover_mean_c$percentage.cover.o<- 0.333
cover_mean_c<-subset(cover_mean_c, cover_mean_c$plot<=48)

# load data from code from Oscar
year_ext <- read.table("data/years_to_extinction_triplets.txt", sep = "\t", header = T)
year_ext_se <- read.table("data/years_to_extinction_triplets_SE.txt", sep = "\t", header = T)

# look at the evolution of cover per species and see if it relates to the count of extinctions
count_lim_sp<-year_ext %>%
  count(splim, name = "count") %>%
  arrange(desc(count))
count_lim_sp$species<-count_lim_sp$splim
count_lim_sp$splim<-NULL

cover_mean_c$delta.june<-cover_mean_c$percentage.cover.x-cover_mean_c$percentage.cover.o
cover_mean_c$delta.august<-cover_mean_c$percentage.cover.y-cover_mean_c$percentage.cover.o
cover_mean_c$delta.mean<-(cover_mean_c$delta.june+cover_mean_c$delta.august)/2

cover_mean_c$delta.june.aug<-cover_mean_c$percentage.cover.x-cover_mean_c$percentage.cover.y

delta_cover_sp<-cover_mean_c %>%
  group_by(species) %>%
  summarize(delta.mean = mean(delta.mean, na.rm = TRUE))

comp_splim_delta<-merge(count_lim_sp,delta_cover_sp, by.y= "species")
comp_splim_delta<-comp_splim_delta %>%
  arrange(desc(count))

ggplot(comp_splim_delta, aes(x = count, y = delta.mean)) +
  geom_point() +
  geom_smooth(method = "lm")
model <- lm(delta.mean ~ count, data = comp_splim_delta)
summary(model)


# some exploration

comp_all_cover<-NULL
for (i in 1:length(sp)) {
  spfocus<-sp[i]
  sp_cover<-subset(cover_mean_c, cover_mean_c$species==spfocus)
  sp_year_ext<-subset(year_ext, grepl(spfocus, year_ext$triplets)==T)
  sp_year_ext$plot<-rownames(sp_year_ext)
  
  comp_sp<-merge(sp_year_ext,sp_cover, by.y = "plot") %>%
    select(plot,splim,year_to_extinction, delta.mean) %>%
    arrange(desc(year_to_extinction))
  
  comp_sp[is.na(comp_sp)] <- 100
  comp_sp$same<-comp_sp$splim == spfocus
  
  comp_all_cover<-rbind(comp_all_cover, comp_sp)
}

ggplot(comp_all_cover, aes(x=year_to_extinction,y=delta.mean, color= same))+
  geom_point()



ggplot(comp_hl, aes(x=year_to_extinction,y=delta.mean, color= same))+
  geom_point()



ra_focus<-subset(cover_mean_c, cover_mean_c$species=="Ra")
ra_year_ext<-subset(year_ext, grepl("Ra", year_ext$triplets)==T)
ra_year_ext$plot<-rownames(ra_year_ext)

comp_ra<-merge(ra_year_ext,ra_focus, by.y = "plot") %>%
  select(plot,splim,year_to_extinction, delta.mean) %>%
  arrange(desc(year_to_extinction)) 

comp_ra[is.na(comp_ra)] <- 100
comp_ra$same<-comp_ra$splim == "Ra"

ggplot(comp_ra, aes(x=year_to_extinction,y=delta.mean, color= same))+
  geom_point()



to_focus<-subset(cover_mean_c, cover_mean_c$species=="To")
to_year_ext<-subset(year_ext, grepl("To", year_ext$triplets)==T)
to_year_ext$plot<-rownames(to_year_ext)

comp_to<-merge(to_year_ext,to_focus, by.y = "plot") %>%
  select(plot,splim,year_to_extinction, delta.mean) %>%
  arrange(desc(year_to_extinction)) 

comp_to[is.na(comp_to)] <- 100
comp_to$same<-comp_to$splim == "To"

ggplot(comp_to, aes(x=year_to_extinction,y=delta.mean, color= same))+
  geom_point()




pt_focus<-subset(cover_mean_c, cover_mean_c$species=="Pt")
pt_year_ext<-subset(year_ext, grepl("Pt", year_ext$triplets)==T)
pt_year_ext$plot<-rownames(pt_year_ext)

comp_pt<-merge(pt_year_ext,pt_focus, by.y = "plot") %>%
  select(plot,splim,year_to_extinction, delta.mean) %>%
  arrange(desc(year_to_extinction)) 

comp_pt[is.na(comp_pt)] <- 100
comp_pt$same<-comp_pt$splim == "Pt"

ggplot(comp_pt, aes(x=year_to_extinction,y=delta.mean, color= same))+
  geom_point()



comp_multiple<-rbind(comp_hl,comp_ra, comp_pt, comp_to)


