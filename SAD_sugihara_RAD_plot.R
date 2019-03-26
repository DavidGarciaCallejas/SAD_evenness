#####################
# calculate rank-abundance plots for SADs obtained with the script "SAD_sugihara_abundances_model"
#####################

# d1 <- readr::read_delim(file = "./results/abundances_niche_model_trophic_guild_DD.csv",delim = ";")
# d2 <- readr::read_delim(file = "./results/abundances_niche_model_trophic_guild_DP.csv",delim = ";")
# d3 <- readr::read_delim(file = "./results/abundances_niche_model_trophic_guild_RF.csv",delim = ";")
# model.data <- bind_rows(d1,d2,d3)
# model.data$resource.distribution <- "skewed"

model.data <- readr::read_delim(file = "./results/abundances_niche_model_trophic_guild_complete.csv",delim = ";")

# data cleaning (modify factor names, etc)
# model.data$trophic.guild <- recode(model.data$trophic.guild,"primary Producers" = "resources",
#                                    "herbivores" = "basal", "omnivores" = "intermediate", "carnivores" = "top")
model.data$trophic.guild <- factor(model.data$trophic.guild,levels = c("resources","basal","intermediate","top"))

#######
# subset richness levels (or other factor)
# model.data <- subset(model.data, richness == 100)
# model.data <- subset(model.data, trophic.guild != "resources")
#######
# some info
# sp.per.guild <- model.data %>% group_by(trophic.guild, connectance, niche.apport, resource.distribution) %>% summarise(num = n())
# num.sp.plot <- ggplot(sp.per.guild,aes(x = trophic.guild, y = num, color = niche.apport)) + 
#   geom_point()+
#   facet_grid(connectance~resource.distribution)+
#   NULL
# num.sp.plot
#######


resource.dist.levels <- unique(model.data$resource.distribution)
richness.levels <- unique(model.data$richness)
connectance.levels <- unique(model.data$connectance)
apportionment.levels <- unique(model.data$niche.apport)
trophic.guilds <- unique(model.data$trophic.guild)
replicates <- max(model.data$replicate)

rank.abundances <- model.data[,c("resource.distribution","richness","connectance","niche.apport","replicate","ID","N_predicted","trophic.guild")]
names(rank.abundances)[which(names(rank.abundances) == "N_predicted")] <- "abundance"

# for averaging
average.abundances <- NULL

# i.res.dist <- 1
# i.richness <- 1
# i.connectance <- 1
# i.apport <- 1
# i.tl <- 1

for(i.res.dist in 1:length(resource.dist.levels)){
  for(i.richness in 1:length(richness.levels)){
    for(i.connectance in 1:length(connectance.levels)){
      for(i.apport in 1:length(apportionment.levels)){
        for(i.tl in 1:length(trophic.guilds)){
          
          my.data <- subset(rank.abundances,richness == richness.levels[i.richness] &
                              resource.distribution == resource.dist.levels[i.res.dist] &
                              connectance == connectance.levels[i.connectance] &
                              niche.apport == apportionment.levels[i.apport] &
                              trophic.guild == trophic.guilds[i.tl])
          
          max.sp <- richness.levels[i.richness]
          avg.abund <- data.frame(abundance = numeric(max.sp))
          
          for(i.rep in 1:max(my.data$replicate)){
            temp.abund <- my.data$abundance[my.data$replicate == i.rep]
            temp.abund <- sort(temp.abund,decreasing = T)
            if(length(temp.abund) < max.sp){
              temp.abund[(length(temp.abund)):max.sp] <- 0
            }
            avg.abund$abundance <- avg.abund$abundance + temp.abund
          }
          avg.abund$resource.distribution <- resource.dist.levels[i.res.dist]
          avg.abund$richness <- richness.levels[i.richness]
          avg.abund$connectance <- connectance.levels[i.connectance]
          avg.abund$niche.apport <- apportionment.levels[i.apport]
          avg.abund$trophic.guild <- trophic.guilds[i.tl]
          
          avg.abund$abundance <- avg.abund$abundance/max(my.data$replicate)
          avg.abund <- subset(avg.abund,abundance > 0)
          
          average.abundances <- bind_rows(average.abundances,avg.abund)
        }
      }
    }
  }
}

average.abundances <- average.abundances[,c("resource.distribution","richness","connectance","niche.apport","trophic.guild","abundance")]

rank.abundances <- average.abundances %>% group_by(resource.distribution,richness,connectance,niche.apport,trophic.guild) %>% mutate(relative.abund = abundance/sum(abundance))
rank.abundances <- rank.abundances %>% group_by(resource.distribution,richness,connectance,niche.apport,trophic.guild) %>% mutate(species.rank = rank(-relative.abund,ties.method = "first"))
rank.abundances <- arrange(rank.abundances,resource.distribution,connectance,niche.apport,trophic.guild,desc(relative.abund))

# rename levels
rank.abundances$niche.apport <- recode(rank.abundances$niche.apport,"random.assortment" = "random assortment",
                                       "random.fraction" = "random fraction",
                                       "dominance.decay" = "dominance decay",
                                       "dominance.preemption" = "dominance preemption")
# reorder levels
rank.abundances$niche.apport <- factor(rank.abundances$niche.apport,levels = c("random assortment","dominance decay","random fraction","dominance preemption"))

#######################
# in the following lines, comment or uncomment depending on whether you want the complete plot or only a subset
rank.abundances.plot <- subset(rank.abundances, resource.distribution == "skewed" & 
                            connectance == 0.2 &
                            trophic.guild != "resources")

my.palette <- c("gray60","#009E73","#E69F00","#D55E00")
# my.palette <- c("gray80","gray60","gray40","gray20")

selected.rad.plot <- ggplot(rank.abundances.plot,aes(x = species.rank,y = relative.abund, group = niche.apport)) +
  
  geom_line(aes(color = niche.apport), size = 1.1) +#(aes(linetype = site.ID)) +
  # geom_line(aes(color = niche.apport), size = 1.1) +#(aes(linetype = site.ID)) +
  
  # geom_point(aes(fill = niche.apport),shape = 21, size = 1.5) +#(aes(shape = site.ID)) + 
  geom_point(aes(fill = niche.apport, shape = niche.apport), size = 2) +#(aes(shape = site.ID)) + 
  
  # facet_grid(trophic.guild~connectance+richness,scales = "free")+
  # facet_grid(fct_rev(trophic.guild)~connectance+resource.distribution,scales = "free")+
  facet_grid(fct_rev(trophic.guild)~.,scales = "free")+
  
  scale_shape_manual(values = c(21,22,23,24))+
  
  scale_color_manual(values = my.palette)+
  scale_fill_manual(values = my.palette)+
  
  xlim(0,20)+
  # scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
  
  xlab("species rank") + ylab("relative abundance") +
  DGC::theme_Publication()+
  theme(legend.title=element_blank())+
  theme(strip.background = element_blank())+#,strip.text.x = element_blank()) +
  #guides(color=FALSE)+#, fill = FALSE)+
  NULL

tiff(filename = "./results/images/sugihara_RAD_c02.tiff", res=600, compression = "lzw", width = 4500, height = 3000, units = "px")
selected.rad.plot
dev.off()

# ggsave(selected.rad.plot, filename = "./results/images/sugihara_RAD_c02.pdf", device = cairo_pdf, width = 22, height = 14, units = "cm")
