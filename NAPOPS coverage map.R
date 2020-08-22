#### BBS route map

library(bbsBayes)
library(ggplot2)
library(sf)
library(tidyverse)

load("coords_raw.rda")

pc = st_as_sf(coords,coords = c("lon","lat"), crs = 4326)

pct <- pc %>% st_transform(crs = 102008)

bcrf <- read_sf(dsn = "C:/Users/smithac/Documents/BBSStrategicReview",
                layer = 'bcrfinalg0_Project2_Dissolve')
bcr <- bcrf %>% st_transform(st_crs(pct))


locat = system.file("maps",
                    package="bbsBayes")
map.file = "BBS_usgs_strata"

strat = sf::read_sf(dsn = locat,
                    layer = map.file)

strat1 <- strat %>% st_transform(st_crs(pct))

#pct1 = pct[1:100,]

tt <- st_intersects(strat1,pct)
tt_df <- data.frame(ST_12 = strat1$ST_12,
                    ncounts = lengths(tt),
                    sqrt_ncounts = sqrt(lengths(tt)),
                    stringsAsFactors = FALSE)
 
## read in BAM data summaries
load("bam_counts.rda")
bam = counts
bam$bcr_prov <- as.character(bam$bcr_prov)
bam$bcr_prov[c(31,41:42)] <- c("PE-14","YT-4","YT-6")
bam$ST_12 = paste0(c(rep("CA-",13),rep("US-",5),
                     rep("CA-",21),rep("US-",1),
                     rep("CA-",2)),bam$bcr_prov)
bam = bam[-which(bam$bcr_prov %in% c("MN-0","ON-0")),]
for(sst in as.character(bam$ST_12)){
  ww = which(tt_df$ST_12 == sst)
  wwb = which(bam$ST_12 == sst)
  if(length(ww) > 0){
  tt_df[ww,"ncounts"] <- tt_df[ww,"ncounts"]+bam[wwb,"count"]
}}

tt_df$nc_cat <- cut(tt_df$ncounts,breaks = c(-1,0,(c(100,200,500,1000,15000))))
tt_df$sqrt_ncounts <- sqrt(tt_df$ncounts)
strat2 <- left_join(strat1,tt_df)

strat3 <- filter(strat2,ncounts < 100)

pdf("NA-POPS_data_count_map.pdf")
mp = ggplot()+
  geom_sf(data = bcr,fill = viridis::cividis(1,begin = 1),colour = grey(0.75))+
  geom_sf(data = strat2,aes(fill = sqrt_ncounts),colour = NA)+
  geom_sf(data = strat3,colour = grey(0.75),fill = "white",alpha = 0.01)+
  scale_color_viridis_c(aesthetics = "fill",direction = -1)+
  theme(legend.position = "bottom")
print(mp)
dev.off()

# strat2$sqrt_density <- strat2$ncounts/log(strat2$AREA_1)
# 
# pdf("NA-POPS_count_density_map.pdf")
# mp = ggplot()+
#   geom_sf(data = bcr,size = 0.05,fill = viridis::cividis(1,begin = 1))+
#   geom_sf(data = strat2,aes(fill = sqrt_density),size = 0.05)+
#   scale_color_viridis_c(aesthetics = "fill",direction = -1)+
#   theme(legend.position = "bottom")
# print(mp)
# dev.off()


