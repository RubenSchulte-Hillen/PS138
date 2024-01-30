
  #######
  
# loading the required 
library(ggplot2)
library(ggmap)
library(readxl)
AF <- read_excel("Autofim Data R.xlsx", 
                 sheet = "Ruben Autofim (3)")
AF
colnames(AF)[12] ="ID"
colnames(AF)[2] = "Date"

AF <- AF[,-8:-22]
AFcor <- subset(AF,!is.na(ID))
vecLat <- AFcor[,4]
vecLon <- AFcor[,6]

vecLat <- vecLat/100000
vecLat

vecLon <- vecLon/100000
vecLon

# creating a sample data.frame with your lat/lon points
lon <- vecLon
lat <- vecLat
df <- as.data.frame(cbind(lon,lat))
df2 <- as.data.frame(cbind(lon,lat,AFcor[,2],AFcor[,12]))

colnames(df2)[3] = "Date"


  
library(rgeos)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(viridisLite)
library(gridExtra)

# Generate shape file of the North Pole
world_Npole <- ne_countries(scale = "medium", returnclass = "sf") %>% st_transform(3574)

# Assuming you have a date column named 'Date' in df2, convert it to a suitable format
# and add it as a new column to PP_sf, fixing missing leading zeros
PP_sf <- df2 %>%
  st_as_sf(coords = c('Position LON', 'Position LAT')) %>%
  st_set_crs(4326) %>%
  mutate(
    DateFixed = as.Date(
      ifelse(
        nchar(Date) == 7, 
        paste0("0", Date),
        Date
      ),
      format = "%d%m%Y"
    )
  )

# Create a color scale based on the DateFixed column
date_scale <- scale_color_gradient(name = "Date", low = "blue", high = "red")

# Basic map with points colored by DateFixed
ggplot() + 
  theme_bw() +
  theme(legend.position = "none") +
  geom_sf(aes(), data = world_Npole) + 
  geom_sf(aes(color = DateFixed), data = PP_sf) +
  coord_sf(xlim = c(-1000000, 1500000), ylim = c(-1500000, 1000000)) +
  date_scale

