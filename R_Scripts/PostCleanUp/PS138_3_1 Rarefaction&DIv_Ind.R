############################################################################################
###  PS138_Arcwatch1###
#https://github.com/matthiaswietz/FRAM_eDNA/blob/main/RarefacDiversity.R
############################################################################################

# This script: rarefaction and alpha-diversity 
setwd("~/PS138")
load("~/PS138/PS138_16S_final.Rdata")
load("ampvis_subsets.Rdata")


iNEXT.bac <- readRDS("iNext_bac_final.RDS")
# Load packages

library(iNEXT)
library(olsrr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tibble)
library(gtools)

############################################################################################
###  RAREFACTION AND COVERAGE  ###
############################################################################################


###################################
"
iNEXT.bac <- iNEXT(
  ASV.bac, datatype="abundance", 
  q=c(0), conf=0.95, nboot=100)
"
###################################
ENV.dc <- ampvis.bac$metadata
ASV.bac <- ampvis.bac$abund
TAX.bac <- ampvis.bac$tax
###################################
## RAREFACTION ##
rarefac.bac <- fortify(iNEXT.bac, type=1) %>%
  right_join(ENV.dc, by=c("Assemblage"="CTD #"))

rarefac.line.bac <- rarefac.bac[which(
  rarefac.bac$Method != "Observed"),]
rarefac.line.bac$Method <- factor(
  rarefac.line.bac$Method,
  c("Interpolated","Extrapolated"),
  c("interpolation","extrapolation"))

###################################

## COVERAGE ##
cover.bac <- fortify(iNEXT.bac, type=2) %>%
  right_join(ENV.dc, by=c("Assemblage"="CTD #"))

cover.line.bac <- cover.bac [which(
  cover.bac$Method != "Observed"),]
cover.line.bac$Method <- factor(
  cover.line.bac$Method,
  c("Interpolated","Extrapolation"),
  c("Interpolation","Extrapolation"))


###################################

## COMBINE + PLOT

coverage <- cover.bac
cover.line <- cover.line.bac
rarefaction <- rarefac.bac
rarefac.line <- rarefac.line.bac
###################################

ENV.dc$depth_category <- factor(
  ENV.dc$depth_category,
  levels = c(
    "MP","Ice-TS","Ice-BS","UIW", "2", "10", "chl-max","25", "50",  "100", "200", "500",
    "1000", "1500", "2000", "3000","20m from bottom", "5m from bottom",  "bottom"))
iNext2 <- iNEXT.bac
mergded_df <- merge(iNEXT.bac$DataInfo, ENV.dc, by.x = "Assemblage", by.y = "CTD #", all.x = TRUE)
iNext2$DataInfo$depth_category <- mergded_df$depth_category

###############################
# Define custom color palette (from NMDS plot)
ice_colors <- colorRampPalette(c("#00BCD4","#0D47A1","#000033"))(3)
euphotic_colors <- colorRampPalette(c("orange","lawngreen","darkgreen","black"))(8) # Gradient from green to black
deep_ocean_colors <- colorRampPalette(c("#FF0FF0","red", "#330000","black"))(8)

custom_palette <- c(ice_colors, euphotic_colors, deep_ocean_colors)


custom_order <- c(
  "MP","Ice-TS","Ice-BS","UIW","2", "10", "chl-max","25", "50", "100", "200", "500",
  "1000", "1500", "2000", "3000","20m from bottom", "5m from bottom",  "bottom")
###################################
#####INext Intern plot######

fusion_df <- merge(iNext2$DataInfo, ENV.dc, by.x = "Assemblage", by.y = "CTD #", all.x = TRUE)
iNext2$DataInfo$depth_category <- fusion_df$depth_category
iNext2$DataInfo$depth_category <- ENV.dc$depth_category[match(iNext2$DataInfo$Assemblage, ENV.dc$`CTD #`)]


ggiNEXT(iNext2, type=1)+#,colour = "depth_category" )+
  # scale_colour_manual(values = custom_colors)+
  labs(x="Sample size", y="Species richness") +
  theme_bw() + 
  theme(legend.position="none")
#scale_colour_manual(aesthetics = "colour",values = custom_colors, breaks = custom_order, labels = custom_order)

# Now create the plot
"ggiNEXT(iNext2, type=1) +
  aes(color = iNext2$DataInfo$depth_category) +
  labs(x="Sample size", y="Species richness") +
  theme_bw() +
  theme(legend.position="right",
        legend.title = element_blank())+
  scale_colour_manual(values = custom_palette,
                      breaks = custom_order,
                      labels = custom_order)+
  guides(color = guide_legend(ncol = 1))"

plot(iNEXT.bac, type=2)#, colour = iNEXT.bac$depth_category) # S3 method
ggiNEXT(iNEXT.bac, type=2)+labs(x="Sample size", y="Completness (coverage based)") +
  theme_bw() + 
  scale_x_continuous(limits = c(0,1e+5)) +
  theme(legend.position="none")  

plot(colSums(ASV.bac))





rarefaction_obs <- subset(rarefaction, Method == "Observed")
rarefaction_rar <- subset(rarefaction, Method == "Rarefaction")
rarefaction_ext <- subset(rarefaction, Method == "Extrapolation")

ggplot(rarefac.line, aes(x=x, y=y, colour=depth_category)) +
  geom_line(lwd=1) +
  scale_x_continuous(limits = c(0,1e+6)) +
  labs(x="Sample size", y="Species richness") +
  theme_bw() + 
  theme(legend.position="bottom")+
  scale_colour_manual(values = custom_palette, breaks = custom_order, labels = custom_order)

##
ggplot(rarefac.line, aes(x=x, y=y, colour=ID_running)) +
  geom_line(lwd=0.1) +
  scale_x_continuous(limits = c(0,1e+5)) +
  labs(x="Sample size", y="Species richness") +
  theme_bw() + 
  theme(legend.position="bottom") 
#scale_colour_manual(values = custom_colors, breaks = custom_order, labels = custom_order)

#################################################
ggplot(cover.line, aes(
  x=x, y=y, colour=depth_category))+ 
  geom_line(
    aes(linetype=Method), 
    lwd = 0.5, data=coverage) +
  #  scale_colour_discrete(guide ="none") +
  scale_x_continuous(
    limits = c(0,1e+6)) +
  scale_y_continuous(
    breaks = seq(0.9,1,0.05), 
    limits = c(0.9,1)) +
  #  facet_grid(locus_tag~mooring_full) +
  labs(x="Sample size", y="Sample coverage") +
  theme_bw(base_size = 12) + 
  theme(legend.position="bottom")+ 
  scale_colour_manual(values = custom_palette, breaks = custom_order, labels = custom_order)

ggplot(rarefaction_rar, aes(
  x=x, y=y, colour=depth_category)) +
  geom_line(aes(
    linetype=Method), 
    lwd=0.5, data=rarefaction_rar) +
  # scale_colour_discrete(guide = "none") +
  scale_x_continuous(limits = c(0,1e+5)) +
  #facet_grid(locus_tag~mooring_full) +
  labs(x="Sample size", y="Species richness") +
  theme_bw() + 
  theme(legend.position="bottom")+
  scale_colour_manual(values = custom_palette, breaks = custom_order, labels = custom_order)

ggplot(rarefaction_ext, aes(
  x=x, y=y, colour=depth_category)) +
  geom_line(aes(
    linetype=Method), 
    lwd=0.5, data=rarefaction_ext) +
  # scale_colour_discrete(guide = "none") +
  scale_x_continuous(limits = c(0,1e+5)) +
  #facet_grid(locus_tag~mooring_full) +
  labs(x="Sample size", y="Species richness") +
  theme_bw() + 
  theme(legend.position="bottom")+
  scale_colour_manual(values = custom_palette, breaks = custom_order, labels = custom_order)



###########################################################################################
###  REFORMATTING ###
############################################################################################

richness <- iNEXT.bac$AsyEst[
  iNEXT.bac$AsyEst$Diversity=="Species richness",] %>%
  arrange(Assemblage) %>%
  right_join(ENV.dc, by=c("Assemblage"="CTD #"))

simpson <- iNEXT.bac$AsyEst[
  iNEXT.bac$AsyEst$Diversity=="Simpson diversity",] %>%
  arrange(Assemblage) %>%
  right_join(ENV.dc, by=c("Assemblage"="CTD #"))

div.bac <- data.frame(
  Event_ID = ENV.dc$`Event-ID`,
  Site = ENV.dc$`CTD #`, 
  locus_tag = ENV.dc$locus_tag,
  richness = richness$Observed,
  simpson = simpson$Observed) 

###################################

###########################################################
## COVERAGE ##
ggplot(coverage, aes(
  x=x, y=y))+ 
  geom_line(
    aes(), 
    lwd = 0.5, data=cover.line) +
  scale_colour_discrete(guide ="none") +
  scale_x_continuous(
    limits = c(0,1e+5)) +
  scale_y_continuous(
    breaks = seq(0.9,1,0.05), 
    limits = c(0.9,1)) +
  # facet_grid(~depth_numeric) +
  labs(x="Sample size", y="Sample coverage") +
  theme_bw(base_size = 12) + 
  theme(legend.position="bottom")

ggplot(coverage, aes(
  x=x, y=y, colour=Assemblage))+ 
  geom_line(
    aes(linetype=gear), 
    lwd = 0.5, data=cover.line) +
  scale_colour_discrete(guide ="none") +
  scale_x_continuous(
    limits = c(0,1e+5)) +
  scale_y_continuous(
    breaks = seq(0.9,1,0.05), 
    limits = c(0.9,1)) +
  # facet_grid(~depth_numeric) +
  labs(x="Sample size", y="Sample coverage") +
  theme_bw(base_size = 12) + 
  theme(legend.position="bottom")


ggplot(rarefaction, aes(
  x=x, y=y, colour=sample)) +
  geom_line(aes(
    linetype=gear), 
    lwd=0.5, data=rarefac.line) +
  #scale_colour_discrete(guide = "none") +
  scale_x_continuous(limits = c(0,2e+5)) +
  facet_grid(~ENV.dc$depth_category) +
  #facet_grid(cols = Assemblage )
  labs(x="Sample size", y="Species richness") +
  theme_bw() + 
  theme(legend.position="bottom")



###########################################################################################
###  REFORMATTING ###
############################################################################################

richness <- iNEXT.bac$AsyEst[
  iNEXT.bac$AsyEst$Diversity=="Species richness",] %>%
  arrange(Assemblage) %>%
  right_join(ENV.dc, by=c("Assemblage"="CTD #"))

simpson <- iNEXT.bac$AsyEst[
  iNEXT.bac$AsyEst$Diversity=="Simpson diversity",] %>%
  arrange(Assemblage) %>%
  right_join(ENV.dc, by=c("Assemblage"="CTD #"))

shannon <- iNEXT.bac$AsyEst[
  iNEXT.bac$AsyEst$Diversity == "Shannon diversity",] %>%
  arrange(Assemblage) %>%
  right_join(ENV.dc, by = c("Assemblage" = "CTD #"))

shannon$Observed <- log(shannon$Observed)
evenness <- shannon$Observed / log(richness$Observed)

div.bac <- data.frame(
  depth_category = ENV.dc$depth_category,
  sample = ENV.dc$'CTD #',
  locus_tag = ENV.dc$locus_tag,
  richness = richness$Observed,
  shannon = shannon$Observed,
  simpson = simpson$Observed,
  evenness = evenness
)

###  MERGE ###

##### inv simpson ####
# Merge iNEXT results with metadata
merged_df <- merge(iNEXT.bac$AsyEst, ENV.dc, by.x = "Assemblage", by.y = "CTD #", all.x = TRUE)
# Calculate Inverse Simpson

# Reshape the data to long format
alpha_div_long <- div.bac %>%
  pivot_longer(cols = c(richness, shannon, simpson, evenness), names_to = "IndexType", values_to = "IndexValue")


Div_Indices <- ggplot(alpha_div_long, aes(x = depth_category, y = IndexValue, fill = depth_category)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.6) +
  facet_wrap(~IndexType, scales = "free_y") +
  labs(x = "Depth Category", y = "Diversity Index Value",
       title = "Diversity Indices of Prokaryotic Communities") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = custom_palette) +
  scale_x_discrete(limits = levels(ENV.dc$depth_category))

Div_Indices






#### no outliers in inv simpson ####
# First, identify and remove outliers for InvSimpson
Q1 <- quantile(alpha_div_long$IndexValue[alpha_div_long$IndexType == "InvSimpson"], 0.25)
Q3 <- quantile(alpha_div_long$IndexValue[alpha_div_long$IndexType == "InvSimpson"], 0.75)
IQR <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

alpha_div_long_no_outliers <- alpha_div_long %>%
  filter(!(IndexType == "InvSimpson" & (IndexValue < lower_bound | IndexValue > upper_bound)))

# Then, create the plot
ggplot(alpha_div_long_no_outliers, aes(x = depth_category, y = IndexValue))+#, fill = custom_palette) +
  geom_boxplot(outlier.shape = 32, alpha = 0.6) +
  facet_wrap(~IndexType, scales = "free_y") +
  labs(x = "Depth Category", y = "Diversity Index Value",
       title = "Diversity Indices of Prokaryotic Communities") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = custom_palette) +
  scale_x_discrete(limits = levels(ENV.dc$depth_category)) +
  geom_boxplot(outlier.shape = 32, fill = NA, color = NA, alpha = 1)
# Reshape the data to long format
AlphaDiv <- div.bac
alpha_div_long <- AlphaDiv %>%
  pivot_longer(cols = c(richness, simpson, InvSimpson, shannon),
               names_to = "IndexType",
               values_to = "IndexValue")

# Create the plot
ggplot(alpha_div_long, aes(x = depth_category, y = IndexValue)) +
  geom_boxplot(outlier.shape = 21) +
  labs(x = "Depth Category", y = "Diversity Index Value",
       title = "Diversity Indices of Prokaryotic Communities") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

######
data <-  readr::read_rds("iNext_bac_final_247.RDS")
amp <-  readr::read_rds("PS138_Ampvis-final.RDS")

richness <- data$AsyEst[
  data$AsyEst$Diversity=="Species richness",] %>%
  arrange(Assemblage) 
simpson <- data$AsyEst[
  data$AsyEst$Diversity=="Simpson diversity",] %>%
  arrange(Assemblage) 
shannon <- data$AsyEst[
  data$AsyEst$Diversity=="Shannon diversity",] %>%
  arrange(Assemblage) 

# compile; calculate evenness
shannon$Observed <- log(shannon$Observed)
evenness <- shannon$Observed / log(richness$Observed)

diversity <- data.frame(
  clip_id = richness$Assemblage,
  richness = richness$Observed,
  simpson = simpson$Observed,
  shannon = shannon$Observed,
  evenness = evenness) %>%
  mutate(evenness = shannon/log(richness)) %>%
  left_join(meta, by=c("clip_id"="CTD #"))

ice_colors <- colorRampPalette(c("#00BCD4","#0D47A1","#000033"))(3)
euphotic_colors <- colorRampPalette(c("yellow","orange3","lawngreen","#003300"))(8) # Gradient from green to black
deep_ocean_colors <- colorRampPalette(c("#FF0FF0","#660660","red", "#330300","black"))(8)

custom_palette <- c(ice_colors, euphotic_colors, deep_ocean_colors)


# Define custom order
custom_order <- c(
  "MP","Ice-TS","Ice-BS","UIW","2", "10", "chl-max","25", "50", "100", "200", "500",
  "1000", "1500", "2000", "3000","20m from bottom", "5m from bottom",  "bottom")

# Create the plot, filtering out unwanted depths
diversity %>%
  filter(depth_category %in% custom_order) %>%  # Add this line to filter
  mutate(depth = factor(depth_category, levels = custom_order)) %>%
  reshape2::melt(id.vars = "depth", measure.vars = c("richness", "simpson","evenness")) %>%
  ggplot(aes(x = depth, y = value, fill = depth)) +
  geom_boxplot(alpha=0.7) +
  facet_grid(variable ~ ., scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = custom_palette)

###############################################################
'
# remove temp-data
rm(richness, simpson, 
   coverage, cover.bac, cover.euk, 
   rarefaction, rarefac.bac, rarefac.euk,
   div.bac, div.euk)

rm(list = ls(pattern =
               ".point.*|.line.*"))
'
# Save data
save(iNEXT.bac, file="iNEXT_bac.RData")
###############################

