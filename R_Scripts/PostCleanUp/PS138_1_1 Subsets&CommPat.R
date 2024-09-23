setwd("~/PS138")
load("~/PS138/PS138_16S_final_species_fixed.Rdata")
library(ampvis2)
library("ape")
library("phyloseq")
library(ggplot2)
library(dplyr)
library(DECIPHER)
library(phangorn)
library(dada2)
library(ggtree)
library(phyloseq)
library(dendextend)
library(phyloseq)
library(SRS)
library(tidyr)
library(microbiome)
library(scico)


colnames(ASV.bac)
colSums(ASV.bac[,1:249])
plot(colSums(ASV.bac[,1:249]),pch=19,col='blue') 

# load sample information and color accordingly:
sam <- ENV
head(sam) # sample names differ from column names in asv data (- instead of .), change:


# add the total read counts to the table:
tc <- as.data.frame(cbind(colnames(ASV.bac[,1:249]),colSums(ASV.bac[,1:249])))
head(tc)
tc$V2 <- as.numeric(tc$V2)
sam$tc <- tc$V2[match(ENV$'CTD #',tc$V1)]
head(sam)
str(sam)

sam.n <- sam[sam$tc>500,]

ASV.n <- ASV.bac[,colnames(ASV.bac) %in% sam.n$'CTD #']
Cmin <- min(colSums(ASV.n))
#Cmin :  11276
SRS_ASV.n <- SRS(ASV.n,Cmin)

#ASV.n unchanged counts
#SRS_ASV.n normelized counts
plot(colSums(ASV.n),pch=19,col='blue')
plot(colSums(SRS_ASV.n),pch=19,col='blue')


depth_classes <- c(2, 10, 25, 50, 75, 100, 200, 500, 1000, 1500, 2000, 3000, 4000, 4500)

# Convert 'depth [m]' to numeric, replacing non-numeric values with NA
sam.n$numeric_depth <- as.numeric(as.character(sam.n$`depth [m]`), na.rm = TRUE)

# Handle NAs by assigning them to a category (adjust as needed)
sam.n$numeric_depth_category <- ifelse(is.na(sam.n$numeric_depth), NA,
                                       depth_classes[findInterval(sam.n$numeric_depth,
                                                                  depth_classes, left.open = TRUE) + 1])
sam.n$numeric_depth_category[!is.na(sam.n$layer)] <- sam.n$layer[!is.na(sam.n$layer)]
sam.n$depth_category <- sam.n$numeric_depth_category

# Print the result
print(sam.n$depth_category) 
#alignment <- AlignSeqs(DNAStringSet(seqs_tree), anchor=NA,verbose=FALSE)
rows_to_change <- sam.n$gear == "OFOBS"
sam.n$depth_category[rows_to_change] <- "bottom"

sam.n$depth_category <- gsub("bottom-5", "5m from bottom", sam.n$depth_category)
sam.n$depth_category <- gsub("meltpond", "MP", sam.n$depth_category)


# Check if the value in numeric_depth_category is "vent-plume"
index <- sam.n$numeric_depth_category == "vent-plume"

# Change the corresponding value in depth_category to "vent-plume"
sam.n$depth_category[index] <- "vent-plume"


sam1 <- sam.n[, c("CTD #", names(sam.n)[names(sam.n) != "CTD #"])]


ampvis.bac_SRS <- data.frame(
  OTU = rownames(SRS_ASV.n),
  SRS_ASV.n, TAX.bac, check.names=F)
ampvis.bac_SRS <- amp_load(ampvis.bac_SRS, sam1)

ampvis.bac_hel <- data.frame(
  OTU = rownames(ASV.bac.hel),
  ASV.bac.hel, TAX.bac, check.names=F)
ampvis.bac_hel <- amp_load(ampvis.bac_hel, sam1)

ampvis.bac <- data.frame(
  OTU = rownames(ASV.n),
  ASV.n, TAX.bac, check.names=F)
ampvis.bac <- amp_load(ampvis.bac, sam1)

depth_c <-c(
  "MP","Ice-TS","Ice-BS","ice-edge","UIW","2", "10", "chl-max","25", "50", "75", "100", "200", "500",
  "1000", "1500", "2000", "3000", "4000","20m from bottom", "5m from bottom",  "bottom","vent-plume")
####################

library(dplyr)



#####################################################
print(unique(sam1$depth_category))

ampvis.bac$metadata$depth_category <- factor(
  ampvis.bac$metadata$depth_category,
  levels = c(
    "MP","Ice-TS","Ice-BS","ice-edge","UIW","2", "10", "chl-max","25", "50", "75", "100", "200", "500",
    "1000", "1500", "2000", "3000", "4000","20m from bottom", "5m from bottom",  "bottom","vent-plume"))

ampvis.bac_hel$metadata$depth_category <- factor(
  ampvis.bac_hel$metadata$depth_category,
  levels = c(
    "MP","Ice-TS","Ice-BS","ice-edge","UIW","2", "10", "chl-max","25", "50", "75", "100", "200", "500",
    "1000", "1500", "2000", "3000", "4000","20m from bottom", "5m from bottom",  "bottom","vent-plume"))

ampvis.bac_SRS$metadata$depth_category <- factor(
  ampvis.bac_SRS$metadata$depth_category,
  levels = c(
    "MP","Ice-TS","Ice-BS","ice-edge","UIW","2", "10", "chl-max","25", "50", "75", "100", "200", "500",
    "1000", "1500", "2000", "3000", "4000","20m from bottom", "5m from bottom",  "bottom","vent-plume"))


#save.image("ampvis_bac.Rdata")
######################################################################
################### plots subsets etc  ###############################

ampvis.bac <- amp_filter_samples(ampvis.bac,
                                 depth_category %in% c("MP","Ice-TS","Ice-BS","UIW","2", "10","25","chl-max","50","100", "200", "500", "1000", "1500", "2000", "3000","20m from bottom","5m from bottom","bottom"),
                                 normalise = FALSE  )  

ampvis.bac_hel <- amp_filter_samples(ampvis.bac_hel,
                                     depth_category %in% c("MP","Ice-TS","Ice-BS","UIW","2", "10","25","chl-max","50","100", "200", "500", "1000", "1500", "2000", "3000","bottom","20m from bottom","5m from bottom"),
                                     normalise = FALSE  )  

ampvis.bac_SRS <- amp_filter_samples(ampvis.bac_SRS,
                                     depth_category %in% c("MP","Ice-TS","Ice-BS","UIW","2", "10","25","chl-max","50","100", "200", "500", "1000", "1500", "2000", "3000","bottom","20m from bottom","5m from bottom"),
                                     normalise = FALSE  )  

#saveRDS(ampvis.bac,"ampvis_bac_247_PS138.RDS")
###############################################
ampvis.bac_nofobs <- amp_filter_samples(ampvis.bac,
                                        gear %in% c("CTD","si_corer_9cm","hand_pump"),
                                        normalise = FALSE )
ampvis_hel_nofobs <- amp_filter_samples(ampvis.bac_hel,
                                        gear %in% c("CTD","si_corer_9cm","hand_pump"),
                                        normalise = FALSE )
ampvis_SRS_nofobs <- amp_filter_samples(ampvis.bac_SRS,
                                        gear %in% c("CTD","si_corer_9cm","hand_pump"),
                                        normalise = FALSE )
###############################################
ice.amp <- amp_filter_samples(ampvis.bac_nofobs,
                              depth_category %in% c("MP","Ice-TS","Ice-BS","UIW"),
                              normalise = FALSE)


photo.amp <- amp_filter_samples(ampvis.bac_nofobs,
                                depth_category %in% c("MP","Ice-TS","Ice-BS","UIW","2", "10","25","chl-max","50"),
                                normalise = FALSE)

surface.amp <- amp_filter_samples(ampvis.bac_nofobs,
                                  depth_category %in% c("2", "10","25","'chl-max'"),
                                  normalise = FALSE)

pelagic.amp <- amp_filter_samples(ampvis.bac_nofobs,
                                  depth_category %in% c("50","100", "200", "500", "1000", "1500", "2000", "3000", "4000"),
                                  normalise = FALSE  )                               

CTD.amp <- amp_filter_samples(ampvis.bac_nofobs,
                              depth_category %in% c("2", "10","25","chl-max","50","100", "200", "500", "1000", "1500", "2000", "3000", "4000","bottom","20m from bottom","5m from bottom"),
                              normalise = FALSE  )                               

bottom.amp <- amp_filter_samples(ampvis.bac_nofobs,
                                 depth_category %in% c("bottom","'20m from bottom'","'5m from bottom'"),
                                 normalise = FALSE )                                


ampvis.bac_vp_ie <- amp_filter_samples(ampvis.bac,
                                       depth_category %in% c(  "MP","Ice-TS","Ice-BS","UIW","2", "10", "chl-max","25", "50", "100", "200", "500",
                                                               "1000", "1500", "2000", "3000","20m from bottom", "5m from bottom","bottom"),
                                       normalise = FALSE  )  

ampvis.bac<- ampvis.bac_vp_ie
ampvis.bac_nofobs <- amp_filter_samples(ampvis.bac,
                                        gear %in% c("CTD","si_corer_9cm","hand_pump"),
                                        normalise = FALSE )

ampvis.bac_CTD <- amp_filter_samples(ampvis.bac_nofobs,
                                     gear %in% c("CTD"),
                                     normalise = FALSE )

ampvis.bac_ICE <- amp_filter_samples(ampvis.bac_nofobs,
                                     gear %in% c("si_corer_9cm","hand_pump"),
                                     normalise = FALSE )

ampvis_depth <- amp_filter_samples(ampvis_hel_nofobs,
                                   depth_category %in% c( "4000","20m from bottom", "5m from bottom","bottom"),
                                   normalise = FALSE  )  

station_amp <- amp_subset_samples(ampvis.bac, 
                                  grepl("^ice station [0-9]+$", sample, ignore.case = TRUE))


transect_amp_1 <- amp_subset_samples(ampvis.bac, 
                                     grepl("Transect", sample, ignore.case = TRUE))

# Subset the original ampvis object to get the samples you want to add
transect_amp_add <- amp_subset_samples(ampvis.bac, 
                                       grepl("^ice station [7-9]+$", sample, ignore.case = TRUE))

existing_samples <- (transect_amp_1$metadata$`CTD #`)
new_samples <- (transect_amp_add$metadata$`CTD #`)

# Combine the sample names
all_transect_samples <- c(existing_samples, new_samples)

# Use amp_subset to create a new object with all samples
transect_amp <- amp_filter_samples(ampvis.bac_nofobs,`CTD #`%in% all_transect_samples)

#saveRDS(ampvis.bac,"PS138_Ampvis-final.rds")
save(ampvis.bac, photo.amp,transect_amp,station_amp,ampvis_depth,ampvis.bac_nofobs,ampvis_SRS_nofobs,ampvis_hel_nofobs,ice.amp, surface.amp, pelagic.amp, ampvis.bac_vp, ampvis.bac_nofobs, file = "ampvis_subsets.Rdata")