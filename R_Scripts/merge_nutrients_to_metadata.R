#merge metadata to CTD_Nutrient_Data and ICE_Nutrient_Data
setwd("~/PS138")
library(dplyr)
library(readr)
library(writexl)
library(data.table)
library(readxl)
library(readxl)
library(readr)



Ice_Nutrients<- read_excel("PS138_Lena_nutrients.xlsx")
CTD_Nutrients<- read_csv("Nutrients_Sinhue_fixed(in).csv")

PS138_metadata_16S<- read_excel("PS138_metadata_16S_CTD_ICE_ISCA_130524.xlsx")

#################

Nutrients_Merged <- bind_rows(CTD_Nutrients, Ice_Nutrients)
Metadata_Nutrients <- merge(PS138_metadata_16S, Nutrients_Merged, by = "Niskin_Id", all.x = TRUE)



# Write the merged data frame to a new Excel file
library(openxlsx)
write.xlsx(Metadata_Nutrients, "PS138_metadata_16S_all_nutrients_fixed.xlsx", rowNames = FALSE)
