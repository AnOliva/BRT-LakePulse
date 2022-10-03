# Aim : combine various environmental variable sources and
# normalize the full dataset
# Performed in R V4.1.2

## working directory
setwd("/2_BRT_preprocessing/results")

## libraries
library(readxl)
library(dplyr)
library(data.table) ## v >= 1.9.6

## loading data
# GIS data: air temperature, precipitations...
GIS_combined = read.csv("/combined_GIS_variables.txt", sep=";")
# LULC
LULC = read.csv("/LakePulse_lulc_simple_qc_v2_20220310.csv", sep=";")
# depth
depth = read.csv("/LakePulse_basic_info_depth_qc_20200901.csv", sep=";")
# optics
a_CDOM = read_excel("/a_CDOM_LakePulse.xlsx") 
secchi = read.csv("/LakePulse_secchi_qc_20211223.csv", sep=";") 
TSS_MSPM_OSPM =  read_excel("/QC_TSS_MSPM_OSPM_2017.2019_update20210407.xlsx", sheet = "QC_TSS_MSPM_OSPM_2017.2019")
# altitude, longitude, latitude, lake area, province, ecozone, ecoprovince, depth
basic_info = read.csv("/LakePulse_basic_info_qc_20220310.csv", sep=";")
# RBR: Wtemp, dissolvedO2, chla, salinity
RBR_Wtemp = read.csv("/LakePulse_rbr_top_bottom_temperature_qc_20201110.csv", sep=";")
RBR_salinity = read.csv("/LakePulse_rbr_top_bottom_salinity_qc_20201110.csv", sep=";")
# chlorophylle filtered
chla_am_pm = read.csv("/LakePulse_chla_qc_20201218.csv", sep=";")
# phosphorus
phosphorus = read.csv("/LakePulse_tp_qc_20201116.csv", sep=";")
nitrogen = read.csv("/LakePulse_tn_qc_20201028.csv", sep=";")
# hydrolakes
hydrolake = read.csv("/LakePulse_hydrolake_20200901.csv", sep=";")
# ecoumene
ecoumene_lakes = read.csv("/ecumenes_lake.csv")
# NOX and SRP
NOX = read.csv("/LakePulse_nox_qc_20210517.csv", sep=";")
SRP = read.csv("/LakePulse_srp_qc_20201124.csv", sep=";")
# geology
#soil = read.csv("E/LakePulse_soil_landscape_20220310.csv", sep=";")
#drainage = readxl::read_excel("/calculateSoilLandscapesOfCanadaDRAIN_final20220407.xlsx")
# slope
slope_0_100_m = read.csv("/slope_0_100_m_LP.csv", sep=";")


### normalize and QC the flags

## CDOM 
a_CDOM = t(a_CDOM)
colnames(a_CDOM) = a_CDOM[1,]
a_CDOM = a_CDOM[-1,]
# Extract absorption at specific wavelength
# suggested a 300: Clark, Catherine D., et al. "A Study of Fecal Coliform Sources at a Coastal Site Using Colored Dissolved Organic Matter (CDOM) as a Water Source Tracer." Marine Pollution Bulletin, vol. 54, no. 9, 2007, pp. 1507-13, doi:10.1016/j.marpolbul.2007.04.019.
# 280 and 355 suggested by: Madonia et al, Chromophoric Dissolved Organic Matter as a Tracerof Fecal Contamination for Bathing Water QualityMonitoring in the Northern Tyrrhenian Sea(Latium, Italy)
# remote sensing satellite estimates around 400 nm
a_CDOM = a_CDOM[,c("280", "355", "400")] 
# transform to dataframe
a_CDOM = as.data.frame(a_CDOM)
a_CDOM$lakepulse_id = rownames(a_CDOM)
# rename column
colnames(a_CDOM) = c("aCDOM280", "aCDOM355", "aCDOM400", "lakepulse_id")
# transform non-numeric values to numeric values
a_CDOM[,c("aCDOM280","aCDOM355", "aCDOM400")] = sapply(a_CDOM[,c("aCDOM280","aCDOM355", "aCDOM400")],as.numeric)

## RBR
# keep data below the surface in Epilimnion as eDNA was sampled
RBR_salinity = RBR_salinity[RBR_salinity$depth == "Tube length", ] 
RBR_Wtemp = RBR_Wtemp[RBR_Wtemp$depth == "Tube length", ] 
# flag QC
# temperature flags: 67 (estimated value), 32 (potentially incorrect), 26 (sample not collected)
# remove -1 flags + transform 32 values into NA because we are not sure of the values and keep the rest
# Wtemp flags
RBR_Wtemp$Wtemp =  ifelse(grepl("-1", RBR_Wtemp$temp_flags, fixed = TRUE), NA, RBR_Wtemp$temp_mean)
# Salinity flags
RBR_salinity$salinity =  ifelse(grepl("-1", RBR_salinity$salinity_flags, fixed = TRUE), NA, RBR_salinity$salinity_mean)

## secchi depth
# spread the table 
secchi.spread = dcast(setDT(secchi), lakepulse_id ~ site_name, value.var = c("flags", "mean_m", "secchi_comment")) 
# some lakes have one value either index or optic sites.
# Choose in priority Index (morning sampling with eDNA).
secchi.spread$secchi = secchi.spread$`mean_m_Index Site`
secchi.spread$flags_secchi = secchi.spread$`flags_Index Site`
# if there is no Index Site secchi, take the optic site.
secchi.spread$flags_secchi[is.na(secchi.spread$secchi)] = secchi.spread$`flags_Optics Site`[is.na(secchi.spread$secchi)]
secchi.spread$secchi[is.na(secchi.spread$secchi)] = as.numeric(secchi.spread$`mean_m_Optics Site`[is.na(secchi.spread$secchi)]) 
# add max depth found during sampling either at index or coring site
secchi.spread = merge(secchi.spread, depth[,c("lakepulse_id", "max_depth_found")], by="lakepulse_id") 
# flags
# 15 Difference between replicates or secchi (disappearance, reappearance) is <= 10%
# 16 Difference between replicates or secchi (disappearance, reappearance) is > 10% and <=20 %
# 17 Difference between replicates or secchi (disappearance, reappearance) is between > 20% and <=50 %
# 25 Difference between depth of disappearance and reappearance less than 10 cm
# 26 Sample or data not collected
# 27 Secchi depth greater than lake depth (secchi was visible at the bottom of the lake) => could be used in classes (low/medium/high clarity or maybe in ratio secchi mean /depth?)
# 40 Attention: please read comments
# calculate ratio secchi/depth
# create new column with 27 flag = depth, calculate ratio to the depth of the lake
# if secchi depth is greater than lake depth, then write lake depth otherwise keep secchi depth.
secchi.spread$secchi_adjusted = ifelse(secchi.spread$flags_secchi == "27",
                                       secchi.spread$max_depth_found, secchi.spread$secchi)
# calculate the ratio secchi/depth
secchi.spread$Zsecchi_ratio = secchi.spread$secchi_adjusted / secchi.spread$max_depth_found

## geographical and sampling method data
basic_info = basic_info[,c("lakepulse_id","latitude", "longitude", 
                           "province", "team", "ecozone", "sampling_date", "area")]
# convert to data
basic_info$sampling_date = as.Date(basic_info$sampling_date)
# extract the year of sampling
basic_info$year = format(as.Date(basic_info$sampling_date, format="%d-%m-%Y"),"%Y")
# create 0-1 columns for each ecozone
basic_info$amount = 1
basic_info$ecozone2 = basic_info$ecozone
basic_info = tidyr::spread(basic_info, key=ecozone2, amount, fill = 0)
# as not many lakes were sampled in the North, combine the ecozones
basic_info$NorthEcoz = rowSums(basic_info[,c("Boreal Cordillera" , "Taiga Cordillera","Taiga Plains")], na.rm=TRUE)
# merge Montane Cordillera and Semi Arid Plateau
basic_info$WestMont = rowSums(basic_info[,c("Montane Cordillera" , "Semi-Arid Plateaux")], na.rm=TRUE)
# ecozone grouped
basic_info$ecozone = ifelse(basic_info$ecozone == "Montane Cordillera" | basic_info$ecozone == "Semi-Arid Plateaux", "WestMont", basic_info$ecozone)
basic_info$ecozone = as.factor(basic_info$ecozone)

# reorder the columns
colnames(basic_info) = c("lakepulse_id","latitude","longitude","province","team","ecozone","sampling_date",
                         "Alake",
                         "year","AtHigh", "AtMar","BorCor","BorPla","BorShi",
                         "MixPla","MonCor", "PaciMar", "Prairies",
                         "SemiAPl", "TaigaCor", "TaigaPla","NorthEcoz", "WestMont")


## TSS 
# -1 (not good data),
# 19 (Only one replicate when two or more expected),
# 24 (Blank subtraction led to a negative value (the value may have been set to 0)),
# 26,
# 80 (Sample fragmented and recovered but may have lost some pieces or particles) 
# transform all the -1 data into NA only
TSS_MSPM_OSPM$TSS = ifelse(TSS_MSPM_OSPM$flag_tss == "-1", NA, TSS_MSPM_OSPM$mean_tss)
## MSPM: -1 (not good data),
# 19 (Only one replicate when two or more expected),
# 24 (Blank subtraction led to a negative value (the value may have been set to 0)),
# 26,
# 77 (Estimated value: Filtrated water quantity incorrect or uncertain / Filtrated water estimated from another sample)
# 80 (Sample fragmented and recovered but may have lost some pieces or particles) 
# transform all the -1 data into NA only
TSS_MSPM_OSPM$MSPM = ifelse(TSS_MSPM_OSPM$flag_mspm == "-1", NA, TSS_MSPM_OSPM$mean_mspm)
## OSPM: -1 (not good data),
# 19 (Only one replicate when two or more expected),
# 24 (Blank subtraction led to a negative value (the value may have been set to 0)),
# 26,
# 80 (Sample fragmented and recovered but may have lost some pieces or particles) 
# transform all the -1 data into NA only
TSS_MSPM_OSPM$OSPM = ifelse(TSS_MSPM_OSPM$flag_ospm == "-1", NA, TSS_MSPM_OSPM$mean_ospm)
## transform non-numeric values to numeric values
TSS_MSPM_OSPM[,c("TSS","MSPM", "OSPM")] = sapply(TSS_MSPM_OSPM[,c("TSS","MSPM", "OSPM")],as.numeric)

## Chla AM
# select chlorophyll at sampling time of eDNA (AM)
chla_am_pm = chla_am_pm[chla_am_pm$sampling_period == "AM",]
# remove -1 flags
chla_am_pm$chla =  ifelse(grepl("-1", chla_am_pm$chla_flags, fixed = TRUE), NA, chla_am_pm$chla_concentration)

## population 
# add watershed area
pop_density = merge(GIS_combined[,c("population_samplingYear"), drop = F], LULC[,c("lakepulse_id", "total_km2_area")],
                    by.x = 0, by.y = "lakepulse_id")
# calculate population density within the watershed
pop_density$pop_wshed = pop_density$population_samplingYear / pop_density$total_km2_area
# rename columns
colnames(pop_density) = c("lakepulse_id","npop","total_km2_area","pop_wshed")

## phosphorus
# select the epilimnion where eDNA was sampled
phosphorus = phosphorus[phosphorus$depth == "Epilimnion",]
# remove -1 flags
phosphorus$TP = ifelse(grepl("-1", phosphorus$tp_flags, fixed = TRUE), NA, phosphorus$tp_concentration)

## nitrogen
# select epilimnion where eDNA was sampled
nitrogen = nitrogen[nitrogen$depth == "Epilimnion",]
# remove -1 flags
nitrogen$TN = ifelse(grepl("-1", nitrogen$tn_flags, fixed = TRUE), NA, nitrogen$tn_concentration)

## NOX
# remove flags -1
NOX = NOX[!grepl("-1", NOX$no3_flags),]
# keep epilimnion data
NOX = NOX[!grepl("Hypolimnion", NOX$depth),]

## SRP
# remove flags -1
SRP = SRP[!grepl("-1", SRP$srp_flags),]
# keep epilimnion data
SRP = SRP[!grepl("Hypolimnion", SRP$depth),]

## SRP and NoX into No3
SRP_NO3 = merge(SRP[,c("lakepulse_id", "srp_concentration")], NOX[,c("lakepulse_id", "no3_concentration")], by = "lakepulse_id")
SRP_NO3$SRPNo3 = SRP_NO3$srp_concentration / SRP_NO3$no3_concentration

## LULC
# potential anthropic sources of contamination
# fraction
LULC$fraction_anthro = LULC$fraction_agriculture + LULC$fraction_pasture + LULC$fraction_urban
# area
LULC$area_km2_anthro = LULC$area_km2_agriculture + LULC$area_km2_pasture + LULC$area_km2_urban
# potential agricultural and livestock sources of contaminatiob
# fraction
LULC$fraction_agripas = LULC$fraction_agriculture + LULC$fraction_pasture
# area
LULC$area_km2_agripas = LULC$area_km2_agriculture + LULC$area_km2_pasture
# ratio Agripas and Urban by Awatershed
LULC$agripasNat = (LULC$area_km2_agripas) / (LULC$area_km2_nat_landscapes)
LULC$urbanNat = (LULC$area_km2_urban) / (LULC$area_km2_nat_landscapes)

## hydrolake
# remove flags -1
hydrolake = hydrolake[!grepl("-1", hydrolake$flags),]
# elevation is -2, restime is -1, change to 0
hydrolake[hydrolake[,c("elevation")] < 0 & ! is.na(hydrolake[,c("elevation")]), ]$elevation = 0
hydrolake[hydrolake[,c("res_time")] < 0 & ! is.na(hydrolake[,c("res_time")]), ]$res_time = 0

## GIS
# rename columns
GIS_combined = GIS_combined[,c("TE.allMonthsAvg","TP.allMonthsAvg",
                               "Temp.allMonthsAvg", "daysWithIce","fireNbr_since2012",
                               "manure_wshed","animal_unit","animal_unit_density"), drop = F]

colnames(GIS_combined) = c("evaporation","precipitation",
                           "Tair", "daysWithIce","fireNbr",
                           "manure_wshed","animal_unit","animal_unit_density")
# transform NA into 0
GIS_combined$manure_wshed = ifelse(is.na(GIS_combined$manure_wshed), 0, GIS_combined$manure_wshed)
GIS_combined$animal_unit = ifelse(is.na(GIS_combined$animal_unit), 0, GIS_combined$animal_unit)
GIS_combined$animal_unit_density = ifelse(is.na(GIS_combined$animal_unit_density), 0, GIS_combined$animal_unit_density)

## slope
colnames(slope_0_100_m) = c("idLakePulse","slope100m")

### combine variables
# arrange in same order
basic_info = arrange(basic_info, lakepulse_id)
GIS_combined = arrange(GIS_combined, rownames(GIS_combined))
chla_am_pm = arrange(chla_am_pm, lakepulse_id)
depth = arrange(depth, lakepulse_id)
LULC = arrange(LULC, lakepulse_id)
RBR_Wtemp = arrange(RBR_Wtemp, lakepulse_id)
RBR_salinity = arrange(RBR_salinity, lakepulse_id)
secchi.spread = arrange(secchi.spread, lakepulse_id)
phosphorus = arrange(phosphorus, lakepulse_id)
nitrogen = arrange(nitrogen, lakepulse_id)
pop_density = arrange(pop_density, lakepulse_id)
hydrolake = arrange(hydrolake, lakepulse_id)
ecoumene_lakes = arrange(ecoumene_lakes, lakepulse)
SRP_NO3 = arrange(SRP_NO3, lakepulse_id)
slope = arrange(slope_0_100_m, idLakePulse)

# cbind
table_combined = cbind(basic_info[,c("lakepulse_id","latitude","longitude","province","team", "ecozone",
                                     "AtHigh", "AtMar","BorCor","BorPla","BorShi",
                                     "MixPla","MonCor", "PaciMar", "Prairies",
                                     "SemiAPl", "TaigaCor", "TaigaPla","NorthEcoz","WestMont",
                                     "year", "Alake")], 
                       
                       LULC[,c("area_km2_forestry", "area_km2_agriculture", "area_km2_nat_landscapes",
                               "area_km2_pasture", "area_km2_urban", "area_km2_grassland",
                               "area_km2_anthro", "area_km2_agripas",
                               "fraction_agriculture", "fraction_nat_landscapes",
                               "fraction_forestry", "fraction_urban", 
                               "fraction_pasture", "fraction_grassland",
                               "fraction_anthro", "fraction_agripas",
                               "agripasNat", "urbanNat",
                               "total_km2_area")],
                       
                       GIS_combined,
                       
                       ecoumene_lakes[,c("ecumene_agr", "ecumene_pop")],
                       
                       secchi.spread[,c("secchi_adjusted"), drop = F],
                       
                       chla_am_pm[,c("chla"), drop = F],
                       
                       phosphorus[,c("TP"), drop = F],
                       
                       pop_density[,c("npop","pop_wshed")],
                       
                       slope[,c("slope100m"), drop = F],
                       
                       RBR_Wtemp[,c("Wtemp"), drop = F],
                       
                       RBR_salinity[,c("salinity"), drop = F]
)
# add TSPM
table_combined = full_join(table_combined, TSS_MSPM_OSPM[,c("lakeid","TSS")], by = c("lakepulse_id" ="lakeid"))
# add CDOM
table_combined = full_join(table_combined, a_CDOM, by ="lakepulse_id")
# add hydrolake
table_combined = full_join(table_combined, hydrolake[,c("lakepulse_id", #"lake_area",
                                                        "depth_avg"#,
                                                        #"dis_avg","res_time","elevation",
                                                        #"slope_100", "wshd_area", "vol_total"
)], by ="lakepulse_id")
# add SRPNO3
table_combined = full_join(table_combined, SRP_NO3[,c("lakepulse_id", "SRPNo3")], by ="lakepulse_id")
# add nitrogen
table_combined = full_join(table_combined, nitrogen[,c("lakepulse_id","TN")], by ="lakepulse_id")



### Transformation of variables

# inverse hyperbolic since function:
# https://stackoverflow.com/questions/33612312/how-to-reverse-inverse-hyperbolic-sine-transformation-in-r
ihs = function(x) {
  y = log(x + sqrt(x^2 + 1))
  return(y)
}
# transformation without min() #"vol_total","slope_100", "wshd_area",
table_combined[,c("secchi_adjusted", "aCDOM280","aCDOM355", "aCDOM400",
                  "depth_avg", "Alake", "TP", "TN","total_km2_area",
                  "chla")] = log10(table_combined[,c("secchi_adjusted", "aCDOM280","aCDOM355", "aCDOM400",
                                                     "depth_avg", "Alake", "TP", "TN","total_km2_area",
                                                     "chla")])
# transformation with min() "res_time","dis_avg", "wind1day", "wind7day", "wind30days"
table_combined[,c("pop_wshed", "animal_unit", "area_km2_anthro",
                  "manure_wshed", "npop", "area_km2_nat_landscapes", "area_km2_pasture",
                  "area_km2_agriculture","area_km2_grassland",
                  "area_km2_urban", "area_km2_forestry", "area_km2_agripas", 
                  "TSS", "urbanNat", "agripasNat"
)] = ihs(table_combined[,c("pop_wshed", "animal_unit", "area_km2_anthro",
                           "manure_wshed", "npop", "area_km2_nat_landscapes", "area_km2_pasture",
                           "area_km2_agriculture","area_km2_grassland",
                           "area_km2_urban", "area_km2_forestry", "area_km2_agripas", 
                           "TSS", "urbanNat", "agripasNat")])

# transformation sqrt"TP.SummerMonthsAvg", "TP.allMonthsAvg","elevation",#"precip_sum7day","precip_sum30day",
table_combined[,c("fraction_anthro", "fraction_agriculture","animal_unit_density",
                  "fraction_pasture","fraction_grassland", "fraction_nat_landscapes",
                  "fraction_urban", "fraction_forestry", "fraction_agripas",
                  "slope100m")] = sqrt(table_combined[,c("fraction_anthro", "fraction_agriculture","animal_unit_density",
                                                           "fraction_pasture","fraction_grassland", "fraction_nat_landscapes",
                                                           "fraction_urban", "fraction_forestry", "fraction_agripas",
                                                           "slope100m")])

## remove column that have occurrences < 5% of the dataset
table_combined = table_combined[, colSums(table_combined != 0, na.rm = T) > (nrow(table_combined)*0.05)]

## remove categorical variables present in more than 90% or less than 5% of the dataset
cat.var = c("AtHigh","AtMar","BorCor","BorPla","BorShi",
            "MixPla","MonCor","PaciMar","Prairies",
            "SemiAPl","TaigaCor","TaigaPla","NorthEcoz",
            "WestMont","ecumene_agr","ecumene_pop")
# select the variables
var.to.delete = c(colnames(table_combined[ , names(table_combined) %in% cat.var][, colSums(table_combined[ , names(table_combined) %in% cat.var] != 0, na.rm = T) < (nrow(table_combined)*0.05), drop = F]),
                  colnames(table_combined[ , names(table_combined) %in% cat.var][, colSums(table_combined[ , names(table_combined) %in% cat.var] != 0, na.rm = T) > (nrow(table_combined)*0.95), drop = F]))
# remove
table_combined = table_combined[ , ! names(table_combined) %in% var.to.delete]

## NO EFFECT ON BRT + scaling continuous variables except 0-1 variables (SHOULD BE PERFORMED ON SELECTED OBS, NOT 664 lakes)
#table_combined[ , ! names(table_combined) %in% cat.var] = scale(table_combined[ , ! names(table_combined) %in% cat.var])

## load upscaled dataset to add idLatLong column
lakeid <- read.csv("/QCed_lakepulse_id_vs_upscaling_idLatLong.csv", sep=";")
table_combined = merge(table_combined, lakeid[,c("idLatLong_upscaling", "lakepulse_id")], by = "lakepulse_id", all.x = T)
# rename the idLatLong column
colnames(table_combined)[colnames(table_combined) == 'idLatLong_upscaling'] <- 'idLatLong'

## export
write.table(table_combined, "variables_combined.txt", row.names = FALSE, sep= ";")



