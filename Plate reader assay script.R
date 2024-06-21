#Load required packages
library(readxl)
library(tidyverse)
#Import datasets

#Names of samples in wells (necessary for later)
##Experiment 27/02/2024
sample_name_mapping_27_02_2024 <- c("A1" = "WT",
                         "C1" = "WT",
                         "E1" = "WT",
                         "G1" = "WT", 
                         "A3" = "clone 1",
                         "C3" = "clone 1",
                         "E3" = "clone 1",
                         "G3" = "clone 1",
                         "A5" = "clone 2",
                         "C5" = "clone 2",
                         "E5" = "clone 2",
                         "G5" = "clone 2",
                         "A7" = "clone 3",
                         "C7" = "clone 3",
                         "E7" = "clone 3",
                         "G7" = "clone 3",
                         "A9" = "clone 4",
                         "C9" = "clone 4",
                         "E9" = "clone 4",
                         "G9" = "clone 4",
                         "A11" = "clone 5",
                         "C11" = "clone 5",
                         "E11" = "clone 5",
                         "G11" = "clone 5",
                         "A12" = "- CTRL",
                         "C12" = "- CTRL",
                         "E12" = "- CTRL",
                         "G12" = "- CTRL") #pass all the names based on the wells that have been used
## Assign specific values to samples
final_plate_reader_merged$Medium[final_plate_reader_merged$Sample %in% c(paste0("A",  seq(from = 1, to = 12, by = 2)), paste0("C",  seq(from = 1, to = 12, by = 2)), "A12", "C12")] <- "MM" #add more variables as done here if required
final_plate_reader_merged$Medium[final_plate_reader_merged$Sample %in% c(paste0("E",  seq(from = 1, to = 12, by = 2)), paste0("G",  seq(from = 1, to = 12, by = 2)), "E12", "G12")] <- "LB"
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% c(paste0("A",  seq(from = 1, to = 12, by = 2)), paste0("E",  seq(from = 1, to = 12, by = 2)), "A12", "E12")] <- "0 mM"
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% c(paste0("C",  seq(from = 1, to = 12, by = 2)), paste0("G",  seq(from = 1, to = 12, by = 2)), "C12", "G12")] <- "1 mM"

#-------------------------------------------------------

##Experiment 07/03/2024
sample_name_mapping_07_03_2024 <- c("A1" = "clone 1",
                         "B1" = "clone 8",
                         "D1" = "clone 1",
                         "E1" = "clone 8", 
                         "A2" = "clone 2",
                         "B2" = "clone 9",
                         "D2" = "clone 2",
                         "E2" = "clone 9",
                         "A3" = "clone 3",
                         "B3" = "clone 10",
                         "D3" = "clone 3",
                         "E3" = "clone 10",
                         "A4" = "clone 4",
                         "B4" = "clone 11",
                         "D4" = "clone 4",
                         "E4" = "clone 11",
                         "A5" = "clone 5",
                         "B5" = "clone 12",
                         "D5" = "clone 5",
                         "E5" = "clone 12",
                         "A6" = "clone 6",
                         "B6" = "clone 13",
                         "D6" = "clone 6",
                         "E6" = "clone 13",
                         "A7" = "clone 7",
                         "B7" = "WT",
                         "D7" = "clone 7",
                         "E7" = "WT") #pass all the names based on the wells that have been used

## Assign specific values to samples
final_plate_reader_merged$Medium[final_plate_reader_merged$Sample %in% c(paste0("A",  seq(from = 1, to = 7, by = 2)), paste0("B",  seq(from = 1, to = 7, by = 2)))] <- "MM" #add more variables as done here if required
final_plate_reader_merged$Medium[final_plate_reader_merged$Sample %in% c(paste0("D",  seq(from = 1, to = 7, by = 2)), paste0("E",  seq(from = 1, to = 7, by = 2)))] <- "MM"
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% c(paste0("A",  seq(from = 1, to = 7, by = 2)), paste0("B",  seq(from = 1, to = 7, by = 2)))] <- "0 mM"
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% c(paste0("D",  seq(from = 1, to = 7, by = 2)), paste0("E",  seq(from = 1, to = 7, by = 2)))] <- "1 mM"

#-------------------------------------------------------

##Experiment 24/04/2024
sample_name_mapping_24_04_2024 <- c("A1" = "BAA002",
                                    "A2" = "BAA013",
                                    "A3" = "BAA014",
                                    "A4" = "BAA015", 
                                    "A5" = "BAA016",
                                    "A6" = "BAA020",
                                    "A7" = "BAA033",
                                    "A8" = "BAA034",
                                    "A9" = "WT",
                                    "A10" = "CLONE 2",
                                    "A11" = "CLONE 13",
                                    "A12" = "- CTRL",
                                    "C1" = "BAA002",
                                    "C2" = "BAA013",
                                    "C3" = "BAA014",
                                    "C4" = "BAA015", 
                                    "C5" = "BAA016",
                                    "C6" = "BAA020",
                                    "C7" = "BAA033",
                                    "C8" = "BAA034",
                                    "C9" = "WT",
                                    "C10" = "CLONE 2",
                                    "C11" = "CLONE 13",
                                    "C12" = "- CTRL") #pass all the names based on the wells that have been used

## Assign specific values to samples
final_plate_reader_merged$Medium[final_plate_reader_merged$Sample %in% c(paste0("A",  seq(from = 1, to = 12)), paste0("C",  seq(from = 1, to = 12)))] <- "MM" #add more variables as done here if required
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% paste0("A",  seq(from = 1, to = 12))] <- "0 mM"
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% paste0("C",  seq(from = 1, to = 12))] <- "10 mM"

#-------------------------------------------------------

##Experiment 14/05/2024
sample_name_mapping_14_05_2024 <- c("C1" = "BAA033",
                                    "C3" = "BAA034",
                                    "C5" = "WT",
                                    "C7" = "- CTRL", 
                                    "E1" = "BAA033",
                                    "E3" = "BAA034",
                                    "E5" = "WT",
                                    "E7" = "- CTRL") #pass all the names based on the wells that have been used

## Assign specific values to samples
final_plate_reader_merged$Medium[final_plate_reader_merged$Sample %in% c(paste0("C",  seq(from = 1, to = 7, by = 2)), paste0("E",  seq(from = 1, to = 7, by = 2)))] <- "MM" #add more variables as done here if required
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% paste0("C",  seq(from = 1, to = 7, by = 2))] <- "0 mM"
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% paste0("E",  seq(from = 1, to = 7, by = 2))] <- "10 mM"

#-------------------------------------------------------

##Experiment 17/05/2024
sample_name_mapping_17_05_2024 <- c("A1" = "BAA046_2",
                                    "A3" = "BAA046_5",
                                    "A5" = "BAA047",
                                    "A7" = "BAA048_2",
                                    "A9" = "BAA048_3",
                                    "A11" = "BAA048_4",
                                    "C1" = "BAA013",
                                    "C3" = "- CTRL",
                                    "E1" = "BAA046_2",
                                    "E3" = "BAA046_5",
                                    "E5" = "BAA047",
                                    "E7" = "BAA048_2",
                                    "E9" = "BAA048_3",
                                    "E11" = "BAA048_4",
                                    "G1" = "BAA013",
                                    "G3" = "- CTRL") #pass all the names based on the wells that have been used

## Assign specific values to samples
final_plate_reader_merged$Medium[final_plate_reader_merged$Sample %in% c(paste0("A",  seq(from = 1, to = 11, by = 2)), c("C1", "C3"), paste0("E",  seq(from = 1, to = 11, by = 2)), c("G1", "G3"))] <- "MM" #add more variables as done here if required
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% c(paste0("A",  seq(from = 1, to = 11, by = 2)), "C1", "C3")] <- "0 mM"
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% c(paste0("E",  seq(from = 1, to = 11, by = 2)), "G1", "G3")] <- "10 mM"

#------------------------------
##Experiment 06/06/2024
sample_name_mapping_06_06_2024 <- c("A1" = "BAA013_old_cryo",
                                    "A2" = "BAA033_2",
                                    "A3" = "BAA033_4",
                                    "A4" = "BAA033_8",
                                    "A5" = "BAA033_1",
                                    "A6" = "BAA033_5",
                                    "A7" = "BAA033_6",
                                    "A8" = "BAA033_7",
                                    "A9" = "BAA033_old",
                                    "A10" = "BAA034",
                                    "A11" = "1S145",
                                    "A12" = "- CTRL",
                                    "B1" = "BAA013",
                                    "B2" = "BAA033_2",
                                    "B3" = "BAA033_4",
                                    "B4" = "BAA033_8",
                                    "B5" = "BAA033_1",
                                    "B6" = "BAA033_5",
                                    "B7" = "BAA033_6",
                                    "B8" = "BAA033_7",
                                    "B9" = "BAA033_old",
                                    "B10" = "BAA034",
                                    "B11" = "1S145",
                                    "B12" = "- CTRL") #pass all the names based on the wells that have been used

## Assign specific values to samples
final_plate_reader_merged$Medium[final_plate_reader_merged$Sample %in% c(paste0("A",  seq(from = 1, to = 12)), paste0("B",  seq(from = 1, to = 12)))] <- "MM" #add more variables as done here if required
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% c(paste0("A",  seq(from = 1, to = 12)))] <- "0 mM"
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% c(paste0("B",  seq(from = 1, to = 12)))] <- "10 mM"

#------------------------------
##Experiment 17/06/2024
sample_name_mapping_17_06_2024 <- c("A1" = "1S145",
                                    "A2" = "BAA013",
                                    "A3" = "BAA024",
                                    "A4" = "BAA033_1",
                                    "A5" = "BAA033_2",
                                    "A6" = "BAA033_3",
                                    "A7" = "BAA034",
                                    "A8" = "BAA049",
                                    "A9" = "- CTRL",
                                    "C1" = "1S145",
                                    "C2" = "BAA013",
                                    "C3" = "BAA024",
                                    "C4" = "BAA033_1",
                                    "C5" = "BAA033_2",
                                    "C6" = "BAA033_3",
                                    "C7" = "BAA034",
                                    "C8" = "BAA049",
                                    "C9" = "- CTRL") #pass all the names based on the wells that have been used

## Assign specific values to samples
final_plate_reader_merged$Medium[final_plate_reader_merged$Sample %in% c(paste0("A",  seq(from = 1, to = 9)), paste0("C",  seq(from = 1, to = 9)))] <- "MM" #add more variables as done here if required
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% c(paste0("A",  seq(from = 1, to = 9)))] <- "0 mM"
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% c(paste0("C",  seq(from = 1, to = 9)))] <- "10 mM"
#==================================

#Experiment 27/02/2024
##OD600 dataset
input_OD <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-02-27 Demetra/Plate_reader_input_OD.xlsx")
##GFP dataset
input_GFP <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-02-27 Demetra/Plate_reader_input_GFP.xlsx")

#------------------------------
#Experiment 07/03/2024
##OD600 dataset
input_OD <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-03-07 Demetra/2024-03-07_Plate_reader_input_OD.xlsx")
##GFP dataset
input_GFP <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-03-07 Demetra/2024-03-07_Plate_reader_input_GFP.xlsx")
#------------------------------

#Experiment 24/04/2024
##OD600 dataset
input_OD <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-04-24 GFP and mCherry/2024-04-24 OD.xlsx")
##GFP dataset
input_GFP <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-04-24 GFP and mCherry/2024-04-24 GFP.xlsx")
##mCherry dataset
input_mCherry <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-04-24 GFP and mCherry/2024-04-24 mCherry.xlsx")

#-------------------------------
#Experiment 14/05/2024
##OD600 dataset
input_OD <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-05-14 BAA033 and BAA034 GFP and mCherry new settings/2024-05-14 OD.xlsx")
##GFP dataset
input_GFP <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-05-14 BAA033 and BAA034 GFP and mCherry new settings/2024-05-14 GFP.xlsx")
##mCherry dataset
input_mCherry <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-05-14 BAA033 and BAA034 GFP and mCherry new settings/2024-05-14 mCherry.xlsx")

#-------------------------------
#Experiment 17/05/2024
##OD600 dataset
input_OD <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-05-17 BAA013 and new strains/2024-05-17 OD.xlsx")
##GFP dataset
input_GFP <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-05-17 BAA013 and new strains/2024-05-17 GFP.xlsx")

#-------------------------------
#Experiment 06/06/2024
##OD600 dataset
input_OD <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-06-06 Alberto and Demetra/2024-06-06 OD.xlsx")
##GFP dataset
input_GFP <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-06-06 Alberto and Demetra/2024-06-06 GFP.xlsx")
##mCherry dataset
input_mCherry <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-06-06 Alberto and Demetra/2024-06-06 mCherry.xlsx")

#-------------------------------
#Experiment 17/06/2024
##OD600 dataset
input_OD <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-06-17 Alberto/2024-06-17 OD.xlsx")
##GFP dataset
input_GFP <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-06-17 Alberto/2024-06-17 GFP.xlsx")
##mCherry dataset
input_mCherry <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-06-17 Alberto/2024-06-17 mCherry.xlsx")
#==================================

#remove date that is created when importing the excel files
input_OD$Time <- gsub("^1899-12-31 ", "", input_OD$Time)
input_GFP$Time <- gsub("^1899-12-31 ", "", input_GFP$Time)

#parse time strings for OD dataset
time_values <- strsplit(as.character(input_OD$Time), ":")
hours <- sapply(time_values, function(x) as.numeric(x[1]))
minutes <- sapply(time_values, function(x) as.numeric(x[2]))
seconds <- sapply(time_values, function(x) as.numeric(x[3]))

#calculate time as fractions of hours for OD dataset
input_OD$Time <- hours + minutes/60 + seconds/3600

#use the same time stamps for the GFP dataset
input_GFP$Time <- input_OD$Time

#remove temperature column
input_OD <- input_OD %>% select(-"T° Read 3:600")
input_GFP <- input_GFP %>% select(-"T° GFP:485,528")

#merge OD and GFP data frames
merged_plate_reader <- merge(input_OD, input_GFP, by = c("Time"), suffixes = c("_OD", "_GFP"))

#rearrange OD dataset
longer_merged_OD <- pivot_longer(merged_plate_reader[, 1:97], 
                                 cols = seq(2,97),
                                 names_to = "Sample_OD",
                                 values_to = "OD_value")

#rearrange GFP dataset
longer_merged_GFP <- pivot_longer(merged_plate_reader[, c(1, 98:193)], 
                              cols = seq(2, 97),
                              names_to = "Sample_GFP",
                              values_to = "GFP_value")


longer_merged_OD$Sample_OD <- substr(longer_merged_OD$Sample_OD, 1, nchar(longer_merged_OD$Sample_OD) - 3) 
longer_merged_OD <- longer_merged_OD %>% dplyr::rename("Sample" = Sample_OD)

longer_merged_GFP$Sample_GFP <- substr(longer_merged_GFP$Sample_GFP, 1, nchar(longer_merged_GFP$Sample_GFP) - 4)
longer_merged_GFP <- longer_merged_GFP %>% dplyr::rename("Sample" = Sample_GFP)

final_plate_reader_merged <- merge(longer_merged_OD, longer_merged_GFP, by = c("Time", "Sample"))

#Add values to samples
## Initialize  new columns with NA values
final_plate_reader_merged$Medium <- NA
final_plate_reader_merged$IPTG <- NA

## Assign specific values to samples (replace with lines under each experiment's date found at the start of the script)
final_plate_reader_merged$Medium[final_plate_reader_merged$Sample %in% c(paste0("A",  seq(from = 1, to = 12, by = 2)), paste0("C",  seq(from = 1, to = 12, by = 2)), "A12", "C12")] <- "MM"
final_plate_reader_merged$Medium[final_plate_reader_merged$Sample %in% c(paste0("E",  seq(from = 1, to = 12, by = 2)), paste0("G",  seq(from = 1, to = 12, by = 2)), "E12", "G12")] <- "LB"
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% c(paste0("A",  seq(from = 1, to = 12, by = 2)), paste0("E",  seq(from = 1, to = 12, by = 2)), "A12", "E12")] <- "0 mM"
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% c(paste0("C",  seq(from = 1, to = 12, by = 2)), paste0("G",  seq(from = 1, to = 12, by = 2)), "C12", "G12")] <- "1 mM"

final_plate_reader_merged <- final_plate_reader_merged %>% drop_na() #Once you have passed in the medium only to the full wells, use this function to remove all other wells that have not been used

#Give sample names to wells
## Define a vector to map old sample names to new ones
sample_name_mapping <- sample_name_mapping_17_05_2024 #choose the one corresponding to the experiment data (see experiment dates at the start of the script)

## Update the Sample column with the new names
final_plate_reader_merged$Sample <- ifelse(final_plate_reader_merged$Sample %in% names(sample_name_mapping), 
                    sample_name_mapping[final_plate_reader_merged$Sample], 
                    final_plate_reader_merged$Sample)

ggplot(data = filter(final_plate_reader_merged, IPTG == "10 mM"),
       mapping = aes(x = Time,
                     group = Sample)) +
  geom_point(aes(y = GFP_value), shape = 1) +
  geom_line(aes(y = GFP_value, col = Sample)) +
  geom_point(aes(y = OD_value*30000)) +
  geom_line(aes(y = OD_value*30000, col = Sample)) +
  scale_y_continuous(sec.axis = sec_axis(~./30000 , name = 'OD value')) +
  labs(y = "GFP value") +
  theme_bw()



#=====================================
#2 fluo channels (e.g. GFP and mCherry)

#remove date that is created when importing the excel files
input_OD$Time <- gsub("^1899-12-31 ", "", input_OD$Time)
input_GFP$Time <- gsub("^1899-12-31 ", "", input_GFP$Time)
input_mCherry$Time <- gsub("^1899-12-31 ", "", input_mCherry$Time)

#parse time strings for OD dataset
time_values <- strsplit(as.character(input_OD$Time), ":")
hours <- sapply(time_values, function(x) as.numeric(x[1]))
minutes <- sapply(time_values, function(x) as.numeric(x[2]))
seconds <- sapply(time_values, function(x) as.numeric(x[3]))

#calculate time as fractions of hours for OD dataset
input_OD$Time <- hours + minutes/60 + seconds/3600

#use the same time stamps for the GFP dataset
input_GFP$Time <- input_OD$Time

#use the same time stamps for the mCherry dataset
input_mCherry$Time <- input_OD$Time

#remove temperature column
input_OD <- input_OD %>% select(-"T° Read 3:600")
input_GFP <- input_GFP %>% select(-"T° GFP_mCherry:485,528")
input_mCherry <- input_mCherry %>% select(-"T° GFP_mCherry:583,613") #NOTE: when using experiments before the one on 14/05/2024, use "T° GFP_mCherry:587,610"

#merge OD, GFP and mCherry data frames
merged_plate_reader_OD_GFP <- merge(input_OD, input_GFP, by = c("Time"), suffixes = c("_OD", "_GFP"))
merged_plate_reader <- merge(merged_plate_reader_OD_GFP, input_mCherry, by = c("Time"), suffixes = c("", "_mCherry"))

#rearrange OD dataset
longer_merged_OD <- pivot_longer(merged_plate_reader[, 1:97], 
                                 cols = seq(2,97),
                                 names_to = "Sample_OD",
                                 values_to = "OD_value")

#rearrange GFP dataset
longer_merged_GFP <- pivot_longer(merged_plate_reader[, c(1, 98:193)], 
                                  cols = seq(2, 97),
                                  names_to = "Sample_GFP",
                                  values_to = "GFP_value")

#rearrange mCherry dataset
longer_merged_mCherry <- pivot_longer(merged_plate_reader[, c(1, 194:289)], 
                                  cols = seq(2, 97),
                                  names_to = "Sample_mCherry",
                                  values_to = "mCherry_value")


#change names of columns and values
longer_merged_OD$Sample_OD <- substr(longer_merged_OD$Sample_OD, 1, nchar(longer_merged_OD$Sample_OD) - 3) 
longer_merged_OD <- longer_merged_OD %>% dplyr::rename("Sample" = Sample_OD)

longer_merged_GFP$Sample_GFP <- substr(longer_merged_GFP$Sample_GFP, 1, nchar(longer_merged_GFP$Sample_GFP) - 4)
longer_merged_GFP <- longer_merged_GFP %>% dplyr::rename("Sample" = Sample_GFP)

longer_merged_mCherry <- longer_merged_mCherry %>% dplyr::rename("Sample" = Sample_mCherry)

final_plate_reader_merged <- merge(longer_merged_OD, longer_merged_GFP, by = c("Time", "Sample"))
final_plate_reader_merged <- merge(final_plate_reader_merged, longer_merged_mCherry, by = c("Time", "Sample"))

#Add values to samples
## Initialize  new columns with NA values
final_plate_reader_merged$Medium <- NA
final_plate_reader_merged$IPTG <- NA

## Assign specific values to samples
final_plate_reader_merged$Medium[final_plate_reader_merged$Sample %in% c(paste0("C",  seq(from = 1, to = 7, by = 2)), paste0("E",  seq(from = 1, to = 7, by = 2)))] <- "MM" #add more variables as done here if required
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% paste0("C",  seq(from = 1, to = 7, by = 2))] <- "0 mM"
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% paste0("E",  seq(from = 1, to = 7, by = 2))] <- "10 mM"


final_plate_reader_merged <- final_plate_reader_merged %>% drop_na() #Once you have passed in the medium only to the full wells, use this function to remove all other wells that have not been used

#Give sample names to wells
## Define a vector to map old sample names to new ones
sample_name_mapping <- sample_name_mapping_17_06_2024

## Update the Sample column with the new names
final_plate_reader_merged$Sample <- ifelse(final_plate_reader_merged$Sample %in% names(sample_name_mapping), 
                                           sample_name_mapping[final_plate_reader_merged$Sample], 
                                           final_plate_reader_merged$Sample)

#Plot with OD and GFP
ggplot(data = filter(final_plate_reader_merged, IPTG == "10 mM"),
       mapping = aes(x = Time,
                     group = Sample)) +
  geom_point(aes(y = GFP_value), shape = 1) +
  geom_line(aes(y = GFP_value, col = Sample)) +
  geom_point(aes(y = OD_value*30000)) +
  geom_line(aes(y = OD_value*30000, col = Sample)) +
  #geom_text(data = . %>% group_by(Sample) %>% slice_tail(n = 1),  # Labels at the end of each curve
            #aes(label = Sample, x = Time, y = OD_value*30000, color = Sample),
            #hjust = 1.2, vjust = 0.5, size = 3) +  # Adjust position and size of the labels
  scale_y_continuous(sec.axis = sec_axis(~./30000 , name = 'OD value')) +
  labs(y = "GFP value") +
  theme_bw() +
  facet_wrap(~Sample)

#Plot with OD and mCherry
ggplot(data = filter(final_plate_reader_merged, IPTG == "10 mM"),
       mapping = aes(x = Time,
                     group = Sample)) +
  geom_point(aes(y = mCherry_value), shape = 1) +
  geom_line(aes(y = mCherry_value, col = Sample)) +
  geom_point(aes(y = OD_value*5000), shape = 16) +
  geom_line(aes(y = OD_value*5000, col = Sample)) +
  #geom_text(data = . %>% group_by(Sample) %>% slice_tail(n = 1),  # Labels at the end of each curve
            #aes(label = Sample, x = Time, y = OD_value/5000, color = Sample),
            #hjust = 1.2, vjust = 0.5, size = 3) +  # Adjust position and size of the labels
  scale_y_continuous(sec.axis = sec_axis(~./5000 , name = 'OD value')) +
  labs(y = "mCherry value") +
  theme_bw() +
  facet_wrap(~Sample)
