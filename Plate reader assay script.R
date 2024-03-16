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
#------------------------------
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
longer_merged_OD <- longer_merged_OD %>% rename("Sample" = Sample_OD)

longer_merged_GFP$Sample_GFP <- substr(longer_merged_GFP$Sample_GFP, 1, nchar(longer_merged_GFP$Sample_GFP) - 4)
longer_merged_GFP <- longer_merged_GFP %>% rename("Sample" = Sample_GFP)

final_plate_reader_merged <- merge(longer_merged_OD, longer_merged_GFP, by = c("Time", "Sample"))

#Add values to samples
## Initialize  new columns with NA values
final_plate_reader_merged$Medium <- NA
final_plate_reader_merged$IPTG <- NA

## Assign specific values to samples
final_plate_reader_merged$Medium[final_plate_reader_merged$Sample %in% c(paste0("A",  seq(from = 1, to = 12, by = 2)), paste0("C",  seq(from = 1, to = 12, by = 2)), "A12", "C12")] <- "MM" #add more variables as done here if required
final_plate_reader_merged$Medium[final_plate_reader_merged$Sample %in% c(paste0("E",  seq(from = 1, to = 12, by = 2)), paste0("G",  seq(from = 1, to = 12, by = 2)), "E12", "G12")] <- "LB"
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% c(paste0("A",  seq(from = 1, to = 12, by = 2)), paste0("E",  seq(from = 1, to = 12, by = 2)), "A12", "E12")] <- "0 mM"
final_plate_reader_merged$IPTG[final_plate_reader_merged$Sample %in% c(paste0("C",  seq(from = 1, to = 12, by = 2)), paste0("G",  seq(from = 1, to = 12, by = 2)), "C12", "G12")] <- "1 mM"

final_plate_reader_merged <- final_plate_reader_merged %>% drop_na() #Once you have passed in the medium only to the full wells, use this function to remove all other wells that have not been used

#Give sample names to wells
## Define a vector to map old sample names to new ones
sample_name_mapping <- sample_name_mapping_27_02_2024

## Update the Sample column with the new names
final_plate_reader_merged$Sample <- ifelse(final_plate_reader_merged$Sample %in% names(sample_name_mapping), 
                    sample_name_mapping[final_plate_reader_merged$Sample], 
                    final_plate_reader_merged$Sample)

ggplot(data = filter(final_plate_reader_merged, Sample != "clone 1" & Sample != "WT", Medium == "MM", IPTG == "1 mM"),
       mapping = aes(x = Time,
                     group = Sample)) +
  geom_point(aes(y = GFP_value, shape = Sample)) +
  geom_line(aes(y = GFP_value), col = 'limegreen') +
  geom_point(aes(y = OD_value*20000, shape = Sample)) +
  geom_line(aes(y = OD_value*20000), col = "black") +
  scale_y_continuous(sec.axis = sec_axis(~./20000 , name = 'OD value')) +
  labs(y = "GFP value") +
  theme_bw()
