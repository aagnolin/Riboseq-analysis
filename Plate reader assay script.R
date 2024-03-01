#Load required packages
library(readxl)
library(tidyverse)
#Import datasets
#------------------------------
#Experiment 01/03/2024
##OD600 dataset
input_OD <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-02-27 Demetra/Plate_reader_input_OD.xlsx")
##GFP dataset
input_GFP <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/Plate reader assays/2024-02-27 Demetra/Plate_reader_input_GFP.xlsx")

#remove date
input_OD$Time <- gsub("^1899-12-31 ", "", input_OD$Time)
input_GFP$Time <- gsub("^1899-12-31 ", "", input_GFP$Time)

# Parse time strings for OD dataset
time_values <- strsplit(as.character(input_OD$Time), ":")
hours <- sapply(time_values, function(x) as.numeric(x[1]))
minutes <- sapply(time_values, function(x) as.numeric(x[2]))
seconds <- sapply(time_values, function(x) as.numeric(x[3]))

# Calculate time difference in hours for OD dataset
input_OD$Time <- hours + minutes/60 + seconds/3600

#use the same times for the GFP dataset
input_GFP$Time <- input_OD$Time

#remove temperature column
input_OD <- input_OD %>% select(-"T° Read 3:600")
input_GFP <- input_GFP %>% select(-"T° GFP:485,528")

#merge OD and GFP data frames
merged_plate_reader <- merge(input_OD, input_GFP, by = c("Time"), suffixes = c("_OD", "_GFP"))

separate(merged_plate_reader, col = Time, into = c("Sample", "Type"), sep = "_")


longer_merged_OD <- pivot_longer(merged_plate_reader[, 1:97], 
                                 cols = seq(2,97),
                                 names_to = "Sample_OD",
                                 values_to = "OD_value")


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
## Initialize a new column with NA values
final_plate_reader_merged$Medium <- NA

## Assign a specific value to samples A1 to A12
final_plate_reader_merged$Medium[final_plate_reader_merged$Sample %in% paste0("A", 1:12)] <- "LB" #add more variables as done here if required

#final_plate_reader_merged <- final_plate_reader_merged %>% drop_na() #Once you have passed in the medium only to the full wells, use this function to remove all other wells that have not been used

#Give sample names to wells
## Define a vector to map old sample names to new ones
sample_name_mapping <- c("A1" = "WT", 
                         "A3" = "clone 1",
                         "A5" = "clone 2",
                         "A7" = "clone 3",
                         "A9" = "clone 4",
                         "A11" = "clone 5",
                         "A12" = "CTRL") #pass all the names based on the wells that have been used

## Update the Sample column with the new names
final_plate_reader_merged$Sample <- ifelse(final_plate_reader_merged$Sample %in% names(sample_name_mapping), 
                    sample_name_mapping[final_plate_reader_merged$Sample], 
                    final_plate_reader_merged$Sample)

ggplot(data = filter(final_plate_reader_merged, Sample == "clone 1"),
       mapping = aes(x = Time, 
                     y = GFP_value)) +
  geom_point(aes(col = Sample)) +
  geom_line(aes(col = Sample))
