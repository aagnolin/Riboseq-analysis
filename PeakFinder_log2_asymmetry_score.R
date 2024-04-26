# Load required libraries
library(dplyr)
library(readr)

# Extract command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 2) {
  stop("Usage: Rscript your_script.R alt_predict_file.csv gene_info_df.csv <subtract/add_value> <subtract/add_value>")
}

# Assign arguments to variables
alt_predict_file <- args[1]
gene_info_df_file <- args[2]

# Now you can use these variables in your script
#load input files
input_data <- read_csv(alt_predict_file)

#load gene info
gene_info_df <- read_csv(gene_info_df_file)

#Choose target positions first half of genes
target_sequences_PeakFinder <- data.frame(gene_info_df$gene_length)
target_sequences_PeakFinder <- gene_info_df %>%
  mutate(target_start = ifelse(Strand == "+", StartPosition, StartPosition),
         target_end = ifelse(Strand == "+", floor(StartPosition + (gene_length/2)), floor(StartPosition + (gene_length/2)))) %>% 
  select(locus_tag, target_start, target_end, gene_length, StartPosition)

#merge input data and target sequences
input_data_merge <- merge(input_data, target_sequences_PeakFinder, by = "locus_tag", all = TRUE)

#Sequences outside the ORFs do not have a locus tag assigned, but this is necessary if we want to analyse peaks in the upstream region of the ORFs
##Sort target_sequences_PeakFinder by target_start for binary search
target_sequences_PeakFinder_sorted <- target_sequences_PeakFinder[order(target_sequences_PeakFinder$target_start), ]

cat("Finding peaks in first half of genes...\n")

##Binary search to assign locus_tag based on position
assign_locus_tag <- function(position) {
  left <- 1
  right <- nrow(target_sequences_PeakFinder_sorted)
  
  while (left <= right) {
    mid <- floor((left + right) / 2)
    
    if (position >= target_sequences_PeakFinder_sorted[mid, "target_start"] && position <= target_sequences_PeakFinder_sorted[mid, "target_end"]) {
      return(target_sequences_PeakFinder_sorted[mid, c("locus_tag", "target_start", "target_end", "StartPosition")])
    } else if (position < target_sequences_PeakFinder_sorted[mid, "target_start"]) {
      right <- mid - 1
    } else {
      left <- mid + 1
    }
  }
  
  return(NA)
}

##Assign locus_tag, target_start, and target_end for rows with NA locus_tag that are in the target position range
for (i in 1:nrow(input_data_merge)) {
  if (is.na(input_data_merge[i, "locus_tag"])) {
    result <- assign_locus_tag(input_data_merge[i, "position"])
    input_data_merge[i, c("locus_tag", "target_start", "target_end", "StartPosition")] <- result
  }
}

#filter data that is between the provided positions and add relative position
input_data_max <- group_by(input_data_merge, locus_tag) %>% filter(position >= target_start & position <= target_end)
input_data_max <- mutate(input_data_max, relative_position = position - StartPosition) %>%  select(- sequence)

#add gene length
input_data_max <- merge(input_data_max, gene_info_df, by = "locus_tag", all = TRUE) %>% select(- "gene_length.x")

#Write output file 1
write_csv(input_data_max, "output_left_half.csv")

#sum Norm_count in target positions for each gene, then divide by the target sequence length
df_Ribo_reads_target_1 <- input_data_max %>% group_by(locus_tag) %>% summarize(Sum_Norm_count_first_half = sum(Norm_count))

#Second side
#Choose target positions second half of genes
target_sequences_PeakFinder <- data.frame(gene_info_df$gene_length)
target_sequences_PeakFinder <- gene_info_df %>%
  mutate(target_start = ifelse(Strand == "+", ceiling(StartPosition + (gene_length/2)), ceiling(StartPosition + (gene_length/2))),
         target_end = ifelse(Strand == "+", EndPosition, EndPosition)) %>%
  select(locus_tag, target_start, target_end, gene_length, StartPosition)

#merge input data and target sequences
input_data_merge <- merge(input_data, target_sequences_PeakFinder, by = "locus_tag", all = TRUE)

#Sequences outside the ORFs do not have a locus tag assigned, but this is necessary if we want to analyse peaks in the upstream region of the ORFs
##Sort target_sequences_PeakFinder by target_start for binary search
target_sequences_PeakFinder_sorted <- target_sequences_PeakFinder[order(target_sequences_PeakFinder$target_start), ]

cat("Finding peaks in second half of genes...\n")

##Binary search to assign locus_tag based on position
assign_locus_tag <- function(position) {
  left <- 1
  right <- nrow(target_sequences_PeakFinder_sorted)
  
  while (left <= right) {
    mid <- floor((left + right) / 2)
    
    if (position >= target_sequences_PeakFinder_sorted[mid, "target_start"] && position <= target_sequences_PeakFinder_sorted[mid, "target_end"]) {
      return(target_sequences_PeakFinder_sorted[mid, c("locus_tag", "target_start", "target_end", "StartPosition")])
    } else if (position < target_sequences_PeakFinder_sorted[mid, "target_start"]) {
      right <- mid - 1
    } else {
      left <- mid + 1
    }
  }
  
  return(NA)
}

##Assign locus_tag, target_start, and target_end for rows with NA locus_tag that are in the target position range
for (i in 1:nrow(input_data_merge)) {
  if (is.na(input_data_merge[i, "locus_tag"])) {
    result <- assign_locus_tag(input_data_merge[i, "position"])
    input_data_merge[i, c("locus_tag", "target_start", "target_end", "StartPosition")] <- result
  }
}

#filter data that is between the provided positions and add relative position
input_data_max <- group_by(input_data_merge, locus_tag) %>% filter(position >= target_start & position <= target_end)
input_data_max <- mutate(input_data_max, relative_position = position - StartPosition) %>%  select(- sequence)

#add gene length
input_data_max <- merge(input_data_max, gene_info_df, by = "locus_tag", all = TRUE) %>% select(- "gene_length.x")

#Write output file 2
write_csv(input_data_max, "output_right_half.csv")

df_Ribo_reads_target_2 <- input_data_max %>% group_by(locus_tag) %>% summarize(Sum_Norm_count_second_half = sum(Norm_count))
Merged_target_ranges <- merge(df_Ribo_reads_target_1, df_Ribo_reads_target_2, by = "locus_tag", all = TRUE)
Merged_target_ranges <- Merged_target_ranges %>% mutate(log2_asymmetry_score = log2(Sum_Norm_count_second_half/Sum_Norm_count_first_half))

#Write output
write_csv(Merged_target_ranges, "output_asymmetry_score.csv")
