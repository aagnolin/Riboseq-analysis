# Define usage function
usage <- function() {
  cat("Usage: Rscript PeakFinder_two_sided.R alt_predict_file.csv gene_info_df.csv <subtract/add_value_1> <subtract/add_value_2> <subtract/add_value_3> <subtract/add_value_4>\n")
  cat("\n")
  cat("Arguments:\n")
  cat("  alt_predict_file.csv: output file from alt_predict_v2.py\n")
  cat("  gene_info_df.csv: file containing gene information\n")
  cat("  <subtract/add_value_1>: Value to subtract or add relative to start position of each gene (first range)\n")
  cat("  <subtract/add_value_2>: Value to subtract or add relative to start position of each gene (first range)\n")
  cat("  <subtract/add_value_3>: Value to subtract or add relative to start position of each gene (second range)\n")
  cat("  <subtract/add_value_4>: Value to subtract or add relative to end position of each gene (second range)\n")
}

# Extract command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 6) {
  usage()
  quit(status = 1, save = "no")
}

# Assign arguments to variables
alt_predict_file <- args[1]
gene_info_df_file <- args[2]
value_1 <- as.numeric(args[3])
value_2 <- as.numeric(args[4])
value_3 <- as.numeric(args[5])
value_4 <- as.numeric(args[6])

#load necessary libraries
library(dplyr)
library(readr)

#load input files
input_data <- read_csv(alt_predict_file)

#load gene info
gene_info_df <- read_csv(gene_info_df_file)

# Assign names to output files
input_file_name <- tools::file_path_sans_ext(basename(alt_predict_file))
output_left_side <- paste0(input_file_name, "_output_left_side.csv")
output_right_side <- paste0(input_file_name, "_output_right_side.csv")
merged_ratios <- paste0(input_file_name, "_merged_ratios.csv")
plot_name <- paste0(input_file_name, "_plot.pdf")
raw_plot <- paste0(input_file_name, "_raw_plot.csv")

#Choose target positions
target_sequences_PeakFinder <- gene_info_df %>%
  mutate(target_start = ifelse(Strand == "+", StartPosition + value_1, StartPosition + value_1),
         target_end = ifelse(Strand == "+", StartPosition + value_2, StartPosition + value_2)) %>%
  select(locus_tag, target_start, target_end, StartPosition)

#merge input data and target sequences
input_data_merge <- merge(input_data, target_sequences_PeakFinder, by = "locus_tag", all = TRUE)

#Sequences outside the ORFs do not have a locus tag assigned, but this is necessary if we want to analyse peaks in the upstream region of the ORFs
##Sort target_sequences_PeakFinder by target_start for binary search
target_sequences_PeakFinder_sorted <- target_sequences_PeakFinder[order(target_sequences_PeakFinder$target_start), ]

cat("Finding peaks in provided position range...\n")

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

#save output for optional bar plot
data_plot_left <- input_data_max

#Write output file 1
write_csv(input_data_max, output_left_side)

#sum Norm_count in target positions for each gene, then divide by the target sequence length
df_Ribo_reads_target_1 <- input_data_max %>% group_by(locus_tag) %>% summarize(Sum_Norm_count_initiation = sum(Norm_count)) %>% mutate(Translation_initiation_Ratio = Sum_Norm_count_initiation/((abs(value_1) + abs(value_2))))

#Second side
#Choose new target positions
target_sequences_PeakFinder <- gene_info_df %>%
  mutate(target_start = ifelse(Strand == "+", StartPosition + value_3, StartPosition + value_3),
         target_end = ifelse(Strand == "+", EndPosition + value_4, EndPosition + value_4)) %>%
  select(locus_tag, target_start, target_end, StartPosition)

#merge input data and target sequences
input_data_merge <- merge(input_data, target_sequences_PeakFinder, by = "locus_tag", all = TRUE)

#Sequences outside the ORFs do not have a locus tag assigned, but this is necessary if we want to analyse peaks in the upstream region of the ORFs
##Sort target_sequences_PeakFinder by target_start for binary search
target_sequences_PeakFinder_sorted <- target_sequences_PeakFinder[order(target_sequences_PeakFinder$target_start), ]

cat("Finding peaks in provided secondary position range...\n")

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
write_csv(input_data_max, output_right_side)

df_Ribo_reads_target_2 <- input_data_max %>% group_by(locus_tag) %>% summarize(Sum_Norm_count_ORF = sum(Norm_count)) %>% merge(gene_info_df, by = "locus_tag") %>% mutate(ORF_translation_Ratio = Sum_Norm_count_ORF/(gene_length - value_2)) 
Merged_ratios <- merge(df_Ribo_reads_target_1, df_Ribo_reads_target_2, by = "locus_tag", all = TRUE) %>% select(c(-"StartPosition", -"EndPosition", -"Sequence", -"Strand"))

#Write final output
write_csv(Merged_ratios, merged_ratios)

#create final plot
library(ggplot2)

combined <- bind_rows(data_plot_left, input_data_max) %>% select(-locus_tag, 
                                                                 -genome, 
                                                                 -position, 
                                                                 -strand, 
                                                                 -gene.x, 
                                                                 -offset, 
                                                                 -in_orf_90, 
                                                                 -count, 
                                                                 -Norm_count, 
                                                                 -target_start, 
                                                                 -target_end,
                                                                 -EndPosition,
                                                                 -Sequence, 
                                                                 -Strand, 
                                                                 -gene.y, 
                                                                 -gene_length.y,
                                                                 -StartPosition.x,
                                                                 -StartPosition.y)

p <- ggplot(data = combined,
            mapping = aes(x = relative_position)) +
  geom_freqpoly(bins = sum(abs(value_1), 100), linewidth = 0.8) +
  geom_bar(alpha = 0.5, width = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", alpha = 0.6, linewidth = 0.7) +
  scale_x_continuous(limits = c(value_1, 100)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(x = "Position") +
  theme(plot.subtitle = element_text(face = "bold"),
        plot.caption = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.key = element_rect(linetype = "solid")) +
  theme(axis.ticks = element_line(linewidth = 0.5)) + 
  theme(panel.background = element_rect(fill = NA)) + 
  theme(axis.line = element_line(linetype = "solid")) + 
  theme(axis.line = element_line(linetype = "blank"), 
        panel.grid.major = element_line(colour = "gray89",
                                        linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"),
        panel.background = element_rect(linetype = "solid"))

#Write raw data for plot
write_csv(combined, raw_plot)

#generate plot
ggsave(plot = p, filename = plot_name, device = "pdf")

cat("Peak distribution plot has been generated\n")
