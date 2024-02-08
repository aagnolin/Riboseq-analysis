#github aagnolin https://github.com/aagnolin

# Check if the user requests help or specifies the method
if (length(commandArgs(trailingOnly = TRUE)) == 0 || any(commandArgs(trailingOnly = TRUE) %in% c("-h", "--help"))) {
  cat("Script Usage:\n")
  cat("--------------\n")
  cat("Rscript script_name.R <input_folder> <output_folder> <codon_table> [method] [filter_threshold]\n")
  cat("\n")
  cat("Arguments:\n")
  cat("<input_folder>     : Path to the input folder containing CSV files.\n")
  cat("<output_folder>    : Path to the output folder where files will be saved.\n")
  cat("<codon_table>      : Path to the codon usage table Excel file.\n")
  cat("[method]           : Method for calculating codon observations ('cumulative' for cumulative method, default is 'single').\n")
  cat("[filter_threshold] : Optional filter threshold for data (only include data with count >= filter_threshold).\n")
  cat("\n")
  cat("Description:\n")
  cat("This script performs two actions on the alt_predict_v2 output files supplied in the input folder:\n")
  cat("it creates a table containing only the highest peaks per each gene and another table with the pause\n")
  cat("codon of all peaks present in the input file.\n")
  cat("Additionally, it uses the output file with the extracted pause codons and an additional excel file\n")
  cat("(supplied by the user) containing the usage per codon values of the studied organism and generates\n")
  cat("plots comparing codon usage vs codon occupancy for each file and saves them as SVG images.\n")
  cat("\n")
  cat("Example Usage:\n")
  cat("--------------\n")
  cat("Rscript script_name.R /path/to/input/folder /path/to/output/folder /path/to/codon/table.xlsx\n")
  cat("Rscript script_name.R /path/to/input/folder /path/to/output/folder /path/to/codon/table.xlsx cumulative\n")
  cat("Rscript script_name.R /path/to/input/folder /path/to/output/folder /path/to/codon/table.xlsx single 100\n")
  cat("\n")
  quit(status = 0)
}

# Parse the optional argument for method
method <- "single"  # Default method is "single"
if (length(commandArgs(trailingOnly = TRUE)) > 3) {
  specified_method <- commandArgs(trailingOnly = TRUE)[4]
  if (specified_method == "cumulative") {
    method <- specified_method
  }
}

# Parse the optional argument for filter_threshold
filter_threshold <- NA
if (length(commandArgs(trailingOnly = TRUE)) > 4) {
  filter_threshold <- as.numeric(commandArgs(trailingOnly = TRUE)[5])
}

# Load required libraries
library(dplyr)
library(magrittr)
library(stringr)
library(ggplot2)
library(readxl)
library(ggpubr)
library(tools)

# Get the input and output folder paths from command-line arguments
input_folder <- commandArgs(trailingOnly = TRUE)[1]
output_folder <- commandArgs(trailingOnly = TRUE)[2]
codon_table <- commandArgs(trailingOnly = TRUE)[3]

# Get the list of input files in the input folder
input_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)

# Process each input file
for (input_file in input_files) {
  # Read the input data file
  alt_predict_input <- read.csv(input_file)
  
  # Perform the desired actions on the data
  ## Filter non-coding RNA out
  alt_predict_input <- filter(alt_predict_input, !grepl("^BSU_", locus_tag))
  # alt_predict_input <- filter(alt_predict_input, !grepl("^trn_", gene), gene != 'ssrA')
  ## Filter for highest peaks per gene
  filtered_alt_predict <- alt_predict_input %>%
    group_by(gene) %>%
    filter(count == max(count))
  
  # Generate the output file path based on the input file
  output_file <- file.path(output_folder, paste0(file_path_sans_ext(basename(input_file)), "_highest_peaks_output.csv"))
  
  # Write the filtered data to the output file
  write.csv(filtered_alt_predict, file = output_file, row.names = FALSE)
  
  # Print a message indicating the completion
  cat("Processed file:", input_file, "\n")
  
  # Generate files with extracted pause codon for each peak
  cat("Extracting pause codon from input files...\n")
  
  # Calculate the coding frame and extract the codons in the A (Pause), P and E sites 
  Alt_predict_modulus <- alt_predict_input %>%
    mutate(modulus = offset %% 3) %>%
    mutate(modulus = ifelse(is.na(modulus), 5, modulus))
  
  output_NA <- Alt_predict_modulus %>%
    filter(modulus == 5) %>%
    mutate(Pause_codon = NA, P_site = NA, E_site = NA)
  
  output_0 <- Alt_predict_modulus %>%
    filter(modulus == 0) %>%
    mutate(Pause_codon = str_sub(sequence, 51, 53), 
           P_site = str_sub(sequence, 48, 50),
           E_site = str_sub(sequence, 45, 47),
           motif_input_sequence = str_sub(sequence, 42, 65))
  
  output_1 <- Alt_predict_modulus %>%
    filter(modulus == 1) %>%
    mutate(Pause_codon = str_sub(sequence, 50, 52), 
           P_site = str_sub(sequence, 47, 49),
           E_site = str_sub(sequence, 44, 46),
           motif_input_sequence = str_sub(sequence, 41, 64))
  
  output_2 <- Alt_predict_modulus %>%
    filter(modulus == 2) %>%
    mutate(Pause_codon = str_sub(sequence, 49, 51), 
           P_site = str_sub(sequence, 46, 48),
           E_site = str_sub(sequence, 43, 45),
           motif_input_sequence = str_sub(sequence, 40, 63))
  
  Super_codonator_output <- bind_rows(output_NA, output_0, output_1, output_2) %>%
    arrange(position) %>%
    select(-modulus)
  
  #Calculate and add Pause score column
  
  Super_codonator_output <- Super_codonator_output %>% group_by(gene) %>% mutate(Pause_score = count/sum(count))
  
  # Apply optional filtering if the filter_threshold is provided
  if (!is.na(filter_threshold)) {
    Super_codonator_output_filtered <- Super_codonator_output %>%
      filter(count >= filter_threshold)
  }
  
  # Generate plot to show the filter_threshold cutoff on the distribution of the data based on count
  if (!is.na(filter_threshold)) {
    
    library(ggbreak)
    max_count <- max(Super_codonator_output$count)
    
    density_plot <- ggplot(data = Super_codonator_output,
                           mapping = aes(
                             x = count
                           )) +
      geom_density(adjust = 5, fill = 'red', alpha = 0.3) +
      geom_vline(xintercept = filter_threshold, color = "black", linetype = 'dashed', linewidth = 0.5)
    
    # Calculate the fraction of peaks remaining after applying the filter_threshold
    total_counts <- nrow(Super_codonator_output)
    filtered_counts <- nrow(Super_codonator_output_filtered)
    percentage_filtered_counts <- (filtered_counts / total_counts) * 100
    
    # Add the filter_threshold and percentage to the plot
    density_plot <- density_plot +
      scale_x_break(c(100, max_count - 50), expand = F, space = 0.1) +
      scale_y_continuous(expand = c(0,0)) +
      annotate("text", x = 10, y = 0, label = paste("Count Threshold:", filter_threshold), vjust = -1, hjust = 0) +
      annotate("text", x = 10, y = 0, label = paste("Filtered:", round(percentage_filtered_counts, 2), "%"), vjust = -2, hjust = 0)+
      labs(y = "density") +
      theme_bw() + labs(colour = "Sample") +
      theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
            legend.key.height= unit(0.05, 'cm'),
            legend.key.width= unit(0.05, 'cm')) +
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
            panel.background = element_rect(linetype = "solid")) + theme(axis.text.x = element_text(vjust = 0.5,
                                                                                                    angle = 45))
    
    
    # Generate the output file path for the plot
    output_file_density_plot <- file.path(output_folder, paste0(file_path_sans_ext(basename(input_file)), "density_plot_cutoff.svg"))
    
    # Set the desired width and height for the SVG plot
    plot_width <- 4.375  # Specify the width in pixels
    plot_height <- 3.0208  # Specify the height in pixels
    
    # Export the plot as an SVG image with the specified size
    svg(output_file_density_plot, width = plot_width, height = plot_height)
    print(density_plot)
    dev.off()
    
    message('Density plot of peaks with selected threshold has been generated')
  }
  
  # Generate the output file path for the pause codon output
  output_pause_codon <- file.path(output_folder, paste0(file_path_sans_ext(basename(input_file)), "_pause_codon_output.csv"))
  
  # Write the filtered data to the output file
  write.csv(Super_codonator_output, file = output_pause_codon, row.names = FALSE)
  
  # Print a message indicating the completion
  cat("Processed file:", input_file, "\n")
  
  # Print a final completion message
  cat("Dataset successfully written to the output folder\n")
  
  # Read the Codon usage table
  Codon_usage_table <- read_excel(codon_table)
  
  # Calculate codon observations method "single" for Pause_codon
  if (method == "single") {
    codon_observations_pause <- table(Super_codonator_output$Pause_codon)
    codon_observations_P_site <- table(Super_codonator_output$P_site)
    codon_observations_E_site <- table(Super_codonator_output$E_site)
  }
  
  # Calculate codon observations method "cumulative" for Pause_codon
  if (method == "cumulative") {
    codon_observations_pause <- table(rep(Super_codonator_output$Pause_codon, Super_codonator_output$count))
    codon_observations_P_site <- table(rep(Super_codonator_output$P_site, Super_codonator_output$count))
    codon_observations_E_site <- table(rep(Super_codonator_output$E_site, Super_codonator_output$count))
  }
  
  # Calculate total number of entries for each site
  total_entries_pause <- sum(codon_observations_pause)
  total_entries_P_site <- sum(codon_observations_P_site)
  total_entries_E_site <- sum(codon_observations_E_site)
  
  # Calculate frequency by dividing by total entries for each site
  codon_frequency_normalized_pause <- codon_observations_pause / total_entries_pause
  codon_frequency_normalized_P_site <- codon_observations_P_site / total_entries_P_site
  codon_frequency_normalized_E_site <- codon_observations_E_site / total_entries_E_site
  
  # Create data frames for each site
  codon_df_pause <- data.frame(
    Codon = names(codon_frequency_normalized_pause),
    Occupancy_Pause = as.vector(codon_frequency_normalized_pause),
    stringsAsFactors = FALSE
  )
  
  codon_df_P_site <- data.frame(
    Codon = names(codon_frequency_normalized_P_site),
    Occupancy_P_site = as.vector(codon_frequency_normalized_P_site),
    stringsAsFactors = FALSE
  )
  
  codon_df_E_site <- data.frame(
    Codon = names(codon_frequency_normalized_E_site),
    Occupancy_E_site = as.vector(codon_frequency_normalized_E_site),
    stringsAsFactors = FALSE
  )
  
  # Merge occupancy data frames with usage data frame
  Merged_occupancy_and_usage <- merge(codon_df_pause, Codon_usage_table, by = "Codon", all = TRUE)
  Merged_occupancy_and_usage <- merge(Merged_occupancy_and_usage, codon_df_P_site, by = "Codon", all = TRUE)
  Merged_occupancy_and_usage <- merge(Merged_occupancy_and_usage, codon_df_E_site, by = "Codon", all = TRUE)
  
  # Generate the output file path for the codon occupancy output
  output_codon_occupancy <- file.path(output_folder, paste0(file_path_sans_ext(basename(input_file)), "_codon_occupancy_output.csv"))
  
  # Write the filtered data to the output file
  write.csv(Merged_occupancy_and_usage, file = output_codon_occupancy, row.names = FALSE)
  
  # Print a message indicating the completion
  cat("Processed file:", input_file, "\n")
  
  # Print a final completion message
  cat("Dataset successfully written to the output folder\n")
  
  # Generate the plot
  codon_occupancy_plot <- ggplot(data = Merged_occupancy_and_usage,
                                 mapping = aes(x = Usage, y = Occupancy_Pause)) +
    geom_point(size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    #stat_cor(label.y = 0.1, label.x = 0.7, method = "pearson") +
    geom_text(aes(label = Codon), vjust = -0.5, alpha = 0.5, size = 3) +
    scale_y_continuous(limits = c(0, 0.12)) +
    theme_bw() +
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
  
  # Generate the output file path for the codon occupancy plot
  output_file_codon_occupancy_plot <- file.path(output_folder, paste0(file_path_sans_ext(basename(input_file)), "_codon_occupancy_plot.svg"))
  
  # Set the desired width and height for the SVG plot
  plot_width <- 3.75  # Specify the width in pixels
  plot_height <- 3.33  # Specify the height in pixels
  
  # Export the plot as an SVG image with the specified size
  svg(output_file_codon_occupancy_plot, width = plot_width, height = plot_height)
  print(codon_occupancy_plot)
  dev.off()
  
}
