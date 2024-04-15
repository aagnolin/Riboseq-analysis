#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(dplyr)

# help section
print_help <- function() {
  cat("Usage: run_analysis.R <input_file_1> <input_file_2> <output_location> <cutoff_value> <include_ggplot> [<x_label>] [<y_label>]\n\n")
  cat("Arguments:\n")
  cat("  <input_file_1>:     Path to the first input CSV file.\n")
  cat("  <input_file_2>:     Path to the second input CSV file.\n")
  cat("  <output_location>:  Path to the output CSV file.\n")
  cat("  <cutoff_value>:     Cutoff value for filtering.\n")
  cat("  <include_ggplot>:   Whether to include the ggplot output (TRUE or FALSE).\n")
  cat("  [<x_label>]:        Optional. X-axis label for the plot.\n")
  cat("  [<y_label>]:        Optional. Y-axis label for the plot.\n\n")
  cat("Example:\n")
  cat("  Rscript run_analysis.R input_file_1.csv input_file_2.csv output.csv 1 TRUE \"Avg normalised reads X\" \"Avg normalised reads Y\"\n")
}

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if help flag is provided
if ("--help" %in% args || "-h" %in% args) {
  print_help()
  quit(status = 0)
}

# Check if the correct number of arguments is provided
if (length(args) < 5 || length(args) > 7) {
  cat("Error: Incorrect number of arguments!\n\n")
  print_help()
  quit(status = 1)
}


# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file_1 <- args[1]
input_file_2 <- args[2]
output_location <- args[3]
cutoff_value <- as.numeric(args[4])
scatter_plot <- as.logical(args[5])
x_label <- args[6]
y_label <- args[7]

# Read input files
riboseq_1 <- read.csv(input_file_1)
riboseq_2 <- read.csv(input_file_2)

# Merge the datasets by position
Merged_same_position <- merge(riboseq_1, riboseq_2, 'position')

# Filter based on cutoff value
Merged_same_position_cutoff <- filter(Merged_same_position, Norm_count.x >= cutoff_value & Norm_count.y >= cutoff_value)

# Perform total least squares regression
data <- cbind(Merged_same_position_cutoff$Norm_count.x, Merged_same_position_cutoff$Norm_count.y)
logData <- log10(data)
meanX <- colMeans(logData)
mcX <- sweep(logData, 2, meanX, FUN="-")
pca <- svd(mcX)
scores <- pca$u %*% diag(pca$d)
loads <- pca$v

# Calculate explained variance
Xfit <- sweep(scores[,1] %*% t(loads[,1]), 2, meanX, FUN="+")
XfitMc <- scores[,1]%*%t(loads[,1])
adjusted_explained_variance <- sum(XfitMc^2)/sum(mcX^2)*100

# Extracting coefficients of the green line
coefficients_tls_line <- coef(lm(Xfit[, 2] ~ Xfit[, 1] + 1))
equation <- paste("y =", round(coefficients_tls_line[2], 2), "* x +", round(coefficients_tls_line[1], 2))

# Output ggplot if requested
if (scatter_plot) {
  gg <- ggplot(data = as.data.frame(logData),
               mapping = aes(x = V1,
                             y = V2)) +
    geom_point(size = 1, alpha = 0.3) +
    geom_abline(intercept = coefficients_tls_line[1], slope = coefficients_tls_line[2], col = "darkgreen", linewidth = 1) +
    annotate("text", x = min(logData[, 1]), y = max(logData[, 2]), label = paste0("Explained variance = ", round(adjusted_explained_variance, 2), "%"), hjust = 0, vjust = 1, size = 4) +
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
  
  # Set custom axis labels
  if (!is.null(x_label)) {
    gg <- gg + labs(x = x_label)
  }
  if (!is.null(y_label)) {
    gg <- gg + labs(y = y_label)
  }
  
  print(gg)
}

# Rename axes and create output dataset
Output_tls_nt_replicates <- as.data.frame(logData) %>% rename(Norm_count.x = V1, Norm_count.y = V2)
Output_tls_nt_replicates <- cbind(Output_tls_nt_replicates, adjusted_explained_variance, equation)
write.csv(Output_tls_nt_replicates, file = output_location, row.names = FALSE)
