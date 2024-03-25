library(tidyverse)
library(Biostrings)

# Read the reference genome in FASTA format
genome <- readDNAStringSet("C:/Users/aagnoli/OneDrive - UvA/B. subtilis 168 annotation and genome files/WT-Prspb-amyI.fa")

# Read the GFF file
gff_lines <- readLines("C:/Users/aagnoli/OneDrive - UvA/B. subtilis 168 annotation and genome files/wt-prspb-amyI-GFF3_updated_amyM_positions.gff")


#Extract information for B. subtilis genes (including AmyM) 
locus_tags <- list()
start_positions <- list()
end_positions <- list()
sequences <- DNAStringSet()
strands <- character()

# Process each line of the GFF file
for (line in gff_lines) {
  # Check if the line contains gene information
  if (grepl("gene\\s", line)) {
    # Split the line by tabs
    line_parts <- unlist(strsplit(line, "\t"))
    
    # Extract the relevant information
    locus_tag <- gsub(".*locus_tag=(\\S+).*", "\\1", line_parts[9])
    
    # Remove semicolon and text after it
    locus_tag <- gsub(";.*", "", locus_tag)
    
    start_pos <- as.numeric(line_parts[4])
    end_pos <- as.numeric(line_parts[5])
    
    # Extract the strand information from the 7th column
    strand <- line_parts[7]
    
    # Check if the LocusTag starts with "B"
    if (substring(locus_tag, 1, 1) == "B") {
      # Append the values to the lists
      locus_tags <- append(locus_tags, locus_tag)
      start_positions <- append(start_positions, start_pos)
      end_positions <- append(end_positions, end_pos)
      strands <- append(strands, strand)
      
      # Extract the DNA sequence for the current locus tag and adjust for the strand
      if (strand == "+") {
        seq <- subseq(genome, start_pos, end_pos)
      } else if (strand == "-") {
        seq <- reverseComplement(subseq(genome, start_pos, end_pos))
      } else {
        # Handle unrecognized strand information
        seq <- DNAString("")
      }
      
      sequences <- c(sequences, seq)
    }
  }
}

# Create a data frame with all combined information for genes
gene_info_df <- data.frame(
  LocusTag = unlist(locus_tags),
  StartPosition = unlist(start_positions),
  EndPosition = unlist(end_positions),
  Sequence = as.character(sequences),
  Strand = unlist(strands)
)


#load input files (ppGpp data, 8h)
##5-3'
Normalized_212_8h_1_53_full <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/ppGpp/Check RBS peak/Normalized_212_8h_1_53_full.csv")
Normalized_212_8h_2_53_full <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/ppGpp/Check RBS peak/Normalized_212_8h_2_53_full.csv")
Normalized_ppGpp_8h_1_53_full <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/ppGpp/Check RBS peak/Normalized_ppGpp_8h_1_53_full.csv")
Normalized_ppGpp_8h_2_53_full <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/ppGpp/Check RBS peak/Normalized_ppGpp_8h_2_53_full.csv")
##3-5'
Normalized_212_8h_1_35_full <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/ppGpp/Check RBS peak/Normalized_212_8h_1_35_full.csv")
Normalized_212_8h_2_35_full <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/ppGpp/Check RBS peak/Normalized_212_8h_2_35_full.csv")
Normalized_ppGpp_8h_1_35_full <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/ppGpp/Check RBS peak/Normalized_ppGpp_8h_1_35_full.csv")
Normalized_ppGpp_8h_2_35_full <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/ppGpp/Check RBS peak/Normalized_ppGpp_8h_2_35_full.csv")

# Identify target sequences containing RBSs using the gene_info_df from the script of all riboseq analyses. Takes strand directionality into consideration
target_sequences <- gene_info_df %>%
  mutate(target_start = ifelse(Strand == "+", StartPosition - 21, EndPosition - 14),  # Adjust for strand direction
         target_end = ifelse(Strand == "+", StartPosition + 14, EndPosition + 21)) %>%
  select(LocusTag, target_start, target_end)

# SEQUENCE FINDER
# Initialize a list to store the subsequences
subseq_list <- list()

# Iterate over each row of target_sequences
for (i in 1:nrow(target_sequences)) {
  # Extract start and end positions for the ith row
  start_pos <- target_sequences[i, "target_start"]
  end_pos <- target_sequences[i, "target_end"]
  
  # Extract subsequence from the genome
  if (gene_info_df$Strand[i] == "+") {
    subseq_list[[i]] <- as.character(subseq(genome, start_pos, end_pos))
  } else if (gene_info_df$Strand[i] == "-") {
    subseq_list[[i]] <- as.character(reverseComplement(subseq(genome, start_pos, end_pos)))
  } else {
    # Handle unrecognized strand information
    subseq_list[[i]] <- NA  # Or any other suitable handling
  }
}

# Combine all subsequences into a single data frame
subseq_df <- data.frame(
  target_sequence = unlist(subseq_list),
  locus_tag = target_sequences$LocusTag
)
#-------------------------------------------------------------------------------
#PEAK FINDER
#select the input dataset
input_data <- Normalized_212_8h_1_53_full %>% dplyr::rename(LocusTag = "locus_tag")
input_data <- merge(input_data, target_sequences, by = "LocusTag", all = TRUE)

#Sequences outside the ORFs do not have a locus tag assigned, but this is necessary if we want to analyse peaks in the upstream region of the ORFs
##Sort target_sequences by target_start for binary search
target_sequences_sorted <- target_sequences[order(target_sequences$target_start), ]

##Binary search to assign locus_tag based on position
assign_locus_tag <- function(position) {
  left <- 1
  right <- nrow(target_sequences_sorted)
  
  while (left <= right) {
    mid <- floor((left + right) / 2)
    
    if (position >= target_sequences_sorted[mid, "target_start"] && position <= target_sequences_sorted[mid, "target_end"]) {
      return(target_sequences_sorted[mid, c("LocusTag", "target_start", "target_end")])
    } else if (position < target_sequences_sorted[mid, "target_start"]) {
      right <- mid - 1
    } else {
      left <- mid + 1
    }
  }
  
  return(NA)
}

##Assign locus_tag, target_start, and target_end for rows with NA locus_tag that are in the target position range
for (i in 1:nrow(input_data)) {
  if (is.na(input_data[i, "LocusTag"])) {
    result <- assign_locus_tag(input_data[i, "position"])
    input_data[i, c("LocusTag", "target_start", "target_end")] <- result
  }
}

input_data <- input_data %>% dplyr::rename(locus_tag = "LocusTag")
#Filter highest peaks on target regions of each locus tag and calculate the position of those peaks relative to the 1st nt of the start codon
#Filtered_input_data <- input_data %>% filter(locus_tag == "BSU13280" | locus_tag == "BSU13280" | locus_tag == "BSU37350" | locus_tag == "BSU28860" | locus_tag == "BSU36660" | locus_tag == "BSU36650" | locus_tag == "BSU03780" | locus_tag == "BSU14180" | locus_tag == "BSU18060" | locus_tag == "BSU28290" | locus_tag == "BSU08760", position >= target_start & position <= target_end)
input_data_max <- group_by(input_data, locus_tag) %>% filter(position >= target_start & position <= target_end, Norm_count == max(Norm_count))
input_data_max <- mutate(input_data_max, relative_position_RBS_peak = target_end - position)



#merge dataset with target sequences for motif search
input_data_max <- merge(input_data_max, subseq_df, by = "locus_tag")

#plot results
ggplot(data = input_data_max,
       mapping = aes(x = relative_position_RBS_peak)) +
  geom_bar() +
  scale_x_reverse()

#write list of sequences for the motif search in kplogo
write_csv(input_data_max, file = "C:/Users/aagnoli/Desktop/RBS_peak_finder_output_2.csv")

