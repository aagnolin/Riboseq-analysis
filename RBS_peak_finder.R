library(tidyverse)
library(Biostrings)

# Read the reference genome in FASTA format
genome <- readDNAStringSet("C:/Users/aagnoli/OneDrive - UvA/B. subtilis 168 annotation and genome files/WT-Prspb-amyI.fa")

# Read the GFF file
gff_lines <- readLines("C:/Users/aagnoli/OneDrive - UvA/B. subtilis 168 annotation and genome files/wt-prspb-amyI-GFF3_updated_amyM_positions.gff")

# Extract information for B. subtilis genes (including AmyM) 
locus_tags <- list()
start_positions <- list()
end_positions <- list()
sequences <- DNAStringSet()
strands <- character()
gene_names <- list()  # Initialize list for gene names

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
    
    # Extract gene name and remove unwanted string after it
    gene_info <- gsub(".*gene=(\\S+);.*", "\\1", line_parts[9])
    gene_name <- gsub(";.*", "", gene_info)
    
    # Check if the LocusTag starts with "B"
    if (substring(locus_tag, 1, 1) == "B") {
      # Append the values to the lists
      locus_tags <- append(locus_tags, locus_tag)
      start_positions <- append(start_positions, start_pos)
      end_positions <- append(end_positions, end_pos)
      strands <- append(strands, strand)
      gene_names <- append(gene_names, gene_name)  # Append gene name
      
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
  locus_tag = unlist(locus_tags),
  StartPosition = unlist(start_positions),
  EndPosition = unlist(end_positions),
  Sequence = as.character(sequences),
  Strand = unlist(strands),
  gene = unlist(gene_names) 
) %>% mutate(gene_length = EndPosition - StartPosition)

#------------------------------------------------

# SEQUENCE FINDER

# Identify target sequences based on position values. Takes strand directionality into consideration
target_sequences <- gene_info_df %>%
  mutate(target_start = ifelse(Strand == "+", StartPosition + 30, EndPosition - 59),
         target_end = ifelse(Strand == "+", StartPosition + 59, EndPosition - 30)) %>%
  select(LocusTag, target_start, target_end)

# Initialize a list to store the subsequences
subseq_list <- list()

# Iterate over each row of target_sequences
for (i in 1:nrow(target_sequences)) {
  # Extract start and end positions for the ith row
  start_pos <- target_sequences[i, "target_start"]
  end_pos <- target_sequences[i, "target_end"]
  
  # Extract locus tag
  locus_tag <- target_sequences$LocusTag[i]
  
  # Extract subsequence from the genome
  if (gene_info_df$Strand[i] == "+") {
    subseq_list[[i]] <- as.character(subseq(genome, start_pos, end_pos))
  } else if (gene_info_df$Strand[i] == "-") {
    subseq_list[[i]] <- as.character(reverseComplement(subseq(genome, start_pos, end_pos)))
  } else {
    subseq_list[[i]] <- ""
  }
}

# Combine all subsequences into a single data frame
subseq_df <- data.frame(
  target_sequence = unlist(subseq_list),
  locus_tag = target_sequences$LocusTag
)

#-------------------------------------------------------------------------------

#PEAK FINDER

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
#All 1st replicates normalized by total reads between them
N_WT_8h_1_35_full <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/ppGpp translation initiation issue/Normalized_212_8h_1_35_full.csv")
N_ppGpp_8h_1_35_full <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/ppGpp translation initiation issue/Normalized_ppGpp_8h_1_35_full.csv")
N_WT_16h_1_35_full <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/ppGpp translation initiation issue/Normalized_16h-1_ribo_35_full.csv")
N_ppGpp_16h_1_35_full <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/ppGpp translation initiation issue/Normalized_ppGpp_16h_1_35_full.csv")

#Mup data
Normalized_filtered_2X_Spin_53_full_14 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/mup RBS finder/Normalized_filtered_2X_Spin_53_full_14.csv")
Normalized_filtered_2X_Spin_Jun_53_full_13 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/mup RBS finder/Normalized_filtered_2X_Spin_Jun_53_full_13.csv")
Normalized_filtered_2X_Spin_Jun_53_full_12 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/mup RBS finder/Normalized_filtered_2X_Spin_Jun_53_full_12.csv")
Normalized_filtered_2X_Spin_Jun_53_full_11 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/mup RBS finder/Normalized_filtered_2X_Spin_Jun_53_full_11.csv")

#FIXATION CONDITIONS
#normalized based on total number of reads
N_Cm_1 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized by total number of reads/Normalized_Cm_1.csv")
N_Cm_2 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized by total number of reads/Normalized_Cm_2.csv")
N_dsp_1 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized by total number of reads/Normalized_dsp_1.csv")
N_dsp_2 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized by total number of reads/Normalized_dsp_2.csv")
N_dsp_for_1 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized by total number of reads/Normalized_dsp+for_1.csv")
N_dsp_for_2 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized by total number of reads/Normalized_dsp+for_2.csv")
N_for_1 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized by total number of reads/Normalized_for_1.csv")
N_for_2 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized by total number of reads/Normalized_for_2.csv")
N_meth_1 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized by total number of reads/Normalized_meth_1.csv")
N_meth_2 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized by total number of reads/Normalized_meth_2.csv")
N_No_1 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized by total number of reads/Normalized_No_1.csv")
N_No_2 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized by total number of reads/Normalized_No_2.csv")
N_Tet_1 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized by total number of reads/Normalized_Tet_1.csv")
N_Tet_2 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized by total number of reads/Normalized_Tet_2.csv")


#Choose target positions (works in the same way as Sequence Finder. Takes strand directionality into consideration
#-------------
###INITIATION
target_sequences_PeakFinder <- gene_info_df %>%
  mutate(target_start = ifelse(Strand == "+", StartPosition - 15, StartPosition - 15),  #Write the numbers of the positions to select the range for sequence extraction (-15 +15)
         target_end = ifelse(Strand == "+", StartPosition + 5, StartPosition + 5)) %>% #I checked the + and - signs and this looks correct (+5 +5)
  select(locus_tag, target_start, target_end, StartPosition)

###ORF
target_sequences_PeakFinder <- gene_info_df %>%
  mutate(target_start = ifelse(Strand == "+", StartPosition + 5, StartPosition + 5),  #Write the numbers of the positions to select the range for sequence extraction (-15 +15)
         target_end = ifelse(Strand == "+", EndPosition, EndPosition)) %>% #I checked the + and - signs and this looks correct (+5 +5)
  select(locus_tag, target_start, target_end, StartPosition)

###CUSTOM (for CDS plots)
target_sequences_PeakFinder <- gene_info_df %>%
  mutate(target_start = ifelse(Strand == "+", StartPosition - 20, StartPosition - 20),  #Write the numbers of the positions to select the range for sequence extraction (-15 +15)
         target_end = ifelse(Strand == "+", StartPosition + 40, StartPosition + 40)) %>% #I checked the + and - signs and this looks correct (+5 +5)
  select(locus_tag, target_start, target_end, StartPosition)

#ASYMMETRY SCORE 5'
target_sequences_PeakFinder <- data.frame(gene_info_df$gene_length)
target_sequences_PeakFinder <- gene_info_df %>%
  mutate(target_start = ifelse(Strand == "+", StartPosition, StartPosition),  #Write the numbers of the positions to select the range for sequence extraction (-15 +15)
         target_end = ifelse(Strand == "+", floor(StartPosition + (gene_length/2)), floor(StartPosition + (gene_length/2)))) %>% #I checked the + and - signs and this looks correct (+5 +5)
  select(locus_tag, target_start, target_end, gene_length, StartPosition)

#ASYMMETRY SCORE 3'
target_sequences_PeakFinder <- data.frame(gene_info_df$gene_length)
target_sequences_PeakFinder <- gene_info_df %>%
  mutate(target_start = ifelse(Strand == "+", ceiling(StartPosition + (gene_length/2)), ceiling(StartPosition + (gene_length/2))),  #Write the numbers of the positions to select the range for sequence extraction (-15 +15)
         target_end = ifelse(Strand == "+", EndPosition, EndPosition)) %>% #I checked the + and - signs and this looks correct (+5 +5)
  select(locus_tag, target_start, target_end, gene_length, StartPosition)

#-------------

#select the input dataset from the ones loaded
input_data <- Normalized_212_8h_1_53_full
#-----
N_Tet_1 <- merge(N_Tet_1, gene_info_df, by = "gene", all = T) %>% select(-StartPosition, - EndPosition, - Sequence, - Strand, - gene_length.y) #if the previous command does not work it means the dataset does not have the locustag, so use this command by changing the name of your starting dataset
input_data <- N_Tet_1
#-----
input_data <- merge(input_data, target_sequences_PeakFinder, by = "locus_tag", all = TRUE)

#Sequences outside the ORFs do not have a locus tag assigned, but this is necessary if we want to analyse peaks in the upstream region of the ORFs
##Sort target_sequences_PeakFinder by target_start for binary search
target_sequences_PeakFinder_sorted <- target_sequences_PeakFinder[order(target_sequences_PeakFinder$target_start), ]

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
for (i in 1:nrow(input_data)) {
  if (is.na(input_data[i, "locus_tag"])) {
    result <- assign_locus_tag(input_data[i, "position"])
    input_data[i, c("locus_tag", "target_start", "target_end", "StartPosition")] <- result
  }
}

#Filter highest peaks on target regions of each locus tag and calculate the position of those peaks relative to the 1st nt of the start codon
input_data_max <- group_by(input_data, locus_tag) %>% filter(position >= target_start & position <= target_end)
input_data_max <- mutate(input_data_max, relative_position_RBS_peak = (position - StartPosition))



#merge dataset with target sequences for motif search
#input_data_max <- merge(input_data_max, subseq_df, by = "locus_tag")

#add gene length
input_data_max <- merge(input_data_max, gene_info_df, by = "locus_tag", all = TRUE) %>% select(- "gene_length.x")

#plot results (also to double check your selection range is good)
ggplot(data = input_data_max,
       mapping = aes(x = relative_position_RBS_peak)) +
  geom_bar()

#Export table to make CDS plots (if using CUSTOM target sequences)
write.csv(input_data_max, "C:/Users/aagnoli/Desktop/ppGpp_16h_CDS_plot_data.csv", row.names = FALSE) 

#sum Norm_count in target positions for each gene, then divide by the target sequence length
df_Ribo_reads_target_1 <- input_data_max %>% group_by(locus_tag) %>% summarize(Sum_Norm_count_initiation = sum(Norm_count)) %>% mutate(Translation_initiation_Ratio = Sum_Norm_count_initiation/20)
#after repeating the previous script with another set of target positions, use this
df_Ribo_reads_target_2 <- input_data_max %>% group_by(locus_tag) %>% summarize(Sum_Norm_count_ORF = sum(Norm_count)) %>% merge(gene_info_df, by = "locus_tag") %>% mutate(ORF_translation_Ratio = Sum_Norm_count_ORF/(gene_length -5)) #change this, it should be divided by gene length -5 
Merged_ratios <- merge(df_Ribo_reads_target_1, df_Ribo_reads_target_2, by = "locus_tag", all = TRUE) %>% select(c(-"StartPosition", -"EndPosition", -"Sequence", -"Strand"))
Merged_ratios <- Merged_ratios %>% mutate(log2_asymmetry_score = log2(Sum_Norm_count_ORF/Sum_Norm_count_initiation))
#Write into excel
write_xlsx(Merged_ratios, "C:/Users/aagnoli/Desktop/Log2_Tet1.xlsx")

#----------------
#write list of sequences for the motif search in kplogo
write_csv(input_data_max, file = "C:/Users/aagnoli/Desktop/RBS_peak_finder_output_2.csv")

