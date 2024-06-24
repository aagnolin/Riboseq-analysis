library(tidyverse)

#HEATMAP A, P and E SITES
Normalized_3_6_MONO_codon_occ <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/June 2024/Peak_picker results/Normalized_mapped_3_6_MONO_35_full_cutoff_1_pause_score_usage_output.csv") %>% drop_na(Codon)
Normalized_3_6_DISOME_codon_occ <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/June 2024/Peak_picker results/Normalized_mapped_3_6_DISOME_35_full_cutoff_1_pause_score_usage_output.csv") %>% drop_na(Codon)
Normalized_10_MONO_codon_occ <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/June 2024/Peak_picker results/Normalized_mapped_10_MONO_35_full_cutoff_1_pause_score_usage_output.csv") %>% drop_na(Codon)
Normalized_10_DISOME_codon_occ <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/June 2024/Peak_picker results/Normalized_mapped_10_DISOME_35_full_cutoff_1_pause_score_usage_output.csv") %>% drop_na(Codon)

Merged_codon_occ_1 <- merge(Normalized_3_6_MONO_codon_occ, Normalized_3_6_DISOME_codon_occ, by = c("Amino Acid", "Codon", "Usage"), suffixes = c("_3_6_Mono", "_3_6_Disome"))
Merged_codon_occ_2 <- merge(Normalized_10_MONO_codon_occ, Normalized_10_DISOME_codon_occ, by = c("Amino Acid", "Codon", "Usage"), suffixes = c("_10_Mono", "_10_Disome"))
Merged_codon_occ_all <- merge(Merged_codon_occ_1, Merged_codon_occ_2, by = c("Amino Acid", "Codon", "Usage")) %>% 
  pivot_longer(cols = seq(4,15),
               names_to = "Sample",
               values_to = "Normalised codon pause score") %>% 
  mutate(Sample = str_remove(Sample, "Normalised_codon_pause_score_")) %>% 
  separate(Sample, into = c("Site", "Other"), sep = "_site_", remove = FALSE) %>%
  select(-4) %>% 
  rename("Ribosomal site" = Site, "Sample" = Other) 

#NEW
ggplot(data = group_by(Merged_codon_occ_all, Sample),
       mapping = aes(x = factor(`Ribosomal site`, levels = c('A', 'P', 'E'), labels = c('A', 'P', 'E')),
                     y = Codon)) +
  geom_tile(mapping = aes(fill = `Normalised codon pause score`), linewidth = 0.5, color = "black") +
  scale_fill_gradient(low = 'white', high = 'blue') +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  #coord_flip() +
  labs(colour = "Sample",
       x = "Ribosomal site") +
  theme(plot.subtitle = element_text(face = "bold"),
        plot.caption = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.key = element_rect(linetype = "solid")) + 
  theme(axis.ticks = element_line(linewidth = 0.5)) + theme(panel.background = element_rect(fill = NA)) + 
  theme(axis.line = element_line(linetype = "solid")) + 
  theme(axis.line = element_line(linetype = "blank"),
        panel.grid.major = element_line(colour = "gray89",
                                        linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"),
        panel.background = element_rect(linetype = "solid")) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  theme(axis.text.x = element_text(vjust = 0.5) + 
          theme(axis.text.x = element_text(size = 8),
                axis.text.y = element_text(size = 8))) + theme(axis.text.x = element_text(angle = 0)) + theme(axis.text.y = element_text(size = 6)) +
  labs(fill = "Normalised \nCodon \nPause \nScore") +
  facet_grid(~factor(Sample, levels = c('3_6_Mono', '3_6_Disome', '10_Mono', '10_Disome'), labels = c('3_6_Mono', '3_6_Disome', '10_Mono', '10_Disome')))
