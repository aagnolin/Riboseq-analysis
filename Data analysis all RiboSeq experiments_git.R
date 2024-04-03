#Data analysis Riboseq experiments
#github aagnolin https://github.com/aagnolin

library(tidyverse)
library(readxl)

#IMPORT DATASETS

##mupirocin sucrose reference
Mup_reference_filtered_SAM_35_full_rep1 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Reference Mup sucrose/Mup_1_filtered_SAM_35_full.csv")
Mup_reference_filtered_SAM_35_full_rep2 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Reference Mup sucrose/Mup_2_filtered_SAM_35_full.csv")
#-----------------------------------------------------

##QUICK RIBOSEQ METHODS

##December 2021 (Chloramphenicol)
filtered_1X_Spin_35_full_December_2021 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/December 2021/alt_predict_v2_results_December_2021/filtered_1X_Spin_Dec_35_full.csv")
filtered_2X_Spin_35_full_December_2021 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/December 2021/alt_predict_v2_results_December_2021/filtered_2X_Spin_Dec_35_full.csv")
filtered_Spin_EtOH_35_full_December_2021 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/December 2021/alt_predict_v2_results_December_2021/filtered_Spin_EtOH_Dec_35_full.csv")
##January 2023 (Mupirocin)
filtered_2X_Spin_35_full_January_2023 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/January-February 2023/alt_predict_v2_results_January_2023/2X_Spin_filtered_SAM_35_full.csv")
filtered_2X_Spin_depleted_35_full_January_2023 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/January-February 2023/alt_predict_v2_results_January_2023/2X_Spin_depleted_filtered_SAM_35_full.csv")
filtered_2X_PEG_35_full_January_2023 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/January-February 2023/alt_predict_v2_results_January_2023/2X_PEG_filtered_SAM_35_full.csv")
filtered_1X_PEG_Spin_35_full_January_2023 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/January-February 2023/alt_predict_v2_results_January_2023/1X_PEG-Spin_filtered_SAM_35_full.csv")
##June 2023 (Mupirocin)
filtered_2X_Spin_35_full_June_2023 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/June 2023/alt_predict_v2_results_June_2023/filtered_2X_Spin_35_full.csv")
filtered_2X_Spin_depl_35_full_June_2023 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/June 2023/alt_predict_v2_results_June_2023/filtered_2X_Spin_depl_35_full.csv")
filtered_2X_Spin_oligo_35_full_June_2023 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/June 2023/alt_predict_v2_results_June_2023/filtered_2X_Spin_oligo_35_full.csv")
filtered_2X_PEG_35_full_June_2023 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/June 2023/alt_predict_v2_results_June_2023/filtered_2X_PEG_35_full.csv")
filtered_2X_PEG_depl_35_full_June_2023 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/June 2023/alt_predict_v2_results_June_2023/filtered_2X_PEG_depl_35_full.csv")
filtered_2X_PEG_oligo_35_full_June_2023 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/June 2023/alt_predict_v2_results_June_2023/filtered_2X_PEG_oligo_35_full.csv")
filtered_Histag_Standard_35_full_June_2023 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/June 2023/alt_predict_v2_results_June_2023/filtered_Histag_Standard_35_full.csv")
filtered_Histag_Wash_35_full_June_2023 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/June 2023/alt_predict_v2_results_June_2023/filtered_Histag_Wash_35_full.csv")
filtered_Histag_Lysate_35_full_June_2023 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/June 2023/alt_predict_v2_results_June_2023/filtered_Histag_Lysate_35_full.csv")
##July 2023 (Mupirocin)
filtered_2X_PEG_July_35_full_July_2023 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/July 2023/alt_predict_v2_results_July_2023/filtered_2X_PEG_July_35_full.csv")
filtered_1X_PEG_Spin_July_35_full_July_2023 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/July 2023/alt_predict_v2_results_July_2023/filtered_1X_PEG-Spin_July_35_full.csv")
filtered_Histag1_July_35_full_July_2023 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/July 2023/alt_predict_v2_results_July_2023/filtered_Histag1_July_35_full.csv")
filtered_Histag2_July_35_full_July_2023 <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/July 2023/alt_predict_v2_results_July_2023/filtered_Histag2_July_35_full.csv")
###Quality control of datasets
Sequence_length_distributions_all_samples <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Sequence Length Distributions - all sequencing files.xlsx")
Sequence_length_distributions_all_samples$`Date of Experiment` <- format(Sequence_length_distributions_all_samples$`Date of Experiment`, format = "%b-%y")
RNA_fractions_all_samples <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Fractions of RNAs - all sequencing files.xlsx")
RNA_fractions_all_samples$`Date of Experiment` <- format(RNA_fractions_all_samples$`Date of Experiment`, format = "%b-%y")

#-----------------------------------------------------

#16, 64 AND 160h LONG FERMENTATION (Chloramphenicol)
X16h_1_ribo <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/Riboseq datasets and plots/16h-1_Ribo_35_full.csv") 
X16h_2_ribo <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/Riboseq datasets and plots/16h-2_ribo_35_full.csv")
X64h_1_ribo <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/Riboseq datasets and plots/64h-1_ribo_35_full.csv")
X64h_2_ribo <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/Riboseq datasets and plots/64h-2_ribo_35_full.csv")
X160h_1_ribo <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/Riboseq datasets and plots/160h-1_ribo_35_full.csv")
X160h_2_ribo<-read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/Riboseq datasets and plots/160h-2_ribo_35_full.csv")

#-----------------------------------------------------

#FIXATION CONDITIONS
Cm_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Cm_1.csv")
Cm_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Cm_2.csv")
dsp_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/dsp_1.csv")
dsp_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/dsp_2.csv")
dsp_form_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/dsp+for_1.csv")
dsp_form_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/dsp+for_2.csv")
formaldehyde_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/for_1.csv")
formaldehyde_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/for_2.csv")
methanol_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/meth_1.csv")
methanol_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/meth_2.csv")
No_treatment_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/No_1.csv")
No_treatment_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/No_2.csv")
Tet_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Tet_1.csv")
Tet_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Tet_2.csv")

#normalized with script (only ORFs)
N_Cm_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Normalized_Cm_1.csv")
N_Cm_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Normalized_Cm_2.csv")
N_dsp_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Normalized_dsp_1.csv")
N_dsp_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Normalized_dsp_2.csv")
N_dsp_form_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Normalized_dsp_form_1.csv")
N_dsp_form_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Normalized_dsp_form_2.csv")
N_formaldehyde_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Normalized_for_1.csv")
N_formaldehyde_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Normalized_for_2.csv")
N_methanol_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Normalized_meth_1.csv")
N_methanol_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Normalized_meth_2.csv")
N_No_treatment_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Normalized_No_1.csv")
N_No_treatment_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Normalized_No_2.csv")
N_Tet_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Normalized_Tet_1.csv")
N_Tet_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Normalized_Tet_2.csv")

#normalized manually (ORFs-normalized but also intergenic regions)
N_Cm_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized_Cm_1.csv")
N_Cm_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized_Cm_2.csv")
N_dsp_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized_dsp_1.csv")
N_dsp_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized_dsp_2.csv")
N_dsp_form_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized_dsp_form_1.csv")
N_dsp_form_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized_dsp_form_2.csv")
N_formaldehyde_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized_formaldehyde_1.csv")
N_formaldehyde_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized_formaldehyde_2.csv")
N_methanol_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized_methanol_1.csv")
N_methanol_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized_methanol_2.csv")
N_No_treatment_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized_No_treatment_1.csv")
N_No_treatment_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized_No_treatment_2.csv")
N_Tet_1 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized_Tet_1.csv")
N_Tet_2 <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized_Tet_2.csv")
#-----------------------------------------------------

#ppGpp (Chloramphenicol)
WT_8h_1_35_full <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/ppGpp/212_8h_1_35_full.csv")
WT_8h_2_35_full <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/ppGpp/212_8h_2_35_full.csv")
WT_16h_1_35_full <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/Riboseq datasets and plots/16h-1_Ribo_35_full.csv") #This is the same dataset as 16h_1 (long fermentation experiment)
WT_16h_2_35_full <- read_csv(file = "C:/Users/aagnoli/OneDrive - UvA/Riboseq datasets and plots/16h-2_Ribo_35_full.csv") #This is the same dataset as 16h_2 (long fermentation experiment)
ppGpp_8h_1_35_full <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/ppGpp/ppGpp_8h_1_35_full.csv")
ppGpp_8h_2_35_full <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/ppGpp/ppGpp_8h_2_35_full.csv")
ppGpp_16h_1_35_full <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/ppGpp/ppGpp_16h_1_35_full.csv")
ppGpp_16h_2_35_full <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/Results 14-02-23 Yaozu ppGpp and codon exchange analysis/ppGpp/ppGpp_16h_2_35_full.csv")

##Merged replicates
Cm_Avg <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Merged replicates same nt/Cm1_Cm2_same_nt.csv") %>% 
  mutate(Avg_Norm_reads = (Norm_count.x + Norm_count.y)/2) %>%
  select(-genome.y,
         -strand.y,
         -gene.y,
         -gene_length.y,
         -offset.y,
         -in_orf_90.y,
         -count.y,
         -sequence.y,
         -Norm_count.x,
         -Norm_count.y
  )
Tet_Avg <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Merged replicates same nt/Tet1_Tet2_same_nt.csv") %>% 
  mutate(Avg_Norm_reads = (Norm_count.x + Norm_count.y)/2) %>%
  select(-genome.y,
         -strand.y,
         -gene.y,
         -gene_length.y,
         -offset.y,
         -in_orf_90.y,
         -count.y,
         -sequence.y,
         -Norm_count.x,
         -Norm_count.y
  )
dsp_Avg <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Merged replicates same nt/dsp1_dsp2_same_nt.csv") %>% 
  mutate(Avg_Norm_reads = (Norm_count.x + Norm_count.y)/2) %>%
  select(-genome.y,
         -strand.y,
         -gene.y,
         -locus_tag,
         -gene_length.y,
         -offset.y,
         -in_orf_90.y,
         -count.y,
         -sequence.y,
         -Norm_count.x,
         -Norm_count.y
  )

dsp_form_Avg <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Merged replicates same nt/dsp_form1_dsp_form2_same_nt.csv") %>% 
  mutate(Avg_Norm_reads = (Norm_count.x + Norm_count.y)/2) %>%
  select(-genome.y,
         -strand.y,
         -gene.y,
         #-locus_tag,
         -gene_length.y,
         -offset.y,
         -in_orf_90.y,
         -count.y,
         -sequence.y,
         -Norm_count.x,
         -Norm_count.y
  )

form_Avg <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Merged replicates same nt/formaldehyde1_formaldehyde2_same_nt.csv") %>% 
  mutate(Avg_Norm_reads = (Norm_count.x + Norm_count.y)/2) %>%
  select(-genome.y,
         -strand.y,
         -gene.y,
         #-locus_tag,
         -gene_length.y,
         -offset.y,
         -in_orf_90.y,
         -count.y,
         -sequence.y,
         -Norm_count.x,
         -Norm_count.y
  )

methanol_Avg <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Merged replicates same nt/methanol1_methanol2_same_nt.csv") %>% 
  mutate(Avg_Norm_reads = (Norm_count.x + Norm_count.y)/2) %>%
  select(-genome.y,
         -strand.y,
         -gene.y,
         #-locus_tag,
         -gene_length.y,
         -offset.y,
         -in_orf_90.y,
         -count.y,
         -sequence.y,
         -Norm_count.x,
         -Norm_count.y
  )

No_treatment_Avg <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Merged replicates same nt/NoTreatment1_NoTreatment2_same_nt.csv") %>% 
  mutate(Avg_Norm_reads = (Norm_count.x + Norm_count.y)/2) %>%
  select(-genome.y,
         -strand.y,
         -gene.y,
         #-locus_tag,
         -gene_length.y,
         -offset.y,
         -in_orf_90.y,
         -count.y,
         -sequence.y,
         -Norm_count.x,
         -Norm_count.y
  )

WT_8h_Avg <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis ppGpp/WT_8h_rep1_vs_WT_8h_rep2.csv") %>% 
  mutate(Avg_Norm_reads = (Norm_count.x + Norm_count.y)/2) %>%
  select(-genome.y,
         -strand.y,
         -gene.y,
         -gene_length.y,
         -offset.y,
         -in_orf_90.y,
         -count.y,
         -sequence.y,
         #-Pause_codon.x,
         #-Pause_codon.y,
         -Norm_count.x,
         -Norm_count.y
)

WT_16h_Avg <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis ppGpp/WT_16h_rep2_vs_WT_16h_rep1.csv") %>% 
  mutate(Avg_Norm_reads = (Norm_count.x + Norm_count.y)/2) %>%
  select(-genome.y,
         -strand.y,
         -gene.y,
         -gene_length.y,
         -offset.y,
         -in_orf_90.y,
         -count.y,
         -sequence.y,
         #-Pause_codon.x,
         #-Pause_codon.y,
         -Norm_count.x,
         -Norm_count.y
  )
   
ppGpp_8h_Avg <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis ppGpp/ppGpp_8h_rep1_vs_ppGpp_8h_rep2.csv") %>% 
  mutate(Avg_Norm_reads = (Norm_count.x + Norm_count.y)/2) %>%
  select(-genome.y,
         -strand.y,
         -gene.y,
         -gene_length.y,
         -offset.y,
         -in_orf_90.y,
         -count.y,
         -sequence.y,
         #-Pause_codon.x,
         #-Pause_codon.y,
         -Norm_count.x,
         -Norm_count.y
  )
    
ppGpp_16h_Avg <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis ppGpp/ppGpp_16h_rep2_vs_ppGpp_16h_rep1.csv") %>% 
  mutate(Avg_Norm_reads = (Norm_count.x + Norm_count.y)/2) %>%
  select(-genome.y,
         -strand.y,
         -gene.y,
         -gene_length.y,
         -offset.y,
         -in_orf_90.y,
         -count.y,
         -sequence.y,
         #-Pause_codon.x,
         #-Pause_codon.y,
         -Norm_count.x,
         -Norm_count.y
  )  

x2X_Spin_Avg <- read_csv("C:/Users/aagnoli/OneDrive - UvA/WP 1 - RiboSeq/Manuscript Spin-PEG/Data peak comparisons with tls/Avg_2X_Spin_Jun_and_Jan.csv") %>% 
  mutate(Avg_Norm_reads = (Norm_count.x + Norm_count.y)/2) %>%
  select(-genome.y,
         -strand.y,
         -gene.y,
         -gene_length.y,
         -offset.y,
         -in_orf_90.y,
         -count.y,
         -sequence.y,
         #-Pause_codon.x,
         #-Pause_codon.y,
         -Norm_count.x,
         -Norm_count.y
  )  
x2X_PEG_Avg <- read_csv("C:/Users/aagnoli/OneDrive - UvA/WP 1 - RiboSeq/Manuscript Spin-PEG/Data peak comparisons with tls/Avg_2X_PEG_Jul_and_Jan.csv") %>% 
  mutate(Avg_Norm_reads = (Norm_count.x + Norm_count.y)/2) %>%
  select(-genome.y,
         -strand.y,
         -gene.y,
         -gene_length.y,
         -offset.y,
         -in_orf_90.y,
         -count.y,
         -sequence.y,
         #-Pause_codon.x,
         #-Pause_codon.y,
         -Norm_count.x,
         -Norm_count.y
  )  
x1X_PEG_Spin_Avg <- read_csv("C:/Users/aagnoli/OneDrive - UvA/WP 1 - RiboSeq/Manuscript Spin-PEG/Data peak comparisons with tls/Avg_1X_PEG-Spin_Jan_and_Jul.csv") %>% 
  mutate(Avg_Norm_reads = (Norm_count.x + Norm_count.y)/2) %>%
  select(-genome.y,
         -strand.y,
         -gene.y,
         -gene_length.y,
         -offset.y,
         -in_orf_90.y,
         -count.y,
         -sequence.y,
         #-Pause_codon.x,
         #-Pause_codon.y,
         -Norm_count.x,
         -Norm_count.y
  )  
Mup_Avg <- read_csv("C:/Users/aagnoli/OneDrive - UvA/WP 1 - RiboSeq/Manuscript Spin-PEG/Data peak comparisons with tls/Avg_Mup_1_and_2.csv") %>% 
  mutate(Avg_Norm_reads = (Norm_count.x + Norm_count.y)/2) %>%
  select(-genome.y,
         -strand.y,
         -gene.y,
         -gene_length.y,
         -offset.y,
         -in_orf_90.y,
         -count.y,
         -sequence.y,
         #-Pause_codon.x,
         #-Pause_codon.y,
         -Norm_count.x,
         -Norm_count.y
  )  
#Normalize these averaged datasets
riboseq_1 <- 
riboseq_2 <- 

sum(riboseq_1$Avg_Norm_reads)
sum(riboseq_2$Avg_Norm_reads)

#Normalize datasets
Normalized_riboseq_1 <- mutate(riboseq_1, Norm_count = riboseq_1$Avg_Norm_reads*(sum(riboseq_1$Avg_Norm_reads)/sum(riboseq_1$Avg_Norm_reads))) #from now on, use the Norm_count column for further analysis
Normalized_riboseq_2 <- mutate(riboseq_2, Norm_count = riboseq_2$Avg_Norm_reads*(sum(riboseq_1$Avg_Norm_reads)/sum(riboseq_2$Avg_Norm_reads))) #from now on, use the Norm_count column for further analysis

#-----------------------------------------------------

#Paste here the name of the 2 files to compare (set file with higher total number of sequences as riboseq_1 for normalization at later steps)
riboseq_1 <- form_Avg
riboseq_2 <- dsp_form_Avg

#===========================
#FURTHER CLEANING OF alt_predict DATASETS
##some tRNAs and rRNAs still remain in the tables after removal with bowtie2. Moreover, small RNAs like ssrA are not removed with Bowtie2.
riboseq_1 <- filter(riboseq_1, !grepl("^BSU_", locus_tag))
riboseq_2 <- filter(riboseq_2, !grepl("^BSU_", locus_tag))
##use these below if datasets do not have the locus_tag column
riboseq_1 <- filter(riboseq_1, is.na(gene) | !grepl("-", gene) & gene != 'ssrA' & gene != "scr" & gene != "rnpB")
riboseq_2 <- filter(riboseq_2, is.na(gene) | !grepl("-", gene) & gene != 'ssrA' & gene != "scr" & gene != "rnpB")
#===========================

#Check which dataset has the highest number of reads
sum(riboseq_1$count)
sum(riboseq_2$count)

#Normalize datasets
Normalized_riboseq_1 <- mutate(riboseq_1, Norm_count = riboseq_1$count*(sum(riboseq_1$count)/sum(riboseq_1$count)))
Normalized_riboseq_2 <- mutate(riboseq_2, Norm_count = riboseq_2$count*(sum(riboseq_1$count)/sum(riboseq_2$count)))

#==========================================================
#GENE COVERAGE COMPARISONS

#Filter for reads in ORFs
riboseq_1_ORFs <- riboseq_1 %>% filter(in_orf_90 == TRUE) #pass count>10 or another cutoff in the filter function if required
sum(riboseq_1_ORFs$count)
riboseq_2_ORFs <- riboseq_2 %>% filter(in_orf_90 == TRUE) #pass count>10 or another cutoff in the filter function if required
sum(riboseq_2_ORFs$count)

#Normalize datasets (only considering ORFs)

Normalized_riboseq_1_ORF <- mutate(riboseq_1_ORFs, Norm_count = riboseq_1_ORFs$count*(sum(riboseq_1_ORFs$count)/sum(riboseq_1_ORFs$count)))
Normalized_riboseq_2_ORF <- mutate(riboseq_2_ORFs, Norm_count = riboseq_2_ORFs$count*(sum(riboseq_1_ORFs$count)/sum(riboseq_2_ORFs$count)))

#Calculate coverage (count/gene length) per gene

Gene_coverage_Normalized_riboseq_1 <- Normalized_riboseq_1_ORF %>% group_by(gene, gene_length) %>% summarise(Norm_tot_counts = sum(Norm_count))
Gene_coverage_Normalized_riboseq_1 <- Gene_coverage_Normalized_riboseq_1 %>% mutate(gene_coverage = Norm_tot_counts/gene_length)
Gene_coverage_Normalized_riboseq_2 <- Normalized_riboseq_2_ORF %>% group_by(gene, gene_length) %>% summarise(Norm_tot_counts = sum(Norm_count))
Gene_coverage_Normalized_riboseq_2 <- Gene_coverage_Normalized_riboseq_2 %>% mutate(gene_coverage = Norm_tot_counts/gene_length)

#merge datasets by gene
Merged_gene_coverage <- merge(Gene_coverage_Normalized_riboseq_1, Gene_coverage_Normalized_riboseq_2, by = 'gene', suffixes = c('1', '2'))
cutoff_gene_coverage <- 1 #change cutoff for count if required
Merged_gene_coverage <- Merged_gene_coverage %>% filter(Norm_tot_counts1 > cutoff_gene_coverage & Norm_tot_counts2 > cutoff_gene_coverage)

#linear model correlation (R-squared)
linear_model_gene_coverage <- lm(Merged_gene_coverage$gene_coverage1 ~ Merged_gene_coverage$gene_coverage2)
r_squared_gene_coverage <- summary(linear_model_gene_coverage)$r.squared
coeff_gene_coverage <-coefficients(linear_model_gene_coverage)           
intercept_gene_coverage <-coeff_gene_coverage[1] 
slope_gene_coverage <- coeff_gene_coverage[2] 

#Plot gene coverage of different samples
ggplot(data = Merged_gene_coverage,
       mapping = aes(
         x = gene_coverage1,
         y = gene_coverage2
       )) +
  geom_point(size = 1, alpha = 0.3) +
  geom_abline(intercept = intercept_gene_coverage, slope = slope_gene_coverage, color="red", linewidth = 0.8) +
  #geom_smooth(method = 'lm', se = FALSE, color = 'red', linewidth = 0.8) +
  theme_bw() + 
  theme(panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank")) + 
  labs(x = "Gene coverage 2X Spin", y = "Gene coverage Sucrose") +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') + 
  annotate("text", x = min(Merged_gene_coverage$gene_coverage1),
           y = max(Merged_gene_coverage$gene_coverage2), 
           label = paste("R-squared =", round(r_squared_gene_coverage, 3)), 
           hjust = 0, vjust = 1, size = 4) +
  labs(colour = "Sample") +
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
        legend.key.height= unit(0.05, 'cm'),
        legend.key.width= unit(0.05, 'cm')) +
  theme(plot.subtitle = element_text(face = "bold"),
        plot.caption = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.key = element_rect(linetype = "solid")) + 
  theme(axis.ticks = element_line(size = 0.5)) + theme(panel.background = element_rect(fill = NA)) + 
  theme(axis.line = element_line(linetype = "solid")) + theme(axis.line = element_line(linetype = "blank"),
                                                              panel.grid.major = element_line(colour = "gray89",
                                                                                              linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"),
                                                              panel.background = element_rect(linetype = "solid"))

#Pearson correlation
cor(Merged_gene_coverage$gene_coverage1, Merged_gene_coverage$gene_coverage2, method = "pearson")

#====================================================================

#COMPARISON OF RIBOSOME PROFILES

#Profiles of amyM - single sample (change data in plot from 1 to 2 to see the first or second dataset)
ggplot(data = Normalized_riboseq_1_ORF  %>% dplyr::filter(gene == "cspB"),
       mapping = aes(
         x = position)) +
  geom_col(aes(y = Norm_count), fill = 'darkgreen', alpha = 0.5) +
  scale_y_continuous(limits = c(0,2000)) +
  theme_bw() + theme(panel.grid.major = element_line(linetype = "blank"),
                     panel.grid.minor = element_line(linetype = "blank"),
                     axis.title = element_text(face = "bold"),
                     axis.text = element_text(face = "bold")) +labs(x = "Position", y = "count")


#Profiles of amyM - two samples

# Create a full range of positions combining both datasets
all_positions <- union(Normalized_riboseq_1$position, Normalized_riboseq_2$position)

# Create a dataframe with all positions
all_data <- data.frame(position = all_positions)

# Left join each dataset separately
all_data <- all_data %>%
  left_join(Normalized_riboseq_1, by = c('position')) %>%
  left_join(Normalized_riboseq_2, by = c('position'),
            suffix = c('_1', '_2'))

# Filter data for the gene of interest
gene_of_interest <- "glnA" # or gene name
filtered_data <- all_data %>% 
  ##filter(locus_tag_1 == gene_of_interest | locus_tag_2 == gene_of_interest) #use this if you use pass the locus tag to the gene_of_interest vector
  filter(gene_1 == gene_of_interest | gene_2 == gene_of_interest) #use this if you use the gene name on the gene_of_interest

# Plot the data
Plot <- ggplot(data = filtered_data,
               mapping = aes(x = position)) +
  geom_col(aes(y = coalesce(Norm_count_1, 0)), fill = 'navy', size = 1, alpha = 0.5) +
  geom_col(aes(y = coalesce(Norm_count_2, 0)), fill = 'orange', size = 1, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold")) +
  labs(x = "Position", y = "Norm_count") +
  scale_y_continuous(limits = c(0,2000))

#Show plot
Plot

# [Interactive plot]
library(plotly)

# Convert the ggplot to an interactive plot using plotly
ggplotly(Plot)

#Shiny app for comparison on all genes

# Load required libraries
library(dplyr)
library(ggplot2)
library(plotly)
library(shiny)
library(shinyjs)

# Install the shinyWidgets package if you haven't already
# install.packages("shinyWidgets")

# Load the shinyWidgets library
library(shinyWidgets)

# Define the end coordinate in the merged data
max_position <- max(max(Normalized_riboseq_1$position), max(Normalized_riboseq_2$position))

# Define the UI for the Shiny app
ui <- fluidPage(
  useShinyjs(),  # Initialize shinyjs
  titlePanel("Gene Expression Data Visualization"),
  sidebarLayout(
    sidebarPanel(
      # Typable text input with a dropdown menu
      pickerInput(
        inputId = "gene_selection",
        label = "Select or Type Gene Name:",
        choices = unique(c(Normalized_riboseq_1$gene, Normalized_riboseq_2$gene)),
        options = list(
          `actions-box` = TRUE,
          `live-search` = TRUE
        )
      ),
      # Input fields for start and end positions (disabled by default)
      numericInput("start_position", "Start Position:", min = 1, max = max_position, value = 1),
      numericInput("end_position", "End Position:", min = 1, max = max_position, value = max_position),
      # Action button to apply position filter
      actionButton("apply_position", "Apply Position Filter")
    ),
    mainPanel(
      # Plot output
      plotlyOutput("gene_plot")
    )
  )
)

# Define the server logic
server <- function(input, output, session) {
  # Enable/disable position inputs based on gene selection
  observe({
    if (!is.null(input$gene_selection)) {
      shinyjs::enable("start_position")
      shinyjs::enable("end_position")
    } else {
      shinyjs::disable("start_position")
      shinyjs::disable("end_position")
    }
  })
  
  # Update start and end positions based on gene selection
  observeEvent(input$gene_selection, {
    gene_of_interest <- input$gene_selection
    
    # Find the start and end positions for the selected gene
    gene_data <- bind_rows(Normalized_riboseq_1, Normalized_riboseq_2) %>%
      filter(gene == gene_of_interest)
    
    if (nrow(gene_data) > 0) {
      start_position <- min(gene_data$position)
      end_position <- max(gene_data$position)
      
      # Update the start and end position inputs
      updateNumericInput(session, "start_position", value = start_position)
      updateNumericInput(session, "end_position", value = end_position)
    }
  })
  
  # Reactive function to filter data based on the selected gene and positions
  filtered_data <- reactive({
    gene_of_interest <- input$gene_selection
    start_position <- input$start_position
    end_position <- input$end_position
    
    # Filter data by gene
    all_positions <- union(Normalized_riboseq_1$position, Normalized_riboseq_2$position)
    all_data <- data.frame(position = all_positions)
    all_data <- all_data %>%
      left_join(Normalized_riboseq_1, by = c('position')) %>%
      left_join(Normalized_riboseq_2, by = c('position'), suffix = c('_1', '_2'))
    filtered_data <- all_data %>% 
      filter(gene_1 == gene_of_interest | gene_2 == gene_of_interest)
    
    # Filter data by position range
    filtered_data <- filtered_data %>%
      filter(position >= start_position & position <= end_position)
    
    return(filtered_data)
  })
  
  # Create an interactive plot based on the selected gene and positions
  output$gene_plot <- renderPlotly({
    gene_data <- filtered_data()
    Plot <- ggplot(data = gene_data,
                   mapping = aes(x = position)) +
      geom_col(aes(y = coalesce(Norm_count_1, 0)), fill = 'navy', size = 1, alpha = 0.5) +
      geom_col(aes(y = coalesce(Norm_count_2, 0)), fill = 'orange', size = 1, alpha = 0.5) +
      theme_bw() +
      theme(panel.grid.major = element_line(linetype = "blank"),
            panel.grid.minor = element_line(linetype = "blank"),
            axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold")) +
      labs(x = "Position", y = "Norm_count")
    
    # Convert the ggplot to a plotly object
    ggplotly(Plot)
  })
  
  # Filter data based on position filter button click
  observeEvent(input$apply_position, {
    # Trigger data filtering when the button is clicked
    filtered_data()
  })
}

# Run the Shiny app
shinyApp(ui, server)

#CLUSTERING METHOD (LUCAS LIGTVOET) - modified to adapt to compare two samples
#NOTE: DIFFERENTLY FROM BEFORE, THE DATA USED TO CLUSTER AND CREATE THE PLOTS IS NOT NORMALIZED

X1<- Normalized_riboseq_1
X2<- Normalized_riboseq_2

#merge all data frames (this only works for 6 samples)
df_list <- list(
  X1[!is.na(X1$gene), c("gene", "position", "Norm_count")], 
  X2[!is.na(X2$gene), c("gene", "position", "Norm_count")]
)

X.all<-df_list %>% reduce(full_join, by=c("gene","position"))
colnames(X.all)<-c("gene","position" , "sample1", "sample2")
X.all[is.na(X.all)]<-0
X.all.20<- X.all[which(X.all$sample1 >= 20 | X.all$sample2 >= 20),]
X.all.20<-as.data.frame(X.all.20)


#Note that this data is not normalized and not log transformed

#CLustering Input
Data<- X.all.20       #type name of dataframe here (If you used code above don't change it)
# Data should be gene first column, position second column and samples starting from column 3
d <- 3 # distance of binning (I used 1 and 3)


#PEAK CLUSTERING#
get.peak.clust <- function(D,graph=TRUE){
  D <- as.matrix(D[,2,drop=FALSE])
  rownames(D) <- paste0("peak",c(1:nrow(D)))
  
  D <- dist(D)
  hc <- hclust(D, method="single", members = NULL)       #nearest neighbour clustering: very suitable in this context
  #if h=1, then all peaks that have a distance of one bp from eachother are combined
  if(graph){plot(hc)}
  return(hc)
}

get.clustered.peaks <- function(x,hc,distance=2){
  library(dplyr)
  
  D <- x[,-c(1,2)]
  C <- data.frame(start=x[,2], end=x[,2])  #1-based system
  
  #Binning procedure#
  cl <- cutree(hc, h = distance)
  # Group by mean using dplyr
  DD <- data.frame(cl=cl,D)
  Counts <- dplyr::group_by(DD,cl) %>% dplyr::summarise(across(everything(), sum))
  
  CC <- data.frame(cl=cl,C)
  Coordinates <- dplyr::group_by(CC,cl) %>% dplyr::summarise(start=min(start),end=max(end))
  return(data.frame(Coordinates[,-1],Counts[,-1]))
}


#selecting genes with a minimal of 2 values

length<-NULL
for (i in 1:length(unique(Data$gene))) {
  
  length<-c(length,
            length(which(Data$gene %in% unique(Data$gene)[i]))
  )}
genes.names<-unique(Data$gene)[which(length >=2)]


result<-NULL
for (i in 1:length(genes.names)) {
  
  
  DF<-Data[which(Data$gene == genes.names[i]),]
  colnames(DF)<- c("gene", "coordinate_position",	"sample1",	"sample2")
  Groups<-c(1,1,2,2,3,3)
  hc <- get.peak.clust(DF,graph = FALSE)
  clust.peaks<-get.clustered.peaks(DF,hc=hc,distance=d)
  clust.peak<-cbind(gene=rep(genes.names[i], nrow(clust.peaks)),clust.peaks)
  result<-rbind(result,clust.peak)
  
}
result<-result[order(result$start, decreasing = FALSE),]
View(result)
Binned.data<-result

#which sample you want to look at (column number)(column 4 is sample 1)
sample_1 <- for (m1 in 4) {
  
  #for changing the gene you want to look at
  for (i1 in "amyI") {
    
    
    par(mfrow=c(1,1))
    
    Distance1<-NULL
    for (k1 in 1:nrow(Binned.data[which(Binned.data$gene == i1),])) {
      Distance1<-c(Distance1,Binned.data[which(Binned.data$gene == i1),"start"][k1+1]-Binned.data[which(Binned.data$gene == i1),"end"][k1])
    }
  }}  

#which sample you want to look at (column number)(column 5 is sample 2)
sample_2 <- for (m2 in 5) {
  
  #for changing the gene you want to look at
  for (i2 in "amyI") {
    
    
    par(mfrow=c(1,1))
    
    Distance2<-NULL
    for (k2 in 1:nrow(Binned.data[which(Binned.data$gene == i2),])) {
      Distance2<-c(Distance2,Binned.data[which(Binned.data$gene == i2),"start"][k2+1]-Binned.data[which(Binned.data$gene == i2),"end"][k2])
    }
  }}  

#Distance

P1 <- barplot(names=Binned.data[which(Binned.data$gene == i1),"start"],
              height=Binned.data[which(Binned.data$gene == i1),m1],
              col = "navy", 
              width = c((Binned.data[which(Binned.data$gene == i1),"end"]-Binned.data[which(Binned.data$gene == i1),"start"])+1),
              space = c(0,((Distance2[1:(length(Distance2)-1)]-1)/mean((Binned.data[which(Binned.data$gene == i1),"end"]-Binned.data[which(Binned.data$gene == i1),"start"])+1))))

P2 <- barplot(names=Binned.data[which(Binned.data$gene == i2),"start"],
              height=Binned.data[which(Binned.data$gene == i2),m2],
              col = adjustcolor("orange", alpha.f = 0.5),
              width = c((Binned.data[which(Binned.data$gene == i2),"end"]-Binned.data[which(Binned.data$gene == i2),"start"])+1),
              space = c(0,((Distance2[1:(length(Distance2)-1)]-1)/mean((Binned.data[which(Binned.data$gene == i2),"end"]-Binned.data[which(Binned.data$gene == i2),"start"])+1))),
              add = TRUE)

#==================================================

#SCATTER PLOTS FOR COMPARISON (nt and codon level)

##comparison scatter plots on nt level --> NOTE: This includes peaks in the entire genome, inside and outside ORFs

Merged_same_position <- merge(riboseq_1, riboseq_2, 'position')

cutoff <- 1 #cutoff of number of reads per peak for plotting. Change if desired
Merged_same_position_cutoff <- filter(Merged_same_position, Avg_Norm_reads.x >= cutoff & Avg_Norm_reads.y >= cutoff)

#linear regression

ggplot(data = Merged_same_position_cutoff,
       mapping = aes(x = Avg_Norm_reads.x,
                     y = Avg_Norm_reads.y)) +
  geom_point() +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  geom_smooth(method = 'lm') +
  geom_abline(col = 'red')

linear_model_merged_same_position <- lm(Merged_same_position_cutoff$Norm_count.x ~ Merged_same_position_cutoff$Norm_count.y)
summary(linear_model_merged_same_position)$r.squared

cor(Merged_same_position_cutoff$Norm_count.x, Merged_same_position_cutoff$Norm_count.y, method = "pearson")

write_csv(Merged_same_position_cutoff, file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Merged replicates same nt/Tet1_Tet2_same_nt.csv")

#total least squares regression (mean centered data)
data = cbind(Merged_same_position_cutoff$Avg_Norm_reads.x, Merged_same_position_cutoff$Avg_Norm_reads.y)
logData = log10(data)
meanX = colMeans(logData)
mcX = sweep(logData,2,meanX,FUN="-")
pca <- svd(mcX)
scores = pca$u%*%diag(pca$d)
loads = pca$v

# calculate the estimated centered data back using only the first component.. and then add the columns means again to have a proper estimate

Xfit = sweep(scores[,1]%*%t(loads[,1]),2,meanX,FUN="+")

# Calculate the explained variance (as some kind of ‘quality’ measure.. 100=perfect fit

sum((Xfit^2)/sum(logData^2)*100)

#Calculate the explained variance (as some kind of ‘quality’ measure.. 100=perfect fit)
explained_variance <- sum((Xfit^2)/sum(logData^2)*100)

#Calculate the explained variance adjusted
XfitMc <- scores[,1]%*%t(loads[,1])
sum(XfitMc**2)/sum(mcX**2)*100

adjusted_explained_variance <- sum(XfitMc**2)/sum(mcX**2)*100

# Extracting coefficients of the green line
coefficients_tls_line <- coef(lm(Xfit[, 2] ~ Xfit[, 1] + 1))
equation <- paste("y =", round(coefficients_tls_line[2], 2), "* x +", round(coefficients_tls_line[1], 2))

ggplot(data = as.data.frame(logData),
       mapping = aes(x = V1,
                     y = V2)) +
  geom_point() +
  geom_abline(intercept = coefficients_tls_line[1], slope = coefficients_tls_line[2], col = "darkgreen", linewidth = 1) +
  #geom_abline(intercept = 0, slope = 1, col = "red", lwd = 1) +
  annotate("text", x = 1, y = max(logData[, 2]), label = paste0("Explained variance = ", round(adjusted_explained_variance, 2), "%"), hjust = 0.05, vjust = 1, size = 4) +
  theme_bw() +
  labs(x = "Counts Sucrose",
       y = "Counts 1X PEG-Spin") +
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

Output_tls_nt_replicates <- as.data.frame(logData) %>% rename(Avg_Norm_counts.x = V1, Avg_Norm_counts.y = V2)
Output_tls_nt_replicates <- cbind(Output_tls_nt_replicates, adjusted_explained_variance, equation)
write.csv(Output_tls_nt_replicates, file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Normalized only ORFs script, no intergenetic/Merged replicates same nt/tls_same_nt/dsp_avg_vs_NoTreatment_avg_tls.csv", row.names = F)


## 1) comparison scatter plots on codon level --> merge codons based on 0,+1,+2 positions (0,+1,+2 being the nucleotides in a codon) --> NOTE: This ONLY includes peaks inside ORFs

Merged_same_position_codon <- merge(Normalized_riboseq_1, Normalized_riboseq_2, 'position') %>%  mutate(codon_position = ceiling(offset.x/3))
Merged_same_position_codon <- Merged_same_position_codon %>% select(-genome.y,
                                                                    -strand.y,
                                                                    -gene.y,
                                                                    #-locus_tag.y,
                                                                    -gene_length.y,
                                                                    -offset.y,
                                                                    -in_orf_90.y,
                                                                    -count.y,
                                                                    -sequence.y
)


###Sum the number of reads in each codon
Merged_codon_summed_reads <- Merged_same_position_codon %>% 
  group_by(gene.x, codon_position) %>% summarise(sum_reads_per_codon.x = sum(Norm_count.x),
                                                 sum_reads_per_codon.y = sum(Norm_count.y)) %>%  
  filter(!is.na(gene.x))

cutoff <- 1 #cutoff of number of reads per peak for plotting. Change if desired
Merged_codon_summed_reads_cutoff <- filter(Merged_codon_summed_reads, sum_reads_per_codon.x >= cutoff & sum_reads_per_codon.y >= cutoff)

###compute linear model for trend line (log10-transformed to make a line for plotting)
linear_model_trendline <- lm(log10(Merged_codon_summed_reads_cutoff$sum_reads_per_codon.x) ~ log10(Merged_codon_summed_reads_cutoff$sum_reads_per_codon.y))
coeff <-coefficients(linear_model_trendline)           
intercept <-coeff[1] 
slope <- coeff[2] 

###plot data with linear regression

ggplot(data = Merged_codon_summed_reads_cutoff,
       mapping = aes(x = sum_reads_per_codon.x,
                     y = sum_reads_per_codon.y)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_abline(intercept = intercept, slope = slope, color="red") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10")

linear_model_Merged_codon_summed_reads <- lm(Merged_codon_summed_reads_cutoff$sum_reads_per_codon.x ~ Merged_codon_summed_reads_cutoff$sum_reads_per_codon.y)
summary(linear_model_Merged_codon_summed_reads)$r.squared

cor(Merged_codon_summed_reads_cutoff$sum_reads_per_codon.x, Merged_codon_summed_reads_cutoff$sum_reads_per_codon.y, method = "pearson")

write_csv(Merged_codon_summed_reads_cutoff, file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Results/Merged by codon/Replicates codon/form2_and_form1.csv")

#total least squares regression (mean centered data)

data_codon = cbind(Merged_codon_summed_reads_cutoff$sum_reads_per_codon.x, Merged_codon_summed_reads_cutoff$sum_reads_per_codon.y)
logData_codon = log10(data_codon)
meanX_codon = colMeans(logData_codon)
mcX_codon = sweep(logData_codon,2,meanX_codon,FUN="-")
pca_codon <- svd(mcX_codon)
scores_codon = pca_codon$u%*%diag(pca_codon$d)
loads_codon = pca_codon$v

# calculate the estimated centered data back using only the first component.. and then add the columns means again to have a proper estimate

Xfit_codon = sweep(scores_codon[,1]%*%t(loads_codon[,1]),2,meanX_codon,FUN="+")

#Calculate the explained variance (as some kind of ‘quality’ measure.. 100=perfect fit)
explained_variance_codon <- sum((Xfit_codon^2)/sum(logData_codon^2)*100)

# Extracting coefficients of the green line
coefficients_tls_line_codon <- coef(lm(Xfit_codon[, 2] ~ Xfit_codon[, 1] + 1))
equation_codon <- paste("y =", round(coefficients_tls_line_codon[2], 2), "* x +", round(coefficients_tls_line_codon[1], 2))

ggplot(data = as.data.frame(logData_codon),
       mapping = aes(x = V1,
                     y = V2)) +
  geom_point() +
  geom_abline(intercept = coefficients_tls_line_codon[1], slope = coefficients_tls_line_codon[2], col = "green", linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, col = "red", lwd = 1) +
  annotate("text", x = 1, y = max(logData_codon[, 2]), label = paste("Green: ", equation_codon), hjust = 1, vjust = 1, size = 4)

#-----------------------
## 2) comparison scatter plots on codon level - Different codon merge method --> merge codons based on -1,0,+1 positions (0,+1,+2 being the nucleotides in a codon) --> NOTE: This ONLY includes peaks inside ORFs

Merged_same_position_codon_floor <- merge(Normalized_riboseq_1, Normalized_riboseq_2, 'position') %>%  mutate(codon_position = floor(offset.x.x/3))
Merged_same_position_codon_floor <- Merged_same_position_codon_floor %>% select(-genome.y,
                                                                                -strand.y,
                                                                                -gene.y,
                                                                                #-locus_tag.y,
                                                                                -gene_length.y,
                                                                                -offset.y,
                                                                                -in_orf_90.y,
                                                                                -count.y,
                                                                                -sequence.y
)
###Sum the number of reads in each codon
Merged_codon_summed_reads_floor <- Merged_same_position_codon_floor %>% 
  group_by(gene.x.x, codon_position) %>% summarise(sum_reads_per_codon.x = sum(Norm_count.x),
                                                 sum_reads_per_codon.y = sum(Norm_count.y)) %>%  
  filter(!is.na(gene.x.x))

cutoff <- 10 #cutoff of number of reads per peak for plotting. Change if desired
Merged_codon_summed_reads_cutoff_floor <- filter(Merged_codon_summed_reads_floor, sum_reads_per_codon.x >= cutoff & sum_reads_per_codon.y >= cutoff)

###compute linear model for trend line (log10-transformed to make a line for plotting)
linear_model_trendline_floor <- lm(log10(Merged_codon_summed_reads_cutoff_floor$sum_reads_per_codon.x) ~ log10(Merged_codon_summed_reads_cutoff_floor$sum_reads_per_codon.y))
coeff_floor <-coefficients(linear_model_trendline_floor)           
intercept_floor <-coeff_floor[1] 
slope_floor <- coeff_floor[2] 

###compute linear model for r-squared coefficient (NOT log-transformed since it is for data comparison between the two samples)
linear_model_Merged_codon_summed_reads_floor <- lm(Merged_codon_summed_reads_cutoff_floor$sum_reads_per_codon.x ~ Merged_codon_summed_reads_cutoff_floor$sum_reads_per_codon.y)
r_squared_Merged_codon_summed_reads_floor <- summary(linear_model_Merged_codon_summed_reads_floor)$r.squared

ggplot(data = Merged_codon_summed_reads_cutoff_floor,
       mapping = aes(x = sum_reads_per_codon.x,
                     y = sum_reads_per_codon.y)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_abline(intercept = intercept_floor, slope = slope_floor, color="red") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  geom_smooth(method = 'lm') +
  geom_abline() +
  annotate("text", x = min(Merged_codon_summed_reads_cutoff_floor$sum_reads_per_codon.x),
           y = max(Merged_codon_summed_reads_cutoff_floor$sum_reads_per_codon.x), 
           label = paste("R-squared =", round(r_squared_Merged_codon_summed_reads_floor, 3)), 
           hjust = 0, vjust = 1, size = 4)

###compute linear model for r-squared coefficient (NOT log-transformed since it is for data comparison between the two samples)
linear_model_Merged_codon_summed_reads_floor <- lm(Merged_codon_summed_reads_cutoff_floor$sum_reads_per_codon.x ~ Merged_codon_summed_reads_cutoff_floor$sum_reads_per_codon.y)
summary(linear_model_Merged_codon_summed_reads_floor)$r.squared

###compute Pearson correlation coefficient (NOT log-transformed since it is for data comparison between the two samples)
cor(Merged_codon_summed_reads_cutoff_floor$sum_reads_per_codon.x, Merged_codon_summed_reads_cutoff_floor$sum_reads_per_codon.y, method = "pearson")

write_csv(Merged_codon_summed_reads_cutoff_floor, file = "C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Global analysis fixation conditions/Results/Averaged replicates with codon comparisons/Codon -1 results/No_treatment_avg_vs_form+dsp_avg_codon-1.csv")

#total least squares regression (mean centered data)

data_codon_floor = cbind(Merged_codon_summed_reads_cutoff_floor$sum_reads_per_codon.x, Merged_codon_summed_reads_cutoff_floor$sum_reads_per_codon.y)
logData_codon_floor = log10(data_codon_floor)
meanX_codon_floor = colMeans(logData_codon_floor)
mcX_codon_floor = sweep(logData_codon_floor,2,meanX_codon_floor,FUN="-")
pca_codon_floor <- svd(mcX_codon_floor)
scores_codon_floor = pca_codon_floor$u%*%diag(pca_codon_floor$d)
loads_codon_floor = pca_codon_floor$v

# calculate the estimated centered data back using only the first component.. and then add the columns means again to have a proper estimate

Xfit_codon_floor = sweep(scores_codon_floor[,1]%*%t(loads_codon_floor[,1]),2,meanX_codon_floor,FUN="+")

# Calculate the explained variance (as some kind of ‘quality’ measure.. 100=perfect fit

sum((Xfit_codon_floor^2)/sum(logData_codon_floor^2)*100)

#Calculate the explained variance (as some kind of ‘quality’ measure.. 100=perfect fit)
explained_variance_floor <- sum((Xfit_codon_floor^2)/sum(logData_codon_floor^2)*100)

#Calculate the explained variance adjusted
XfitMc <- scores_codon_floor[,1]%*%t(loads_codon_floor[,1])
sum(XfitMc**2)/sum(mcX_codon_floor**2)*100

adjusted_explained_variance_floor <- sum(XfitMc**2)/sum(mcX_codon_floor**2)*100
# Extracting coefficients of the green line
coefficients_tls_line_codon_floor <- coef(lm(Xfit_codon_floor[, 2] ~ Xfit_codon_floor[, 1] + 1))
equation_codon_floor <- paste("y =", round(coefficients_tls_line_codon_floor[2], 2), "* x +", round(coefficients_tls_line_codon_floor[1], 2))

ggplot(data = as.data.frame(logData_codon_floor),
       mapping = aes(x = V1,
                     y = V2)) +
  geom_point(size = 1, alpha = 0.3) +
  geom_abline(intercept = coefficients_tls_line_codon_floor[1], slope = coefficients_tls_line_codon_floor[2], col = "darkgreen", linewidth = 1) +
  #geom_abline(intercept = 0, slope = 1, col = "red", lwd = 1) + #45 degrees line
  annotate("text", x = 1, y = max(logData_codon_floor[, 2]), label = paste0("Explained variance = ", round(adjusted_explained_variance_floor, 2), "%"), hjust = 0.05, vjust = 1, size = 4) +
  theme_bw() +
  labs(x = "Counts Sucrose",
       y = "Counts 1X PEG-Spin") +
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

#==================================================

#SEQUENCE LENGTH DISTRIBUTION

ggplot(data = filter(Sequence_length_distributions_all_samples, Sample == "Histag Lysate", `Date of Experiment` == "Jun-23"),
       mapping = aes(
         x = `Length (nt)`,
         y = Frequency
       )) +
  geom_col(fill = 'purple', colour = 'black', linewidth = 1.3, position = 'dodge') +
  scale_x_continuous(n.breaks = 18) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 23)) +
  theme_bw() +
  labs(colour = "Sample") +
  theme(legend.position = 'none') +
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
        legend.key.height= unit(0.05, 'cm'),
        legend.key.width= unit(0.05, 'cm')) +
  theme(plot.subtitle = element_text(face = "bold"),
        plot.caption = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.key = element_rect(linetype = "solid")) + 
  theme(axis.ticks = element_line(size = 0.5)) + theme(panel.background = element_rect(fill = NA)) + 
  theme(axis.line = element_line(linetype = "solid")) + theme(axis.line = element_line(linetype = "blank"),
                                                              panel.grid.major = element_line(colour = "gray89",
                                                                                              linetype = "blank"), panel.grid.minor = element_line(linetype = "blank"),
                                                              panel.background = element_rect(linetype = "solid"))
#==================================================

#RNA FRACTIONS

ggplot(data = filter(RNA_fractions_all_samples, Sample == "2X PEG", `Date of Experiment` == "Jan-23"), 
       mapping = aes(x = "", y = Percentage, fill = Fraction)
) + 
  geom_bar(width = 1, stat = "identity", col = 'black', linewidth = 0.5) +
  geom_text(aes(label = paste(Percentage, "%")),
            position = position_stack(vjust = 0.5),
            colour = 'white') +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  coord_polar("y") +
  labs(x = NULL, y = NULL)


#=======================================

#PLOT ILE CODONS ON PEAKS COMPARISON OF TWO SAMPLES ON SINGLE GENES

#BiocManager::install("Biostrings")

library(Biostrings)

# Read the reference genome in FASTA format
genome <- readDNAStringSet("C:/Users/aagnoli/OneDrive - UvA/B. subtilis 168 annotation and genome files/WT-Prspb-amyI.fa")

# Read the GFF file
gff_lines <- readLines("C:/Users/aagnoli/OneDrive - UvA/B. subtilis 168 annotation and genome files/wt-prspb-amyI-GFF3_updated_amyM_positions.gff")

# Initialize empty lists to store locus_tags, start positions, end positions, and sequences
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

# Create a data frame
gene_info_df <- data.frame(
  LocusTag = unlist(locus_tags),
  StartPosition = unlist(start_positions),
  EndPosition = unlist(end_positions),
  Sequence = as.character(sequences),
  Strand = unlist(strands)
)

#If you are using the old gff file (old start and end positions of amyI gene) use the lines below

# Assign true coordinates and sequence for the start of AmyI as the gff file includes the 5' region of AmyI
#gene_info_df_amyI <- filter(gene_info_df, LocusTag == "BSU41070")
#New_StartPosition <- 4216116
#New_EndPosition <- 4218275
#New_Sequence <- substr(gene_info_df_amyI$Sequence, start = 461, stop = 2620)
#gene_info_df[gene_info_df$LocusTag == "BSU41070", c("StartPosition", "EndPosition", "Sequence")] <- c(as.numeric(New_StartPosition), as.numeric(New_EndPosition), New_Sequence)

# Define the Isoleucine codons
isoleucine_codons <- c("ATT", "ATC", "ATA")

# Create a list to store positions of the first nucleotide of Isoleucine codons for each locus tag
all_positions <- list()

# Iterate through each sequence in the data frame
for (i in 1:nrow(gene_info_df)) {
  locus_tag <- gene_info_df$LocusTag[i]
  sequence <- gene_info_df$Sequence[i]
  
  # Create a vector to store positions of the first nucleotide of Isoleucine codons for the current sequence
  isoleucine_positions <- numeric()
  
  # Iterate through the sequence and find positions of the first nucleotide of Isoleucine codons
  for (i in seq(1, nchar(sequence), by = 3)) {
    codon <- substr(sequence, i, i + 2)
    if (codon %in% isoleucine_codons) {
      isoleucine_positions <- c(isoleucine_positions, i)
    }
  }
  
  # Append the list of positions to the all_positions list, using the locus tag as the list name
  all_positions[[locus_tag]] <- isoleucine_positions
}

# Merge the positions with the data frame using locus tags
gene_info_df$IsoleucinePositions <- sapply(gene_info_df$LocusTag, function(tag) all_positions[[tag]])

# Create a new column to store the converted positions
gene_info_df$ConvertedPositions <- NA

# Initialize an empty list to store the converted positions
converted_positions_list <- vector("list", nrow(gene_info_df))

#If you are using the old gff file (old start and end positions of amyI gene) use the lines below

#Make sure the StartPosition and EndPosition variables in the gene_info_df data frame are numerical
#gene_info_df$StartPosition <- as.numeric(gene_info_df$StartPosition)
#gene_info_df$EndPosition <- as.numeric(gene_info_df$EndPosition)

# Iterate through the rows and calculate the converted positions
for (i in 1:nrow(gene_info_df)) {
  start_pos <- gene_info_df$StartPosition[i]
  isoleucine_positions <- as.numeric(unlist(gene_info_df$IsoleucinePositions[i]))
  
  # Calculate the converted positions and store as a list
  converted_positions <- start_pos + isoleucine_positions - 1
  converted_positions_list[[i]] <- converted_positions
  
  # Extract isoleucine codons and store as a list
  sequence <- gene_info_df$Sequence[i]
  codons <- lapply(isoleucine_positions, function(pos) substr(sequence, pos, pos + 2))
  gene_info_df$IsoleucineCodons[i] <- list(codons)
}

#select gene for plotting
gene_of_interest <- "BSU41070" #use locus_tag

# Add the list to the data frame
gene_info_df$ConvertedPositions <- converted_positions_list

filtered_gene_info_df <- filter(gene_info_df, LocusTag == gene_of_interest) %>% select(-6)

vlines <- cbind(unlist(filtered_gene_info_df$ConvertedPositions), unlist(filtered_gene_info_df$IsoleucineCodons))
# Create a data frame for vlines
vlines_df <- data.frame(
  Position = vlines[, 1],
  IsoleucineCodon = as.factor(vlines[, 2])
)

vlines_df$Position <- as.numeric(vlines_df$Position)

# Create a full range of positions combining both datasets
all_positions <- union(Normalized_riboseq_1$position, Normalized_riboseq_2$position)

# Create a dataframe with all positions
all_data <- data.frame(position = all_positions)

# Left join each dataset separately
all_data <- all_data %>%
  left_join(Normalized_riboseq_1, by = c('position')) %>%
  left_join(Normalized_riboseq_2, by = c('position'),
            suffix = c('_1', '_2'))

# Filter data for the gene of interest
filtered_data <- all_data %>% 
  filter(locus_tag_1 == gene_of_interest | locus_tag_2 == gene_of_interest) #change gene1 and 2 with locus tag if you use it above

# Plot the data
Plot <- ggplot(data = filter(filtered_data, position >= 4216116, position <= 4218275), mapping = aes(x = position)) + #if you are plotting amyI profile, use this line: data = filter(filtered_data, position >= 4216116, position <= 4218275)
  geom_col(aes(y = coalesce(Norm_count_1, 0)), fill = 'orange', size = 1, alpha = 0.5) +
  geom_col(aes(y = coalesce(Norm_count_2, 0)), fill = 'navy', size = 1, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold")) +
  labs(x = "Position", y = "Norm_count") +
  annotate("text", x = Inf, y = Inf, label = paste("Gene is on", ifelse(filtered_gene_info_df$Strand == "-", "-", "+"), "strand"), color = "black", hjust = 1.03, vjust = 1.2, size = 3) + #indicate if gene coding frame is sense or antisense 
  theme(panel.background = element_rect(fill = "gray48"))

#Visualize only comparison plot
Plot  

# Plot the data with arrows colored by codon
codon_colors <- c("magenta", "limegreen", "cyan")

# Create a custom color palette for the three codons
custom_color_palette <- c("ATA" = codon_colors[1],
                          "ATC" = codon_colors[2],
                          "ATT" = codon_colors[3])

Plot +  
  #scale_y_continuous(limits = c(0,500)) +
  geom_segment(data = vlines_df, aes(x = Position, y = 0,
                                     xend = Position, yend = 0.1,
                                     color = IsoleucineCodon),
               arrow = arrow(length = unit(0.05, "inches"), type = 'closed'), alpha = 0.5, lineend = 'butt', linejoin = 'mitre') +
  scale_color_manual(values = custom_color_palette) +
  guides(color = guide_legend(title = "Isoleucine Codon")) +
  theme(legend.position = c(0.05, 0.95),  # Adjust these values to position the legend
        legend.justification = c(0, 1),
        legend.box.just = "left") + theme(legend.title = element_text(face = "bold"),
    legend.position = "top", legend.direction = "horizontal")

#=====================================

#tRNA STUFF

tRNA_standard <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Test for tRNA fractions Histag ST and WASH June23/Filtered_tRNA_Standard_35_full.csv")
tRNA_wash <- read_csv("C:/Users/aagnoli/OneDrive - UvA/RNA sequencing data/Test for tRNA fractions Histag ST and WASH June23/Filtered_tRNA_Wash_35_full.csv")

sum(tRNA_standard$count)
sum(tRNA_wash$count)

Norm_tRNA_standard <- mutate(tRNA_standard, Norm_count = tRNA_standard$count*(sum(tRNA_standard$count)/sum(tRNA_standard$count)))
Norm_tRNA_wash <- mutate(tRNA_wash, Norm_count = tRNA_wash$count*(sum(tRNA_standard$count)/sum(tRNA_wash$count)))

Norm_tRNA_standard <- filter(Norm_tRNA_standard, grepl('BSU_tRNA', Norm_tRNA_standard$locus_tag))
sum_per_gene_standard <- group_by(Norm_tRNA_standard, gene) %>% summarize(tot_count = sum(Norm_count))

Norm_tRNA_wash <- filter(Norm_tRNA_wash, grepl('BSU_tRNA', Norm_tRNA_wash$locus_tag))
sum_per_gene_wash <- group_by(Norm_tRNA_wash, gene) %>% summarize(tot_count = sum(Norm_count))

ggplot() +
  geom_col(data = sum_per_gene_standard,
           mapping = aes(
             x = gene,
             y = tot_count)) +
  geom_col(data = sum_per_gene_wash,
           mapping = aes(
             x = gene,
             y = tot_count), fill = 'yellow', alpha = 0.5, col = 'black') +
  theme(axis.text.x = element_text(vjust = 0.6, angle = 90))

#==========================================================================

##FOLD-CHANGE ON GENE COVERAGES

#Here I use 4 datasets; 2 

riboseq_1 <- filtered_Histag_Standard_35_full_June_2023
riboseq_2 <- filtered_2X_Spin_35_full_June_2023
riboseq_3 <- filtered_Histag2_July_35_full_July_2023
riboseq_4 <- filtered_2X_Spin_35_full_January_2023

#Normalize datasets (only considering ORFs)

riboseq_1_ORFs <- riboseq_1 %>% filter(in_orf_90 == TRUE)
sum(riboseq_1_ORFs$count)
riboseq_2_ORFs <- riboseq_2 %>% filter(in_orf_90 == TRUE)
sum(riboseq_2_ORFs$count)
riboseq_3_ORFs <- riboseq_3 %>% filter(in_orf_90 == TRUE)
sum(riboseq_3_ORFs$count)
riboseq_4_ORFs <- riboseq_4 %>% filter(in_orf_90 == TRUE)
sum(riboseq_4_ORFs$count)

Normalized_riboseq_1 <- mutate(riboseq_1_ORFs, Norm_count = riboseq_1_ORFs$count*(sum(riboseq_1_ORFs$count)/sum(riboseq_1_ORFs$count)))
Normalized_riboseq_2 <- mutate(riboseq_2_ORFs, Norm_count = riboseq_2_ORFs$count*(sum(riboseq_1_ORFs$count)/sum(riboseq_2_ORFs$count)))
Normalized_riboseq_3 <- mutate(riboseq_3_ORFs, Norm_count = riboseq_3_ORFs$count*(sum(riboseq_1_ORFs$count)/sum(riboseq_3_ORFs$count)))
Normalized_riboseq_4 <- mutate(riboseq_4_ORFs, Norm_count = riboseq_4_ORFs$count*(sum(riboseq_1_ORFs$count)/sum(riboseq_4_ORFs$count)))

#Calculate Gene coverage

Gene_coverage_Normalized_riboseq_1 <- Normalized_riboseq_1 %>% group_by(locus_tag, gene_length, gene) %>% summarise(Norm_tot_counts = sum(Norm_count))
Gene_coverage_Normalized_riboseq_1 <- Gene_coverage_Normalized_riboseq_1 %>% mutate(gene_coverage = Norm_tot_counts/gene_length)
Gene_coverage_Normalized_riboseq_2 <- Normalized_riboseq_2 %>% group_by(locus_tag, gene_length, gene) %>% summarise(Norm_tot_counts = sum(Norm_count))
Gene_coverage_Normalized_riboseq_2 <- Gene_coverage_Normalized_riboseq_2 %>% mutate(gene_coverage = Norm_tot_counts/gene_length)
Gene_coverage_Normalized_riboseq_3 <- Normalized_riboseq_3 %>% group_by(locus_tag, gene_length, gene) %>% summarise(Norm_tot_counts = sum(Norm_count))
Gene_coverage_Normalized_riboseq_3 <- Gene_coverage_Normalized_riboseq_3 %>% mutate(gene_coverage = Norm_tot_counts/gene_length)
Gene_coverage_Normalized_riboseq_4 <- Normalized_riboseq_4 %>% group_by(locus_tag, gene_length, gene) %>% summarise(Norm_tot_counts = sum(Norm_count))
Gene_coverage_Normalized_riboseq_4 <- Gene_coverage_Normalized_riboseq_4 %>% mutate(gene_coverage = Norm_tot_counts/gene_length)

Merged_riboseq1_and_riboseq2 <- merge(Gene_coverage_Normalized_riboseq_1, Gene_coverage_Normalized_riboseq_2, by = 'locus_tag', all = T, suffixes = c('1', '2'))
Merged_riboseq3_and_riboseq4 <- merge(Gene_coverage_Normalized_riboseq_3, Gene_coverage_Normalized_riboseq_4, by = 'locus_tag', all = T, suffixes = c('3','4'))
Merged_riboseq <- merge(Merged_riboseq1_and_riboseq2, Merged_riboseq3_and_riboseq4, by = 'locus_tag')

Avg_gene_coverage <- Merged_riboseq %>% rowwise() %>% mutate(average_gene_coverage1_2 = mean(c(gene_coverage1, gene_coverage2), na.rm = TRUE), average_gene_coverage3_4 = mean(c(gene_coverage3, gene_coverage4), na.rm = TRUE))

Merged_Avg_gene_coverage <- Avg_gene_coverage %>%
  rowwise() %>%
  mutate(FC = ifelse(
    is.na(average_gene_coverage1_2) | is.na(average_gene_coverage3_4),
    NA,
    log2(average_gene_coverage1_2) - log2(average_gene_coverage3_4)
  ))

Merged_Avg_gene_coverage$FC <- ifelse(
  is.na(Avg_gene_coverage$average_gene_coverage1_2) | is.na(Avg_gene_coverage$average_gene_coverage3_4),
  NA, # Set the fold change to NA if either of the values is NA
  log2(Avg_gene_coverage$average_gene_coverage1_2) - log2(Avg_gene_coverage$average_gene_coverage3_4)
)

Significant_genes <- filter(Merged_Avg_gene_coverage, FC < -1.5 | FC > 1.5)

###spin depl vs Histag standard
Extra_Normalized_2X_Spin_depl <- mutate(X2X_Spin_depl_ORFs, Norm_count = X2X_Spin_depl_ORFs$count*(sum(X2X_Spin_depl_ORFs$count)/sum(X2X_Spin_depl_ORFs$count)))
Extra_Normalized_Histag_Standard <- mutate(Histag_Standard_ORFs, Norm_count = Histag_Standard_ORFs$count*(sum(X2X_Spin_depl_ORFs$count)/sum(Histag_Standard$count)))

Extra_Gene_coverage_2X_Spin_depl <- Extra_Normalized_2X_Spin_depl %>% group_by(gene, gene_length) %>% summarise(Norm_tot_counts = sum(Norm_count))
Extra_Gene_coverage_2X_Spin_depl <- Extra_Gene_coverage_2X_Spin_depl %>% mutate(gene_coverage = Norm_tot_counts/gene_length)
Extra_Gene_coverage_Histag_Standard <- Extra_Normalized_Histag_Standard %>% group_by(gene, gene_length) %>% summarise(Norm_tot_counts = sum(Norm_count))
Extra_Gene_coverage_Histag_Standard <- Extra_Gene_coverage_Histag_Standard %>% mutate(gene_coverage = Norm_tot_counts/gene_length)

Merged_gene_coverage_2x_Spin_depl_and_Histag_Standard_lol <- merge(Extra_Gene_coverage_2X_Spin_depl, Extra_Gene_coverage_Histag_Standard, by = 'gene', suffixes = c("_Spin_depl", "_Histag"))

Merged_gene_coverage_2x_Spin_depl_and_Histag_Standard_lol$fold_change <- log2(Merged_gene_coverage_2x_Spin_and_Histag_Standard_lol$gene_coverage_Spin / Merged_gene_coverage_2x_Spin_and_Histag_Standard_lol$gene_coverage_Histag)

###spin vs Histag standard
Extra_Normalized_2X_Spin <- mutate(X2X_Spin_ORFs, Norm_count = X2X_Spin_ORFs$count*(sum(Histag_Standard_ORFs$count)/sum(X2X_Spin_depl_ORFs$count)))
Extra_Normalized_Histag_Standard <- mutate(Histag_Standard_ORFs, Norm_count = Histag_Standard_ORFs$count*(sum(Histag_Standard_ORFs$count)/sum(Histag_Standard$count)))

Extra_Gene_coverage_2X_Spin <- Extra_Normalized_2X_Spin %>% group_by(gene, gene_length) %>% summarise(Norm_tot_counts = sum(Norm_count))
Extra_Gene_coverage_2X_Spin <- Extra_Gene_coverage_2X_Spin %>% mutate(gene_coverage = Norm_tot_counts/gene_length)
Extra_Gene_coverage_Histag_Standard <- Extra_Normalized_Histag_Standard %>% group_by(gene, gene_length) %>% summarise(Norm_tot_counts = sum(Norm_count))
Extra_Gene_coverage_Histag_Standard <- Extra_Gene_coverage_Histag_Standard %>% mutate(gene_coverage = Norm_tot_counts/gene_length)

Merged_gene_coverage_2x_Spin_and_Histag_Standard <- merge(Extra_Gene_coverage_2X_Spin, Extra_Gene_coverage_Histag_Standard, by = 'gene', suffixes = c("_Spin", "_Histag"))

Merged_gene_coverage_2x_Spin_and_Histag_Standard$fold_change <- log2(Merged_gene_coverage_2x_Spin_and_Histag_Standard$gene_coverage_Spin / Merged_gene_coverage_2x_Spin_and_Histag_Standard$gene_coverage_Histag)

#====================================================================

#PCA

# Load and preprocess individual datasets
dataset1 <- filtered_2X_Spin_35_full_January_2023  # Load your data from each experiment
dataset2 <- filtered_2X_Spin_35_full_June_2023
dataset3 <- filtered_2X_PEG_35_full_January_2023
# ... other datasets if required

# Combine datasets
combined_data <- bind_rows(
  
  data.frame(Method = "Standard", dataset1),
  data.frame(Method = "Wash", dataset2),
  data.frame(Method = "Lysate", dataset3)
  
  # ... bind other datasets if required
)

# Perform PCA on combined data
pca_result <- prcomp(combined_data[, c(-1,-2,-4,-5,-6,-7,-8,-9,-11)], scale = TRUE)  # Exclude the columns with non-numerical values

# Create a data frame for plotting
pca_df <- data.frame(Method = combined_data$Method, PCs = pca_result$x)

# Rename the PCA columns to match the column names in pca_result
colnames(pca_df)[2:ncol(pca_df)] <- colnames(pca_result$x)

# Create a PCA plot using ggplot2
ggplot(pca_df, aes(x = PC1, y = PC2, color = Method)) +
  geom_point() +
  labs(title = "PCA Plot by Method", x = "Principal Component 1", y = "Principal Component 2")


#=============================
#Create fasta file with small RNAs from B. subtilis for bowtie2 depletion (needs loading of gene_info_df in script above)


library(Biostrings)

# Specify the paths to your GFF and FASTA reference files
genome <- readDNAStringSet("C:/Users/aagnoli/OneDrive - UvA/B. subtilis 168 annotation and genome files/WT-Prspb-amyI.fa")

# Read the GFF file
gff_lines <- readLines("C:/Users/aagnoli/OneDrive - UvA/B. subtilis 168 annotation and genome files/wt-prspb-amyI-GFF3.gff")

# Filter gene_info_df to select entries starting with "BSU_misc"
filtered_gene_info_df <- gene_info_df[grep("^BSU_misc", gene_info_df$LocusTag), ]

# Create a character vector to store the FASTA entries
fasta_entries <- character()

# Loop through the filtered data and create FASTA entries
for (i in 1:nrow(filtered_gene_info_df)) {
  locus_tag <- filtered_gene_info_df$LocusTag[i]
  sequence <- filtered_gene_info_df$Sequence[i]
  fasta_entry <- paste(">", locus_tag, "\n", sequence, sep = "")
  fasta_entries <- c(fasta_entries, fasta_entry)
}

# Write the FASTA entries to a file
writeLines(fasta_entries, "C:/Users/aagnoli/OneDrive - UvA/B. subtilis 168 annotation and genome files/BSU_misc_genes.fasta")
