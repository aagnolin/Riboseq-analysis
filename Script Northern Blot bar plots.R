library(readxl)
library(tidyverse)
NB_quantification_2024_02_12 <- read_excel("C:/Users/aagnoli/OneDrive - UvA/Gels/Northern Blots/2024-02-12 Northern Blot/2024-02-12 Quantification Northern Blot ImageStudio.xlsx")
NB_quantification_2024_02_20 <- read_excel("C:/Users/aagnoli/OneDrive - UvA/Gels/Northern Blots/2024-02-20 Northern Blot/2024-02-20 Quantification Northern Blot ImageStudio.xlsx")
NB_quantification_2024_03_11 <- read_excel("C:/Users/aagnoli/OneDrive - UvA/Gels/Northern Blots/2024-03-11 Northern Blot/2024-03-11 Quantification Northern Blot ImageStudio.xlsx")
NB_quantification_2024_03_14 <- read_excel("C:/Users/aagnoli/OneDrive - UvA/Gels/Northern Blots/2024-03-14 Northern Blot/2024-03-14 Quantification Northern Blot ImageStudio.xlsx")

Merged_NB <- rbind(NB_quantification_2024_02_12, NB_quantification_2024_02_20, NB_quantification_2024_03_11, NB_quantification_2024_03_14)

data_NB <- Merged_NB


#boxplot 13 and 20 duplo
ggplot(data = filter(Merged_NB),
       mapping = aes(x = Name,
                     y = Signal)) +
  geom_boxplot(fill = "grey", width = 0.5, linewidth = 1) +
  geom_point(size = 2) +
  theme_bw()
  


#barplot 2nd experiment
ggplot(data = filter(Merged_NB, `Image Name` == "0007487_02"),
       mapping = aes(x = Name,
                     y = Signal)) +
  geom_col(col = 'black', linewidth = 1, width = 0.5) +
  theme_bw()

df <- data_NB %>% select(Name, Signal) %>% filter(Name == 'BAA013' | Name == "BAA020")
t.test(data = df, Signal ~ Name)
