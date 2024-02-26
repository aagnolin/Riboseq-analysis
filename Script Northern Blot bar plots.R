library(readxl)
NB_quantification_2024_02_12 <- read_excel("C:/Users/aagnoli/OneDrive - UvA/Gels/Northern Blots/2024-02-12 Northern Blot/2024-02-12 Quantification Northern Blot ImageStudio.xlsx")
NB_quantification_2024_02_20 <- read_excel("C:/Users/aagnoli/OneDrive - UvA/Gels/Northern Blots/2024-02-20 Northern Blot/2024-02-20 Quantification Northern Blot ImageStudio.xlsx")

data_NB <- NB_quantification_2024_02_20 

ggplot(data = data_NB,
       mapping = aes(x = Name,
                     y = Signal)) +
  geom_col(col = 'black', linewidth = 1, width = 0.5) +
  theme_bw()
