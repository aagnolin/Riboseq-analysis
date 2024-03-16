library(readxl)
library(tidyverse)

#import data
Assay_2024_03_15 <- read_xlsx("C:/Users/aagnoli/OneDrive - UvA/WP 4 - Translation reporter system/High throughput screening experiment/Assay 2024-03-15.xlsx")

#choose dataset
Amylase_data <- Assay_2024_03_15

#plot data
Amylase_data_summary <- Amylase_data %>%
  group_by(`Time (h)`, Sample) %>%
  summarize(OD600_mean = mean(OD600),
            OD600_sd = sd(OD600),
            Amylase_mean = mean(`Amylase activity`),
            Amylase_sd = sd(`Amylase activity`))

# Plot the data with error bars
ggplot(data = Amylase_data_summary,
       mapping = aes(x = `Time (h)`,
                     group = Sample,
                     shape = Sample)) +
  geom_point(aes(y = OD600_mean), colour = 'blue') +
  geom_line(aes(y = OD600_mean), colour = 'blue') +
  geom_errorbar(aes(ymin = OD600_mean - OD600_sd, ymax = OD600_mean + OD600_sd), width = 0.5, colour = 'blue') +
  geom_point(aes(y = Amylase_mean/4), colour = 'darkred') +
  geom_line(aes(y = Amylase_mean/4), linetype = 'dashed', colour = 'darkred') +
  geom_errorbar(aes(ymin = Amylase_mean/4 - Amylase_sd/4, ymax = Amylase_mean/4 + Amylase_sd/4), width = 0.5, alpha = 0.4, colour = 'darkred') +
  scale_y_continuous(sec.axis = sec_axis(~.*4 , name = 'Amylase activity')) +
  labs(x = "Time (h)",
       y = "OD600") +
  theme_bw() +
  theme(axis.text.y.left = element_text(color = "blue")) +
  theme(axis.text.y.right = element_text(color = "darkred"))

#boxplot
ggplot(data = filter(Amylase_data)) +
  geom_boxplot(aes(x = Sample,
               y = `Amylase activity`)) +
  facet_wrap(~`Time (h)`)

#t test amylase activity by strain
df1 <- Amylase_data %>% filter(`Time (h)` != 0) %>% select(`Amylase activity`, Sample)
t.test(data = df1, `Amylase activity` ~ Sample)
