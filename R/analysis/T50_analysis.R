# Analysis of T50 data
# by Camille Sicangco
# Created 14 Nov 2024

T50_data = read.csv("https://hie-pub.westernsydney.edu.au/07fcd9ac-e132-11e7-b842-525400daae48/WTC_TEMP-PARRA_HEATWAVE-FLUX-PACKAGE_L1/data/WTC_TEMP-PARRA_CM_T50-CI_20161019-20161117_L1.csv")

View(T50_data)

# ANOVA to test for differences based on heatwave and temperature treatment
T50_anova = aov(data = dplyr::filter(T50_data, Date == "2016-11-01"), formula = T50_mean ~ T_treatment * HW_treatment)
summary(T50_anova) # Only heatwave effect is significant


# Calculate average T50 for control and heatwave treatments
T50_data1 = T50_data %>% 
  filter(Date == "2016-11-01") %>% 
  group_by(HW_treatment) %>% 
  summarise(mean(T50_mean)) %>% 
  ungroup()
T50_data1
