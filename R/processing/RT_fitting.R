# Fitting of R-T curves
# by Camille Sicangco
# Created 22 Nov 2024

# Load respiration curves
RT_raw_data = read.csv("data/in/raw/WTC_TEMP-PARRA_CM_GX-RvT_20161021-20161109_L1.csv")

# Clean data
RT_data = RT_raw_data %>% 
  mutate(R = -1 * Photo, Tleaf2 = Tleaf^2, ln_R = log(R),
         HWtrt = ifelse(chamber %in% c("C03", "C04", "C07", "C09", "C10", "C12"),
                        "HW", "C")) %>% 
  filter(leaf %in% c("waltz", "waltz3") 
         & Time == "pre" & !(chamber == "C10" & Date == "2016-10-21"))

# Plot R-T curves
RT_data %>% 
  ggplot(aes(x = Tleaf, y = R)) +
  geom_point() +
  theme_classic() +
  facet_wrap(vars(chamber))

# Fit parameters
lmer.fit = lmer(
  ln_R ~ Tleaf + Tleaf2 + (1 | chamber),
  RT_data)
summary(lmer.fit)
a = -2.398
b = coef(lmer.fit)[[1]][1,2]
c = coef(lmer.fit)[[1]][1,3]
R25 = exp(a + 25*b + 25^2*c)
