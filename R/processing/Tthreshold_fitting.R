# Fitting of FvFm-T curves
# by Camille Sicangco
# Created 22 Nov 2024

T50_raw_data = read.csv("data/in/raw/WTC_TEMP-PARRA_CM_T50_20161019-20161117_L0.csv")

# Plot Fv/Fm measured at start of the heatwave 
T50_raw_data %>% 
  filter(Time == "ON",
         Date %in% c("2016-11-01", "2016-11-02")) %>% 
  ggplot(aes(x = T_treat, y = Fv.Fm, color = as.factor(Replicate))) +
    geom_point() +
  facet_wrap(vars(chamber)) +
  theme_classic()

# Remove outliers
T50_data = T50_raw_data %>% 
  filter(Date %in% c("2016-11-01", "2016-11-02"),
         !(chamber == "C12" & Replicate == 2 & T_treat %in% c(44, 46)),
         !(chamber == "C09" & Replicate == 3 & T_treat == 44)
         )


T50_data.l = split(T50_data, ~chamber)

# Fit Tcrit, T50
fl_fits = sapply(T50_data.l, function(df) {
  max = max(df$Fv.Fm[df$Time == "ON" & df$T_treat == 24])
  
  T50.fit = fitcond(dfr = filter(df, Time == "ON"),
                    varnames = c(K = "Fv.Fm", WP = "T_treat"),
                    Kmax = max,
                    x = 50)
  Tcrit.fit = fitcond(dfr = filter(df, Time == "ON"),
                      varnames = c(K = "Fv.Fm", WP = "T_treat"),
                      Kmax = max,
                      x = 5)
  
  T50 = coef(T50.fit)[2,1]
  Tcrit = coef(Tcrit.fit)[2,1]
  out = c(df$HW_treatment[1], Tcrit, T50)
  names(out) = c("HW_treatment", "Tcrit", "T50")
  return(out)
})
fl_fits = data.frame(t(fl_fits)) 
fl_fits[2:3] = as.numeric(unlist(fl_fits[2:3]))
t.test(T50 ~ HW_treatment, fl_fits)
t.test(Tcrit ~ HW_treatment, fl_fits)
fl_fits %>% 
  group_by(HW_treatment) %>% 
  dplyr::summarise(across(1:2, mean))
