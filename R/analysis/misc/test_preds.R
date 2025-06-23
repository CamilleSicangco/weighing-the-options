# Test ability to predict transpiration when forced with predawn and midday water potentials
# by Camille Sicangco
# Created 10 Oct 2024

# Hydraulics
P50 = 4.31
P88 = 6.64
Weibull = fit_Weibull(P50, P88)
b = Weibull[1,1]
c = Weibull[1,2]
Pcrit = calc_Pcrit(b, c)

# For clean flux data #######################################################
# Calculate predicted Pleaf and E values
P = lapply(-WTC4_LWP$Ps, Ps_to_Pcrit, Pcrit = Pcrit)
E_vecs = mapply(trans_from_vc, P = P, 
                kmax_25 = 1.5, 
                Tair = WTC4_LWP$Tair, 
                b = b, c = c, constant_kmax = TRUE,
                SIMPLIFY = FALSE)
preds = sapply(1:nrow(WTC4_LWP), function(i) {
  Pleaf.i = -WTC4_LWP$Pleaf[i]
  j = which.min(abs(Pleaf.i - P[[i]]))
  Pleaf.pred = P[[i]][j]
  E.pred = E_vecs[[i]][j]
  AUC = E.pred / WTC4_LWP$kmax[i]
  kmax.alt = WTC4_LWP$E[i] / AUC
  return(data.frame(chamber = WTC4_LWP$chamber[i],HWtrt = WTC4_LWP$HWtrt[i], 
                    T_treatment = WTC4_LWP$T_treatment[i],
                    datetime = WTC4_LWP$DateTime_hr[i], 
                    Pleaf.pred, E.pred, Pleaf = Pleaf.i, E = WTC4_LWP$E[i], 
                    Tair = WTC4_LWP$Tair[i],
                    kmax.alt, kmax = WTC4_LWP$kmax[i]
                    ))
})

# Format data frame with predictions
preds = data.frame(t(preds)) 
preds[4:ncol(preds)] <- sapply(preds[4:ncol(preds)],as.numeric)
preds$datetime = as_datetime(preds$datetime)
preds[1:3] <- sapply(preds[1:3],as.character)
preds = preds %>% 
  mutate(Pleaf.diff = Pleaf.pred - Pleaf, E.diff = E.pred - E)

# For extended flux data #######################################################
# Calculate predicted Pleaf and E values
P = lapply(-WTC4_LWP.ext$Ps, Ps_to_Pcrit, Pcrit = Pcrit)
E_vecs = mapply(trans_from_vc, P = P, 
                kmax_25 = WTC4_LWP.ext$kmax, 
                Tair = WTC4_LWP.ext$Tair, 
                b = b, c = c, constant_kmax = TRUE,
                SIMPLIFY = FALSE)
preds.ext = sapply(1:nrow(WTC4_LWP.ext), function(i) {
  Pleaf.i = -WTC4_LWP.ext$Pleaf[i]
  j = which.min(abs(Pleaf.i - P[[i]]))
  Pleaf.pred = P[[i]][j]
  E.pred = E_vecs[[i]][j]
  AUC = E.pred / WTC4_LWP.ext$kmax[i]
  kmax.alt = WTC4_LWP.ext$E[i] / AUC
  return(data.frame(chamber = WTC4_LWP.ext$chamber[i],HWtrt = WTC4_LWP.ext$HWtrt[i], 
                    T_treatment = WTC4_LWP.ext$T_treatment[i],
                    datetime = WTC4_LWP.ext$DateTime_hr[i], 
                    Pleaf.pred, E.pred, Pleaf = Pleaf.i, E = WTC4_LWP.ext$E[i], 
                    Tair = WTC4_LWP.ext$Tair[i],
                    kmax.alt, kmax = WTC4_LWP.ext$kmax[i],
                    AUC = AUC
  ))
})

# Format data frame with predictions
preds.ext = data.frame(t(preds.ext)) 
preds.ext[5:ncol(preds.ext)] <- sapply(preds.ext[5:ncol(preds.ext)],as.numeric)
preds.ext$datetime = ymd_hms(preds.ext$datetime)
preds.ext[1:3] <- sapply(preds.ext[1:3],as.character)
preds.ext = preds.ext %>% 
  mutate(Pleaf.diff = Pleaf.pred - Pleaf, E.diff = E.pred - E)

# Plot predictions minus observations versus Tair
# Problem is kmax value
preds %>% 
  #filter(chamber == "C06") %>% 
  ggplot(aes(x = Tair, y = E.diff, color = HWtrt)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  theme_classic() +
  scale_color_manual(values = c("lightblue3", "orange"))
preds %>% 
  ggplot(aes(x = Tair, y = Pleaf.diff, color = HWtrt)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  theme_classic() +
  scale_color_manual(values = c("lightblue3", "orange"))

# Comparison of back-calculated kmax
# Off, but notably original is based on pre-heatwave data
preds.ext %>% 
  filter(datetime < as_datetime("2016-11-01")) %>% 
  ggplot(aes(x = chamber, y = kmax.alt, fill = HWtrt)) +
  geom_violin(width = 1.4) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = c("lightblue3", "orange")) +
  geom_point(aes(x = chamber, y = kmax), shape = 4, size = 7, color = "darkgreen") +
  theme_classic()

plot(fluxes_df$FluxH2O / fluxes_df$LeafArea)
plot(fluxes_v3$FluxH2O / fluxes_v3$LeafArea, ylim = c(-0.001, 0.004))

