# Processing of the WTC4 dataset
# by Camille Sicangco
# Created 27 September 2024

# Load and clean SWC data ######################################################
input_folder <- "data/in/raw/SWC"
SWC_df = clean_SWC(input_folder)
#write.csv(SWC_df, "data/in/SWC_clean.csv", row.names = FALSE)

# Aggregate to hourly averages
SWC_df$DateTime_hr <- HIEv::nearestTimeStep(SWC_df$DateTime_hr,nminutes=30,"floor")
#dat.hr <- doBy::summaryBy(.~DateTime_hr+chamber+T_treatment+HWtrt+combotrt,FUN=mean,keep.names=T,data=SWC_df)

# Compare temperature measurements #############################################

# Load data
WTC4_temp_data = read.csv("http://hie-pub.westernsydney.edu.au/07fcd9ac-e132-11e7-b842-525400daae48/WTC_TEMP-PARRA_HEATWAVE-FLUX-PACKAGE_L1/data/WTC_TEMP-PARRA_CM_TEMPERATURES-COMBINED_20161010-20161123_L1.csv")

# Histogram of thermocouple, radiometer, and air temperatures
WTC4_temp_data %>% 
  ggplot() +
  geom_histogram(aes(Tleaf), fill = "deeppink", alpha = 0.5) +
  geom_histogram(aes(TargTempC_Avg), fill = "lightblue3", alpha = 0.5) +
  geom_histogram(aes(Tair_al), fill = "grey", alpha = 0.5) +
  theme_classic()

# Tleaf vs Tair for thermocouple and radiometer data
WTC4_temp_data %>% 
  ggplot() +
  geom_point(aes(x = Tair_al, y = Tleaf), color = "pink", shape = 1) +
  geom_point(aes(x = Tair_al, y = TargTempC_Avg), color = "lightblue3", shape = 1) +
  theme_classic() +
  geom_abline(intercept = 0, slope = 1)

# Thermocouple minus radiometer data
## TC temps higher
WTC4_temp_data %>% 
  ggplot() +
  geom_point(aes(x = TargTempC_Avg, y = Tleaf - TargTempC_Avg), color = "pink", shape = 1) +
  #geom_point(aes(x = Tair_al, y = ), color = "lightblue3", shape = 1) +
  theme_classic() +
  geom_abline(intercept = 0, slope = 0)
hist(WTC4_temp_data$Tleaf - WTC4_temp_data$TargTempC_Avg)

WTC4_temp_data$DateTime <- as.POSIXct(WTC4_temp_data$DateTime,format="%Y-%m-%d %T",tz="GMT")

temp_data = WTC4_temp_data %>% 
  filter(DateTime >= min(fluxes_df$DateTime_hr) & DateTime <= max(fluxes_df$DateTime_hr)) %>% 
  mutate(DateTime_hr = floor_date(DateTime, unit = "30 minutes")) %>% 
  group_by(chamber, DateTime_hr) %>% 
  summarise(across(c("Tleaf", "TargTempC_Avg"), mean)) %>% 
  ungroup() %>% 
  as.data.frame()


# Clean flux data ##############################################################
# Load flux data
fluxes_df = read.csv(
  "http://hie-pub.westernsydney.edu.au/07fcd9ac-e132-11e7-b842-525400daae48/WTC_TEMP-PARRA_HEATWAVE-FLUX-PACKAGE_L1/data/WTC_TEMP-PARRA_CM_WTCFLUX-CANOPYTEMP_20161029-20161115_L0.csv"
) %>% 
  rename(Tair = Tair_al, PPFD = PAR, Tcan = TargTempC_Avg, A = Photo, E = Trans) %>% 
  select(-PPFD_Avg)

fluxes_df$DateTime_hr <- as.POSIXct(fluxes_df$DateTime_hr,format="%Y-%m-%d %T",tz="GMT")

# Add thermocouple data
temp_data = WTC4_temp_data %>% 
  filter(DateTime >= min(fluxes_df$DateTime_hr) & DateTime <= max(fluxes_df$DateTime_hr)) %>% 
  mutate(DateTime_hr = floor_date(DateTime, unit = "30 minutes")) %>% 
  group_by(chamber, DateTime_hr) %>% 
  summarise(across(c("Tleaf", "TargTempC_Avg"), mean)) %>% 
  ungroup() %>% 
  as.data.frame()
fluxes_df = left_join(fluxes_df, temp_data, by = c("chamber", "DateTime_hr"))

# Test that radiometer temperatures match - yes
plot(fluxes_df$Tcan - fluxes_df$TargTempC_Avg)
fluxes_df = fluxes_df %>% dplyr::select(!TargTempC_Avg)

fluxes_df$Tdiff <- with(fluxes_df,Tcan-Tair)

# Calculate conductance as the simple ratio between Trans and VPD. But use leaf to air VPD. Leaf to air VPD tends to be greater than air VPD, as leaves tend to be warmer than air. Plot the relationships between photosynthesis, conductance, and leaftoairVPD. Note the shift in points at extreme temperatures and VPD in teh heatwave treatment, consistent with lower photosynthesis than would be expected given the measured gs.
fluxes_df$Dleaf <- VPDairToLeaf(VPD=fluxes_df$VPD,Tair=fluxes_df$Tair,Tleaf=fluxes_df$Tcan)
fluxes_df$gs <- fluxes_df$E/fluxes_df$Dleaf/10
fluxes_df$Atogs <- with(fluxes_df,A/(gs*1000))

# Aggregate to hourly averages
fluxes_df$DateTime_hr <- HIEv::nearestTimeStep(fluxes_df$DateTime_hr,nminutes=30,"floor")
dat.hr <- doBy::summaryBy(.~DateTime_hr+chamber+T_treatment+HWtrt+combotrt,FUN=mean,keep.names=T,data=fluxes_df)


# Combine flux and SWC data 
fluxes_df$sw = match_data(fluxes_df, SWC_df, "sw")
fluxes_df$sw0 = match_data(fluxes_df, SWC_df, "sw0")

# Subset to just the heatwave time periods
starttime <- as.POSIXct("2016-10-20 00:00:00",format="%Y-%m-%d %T",tz="GMT")
endtime <- as.POSIXct("2016-11-11 20:00:00",format="%Y-%m-%d %T",tz="GMT")
WTC4_data <- subset(fluxes_df,DateTime_hr>starttime & DateTime_hr<endtime) 


# Estimate soil water potential ################################################
WP_df = read.csv(
  "http://hie-pub.westernsydney.edu.au/07fcd9ac-e132-11e7-b842-525400daae48/WTC_TEMP-PARRA_HEATWAVE-FLUX-PACKAGE_L1/data/WTC_TEMP-PARRA_CM_WATERPOTENTIAL-HEATWAVE_20161019-20161107_L0.csv"
)

# Filter predawn data
predawn_df = WP_df %>% filter(Timing == "pre-dawn")
predawn_df$Date = as.Date(predawn_df$Date)
predawn_df$day = as.numeric(as.Date(predawn_df$Date) - as.Date(starttime))
  
# Filter out pre-drought measure - don't do, causes most fits to have positive slope
#predawn_df = predawn_df %>% filter(Date != "2016-10-19")

# Fit linear regression through predawn data to estimate soil water potential
predawn_fits = plyr::ddply(predawn_df, "chamber", function(x) {
  model = lm(LWP ~ day, data = x)
  coef(model)
})

# Apply regressions to estimate Ps over the course of the heatwave
predawn_est = data.frame(date = rep(seq(starttime, endtime, by = "days"), times = 12),
                         chamber = rep(unique(WP_df$chamber), each = 23),
                         Ps = as.numeric(NA),
                         stringsAsFactors = FALSE
                         )
predawn_est$day = as.numeric(as.Date(predawn_est$date) - as.Date(starttime))

predawn_est$Ps = sapply(1:nrow(predawn_est), function(i) {
  intercept = predawn_fits$`(Intercept)`[predawn_fits$chamber == predawn_est$chamber[i]] 
  slope = predawn_fits$day[predawn_fits$chamber == predawn_est$chamber[i]] 
  Ps = intercept + slope * predawn_est$day[i]
  return(Ps)
})

## Fig S2 ---------------------------
predawn_LWP =  
  predawn_df %>% 
  filter(tissue == "leaf") %>% 
  rename(Ps = LWP) %>% 
  mutate(date = as.POSIXct(Date)) %>% 
  as.data.frame()

FigS2 = ggplot(NULL, aes(x = date, y = Ps, color = chamber)) +
  geom_point(data = predawn_LWP) +
  geom_line(data = predawn_est) +
  theme_classic() +
  guides(color = guide_legend(title = "Chamber")) +
  xlab("Date") +
  ylab(expression("Predawn  "*psi[leaf]*" (MPa)")) +
  theme(text = element_text(size = 14))
ggsave("figs/FigS2_Pleaf_timeseries.tiff", FigS2, height = 6, width = 10)

# Combine with other WTC data
WTC4_data$Ps = sapply(1:nrow(WTC4_data), function(i) {
  date = as.Date(WTC4_data$DateTime_hr[i])
  j = which(as.Date(predawn_est$date) == date & predawn_est$chamber == WTC4_data$chamber[i])
  Ps = predawn_est$Ps[j] * -1
  return(Ps)
})

WTC4_data  %>%  ggplot(aes(x = DateTime_hr, y = -Ps, color = chamber)) + geom_point() + theme_classic()
# Chambers 7, 9, 11 have increasing trend in predawn LWP, but decreasing SWC
WTC4_data %>% ggplot(aes(x = DateTime_hr, y = sw, color = chamber)) + geom_line() + theme_classic()
WTC4_data %>% 
  filter(chamber %in% c("C02", "C07", "C09", "C11")) %>% 
  ggplot(aes(x = DateTime_hr, y = sw0, color = chamber)) + geom_line() + theme_classic()
predawn_df %>% 
  filter(chamber %in% c("C02", "C07", "C09", "C11")) %>% 
  ggplot(aes(x = Date, y = LWP, color = chamber)) + geom_point() + theme_classic()

# Hold constant at lowest measured value
WTC4_data$Ps[WTC4_data$chamber == "C02"] = max(WTC4_data$Ps[WTC4_data$chamber == "C02"])
WTC4_data$Ps[WTC4_data$chamber == "C07"] = max(WTC4_data$Ps[WTC4_data$chamber == "C07"])
WTC4_data$Ps[WTC4_data$chamber == "C09"] = max(WTC4_data$Ps[WTC4_data$chamber == "C09"])
WTC4_data$Ps[WTC4_data$chamber == "C11"] = max(WTC4_data$Ps[WTC4_data$chamber == "C11"])

# Write csv of final data frame for analysis
write.csv(WTC4_data, "data/in/WTC4_data.csv", row.names = FALSE)
WTC4_data = read.csv("data/in/WTC4_data.csv")

# Fit kmax
WP_kmax = WP_df %>% filter(Phase == "baseline", tissue == "leaf") %>%
  select(chamber, Date, Timing, LWP) %>%
  group_by(chamber, Date, Timing) %>%
  summarise(LWP = mean(LWP)) %>%
  ungroup() %>%
  pivot_wider(names_from = Timing, values_from = LWP) %>%
  rename(P_MD = midday, P_PD = "pre-dawn")

LeafArea_df = fluxes_df[1:12,] %>% select(chamber, LeafArea)

E_kmax = fluxes_v3 %>% 
  filter(linktime >= "2016-10-19 11:00:00" & linktime <= "2016-10-19 14:00:00") %>% 
    left_join(LeafArea_df, by = "chamber") %>% 
  select(chamber, FluxH2O, LeafArea, Tair) %>%
  mutate(E = FluxH2O / LeafArea * 10^3) %>%
  group_by(chamber) %>%
  slice_max(E)

kmax_df = WP_kmax %>%
  mutate(E = E_kmax$E, Tair = E_kmax$Tair) %>%
  mutate(kmax = E / (P_PD - P_MD)) %>%
  select(chamber, kmax)

write.csv(kmax_df, "data/in/kmax_values.csv", row.names = FALSE)