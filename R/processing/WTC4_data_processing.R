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


# Clean flux data ##############################################################
# Load flux data
fluxes_df = read.csv(
  "http://hie-pub.westernsydney.edu.au/07fcd9ac-e132-11e7-b842-525400daae48/WTC_TEMP-PARRA_HEATWAVE-FLUX-PACKAGE_L1/data/WTC_TEMP-PARRA_CM_WTCFLUX-CANOPYTEMP_20161029-20161115_L0.csv"
)

fluxes_df$DateTime_hr <- as.POSIXct(fluxes_df$DateTime_hr,format="%Y-%m-%d %T",tz="GMT")
fluxes_df$Tdiff <- with(fluxes_df,TargTempC_Avg-Tair_al)

# Calculate conductance as the simple ratio between Trans and VPD. But use leaf to air VPD. Leaf to air VPD tends to be greater than air VPD, as leaves tend to be warmer than air. Plot the relationships between photosynthesis, conductance, and leaftoairVPD. Note the shift in points at extreme temperatures and VPD in teh heatwave treatment, consistent with lower photosynthesis than would be expected given the measured gs.
fluxes_df$Dleaf <- VPDairToLeaf(VPD=fluxes_df$VPD,Tair=fluxes_df$Tair_al,Tleaf=fluxes_df$TargTempC_Avg)
fluxes_df$gs <- fluxes_df$Trans/fluxes_df$Dleaf/10
fluxes_df$Atogs <- with(fluxes_df,Photo/(gs*1000))

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

# Combine with other WTC data
WTC4_data$Ps = sapply(1:nrow(WTC4_data), function(i) {
  date = as.Date(WTC4_data$DateTime_hr[i])
  j = which(as.Date(predawn_est$date) == date & predawn_est$chamber == WTC4_data$chamber[i])
  Ps = predawn_est$Ps[j] * -1
  return(Ps)
})

# Chambers 7, 9, 11 have increasing trend in predawn LWP, but decreasing SWC
WTC4_data %>% ggplot(aes(x = DateTime_hr, y = Ps, color = chamber)) + geom_line() + theme_classic()
WTC4_data %>% 
  filter(chamber %in% c("C07", "C09", "C11")) %>% 
  ggplot(aes(x = DateTime_hr, y = sw0, color = chamber)) + geom_line() + theme_classic()
predawn_df %>% 
  filter(chamber %in% c("C07", "C09", "C11")) %>% 
  ggplot(aes(x = Date, y = LWP, color = chamber)) + geom_point() + theme_classic()

# Write csv of final data frame for analysis
write.csv(WTC4_data, "data/in/WTC4_data.csv", row.names = FALSE)

# Alternative estimation of Ps from soil water retention curve
#theta_sat = 0.185 # m3 m-3
#bch = 4
#Psie = - 0.003 # MPa
#WTC4_data$Ps = Psie * (WTC4_data$sw0 / theta_sat)**(-bch)

# Fit kmax for Sperry and Sicangco models
## Fit kmax ####################################################################
WP_df = read.csv(
  "http://hie-pub.westernsydney.edu.au/07fcd9ac-e132-11e7-b842-525400daae48/WTC_TEMP-PARRA_HEATWAVE-FLUX-PACKAGE_L1/data/WTC_TEMP-PARRA_CM_WATERPOTENTIAL-HEATWAVE_20161019-20161107_L0.csv"
)
fluxes_v3 = read.csv("data/in/raw/WTC_TEMP-PARRA_WTCFLUX-ALLVARS_20160228-20161123_L0.csv")

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
  select(chamber, FluxH2O, LeafArea) %>%
  mutate(E = FluxH2O / LeafArea * 10^3) %>%
  group_by(chamber) %>%
  summarize(E = max(E))

kmax_df = WP_kmax %>%
  mutate(E = E_kmax$E) %>%
  mutate(kmax = E / (P_PD - P_MD)) %>%
  select(chamber, kmax)
write.csv(kmax_df, "data/in/kmax_values.csv", row.names = FALSE)
