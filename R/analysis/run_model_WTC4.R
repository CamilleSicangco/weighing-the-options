# Model testing with the WTC4 dataset
# by Camille Sicangco
# Created 30 Sept 2024

# Load functions
source("R/functions/analysis_functions.R")

# Load data frame
#WTC4_data = read.csv("data/in/WTC4_data.csv")
#WTC4_data$DateTime_hr <- as.POSIXct(WTC4_data$DateTime_hr,format="%Y-%m-%d %T",tz="GMT")

# Control data #################################################################

# Subset the model predictions for the ambient treatment  
control <- subset(WTC4_data,HWtrt=="C" & PAR > 500)

## Fit the Medlyn model ########################################################

#- fit g1 and g0 for the "control" dataset
gs_fits <- nls(gs ~ g0+1.6*(1+g1/sqrt(VPD))*(Photo/400),
               start=list(g0=0.002,g1=4),
               data=subset(control,PAR>500),
               algorithm="port",
               lower=c(0,0),upper=c(0.003,10))

#g0 = unname(coef(gs_fits)[1])
g1 = unname(coef(gs_fits)[2])
#g0=0
g1 = 2.9
g0=0.003

# Fit the control treatment (i.e., not the heatwave trees)
pred1.c <- PhotosynEB(Tair=control$Tair_al,VPD=control$VPD,Wind=8,Wleaf=0.01,StomatalRatio=1,
                    LeafAbs=0.86,
                    PPFD=control$PAR,g1=g1,g0=g0,
                    Vcmax=34,EaV=51780,EdVC=2e5,delsC=640,
                    Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633)
table(pred1.c$failed) # Some energy balance calculations failed
pred1.c$Tdiff <- with(pred1.c,Tleaf-Tair)

## Fit Sperry and Sicangco models ##############################################

# Add kmax values to the data frame
kmax_df = read.csv("data/in/kmax_values.csv")
control = left_join(control, kmax_df, by = "chamber")
#test_s = control[1:20,] # Subset of dataset for testing

# Generate predictions
pred2.c = make_pred(control) # Really slow! Try to optimize

# Check if either model predicts negative gs values
table(pred2.c$gs[pred2.c$Model == "Sicangco"] >= 0) # 158/2341 have negative gs
table(pred2.c$gs[pred2.c$Model == "Sperry"] >= 0) # All predict positive gs
summary(control$PAR[which(pred2.c$gs < 0)]) # Happens for low PAR (< 52.5)

# Save control observations and predictions
save(control, pred1.c, pred2.c, file = "data/out/control_runs.Rdata")
load("data/out/control_runs.Rdata")

## Plot results ################################################################

# Plot Tleaf predictions minus observations for all models
Tdiff.Sic = subset(pred2, Model == "Sicangco")$Tleaf - control$TargTempC_Avg
Tdiff.Sp = subset(pred2, Model == "Sperry")$Tleaf - control$TargTempC_Avg
Tdiff.Med = pred1$Tleaf - control$TargTempC_Avg

data.frame(DateTime_hr = control$DateTime_hr, 
           Tdiff = c(Tdiff.Sic, Tdiff.Sp, Tdiff.Med), 
           model = rep(c("Sicangco", "Sperry", "Medlyn"), each = nrow(control))) %>% 
  ggplot(aes(x = Tdiff, color = model)) +
  geom_density(alpha = 0.2, size = 1) +
  #geom_histogram(aes(y=..density..), alpha=0.5, 
  #               position="identity") +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed")

# Plot transpiration versus canopy temperature 
# E too high for Sperry and Sicangco models, especially as Tleaf increases
# Notably this is exaggerated for chamber 11, which has an increasing trend in SWP
data.frame(E = c(control$Trans, pred1$ELEAF, pred2$E),
           Tleaf = c(control$TargTempC_Avg, pred1$Tleaf, pred2$Tleaf),
           model = rep(c("observed", "Medlyn", "Sicangco", "Sperry"), each = nrow(control)),
           chamber = rep(control$chamber, times = 4)) %>% 
  filter(model =="Sperry") %>% 
  ggplot(aes(x = Tleaf, y = E, color = chamber)) +
  geom_point(alpha = .5, size = 1) +
  theme_classic() #+
  scale_color_manual(
    values = c("observed" = "black", "Medlyn" = "hotpink2", 
               "Sicangco" = "darkorange", "Sperry" = "blue"))
pred2 %>% ggplot(aes(x = Tleaf, y = E, color = Model)) + geom_point() + theme_classic()

# Heatwave data ################################################################

# Subset the model predictions for the ambient treatment  
heatwave <- subset(WTC4_data,HWtrt=="HW" & PAR > 500)

## Fit the Medlyn model ########################################################

# Fit the heatwave treatment (i.e., not the heatwave trees)
pred1.hw <- PhotosynEB(Tair=heatwave$Tair_al,VPD=heatwave$VPD,Wind=8,Wleaf=0.01,StomatalRatio=1,
                    LeafAbs=0.86,
                    PPFD=heatwave$PAR,g1=g1,g0=g0,
                    Vcmax=34,EaV=51780,EdVC=2e5,delsC=640,
                    Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633)
table(pred1.hw$failed) # Some energy balance calculations failed
pred1.hw$Tdiff <- with(pred1.hw,Tleaf-Tair)

## Fit Sperry and Sicangco models ##############################################

# Add kmax values to the data frame
heatwave = left_join(heatwave, kmax_df, by = "chamber")
#test_s = control[1:20,] # Subset of dataset for testing

# Generate predictions
pred2.hw = make_pred(heatwave) # Really slow! Try to optimize

# Check if either model predicts negative gs values
table(pred2$gs[pred2$Model == "Sicangco"] >= 0) # No negative gs
table(pred2$gs[pred2$Model == "Sperry"] >= 0) # No negative gs

# Save heatwave observations and predictions
save(heatwave, pred1.hw, pred2.hw, file = "data/out/heatwave_runs.Rdata")
load("data/out/heatwave_runs.Rdata")

## Plot results ################################################################ 

# Plot transpiration versus canopy temperature 
# E too high for Sperry and Sicangco models, but trend is good
# Notably this is exaggerated for chambers 7 and 9, which has an "increasing" trend in SWP
data.frame(E = c(heatwave$Trans, pred1.hw$ELEAF, pred2.hw$E),
           Tleaf = c(heatwave$TargTempC_Avg, pred1.hw$Tleaf, pred2.hw$Tleaf),
           model = rep(c("observed", "Medlyn", "Sicangco", "Sperry"), each = nrow(heatwave)),
           chamber = rep(heatwave$chamber, times = 4)) %>% 
  #filter(model =="Sperry") %>% 
  ggplot(aes(x = Tleaf, y = E, color = model)) +
  geom_point(alpha = .5, size = 1) +
  theme_classic() +
scale_color_manual(
  values = c("observed" = "black", "Medlyn" = "hotpink2", 
             "Sicangco" = "darkorange", "Sperry" = "blue"))
pred2.hw %>% ggplot(aes(x = Tleaf, y = E, color = Model)) + geom_point() + theme_classic()
