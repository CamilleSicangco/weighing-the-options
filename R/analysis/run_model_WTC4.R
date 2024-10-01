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

# Test leaf temperature predictions only
Tleaves.c = try(calc_Tleaf(control$Tair_al, VPD = control$VPD, E = control$Trans, 
           PPFD = control$PAR, Wind = 8, Wleaf = 0.01, LeafAbs = 0.86))
Tleaves.c[grep("Error", Tleaves.c)] = NA
Tleaves.c = as.numeric(Tleaves.c)

# Tends to over-predict leaf temperature slightly
Tdiff.c = Tleaves.c - control$TargTempC_Avg
summary(Tdiff.c)
hist(Tdiff.c)
plot(x = control$TargTempC_Avg, y = Tdiff.c)
abline(h = 0, col = "blue")

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
Tdiff.Sic = subset(pred2.c, Model == "Sicangco")$Tleaf - control$TargTempC_Avg
Tdiff.Sp = subset(pred2.c, Model == "Sperry")$Tleaf - control$TargTempC_Avg
Tdiff.Med = pred1.c$Tleaf - control$TargTempC_Avg

# Plot Tdiff by model
data.frame(DateTime_hr = control$DateTime_hr, 
           Tdiff = c(Tdiff.Sic, Tdiff.Sp, Tdiff.Med), 
           model = rep(c("Sicangco", "Sperry", "Medlyn"), each = nrow(control))) %>% 
  ggplot(aes(x = Tdiff, color = model)) +
  geom_density(alpha = 0.2, size = 1) +
  #geom_histogram(aes(y=..density..), alpha=0.5, 
  #               position="identity") +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed")

# Create summary data frame
out.c = data.frame(chamber = rep(control$chamber, times = 4),
                   model = rep(c("observed", "Medlyn", "Sicangco", "Sperry"), each = nrow(control)),
                   E = c(control$Trans, pred1.c$ELEAF, pred2.c$E),
                   A = c(control$Photo, pred1.c$ALEAF, pred2.c$A),
                   gs = c(control$gs, pred1.c$gw, pred2.c$gs),
                   Dleaf = c(control$Dleaf, pred1.c$VPDleaf, pred2.c$Dleaf),
                   Tleaf = c(control$TargTempC_Avg, pred1.c$Tleaf, pred2.c$Tleaf))

out.c$model = factor(out.c$model, levels = c("Sicangco", "Sperry", "Medlyn", "observed"))
# Plot transpiration versus canopy temperature 
# E too high for Sperry and Sicangco models, especially as Tleaf increases
# Notably this is exaggerated for chamber 11, which has an increasing trend in SWP
EvT.c = out.c %>% 
  #filter(model =="Sperry") %>% 
  ggplot(aes(x = Tleaf, y = E, color = model)) +
  geom_point(alpha = .5, size = 1) +
  theme_classic() +
  scale_color_manual(
    values = c("observed" = "black", "Medlyn" = "#D81B60", 
               "Sicangco" = "#FFC107", "Sperry" = "#1E88E5"))  +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  ylab(expression("E"*" (mmol m"^-2*"s"^-1*")")) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) #+
  #scale_shape_manual(values = c("observed" = 16, "Medlyn" = 3, "Sicangco" = 3, "Sperry" = 3))

  
# Plot photosynthesis versus canopy temperature 
AvT.c = out.c %>% 
  #filter(model =="Sperry") %>% 
  ggplot(aes(x = Tleaf, y = A, color = model)) +
  geom_point(alpha = .5, size = 1) +
  theme_classic() +
scale_color_manual(
  values = c("observed" = "black", "Medlyn" = "#D81B60", 
             "Sicangco" = "#FFC107", "Sperry" = "#1E88E5")) +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  ylab(expression("A (" * mu * "mol m"^-2*"s"^-1*")")) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))# +
  #scale_shape_manual(values = c("observed" = 16, "Medlyn" = 3, "Sicangco" = 3, "Sperry" = 3))

  
# Heatwave data ################################################################

# Subset the model predictions for the ambient treatment  
heatwave <- subset(WTC4_data,HWtrt=="HW" & PAR > 500)

# Test leaf temperature predictions only
Tleaves.hw = try(calc_Tleaf(heatwave$Tair_al, VPD = heatwave$VPD, E = heatwave$Trans, 
                           PPFD = heatwave$PAR, Wind = 8, Wleaf = 0.01, LeafAbs = 0.86))
Tleaves.hw[grep("Error", Tleaves.hw)] = NA
Tleaves.hw = as.numeric(Tleaves.hw)

# Tends to over-predict leaf temperature slightly
Tdiff.hw = Tleaves.hw - heatwave$TargTempC_Avg
summary(Tdiff.hw)
hist(Tdiff.hw)
plot(x = heatwave$TargTempC_Avg, y = Tdiff.hw)
abline(h = 0, col = "blue")

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

# Plot Tleaf predictions minus observations for all models
Tdiff.Sic.hw = subset(pred2.hw, Model == "Sicangco")$Tleaf - heatwave$TargTempC_Avg
Tdiff.Sp.hw = subset(pred2.hw, Model == "Sperry")$Tleaf - heatwave$TargTempC_Avg
Tdiff.Med.hw = pred1.hw$Tleaf - heatwave$TargTempC_Avg

# Plot Tdiff by model
data.frame(DateTime_hr = heatwave$DateTime_hr, 
           Tdiff = c(Tdiff.Sic.hw, Tdiff.Sp.hw, Tdiff.Med.hw), 
           model = rep(c("Sicangco", "Sperry", "Medlyn"), each = nrow(heatwave))) %>% 
  ggplot(aes(x = Tdiff, color = model)) +
  geom_density(alpha = 0.2, size = 1) +
  #geom_histogram(aes(y=..density..), alpha=0.5, 
  #               position="identity") +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed")

# Create summary data frame
out.hw = data.frame(chamber = rep(heatwave$chamber, times = 4),
                    model = rep(c("observed", "Medlyn", "Sicangco", "Sperry"), each = nrow(heatwave)),
                    E = c(heatwave$Trans, pred1.hw$ELEAF, pred2.hw$E),
                    A = c(heatwave$Photo, pred1.hw$ALEAF, pred2.hw$A),
                    gs = c(heatwave$gs, pred1.hw$gw, pred2.hw$gs),
                    Dleaf = c(heatwave$Dleaf, pred1.hw$VPDleaf, pred2.hw$Dleaf),
                    Tleaf = c(heatwave$TargTempC_Avg, pred1.hw$Tleaf, pred2.hw$Tleaf),
                    Tair = c(heatwave$Tair_al, pred1.hw$Tair, pred2.hw$Tair)) 
out.hw$model = factor(out.hw$model, levels = c("Sicangco", "Sperry", "Medlyn", "observed"))
# Plot transpiration versus canopy temperature 
# E too high for Sperry and Sicangco models, but trend is good
# Notably this is exaggerated for chambers 7 and 9, which has an "increasing" trend in SWP
EvT.hw = out.hw %>% 
  #filter(model =="Sperry") %>% 
  ggplot(aes(x = Tleaf, y = E, color = model)) +
  geom_point(alpha = .5, size = 1) +
  theme_classic() +
  scale_color_manual(
    values = c("observed" = "black", "Medlyn" = "#D81B60", 
               "Sicangco" = "#FFC107", "Sperry" = "#1E88E5")) +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  ylab(expression("E"*" (mmol m"^-2*"s"^-1*")")) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) #+
  #scale_shape_manual(values = c("observed" = 16, "Medlyn" = 3, "Sicangco" = 3, "Sperry" = 3))


# Plot A versus canopy temperature 
AvT.hw = out.hw %>% ggplot(aes(x = Tleaf, y = A, color = model)) +
  geom_point(alpha = .5, size = 1) +
  theme_classic() +
  scale_color_manual(
    values = c("observed" = "black", "Medlyn" = "#D81B60", 
               "Sicangco" = "#FFC107", "Sperry" = "#1E88E5")) +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  ylab(expression("A (" * mu * "mol m"^-2*"s"^-1*")")) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) #+
  #scale_shape_manual(values = c("observed" = 16, "Medlyn" = 3, "Sicangco" = 3, "Sperry" = 3))

# Combine all plots (i.e. recreate Drake et al. Fig 5 with models)
AEvT.plt = ggarrange(AvT.c + ylim(-3, 13) + ggtitle("Control") + theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)), 
          AvT.hw + ylim(-3, 13) + ggtitle("Heatwave") + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)), 
          EvT.c + ylim(-1, 11) + theme(axis.title.x = element_blank()), 
          EvT.hw + ylim(-1,11) + theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
          nrow = 2, ncol = 2, common.legend = TRUE, legend = "right")
AEvT.plt = annotate_figure(AEvT.plt,
                bottom = textGrob(expression("T"[canopy]*" (\u00B0C)"), gp = gpar(cex = 1)))
ggsave(plot = AEvT.plt, filename = "figs/AEvT_plt.pdf")
