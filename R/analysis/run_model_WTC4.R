# Model testing with the WTC4 dataset
# by Camille Sicangco
# Created 30 Sept 2024

# Load functions
source("R/functions/analysis_functions.R")

# Load data framehttp://127.0.0.1:29759/graphics/plot_zoom_png?width=515&height=769
#WTC4_data = read.csv("data/in/WTC4_data.csv")
#WTC4_data$DateTime_hr <- as.POSIXct(WTC4_data$DateTime_hr,format="%Y-%m-%d %T",tz="GMT")


# Control data #################################################################

# Subset the model predictions for the ambient treatment  
control <- subset(WTC4_data,HWtrt=="C" & PPFD > 500)

# Test energy balance by prescribing E
pred_prescE.c = prescribedE_pred(df = filter(control, !is.na(gs)))

diff_long.c = pred_prescE.c %>% 
  select(ends_with(".diff")) %>%
  pivot_longer(cols = ends_with(".diff"), names_to = "var", values_to = "diff",
               names_pattern = "(.*).diff")
obs_long.c = pred_prescE.c %>% 
  select(ends_with(".obs")) %>%
  select(!E.obs) %>% 
  pivot_longer(cols = ends_with("obs"), names_to = "var", values_to = "obs",
               names_pattern = "(.*).obs")
prescE_long.c = cbind(diff_long.c, obs_long.c[2])

prescE_long.c %>% 
  ggplot(aes(x = obs, y = diff)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, color = "blue") +
  facet_wrap(vars(var), scales = "free") +
  theme_classic() +
  ggtitle("Control")

## Fit GAM #####################################################################
gam_E.c = gam(E ~ s(Tcan), data = control)
gam_A.c = gam(A ~ s(Tcan), data = control)


## Fit the Medlyn model ########################################################

#- fit g1 and g0 for the "control" dataset
gs_fits <- nls(gs ~ g0+1.6*(1+g1/sqrt(VPD))*(Photo/400),
               start=list(g0=0.002,g1=4),
               data=subset(control,PPFD>500),
               algorithm="port",
               lower=c(0,0),upper=c(0.003,10))

#g0 = unname(coef(gs_fits)[1])
g1 = unname(coef(gs_fits)[2])
#g0=0
g1 = 2.9
g0=0.003

# Fit the control treatment (i.e., not the heatwave trees)
Rd0 = 0.92
TrefR = 25
Rd = Rd0 * exp(0.1012 * (control$Tcan - TrefR) - 5e-04 * (control$Tcan^2 - TrefR^2))
pred1.c <- PhotosynEB(Tair=control$Tair,VPD=control$VPD,Wind=8,Wleaf=0.01,StomatalRatio=1,
                    LeafAbs=0.5,
                    PPFD=control$PPFD,g1=g1,g0=g0,
                    Vcmax=34,EaV=51780,EdVC=2e5,delsC=640,
                    Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633)
table(pred1.c$failed) # Some energy balance calculations failed
pred1.c$Tdiff <- with(pred1.c,Tleaf-Tair)

## Fit Sperry and Sicangco models ##############################################

# Add kmax values to the data frame
#kmax_df = read.csv("data/in/kmax_values.csv") 
#control = left_join(control, kmax_df, by = "chamber")
control = control %>% mutate(kmax = 0.7)
#test_s = control[1:20,] # Subset of dataset for testing

# Generate predictions
pred2.c = make_pred(df = control, models = "final", LeafAbs = 0.5, 
                    Tcrit_hw = 42.1, T50_hw = 48.4) # Really slow! Try to optimize
#pred2.c = make_pred(df = control, models = "final", Tcrit = 46, T50 = 48) # Really slow! Try to optimize

# Check if either model predicts negative gs values
table(pred2.c$gs[pred2.c$Model == "Sicangco"] >= 0) # 158/2341 have negative gs
table(pred2.c$gs[pred2.c$Model == "Sperry"] >= 0) # All predict positive gs
summary(control$PPFD[which(pred2.c$gs < 0)]) # Happens for low PPFD (< 52.5)

### Fit intermediate models, i.e. Sperry with either CGnet or TC ###############

# Generate predictions
pred3.c = make_pred(df = control, models = "intermediate") # Really slow! Try to optimize
# Check if either model predicts negative gs values
table(pred3.c$gs[pred3.c$Model == "Sperry + CGnet"] >= 0) # 158/2341 have negative gs
table(pred3.c$gs[pred3.c$Model == "Sperry + TC"] >= 0) # All predict positive gs
summary(control$PPFD[which(pred3.c$gs < 0)]) # Happens for low PPFD (< 52.5)

# Save control observations and predictions
save(control, pred1.c, pred2.c, file = "data/out/control_runs_kmaxpt7_fittedTcritT50.Rdata")
load("data/out/control_runs_kmax1.Rdata")

## Plot results ################################################################

# Plot Tleaf predictions minus observations for all models
Tdiff.Sic = subset(pred2.c, Model == "Sicangco")$Tleaf - control$Tcan
Tdiff.Sp = subset(pred2.c, Model == "Sperry")$Tleaf - control$Tcan
Tdiff.Med = pred1.c$Tleaf - control$Tcan

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
out.c = data.frame(datetime = rep(control$DateTime_hr, times = 4),
                   chamber = rep(control$chamber, times = 4),
                   model = rep(c("observed", "Medlyn", "Sicangco", "Sperry"), each = nrow(control)),
                   E = c(control$E, pred1.c$ELEAF, pred2.c$E),
                   A = c(control$A, pred1.c$ALEAF, pred2.c$A),
                   gs = c(control$gs, pred1.c$gw, pred2.c$gs),
                   Dleaf = c(control$Dleaf, pred1.c$VPDleaf, pred2.c$Dleaf),
                   Tleaf = c(control$Tcan, pred1.c$Tleaf, pred2.c$Tleaf))

out.c$model = factor(out.c$model, levels = c("Sicangco", "Sperry", "Medlyn", "observed"))
# Plot transpiration versus canopy temperature 
# E too high for Sperry and Sicangco models, especially as Tleaf increases
# Notably this is exaggerated for chamber 11, which has an increasing trend in SWP
EvT.c = plotGAM(gam_E.c, smooth.c = "Tcan") +
  geom_point(data = out.c, aes(x = Tleaf, y = E, color = model), alpha = .5, size = .5) +
  theme_classic() +
  scale_color_manual(
    values = c("observed" = "black", "Medlyn" = "#D81B60", 
               "Sicangco" = "#FFC107", "Sperry" = "#1E88E5")) +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  ylab(expression("E"*" (mmol m"^-2*"s"^-1*")")) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)),
         linetype = "none") +
  theme(plot.title = element_blank())

# Plot A versus canopy temperature 
AvT.c = plotGAM(gam_A.c, smooth.c = "Tcan") +
  geom_point(data = out.c, aes(x = Tleaf, y = A, color = model), alpha = .5, size = .5) +
  theme_classic() +
  scale_color_manual(
    values = c("observed" = "black", "Medlyn" = "#D81B60", 
               "Sicangco" = "#FFC107", "Sperry" = "#1E88E5")) +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  ylab(expression("A (" * mu * "mol m"^-2*"s"^-1*")")) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)),
         linetype = "none") +
  theme(plot.title = element_blank())

# Plot Tleaf predictions vs observations
Tleaf_pred_obs.plt = out.c %>% 
  pivot_wider(names_from = model, values_from = Tleaf, id_cols = c(chamber, datetime)) %>% 
  pivot_longer(cols = Medlyn:Sperry, names_to = "model", values_to = "Tleaf_pred") %>% 
  rename(Tleaf = observed) %>% 
  ggplot() +
  geom_point(aes(x = Tleaf, y = Tleaf_pred, color = model), alpha = .5) +
  theme_classic() +
  scale_color_manual(
    values = c("Medlyn" = "#D81B60", 
               "Sicangco" = "#FFC107", "Sperry" = "#1E88E5")) +
  xlab(expression("observed T"[leaf]*" (\u00B0C)")) +
  ylab(expression("predicted T"[leaf]*" (\u00B0C)")) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)),
         linetype = "none") +
  theme(plot.title = element_blank()) + 
  geom_abline(slope = 1)

out.c %>% 
  pivot_wider(names_from = model, values_from = Tleaf, id_cols = c(chamber, datetime)) %>% 
  pivot_longer(cols = Medlyn:Sperry, names_to = "model", values_to = "Tleaf_pred") %>% 
  rename(Tleaf = observed) %>% 
  ggplot() +
  geom_point(aes(x = Tleaf, y = Tleaf_pred - Tleaf, color = model), alpha = .5) +
  theme_classic() +
  scale_color_manual(
    values = c("Medlyn" = "#D81B60", 
               "Sicangco" = "#FFC107", "Sperry" = "#1E88E5")) +
  xlab(expression("observed T"[leaf]*" (\u00B0C)")) +
  ylab(expression("predicted - observed T"[leaf]*" (\u00B0C)")) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)),
         linetype = "none") +
  theme(plot.title = element_blank()) + 
  geom_hline(yintercept = 0)


# Heatwave data ################################################################

# Subset the model predictions for the ambient treatment  
heatwave <- subset(WTC4_data,HWtrt=="HW" & PPFD > 500)
pred_prescE.hw = prescribedE_pred(df = filter(heatwave, !is.na(gs)))

diff_long.hw = pred_prescE.hw %>% 
  select(ends_with(".diff")) %>%
  pivot_longer(cols = ends_with(".diff"), names_to = "var", values_to = "diff",
               names_pattern = "(.*).diff")
obs_long.hw = pred_prescE.hw %>% 
  select(ends_with(".obs")) %>%
  select(!E.obs) %>% 
  pivot_longer(cols = ends_with("obs"), names_to = "var", values_to = "obs",
               names_pattern = "(.*).obs")
prescE_long.hw = cbind(diff_long.hw, obs_long.hw[2])

prescE_long.hw %>% 
  ggplot(aes(x = obs, y = diff)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, color = "blue") +
  facet_wrap(vars(var), scales = "free") +
  theme_classic() + ggtitle("Heatwave")

# Test leaf temperature predictions only
Tleaves.hw = try(calc_Tleaf(heatwave$Tair, VPD = heatwave$VPD, E = heatwave$E, 
                           PPFD = heatwave$PPFD, Wind = 8, Wleaf = 0.01, LeafAbs = 0.86))
Tleaves.hw[grep("Error", Tleaves.hw)] = NA
Tleaves.hw = as.numeric(Tleaves.hw)

# Tends to over-predict leaf temperature slightly
Tdiff.hw = Tleaves.hw - heatwave$Tcan
summary(Tdiff.hw)
hist(Tdiff.hw)
plot(x = heatwave$Tcan, y = Tdiff.hw)
abline(h = 0, col = "blue")

## Fit GAM #####################################################################
gam_E.hw = gam(E ~ s(Tcan), data = heatwave)
gam_A.hw = gam(A ~ s(Tcan), data = heatwave)


## Fit the Medlyn model ########################################################

# Fit the heatwave treatment (i.e., not the heatwave trees)
pred1.hw <- PhotosynEB(Tair=heatwave$Tair,VPD=heatwave$VPD,Wind=8,Wleaf=0.01,StomatalRatio=1,
                    LeafAbs=0.5,
                    PPFD=heatwave$PPFD,g1=g1,g0=g0,
                    Vcmax=34,EaV=51780,EdVC=2e5,delsC=640,
                    Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633, Ca = 420)

table(pred1.hw$failed) # Some energy balance calculations failed
pred1.hw$Tdiff <- with(pred1.hw,Tleaf-Tair)

## Fit Sperry and Sicangco models ##############################################

# Add kmax values to the data frame
#heatwave = left_join(heatwave, kmax_df, by = "chamber") #%>% rename(Tair = Tair, PPFD = PPFD)
heatwave = heatwave %>% mutate(kmax = 0.7)

#test_s = control[1:20,] # Subset of dataset for testing

# Generate predictions
pred2.hw = make_pred(df = heatwave, models = "final", 
                     Tcrit_hw = 43.4, T50_hw = 49.6, LeafAbs = 0.5) # Really slow! Try to optimize
#pred2.hw = make_pred(df = heatwave, models = "final", Tcrit = 46, T50 = 48) # Really slow! Try to optimize

plot(pred2.hw$P[pred2.hw$Model == "Sperry"])
plot(heatwave$Ps)

# Check if either model predicts negative gs values
table(pred2$gs[pred2$Model == "Sicangco"] >= 0) # No negative gs
table(pred2$gs[pred2$Model == "Sperry"] >= 0) # No negative gs

### Fit intermediate models, i.e. Sperry with either CGnet or TC ###############

# Generate predictions
pred3.hw = make_pred(df = heatwave, models = "intermediate") # Really slow! Try to optimize
# Check if either model predicts negative gs values
table(pred3.hw$gs[pred3.hw$Model == "Sperry + CGnet"] >= 0) # 158/2341 have negative gs
table(pred3.hw$gs[pred3.hw$Model == "Sperry + TC"] >= 0) # All predict positive gs
summary(control$PPFD[which(pred3.hw$gs < 0)]) # Happens for low PPFD (< 52.5)


# Save heatwave observations and predictions
save(heatwave, pred1.hw, pred2.hw, file = "data/out/heatwave_runs_kmaxpt7_fittedTcritT50.Rdata")
load("data/out/heatwave_runs_kmaxpt8.Rdata")

## Plot results ################################################################ 

# Plot Tleaf predictions minus observations for all models
Tdiff.Sic.hw = subset(pred2.hw, Model == "Sicangco")$Tleaf - heatwave$Tcan
Tdiff.Sp.hw = subset(pred2.hw, Model == "Sperry")$Tleaf - heatwave$Tcan
Tdiff.Med.hw = pred1.hw$Tleaf - heatwave$Tcan

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
out.hw = data.frame(datetime = rep(heatwave$DateTime_hr, times = 4),
                    chamber = rep(heatwave$chamber, times = 4),
                    model = rep(c("observed", "Medlyn", "Sicangco", "Sperry"), each = nrow(heatwave)),
                    E = c(heatwave$E, pred1.hw$ELEAF, pred2.hw$E),
                    A = c(heatwave$A, pred1.hw$ALEAF, pred2.hw$A),
                    gs = c(heatwave$gs, pred1.hw$gw, pred2.hw$gs),
                    Dleaf = c(heatwave$Dleaf, pred1.hw$VPDleaf, pred2.hw$Dleaf),
                    Tleaf = c(heatwave$Tcan, pred1.hw$Tleaf, pred2.hw$Tleaf),
                    Tair = c(heatwave$Tair, pred1.hw$Tair, pred2.hw$Tair),
                    VPD = rep(heatwave$VPD, times = 4),
                    PPFD = rep(heatwave$PPFD, times = 4)) 
out.hw$model = factor(out.hw$model, levels = c("Sicangco", "Sperry", "Medlyn", "observed"))

pred3.hw.reformat = pred3.hw %>% 
  select(!P) %>% 
  rename(model = Model) %>% 
  mutate(chamber = rep(heatwave$chamber, times = 2), .before = 1)
out.hw.all = bind_rows(out.hw, pred3.hw.reformat)


# Plot transpiration versus canopy temperature 
# E too high for Sperry and Sicangco models, but trend is good
# Notably this is exaggerated for chambers 7 and 9, which has an "increasing" trend in SWP
EvT.hw = plotGAM(gam_E.hw, smooth.c = "Tcan") +
  geom_point(data = out.hw, aes(x = Tleaf, y = E, color = model), alpha = .5, size = .5) +
  theme_classic() +
  scale_color_manual(
    values = c("observed" = "black", "Medlyn" = "#D81B60", 
               "Sicangco" = "#FFC107", "Sperry" = "#1E88E5")) +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  ylab(expression("E"*" (mmol m"^-2*"s"^-1*")")) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)),
         linetype = "none") +
  theme(plot.title = element_blank())

EvT.hw.all = plotGAM(gam_E.hw, smooth.c = "Tcan") +
  geom_point(data = out.hw.all, aes(x = Tleaf, y = E, color = model), alpha = .5, size = .5) +
  theme_classic() +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  ylab(expression("E"*" (mmol m"^-2*"s"^-1*")")) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)),
         linetype = "none") +
  theme(plot.title = element_blank()) +
  scale_color_manual(
    values = c("observed" = "black", "Medlyn" = "#D81B60", 
               "Sicangco" = "#FFC107", "Sperry" = "#1E88E5",
               "Sperry + CGnet" = "#AA3377", "Sperry + TC" = "#228833"))
ggsave(plot = EvT.hw.all, filename = "figs/EvT_hw_all.pdf", width = 8, height = 7)


# Plot A versus canopy temperature 
AvT.hw = plotGAM(gam_A.hw, smooth.c = "Tcan") +
  geom_point(data = out.hw, aes(x = Tleaf, y = A, color = model), alpha = .5, size = .5) +
  theme_classic() +
  scale_color_manual(
    values = c("observed" = "black", "Medlyn" = "#D81B60", 
               "Sicangco" = "#FFC107", "Sperry" = "#1E88E5")) +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  ylab(expression("A (" * mu * "mol m"^-2*"s"^-1*")")) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)),
         linetype = "none") +
  theme(plot.title = element_blank())

# Combine all plots (i.e. recreate Drake et al. Fig 5 with models)
AEvT.plt = ggarrange(AvT.c + ylim(-2.5, 13) + ggtitle("Control") + theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)), 
          AvT.hw + ylim(-2.5, 13) + ggtitle("Heatwave") + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)), 
          EvT.c + ylim(-1, 4) + theme(axis.title.x = element_blank()), 
          EvT.hw + ylim(-1,4) + theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
          nrow = 2, ncol = 2, common.legend = TRUE, legend = "right")
AEvT.plt = annotate_figure(AEvT.plt,
                bottom = text_grob(expression("T"[canopy]*" (\u00B0C)")))
AEvT.plt
ggsave(plot = AEvT.plt, filename = "figs/AEvT_plt_kmaxpt7_fittedTcritT50.pdf", width = 8, height = 7)


pred2.hw %>% 
  ggplot(aes(x = Tleaf, y = -P, color = Model)) +
  geom_point() +
  theme_classic()

Tleaf_pred_obs.plt = out.hw %>% 
  pivot_wider(names_from = model, values_from = Tleaf, id_cols = c(chamber, datetime)) %>% 
  pivot_longer(cols = Medlyn:Sperry, names_to = "model", values_to = "Tleaf_pred") %>% 
  rename(Tleaf = observed) %>% 
  ggplot() +
  geom_point(aes(x = Tleaf, y = Tleaf_pred, color = model), alpha = .5) +
  theme_classic() +
  scale_color_manual(
    values = c("Medlyn" = "#D81B60", 
               "Sicangco" = "#FFC107", "Sperry" = "#1E88E5")) +
  xlab(expression("observed T"[leaf]*" (\u00B0C)")) +
  ylab(expression("predicted T"[leaf]*" (\u00B0C)")) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)),
         linetype = "none") +
  theme(plot.title = element_blank()) + 
  geom_abline(slope = 1)

out.hw %>% 
  pivot_wider(names_from = model, values_from = Tleaf, id_cols = c(chamber, datetime)) %>% 
  pivot_longer(cols = Medlyn:Sperry, names_to = "model", values_to = "Tleaf_pred") %>% 
  rename(Tleaf = observed) %>% 
  ggplot() +
  geom_point(aes(x = Tleaf, y = Tleaf_pred - Tleaf, color = model), alpha = .5) +
  theme_classic() +
  scale_color_manual(
    values = c("Medlyn" = "#D81B60", 
               "Sicangco" = "#FFC107", "Sperry" = "#1E88E5")) +
  xlab(expression("observed T"[leaf]*" (\u00B0C)")) +
  ylab(expression("predicted - observed T"[leaf]*" (\u00B0C)")) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)),
         linetype = "none") +
  theme(plot.title = element_blank()) + 
  geom_hline(yintercept = 0)

bind_rows(out.c, out.hw) %>% 
  #filter(model == "Medlyn") %>% 
  pivot_wider(names_from = model, values_from = Tleaf, id_cols = c(chamber, datetime)) %>% 
  #pivot_longer(cols = Medlyn:Sperry, names_to = "model", values_to = "Tleaf_pred") %>% 
  rename(Tleaf = observed, Tleaf_pred = Medlyn) %>% 
  ggplot() +
  geom_point(aes(x = Tleaf, y = Tleaf_pred - Tleaf), alpha = .5) +
  theme_classic() +
  scale_color_manual(
    values = c("Medlyn" = "#D81B60", 
               "Sicangco" = "#FFC107", "Sperry" = "#1E88E5")) +
  xlab(expression("observed T"[leaf]*" (\u00B0C)")) +
  ylab(expression("predicted - observed T"[leaf]*" (\u00B0C)")) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)),
         linetype = "none") +
  theme(plot.title = element_blank()) + 
  geom_hline(yintercept = 0)


heatwave %>% 
  ggplot(aes(x = Tair, y = Tleaf)) +
  geom_point() +
  theme_classic()+ 
  geom_abline(slope = 1)

out.hw %>% 
  ggplot(aes(x = model, y = Tleaf, fill = model)) +
  geom_boxplot() +
  geom_violin() +
  scale_color_manual(
    values = c("observed" = "darkgrey", "Medlyn" = "#D81B60", 
               "Sicangco" = "#FFC107", "Sperry" = "#1E88E5")) +
  theme_classic() +
  geom_hline(yintercept = 46, linetype = "dashed") #+
  annotate("text", y=Tcrit, label="Tcrit", vjust = -0.5)
  