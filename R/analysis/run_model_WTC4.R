# Model testing with the WTC4 dataset
# by Camille Sicangco
# Created 30 Sept 2024

# Load data framehttp://127.0.0.1:29759/graphics/plot_zoom_png?width=515&height=769
WTC4_data = read.csv("data/in/WTC4_data.csv")
WTC4_data$DateTime_hr <- as.POSIXct(WTC4_data$DateTime_hr,format="%Y-%m-%d %T",tz="GMT")

# Plot A-Tleaf response without filtering
# Note that very few values of A < 0 occur at high temperatures
WTC4_data %>% #filter(HWtrt == "HW") %>% 
  ggplot(aes(x = Tleaf, y = A, color = HWtrt)) +
  geom_point() +
  theme_classic() +
  geom_hline(yintercept = 0)

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
                    Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633, Rd = 0)
pred1.c$P = get_Pleaves_Medlyn(pred1.c, P50 = P50, P88 = P88, 
                                Ps = control$Ps, kmax_25 = 0.5, constant_kmax = TRUE)

table(pred1.c$failed) # Some energy balance calculations failed
pred1.c$Tdiff <- with(pred1.c,Tleaf-Tair)

## Fit Sperry and Sicangco models ##############################################

# Generate predictions
pred2.c = make_pred(df = control,
                    Tcrit = 43.4, T50 = 49.6, LeafAbs = 0.5, P50 = 4.07, P88 = 5.50,
                    Wind = 8, Wleaf = 0.025,
                    Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
                    Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92,
                    kmax_25 = 0.5) # Really slow! Try to optimize

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
save(control, pred1.c, pred2.c, file = "data/out/control_runs_kmaxpt5_new.Rdata")
load("data/out/control_runs_kmaxpt7_fittedTcritT50_JPhydraulics.Rdata")

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
out.c = data.frame(datetime = rep(control$DateTime_hr, times = 6),
                   chamber = rep(control$chamber, times = 6),
                   Model = rep(c("observed", "Medlyn", "Sperry", 
                                 "Sperry + varkmax", "Sperry + CGnet", 
                                 "Sperry + CGnet + TC"),
                               each = nrow(control)),
                   E = c(control$E, pred1.c$ELEAF, pred2.c$E),
                   A = c(control$A, pred1.c$ALEAF, pred2.c$A),
                   gs = c(control$gs, pred1.c$gw, pred2.c$gs),
                   Dleaf = c(control$Dleaf, pred1.c$VPDleaf, pred2.c$Dleaf),
                   Tleaf = c(control$Tcan, pred1.c$Tleaf, pred2.c$Tleaf))

out.c$Model = factor(out.c$Model, levels = c("observed", "Medlyn", "Sperry", 
                                             "Sperry + varkmax", "Sperry + CGnet", 
                                             "Sperry + CGnet + TC"))
# Plot transpiration versus canopy temperature 
# E too high for Sperry and Sicangco models, especially as Tleaf increases
# Notably this is exaggerated for chamber 11, which has an increasing trend in SWP
EvT.c = plotGAM(gam_E.c, smooth.c = "Tcan") +
  geom_point(data = out.c, aes(x = Tleaf, y = E, color = Model), shape = 1, size = .5) +
  theme_classic() +
  scale_color_manual(
    values = palette,
    labels = c("Observations",
               "USO",
               "ProfitMax",
               expression("ProfitMax"[k[max](T)]),
               expression("ProfitMax"[net]),
               expression("ProfitMax"[TC])
    )) +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  ylab(expression("E"*" (mmol m"^-2*"s"^-1*")")) + 
  guides(color = guide_legend(override.aes = list(shape = 19, size = 2)),
         linetype = "none") +
  theme(plot.title = element_blank())
EvT.c

# Plot A versus canopy temperature 
AvT.c = plotGAM(gam_A.c, smooth.c = "Tcan") +
  geom_point(data = out.c, aes(x = Tleaf, y = A, color = Model), shape = 1, size = .5) +
  theme_classic() +
  scale_color_manual(
    values = palette,
    labels = c("Observations",
               "USO",
               "ProfitMax",
               expression("ProfitMax"[k[max](T)]),
               expression("ProfitMax"[net]),
               expression("ProfitMax"[TC])
    )) +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  ylab(expression("A (" * mu * "mol m"^-2*"s"^-1*")")) + 
  guides(color = guide_legend(override.aes = list(shape = 19, size = 2)),
         linetype = "none") #+
  #theme(plot.title = element_blank())
AvT.c

# Plot Tleaf predictions vs observations
Tleaf_pred_obs.plt = out.c %>% 
  pivot_wider(names_from = model, values_from = Tleaf, id_cols = c(chamber, datetime)) %>% 
  pivot_longer(cols = Medlyn:Sperry, names_to = "model", values_to = "Tleaf_pred") %>% 
  rename(Tleaf = observed) %>% 
  ggplot() +
  geom_point(aes(x = Tleaf, y = Tleaf_pred, color = model), shape = 1) +
  theme_classic() +
  scale_color_manual(
    values = c("Medlyn" = "#D81B60", 
               "Sicangco" = "#FFC107", "Sperry" = "#1E88E5")) +
  xlab(expression("observed T"[leaf]*" (\u00B0C)")) +
  ylab(expression("predicted T"[leaf]*" (\u00B0C)")) + 
  guides(color = guide_legend(override.aes = list(shape = 19, size = 2)),
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
                           PPFD = heatwave$PPFD, Wind = 4, Wleaf = 0.04, LeafAbs = 0.5))
Tleaves.hw[grep("Error", Tleaves.hw)] = NA
Tleaves.hw = as.numeric(Tleaves.hw)

# Tends to over-predict leaf temperature slightly
Tdiff.hw = Tleaves.hw - heatwave$Tleaf
summary(Tdiff.hw)
hist(Tdiff.hw)
plot(x = heatwave$Tleaf, y = Tdiff.hw, col = "deeppink")
points(x = heatwave$Tcan, y = Tleaves.hw - heatwave$Tcan)
abline(h = 0, col = "blue")

## Fit GAM #####################################################################
gam_E.hw = gam(E ~ s(Tcan), data = heatwave)
gam_A.hw = gam(A ~ s(Tcan), data = heatwave)


## Fit the Medlyn model ########################################################

Rd = Rd0 * exp(0.1012 * (heatwave$Tcan - TrefR) - 5e-04 * (heatwave$Tcan^2 - TrefR^2))
# Fit the heatwave treatment (i.e., not the heatwave trees)
pred1.hw <- PhotosynEB(Tair=heatwave$Tair,VPD=heatwave$VPD,
                       Wind=8,Wleaf=0.02,StomatalRatio=1,LeafAbs=0.5,
                       PPFD=heatwave$PPFD,g1=g1,g0=g0,
                       Vcmax=36.4,EaV=62307,EdVC=2e5,delsC=639,
                       Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Ca = 420,
                       Rd0 = 0.92)

get_Pleaves_Medlyn = Vectorize(get_Pleaf_Medlyn, "Ps")
pred1.hw$P = get_Pleaves_Medlyn(pred1.hw, P50 = P50, P88 = P88, 
                                  Ps = heatwave$Ps, kmax_25 = 0.5, constant_kmax = TRUE)

table(pred1.hw$failed) # Some energy balance calculations failed
pred1.hw$Tdiff <- with(pred1.hw,Tleaf-Tair)

## Fit Sperry and Sicangco models ##############################################

# Add kmax values to the data frame
#heatwave = left_join(heatwave, kmax_df, by = "chamber") #%>% rename(Tair = Tair, PPFD = PPFD)
#heatwave = heatwave %>% mutate(kmax = 0.7)

heatwave_subset = heatwave %>% 
  filter(Tleaf >= 33) %>% 
  mutate(DateTime_rd = floor_date(DateTime_hr, unit = "1 hour")) %>% 
  group_by(DateTime_rd, chamber, T_treatment, HWtrt, combotrt) %>% 
  summarise(across(PPFD:Ps, mean)) %>% 
  ungroup() %>% 
  rename(DateTime_hr = DateTime_rd)

#test_s = control[1:20,] # Subset of dataset for testing

# Generate predictions
pred2.hw = make_pred(df = heatwave, 
                     Tcrit = 43.4, T50 = 49.6, LeafAbs = 0.5, P50 = 4.07, P88 = 5.50,
                     Wind = 8, Wleaf = 0.025,
                     Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
                     Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92,
                     kmax_25 = 0.5) # Really slow! Try to optimize
#pred2.hw = make_pred(df = heatwave, models = "final", Tcrit = 46, T50 = 48) # Really slow! Try to optimize

plot(pred2.hw$Tleaf, pred2.hw$E)

plot(pred2.hw$P[pred2.hw$Model == "Sperry"])
plot(heatwave$Ps)

# Check if either model predicts negative gs values
table(pred2$gs[pred2$Model == "Sicangco"] >= 0) # No negative gs
table(pred2$gs[pred2$Model == "Sperry"] >= 0) # No negative gs


### Fit intermediate models, i.e. Sperry with either CGnet or TC ###############

# Generate predictions
pred3.hw = make_pred(df = heatwave, models = "intermediate", 
                     Tcrit_hw = 43.4, T50_hw = 49.6, LeafAbs = 0.5, P50 = 4.07, P88 = 5.50) # Really slow! Try to optimize
# Check if either model predicts negative gs values
table(pred3.hw$gs[pred3.hw$Model == "Sperry + CGnet"] >= 0) # 158/2341 have negative gs
table(pred3.hw$gs[pred3.hw$Model == "Sperry + TC"] >= 0) # All predict positive gs
summary(control$PPFD[which(pred3.hw$gs < 0)]) # Happens for low PPFD (< 52.5)


# Save heatwave observations and predictions
save(heatwave, pred1.hw, pred2.hw, file = "data/out/heatwave_runs_varkmaxpt5_new.Rdata")
load("data/out/heatwave_runs_varkmaxpt5_fittedTcritT50_JPhydraulics_Wind8_Wleafpt02_DKparams_CGnetorig_5models.Rdata")

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
out.hw = data.frame(datetime = rep(heatwave$DateTime_hr, times = 6),
                    chamber = rep(heatwave$chamber, times = 6),
                    Model = rep(c("observed", "Medlyn", "Sperry", 
                                  "Sperry + varkmax", "Sperry + CGnet", 
                                  "Sperry + CGnet + TC"), each = nrow(heatwave)),
                    E = c(heatwave$E, pred1.hw$ELEAF, pred2.hw$E),
                    A = c(heatwave$A, pred1.hw$ALEAF, pred2.hw$A),
                    gs = c(heatwave$gs, pred1.hw$gw, pred2.hw$gs),
                    Dleaf = c(heatwave$Dleaf, pred1.hw$VPDleaf, pred2.hw$Dleaf),
                    Tleaf = c(heatwave$Tcan, pred1.hw$Tleaf, pred2.hw$Tleaf),
                    Pleaf = c(rep(NA, nrow(heatwave)), rep(NA, nrow(heatwave)), pred2.hw$P),
                    Tair = c(heatwave$Tair, pred1.hw$Tair, pred2.hw$Tair),
                    VPD = rep(heatwave$VPD, times = 6),
                    PPFD = rep(heatwave$PPFD, times = 6),
                    Ps = rep(heatwave$Ps, times = 6)) 
out.hw$Model = factor(out.hw$Model, levels = c("observed", "Medlyn", "Sperry", "Sperry + varkmax", "Sperry + CGnet", 
                                               "Sperry + CGnet + TC" 
                                               ))

# Plot transpiration versus canopy temperature 
# E too high for Sperry and Sicangco models, but trend is good
# Notably this is exaggerated for chambers 7 and 9, which has an "increasing" trend in SWP
palette = c(Sperry = "#88CCEE", "Sperry + CGnet" = "#332288", "Sperry + CGnet + TC" = "#DDCC77",
            Medlyn = "#CC6677", observed = "black")
palette = c(Sperry = "#1E88E5", "Sperry + varkmax" = "#332288", "Sperry + CGnet" = "#e7d87d", "Sperry + CGnet + TC" = "orange", 
            Medlyn = "#D81B60", observed = "grey50")
EvT.hw = plotGAM(gam_E.hw, smooth.c = "Tcan") +
  geom_point(data = out.hw, aes(x = Tleaf, y = E, color = Model), shape = 1, size = .5) +
  theme_classic() +
  scale_color_manual(
    values = palette,
    labels = c("Observations",
               "USO",
               "ProfitMax",
               expression("ProfitMax"[k[max](T)]),
               expression("ProfitMax"[net]),
               expression("ProfitMax"[TC])
    )) +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  ylab(expression("E"*" (mmol m"^-2*"s"^-1*")")) + 
  guides(color = guide_legend(override.aes = list(shape = 19, size = 2)),
         linetype = "none") +
  theme(plot.title = element_blank())
EvT.hw

# Plot A versus canopy temperature 
AvT.hw = plotGAM(gam_A.hw, smooth.c = "Tcan") +
  geom_point(data = out.hw, aes(x = Tleaf, y = A, color = Model), shape = 1, size = .5) +
  theme_classic() +
  scale_color_manual(
    values = palette,
    labels = c("Observations",
               "USO",
               "ProfitMax",
               expression("ProfitMax"[k[max](T)]),
               expression("ProfitMax"[net]),
               expression("ProfitMax"[TC])
               )) +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  ylab(expression("A (" * mu * "mol m"^-2*"s"^-1*")")) + 
  guides(color = guide_legend(override.aes = list(shape = 19, size = 2)),
         linetype = "none") #+
  #theme(plot.title = element_blank())
AvT.hw
  # Combine all plots (i.e. recreate Drake et al. Fig 5 with models)
AEvT.plt = ggarrange(AvT.c + ylim(-2.5, 13) + 
                       labs(title = "Control",
                            subtitle = expression(bold("(a)"))) + theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)), 
          AvT.hw + ylim(-2.5, 13) +
          labs(title = "Heatwave",
               subtitle = expression(bold("(b)"))) + ggtitle("Heatwave") + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)), 
          EvT.c + ylim(-1, 4) + ggtitle("(c)") + theme(axis.title.x = element_blank(), plot.title = element_text(face = "bold", hjust = 0, vjust = 5)), 
          EvT.hw + ylim(-1,4) + ggtitle("(d)") + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(face = "bold", hjust = 0, vjust = 5)),
          nrow = 2, ncol = 2, common.legend = TRUE, legend = "right")
AEvT.plt = 
  annotate_figure(AEvT.plt,
                bottom = text_grob(expression("T"[canopy]*" (\u00B0C)"), hjust = 1))


AEvT.plt
ggsave(plot = AEvT.plt, filename = "figs/AEvT_plt_new.pdf", width = 8, height = 7)

# gs vs Tleaf
ggplot() +
geom_point(data = out.hw, aes(x = Tleaf, y = gs, color = Model), shape = 1, size = .5) +
  theme_classic() +
  scale_color_manual(
    values = palette,
    labels = c("Observations",
               "USO",
               "ProfitMax",
               expression("ProfitMax"[k[max](T)]),
               expression("ProfitMax"[net]),
               expression("ProfitMax"[TC])
    )) +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  #ylab(expression("E"*" (mmol m"^-2*"s"^-1*")")) + 
  guides(color = guide_legend(override.aes = list(shape = 19, size = 2)),
         linetype = "none") +
  theme(plot.title = element_blank())


pred2.hw %>% 
  ggplot(aes(x = Tleaf, y = -P, color = Model)) +
  geom_point() +
  theme_classic()

Tleaf_pred_obs.plt = 
  bind_rows(out.c, out.hw) %>% 
  pivot_wider(names_from = Model, values_from = Tleaf, id_cols = c(chamber, datetime)) %>% 
  pivot_longer(cols = Medlyn:"Sperry + CGnet + TC", names_to = "Model", values_to = "Tleaf_pred") %>% 
  rename(Tleaf = observed) %>% 
  ggplot() +
  geom_point(aes(x = Tleaf, y = Tleaf_pred, color = Model), shape = 1) +
  theme_classic() +
  scale_color_manual(values = palette[-6],
                     labels = c("USO",
                                "ProfitMax",
                                expression("ProfitMax"[k[max](T)]),
                                expression("ProfitMax"[net]),
                                expression("ProfitMax"[TC])
                     )) +
  xlab(expression("observed T"[canopy]*" (\u00B0C) from radiometer")) +
  ylab(expression("predicted T"[canopy]*" (\u00B0C)")) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)),
         linetype = "none") +
  theme(plot.title = element_blank(),axis.title = element_text(size = 14)) + 
  geom_abline(slope = 1) +
  guides(color = guide_legend(override.aes = list(shape = 19, size = 2)))
Tleaf_pred_obs.plt
ggsave(plot = Tleaf_pred_obs.plt, filename = "figs/Tleaf_pred_vs_obs_all.pdf",
       height = 7, width = 10)

out.hw %>% 
  pivot_wider(names_from = model, values_from = Tleaf, id_cols = c(chamber, datetime)) %>% 
  pivot_longer(cols = Medlyn:Sperry, names_to = "model", values_to = "Tleaf_pred") %>% 
  rename(Tleaf = observed) %>% 
  ggplot() +
  geom_point(aes(x = Tleaf, y = Tleaf_pred - Tleaf, color = model), shape = 1) +
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
  geom_hline(yintercept = 46, linetype = "dashed") +
  annotate("text", y=Tcrit, label="Tcrit", vjust = -0.5)

out.hw %>% 
  filter(Model != "observed") %>% 
  ggplot(aes(x = Ps, y = Pleaf, color = Model)) +
  geom_point(shape = 1)+
  scale_color_manual(values = palette[-6],
                     labels = c("USO",
                                "ProfitMax",
                                expression("ProfitMax"[k[max](T)]),
                                expression("ProfitMax"[net]),
                                expression("ProfitMax"[TC])
                     )) +
  theme_classic() +
  #xlab(expression("T"[leaf]*" (\u00B0C)")) +
  ylab(expression(Psi[leaf]*" (-MPa)")) + 
  guides(color = guide_legend(override.aes = list(shape = 19, size = 2)),
         linetype = "none")
ggsave("figs/PleafvT_kmaxpt7_fittedTcritT50_JPhydraulics.pdf")

# Calculate TSM and HSM ########################################################
# TSM
TSM_summary = out.hw %>%
  filter(model != "observed") %>% 
  rename(Model = model) %>% 
  bind_rows(pred3.hw) %>% 
  group_by(Model) %>% 
  summarise(Tleaf_max = max(Tleaf)) %>% 
  mutate(TSM_Tcrit = 43.4 - Tleaf_max,
         TSM_T50 = 49.9 - Tleaf_max) %>% 
  arrange(factor(Model, levels = c("Medlyn", "Sperry", "Sperry + TC", 
                                   "Sperry + CGnet", "Sicangco"))) %>% 
  as.data.frame()
write.csv(TSM_summary, "data/out/TSM_summary.csv", row.names = FALSE)

# HSM
HSM_summary = out.hw %>%
  filter(model != "observed") %>% 
  rename(Model = model, P = Pleaf) %>% 
  bind_rows(pred3.hw) %>% 
  group_by(Model) %>% 
  summarise(Pleaf_min = max(P)) %>% 
  mutate(HSM_P50 = 4.07 - Pleaf_min,
         HSM_P88 = 5.50 - Pleaf_min) %>% 
  arrange(factor(Model, levels = c("Medlyn", "Sperry", "Sperry + TC", 
                                   "Sperry + CGnet", "Sicangco")))
write.csv(HSM_summary, "data/out/HSM_summary.csv", row.names = FALSE)


# Model evaluation #############################################################
eval_model = function(fit_group = c("low", "medium", "high", "none"),
                      treatment = c("control", "heatwave"),
                      model = c("Medlyn", "Sperry", "Sperry + varkmax", 
                                "Sperry + CGnet", "Sperry + CGnet + TC"),
                      variable = c("E", "A", "gs", "Dleaf", "Tleaf")
)
{
  out.c_wide = out.c %>% 
    pivot_wider(names_from = Model, values_from = c(E, A, gs, Dleaf, Tleaf))
  out.hw_wide = out.hw %>% select(datetime:Tleaf) %>% 
    pivot_wider(names_from = Model, values_from = c(E, A, gs, Dleaf, Tleaf))
  
  df = if (treatment == "control") {out.c_wide} else {out.hw_wide}
  
  df = if (fit_group == "low") {
    df %>% filter(Tleaf_observed < 25)
  } else if (fit_group == "medium" & treatment == "control") {
    df %>% filter(Tleaf_observed >= 25)
  } else if (fit_group == "medium" & treatment == "heatwave") {
    df %>% filter(Tleaf_observed >= 25 & Tleaf_observed < 35)
  } else if (fit_group == "high") {
    df %>% filter(Tleaf_observed >= 35)
  } else if (fit_group == "none") {
    df
  }
  
  # Create vector with observations
  obs = df[[paste0(variable, "_observed")]]
  
  # Create vector with predictions
  pred = df[[paste0(variable, "_", model)]]
  
  mae = mae_vec(obs, pred)
  r2 = rsq_vec(obs, pred)
  
  eval_metrics = c("MAE" = mae, "R2" = r2)
  
  return(eval_metrics)
}


eval_model(fit_group = "none",treatment = "heatwave", model = "Sperry", variable = "E")

## Data binned by temperature --------------------------------------------------
model_eval_df = data.frame(treatment = c(rep("control", 30), rep("heatwave", 45)),
                           fit_group = c(rep(c("low", "medium"), each = 15, times = 2), rep("high", 15)),
                           model = rep(c("Medlyn", "Sicangco", "Sperry"), times = 25), 
                           variable = rep(c("E", "A", "gs", "Dleaf", "Tleaf"), each = 3, times = 5)
)
model_evals = sapply(1:nrow(model_eval_df), 
                     function(i) eval_model(model_eval_df$fit_group[i],
                                            model_eval_df$treatment[i],
                                            model_eval_df$model[i], 
                                            model_eval_df$variable[i]))
model_evals = as_data_frame(t(model_evals))
model_eval_df = cbind(model_eval_df, model_evals)

model_eval_df$fit_group = factor(model_eval_df$fit_group, 
                                 levels = c("low", "medium", "high"))

model_eval_df %>% 
  filter(treatment == "heatwave") %>% 
  ggplot(aes(x = model, y = MAE)) +
  geom_col(fill = "darkorange") +
  facet_wrap(fit_group ~ variable, scales = "free", ncol = 5) +
  theme_classic() +
  ylab("MAE")

## Unbinned data----------------------------------------------------------------
model_eval_df_unbinned = data.frame(treatment = rep(c("control", "heatwave"), each = 25),
                                    fit_group = rep("none", 50),
                                    model = rep(c("Medlyn", "Sperry", "Sperry + varkmax", 
                                                  "Sperry + CGnet", "Sperry + CGnet + TC"), 
                                                times = 10), 
                                    variable = rep(c("E", "A", "gs", "Dleaf", "Tleaf"), each = 5, times = 2)
)
model_evals_unbinned = sapply(1:nrow(model_eval_df_unbinned), 
                     function(i) eval_model(model_eval_df_unbinned$fit_group[i],
                                            model_eval_df_unbinned$treatment[i],
                                            model_eval_df_unbinned$model[i], 
                                            model_eval_df_unbinned$variable[i]))
model_evals_unbinned = as_data_frame(t(model_evals_unbinned))
model_eval_df_unbinned = cbind(model_eval_df_unbinned, model_evals_unbinned)

model_eval_df_unbinned$model = gsub("Sperry", "ProfitMax",
                                    gsub("Medlyn", "USO",
                                         gsub("Sicangco", "Si", model_eval_df_unbinned$model)))

MAE_plt.hw = model_eval_df_unbinned %>% 
  filter(treatment == "heatwave") %>% 
  ggplot(aes(x = model, y = MAE)) +
  geom_col(fill = "darkorange") +
  facet_wrap(vars(variable), scales = "free", ncol = 5) +
  theme_classic() +
  ylab("MAE")

R2_plt.hw = model_eval_df_unbinned %>% 
  filter(treatment == "heatwave") %>% 
  ggplot(aes(x = model, y = R2)) +
  geom_col(fill = "darkorange") +
  facet_wrap(vars(variable), scales = "free", ncol = 5) +
  theme_classic() +
  ylab("R2")

MAE_plt.c = model_eval_df_unbinned %>% 
  filter(treatment == "control") %>% 
  ggplot(aes(x = model, y = MAE)) +
  geom_col(fill = "lightblue") +
  facet_wrap(vars(variable), scales = "free", ncol = 5) +
  theme_classic() +
  ylab("MAE")

R2_plt.c = model_eval_df_unbinned %>% 
  filter(treatment == "control") %>% 
  ggplot(aes(x = model, y = R2)) +
  geom_col(fill = "lightblue") +
  facet_wrap(vars(variable), scales = "free", ncol = 5) +
  theme_classic() +
  ylab("R2")

ggarrange(MAE_plt.c + ggtitle("Control"), 
          MAE_plt.hw + ggtitle("Heatwave"), 
          R2_plt.c, 
          R2_plt.hw)

out.hw %>% 
ggplot() +
  geom_point(aes(x = Tleaf, y = gs, color = model), alpha = .5, size = .5) +
  theme_classic() +
  scale_color_manual(
    values = c("observed" = "black", "Medlyn" = "#D81B60", 
               "Sicangco" = "#FFC107", "Sperry" = "#1E88E5")) +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  #ylab(expression("A (" * mu * "mol m"^-2*"s"^-1*")")) + 
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)),
         linetype = "none") +
  theme(plot.title = element_blank())


# Testing for A behaviour at high temps ########################################
test = out.hw %>% filter(Tleaf > 41, model == "Sperry")
gs_Med = out.hw$gs[out.hw$Tleaf > 41 & out.hw$model == "Medlyn"]
test_A_Sp = calc_A(Tleaf = test$Tleaf, g_w = test$gs, VPD = test$VPD, net = TRUE, netOrig = TRUE,
           PPFD = test$PPFD, Wind = 8, 
           Wleaf = 0.01, LeafAbs = 0.5,
           Vcmax=34,EaV=51780,EdVC=2e5,delsC=640,
           Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633)
test_A_M = calc_A(Tleaf = test$Tleaf, g_w = 0.003, VPD = test$VPD, net = TRUE, netOrig = TRUE,
                   PPFD = test$PPFD, Wind = 8, 
                   Wleaf = 0.01, LeafAbs = 0.5,
                   Vcmax=34,EaV=51780,EdVC=2e5,delsC=640,
                   Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633)
plot(test$A)
points(test_A_Sp, col = "blue")
points(test_A_M, col = "deeppink")
points(out.hw$A[out.hw$Tleaf > 41 & out.hw$model == "Medlyn"], col = "deeppink3")
