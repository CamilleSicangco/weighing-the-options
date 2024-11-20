# Model testing over a temperature range
# by Camille Sicangco

# Function to make predictions with Sperry and Sicangco models
make_pred.1 = function(Tair, model) {
  Wind = 8
  Wleaf = 0.01
  LeafAbs = 0.86
  Vcmax=34
  EaV=51780
  EdVC=2e5
  delsC=640
  Jmax = 60
  EaJ=21640
  EdVJ=2e5
  delsJ=633
  Ps = 2
  VPD = 1.5
  PPFD = 700
  
  # Thermal damage
  Tcrit = 48.574
  T50 = 50.17
  #TSF = 4
  #Tcrit = Tcrit - TSF
  #T50 = T50 - TSF
  
  # Hydraulics
  P50 = 4.31
  P88 = 6.64
  Weibull = fit_Weibull(P50, P88)
  b = Weibull[1,1]
  c = Weibull[1,2]
  kmax_25 = 1.5
  Pcrit = calc_Pcrit(b, c)
  P = Ps_to_Pcrit(Ps, Pcrit)
  E_vec = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax = TRUE)
  
  
  # Calculate Amax values
  #Tair_range = seq(30,50, by = 0.5)
  #Amax_net = Amax_overT(P, b, c, kmax_25, Tair_range = Tair_range, constant_kmax = TRUE, net = TRUE, netOrig = FALSE)[1,2]
  #Amax_gross = Amax_overT(P, b, c, kmax_25, Tair_range = Tair_range, constant_kmax = TRUE)[1,2]
  
  # Compute costs and gains
  cost_gain = calc_costgain(P, b, c, kmax_25 = kmax_25, #Amax_net = Amax_net, Amax_gross = Amax_gross, 
                            Tair = Tair, PPFD = PPFD, 
                            VPD = VPD, Tcrit = Tcrit, T50 = T50,
                            #Wind = Wind, Wleaf = Wleaf, 
                            #LeafAbs = LeafAbs,
                            #Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
                            #Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
                            constant_kmax = TRUE
  )
  #ggplot(cost_gain, aes(x = P, y = cost_gain, color = ID)) +geom_line() +theme_classic()
  # Compute marginal costs and gains
  marg_df = marginal_gaincost(cost_gain)
  
  # Solve optimization
  i = if (model == "Sperry") {
    which.min(abs(marg_df$CG_gross - marg_df$HC))
  } else if (model == "Sicangco") {
    which.min(abs(marg_df$CG_net - marg_df$CC))
  } else {
    stop()
  }
  
  P = marg_df$P[i]
  E = E_vec[i]
  Tleaf = calc_Tleaf(Tair = Tair, E = E, VPD = VPD#, 
                     #Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs
                     )
  Dleaf = VPDairToLeaf(Tleaf = Tleaf, Tair = Tair, VPD = VPD)
  gs = calc_gw(E = E, D_leaf = Dleaf)
  A = if (model == "Sperry") {
    calc_A(Tair = Tair, E = E, VPD = VPD, net = FALSE,
           #Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
           #Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
           #Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ
           )
  } else if (model == "Sicangco") {
    calc_A(Tair = Tair, E = E, VPD = VPD, net = TRUE#, 
           #Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
           #Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
           #Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ
           )
  } else {
    stop()
  }
  
  out = c(model, Tair, P, E, Tleaf, Dleaf, gs, A)
  return(out)
}

# Environment
Tair_vec = seq(30,55, by = 1)
PPFD = 700
VPD = 1.5
kmax = 1.5
Ps = 2
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, kmax = kmax,
                         Ps = Ps)
out = make_pred(Tair_sim.df)

# Make predictions (Sperry and Sicangco)
#make_pred(Tair = 50, "Sperry") # test for one Tair
models = c("Sperry", "Sicangco")
sims_vec = data.frame(Tair = rep(Tair_vec, each = 2), 
                      model = rep(models, times = length(Tair_vec)))
out.1 = mapply(make_pred.1, sims_vec$Tair, sims_vec$model) # multiple Tair's
out.1 = data.frame(t(out.1))
names(out.1) = c("Model", "Tair", "P", "E", "Tleaf", "Dleaf", "gs", "A")
out.1 = out.1 %>% mutate(across(Tair:A, as.numeric))

# Make Medlyn predictions
Medlyn_preds = plantecophys::PhotosynEB(Tair=Tair_vec,
                          VPD=VPD,
                          Wind=8,Wleaf=0.01,StomatalRatio=1, LeafAbs=0.86,
                          PPFD=PPFD,g1 = 2.9,g0=0.003,
                          Vcmax=34,EaV=51780,EdVC=2e5,delsC=640,
                          Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633)
Medlyn_preds_sim = Medlyn_preds %>%
  rename(E = ELEAF, Dleaf = VPDleaf, gs = GS, A = ALEAF) %>%
  mutate(Model = "Medlyn", .before = 1) %>%
  select(Model, Tair, E, Tleaf, Dleaf, gs, A)

# Combine all predictions
out_all = bind_rows(out.1, Medlyn_preds_sim)

# Plot predictions #############################################################

# gs vs Tair
out_all %>% ggplot(aes(x = Tair, y = gs, color = Model)) +
  geom_point(size = 3) + geom_line(size = 1) +
  theme_classic() +
  ylab(expression("g"[s]*" (mol m"^-2*"s"^-1*")")) +
  xlab(expression("T"[air]*" (\u00B0C)")) +
  scale_color_manual(values = c("Sicangco" = "#FFC107", "Sperry" = "#1E88E5", Medlyn = "#D81B60"))+
  theme(text = element_text(size = 20),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

# Tleaf vs Tair
r = 2/(T50 - Tcrit)
T95 = log(1/0.95-1)/-r + T50

out_all %>% ggplot(aes(x = Tair, y = Tleaf, color = Model)) +
  geom_point(size = 3) + geom_line(size = 1) +
  theme_classic() +
  ylab(expression("T"[leaf]*" (\u00B0C)")) +
  xlab(expression("T"[air]*" (\u00B0C)")) +
  scale_color_manual(values = c("Sicangco" = "#FFC107", "Sperry" = "#1E88E5", Medlyn = "#D81B60"))+
  theme(text = element_text(size = 20),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept = Tcrit, linetype = "dashed") +
  annotate("text", x=40, y=Tcrit, label="Tcrit", vjust = -0.5) +
  geom_hline(yintercept = T50, linetype = "dotted", col = "darkorange2") +
  annotate("text", x=40, y=T50, label="T50", vjust = -0.5, col = "darkorange2") +
  geom_hline(yintercept = T95, linetype = "dotdash", col = "red") +
  annotate("text", x=40, y=T95, label="T95", vjust = -0.5, col = "red") +
  geom_abline(slope = 1, intercept = 0, col = "grey") +
  annotate("text", x = 40, y = 40, label = "1:1 line", col = "grey", vjust = 3)



# Pleaf vs Tair
out %>% ggplot(aes(x = Tair, y = P, color = Model)) +
  geom_point(size = 3) + geom_line(size = 1) +
  theme_classic() +
  ylab(expression("\u03A8"[leaf]*" (-MPa)")) +
  xlab(expression("T"[air]*" (\u00B0C)")) +
  scale_color_manual(values = c("Sicangco" = "#FFC107", "Sperry" = "#1E88E5"))+
  theme(text = element_text(size = 20),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  geom_hline(yintercept = P50, linetype = "dashed") +
  annotate("text", x=40, y=P50, label="P50", vjust = -0.5) +
  geom_hline(yintercept = P88, linetype = "dotted", col = "darkorange2") +
  annotate("text", x=40, y=P88, label="P88", vjust = -0.5, col = "darkorange2") +
  geom_hline(yintercept = Pcrit, linetype = "dotdash", col = "red") +
  annotate("text", x = 40, y = Pcrit, label = "Pcrit", vjust = -0.5, col = "red")

# All physio vars
out_l = out_all %>% pivot_longer(
  cols = P:A,
  names_to = "var",
  values_to = "pred"
)

out_l %>% 
  ggplot(aes(x = Tair, y = pred, color = Model)) + 
  geom_point() + geom_line() +
  facet_wrap(vars(var), scales = "free") + 
  theme_classic() +
  scale_color_manual(values = c("Sicangco" = "#FFC107", "Sperry" = "#1E88E5", Medlyn = "#D81B60")) #+
  #ggtitle(" - Default photosyn and leaf params + constant kmax")
