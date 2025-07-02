# ProfitMax comparison: Agross vs Anet with constrained Ci
# Created 23 June 2025
# by Camille Sicangco

# Environment
Tair_vec = seq(20,60, by = 1)

# Params
Tcrit = 43.4
T50 = 49.6
P50 = 4.07
P88 = 5.50
Wind = 8
Wleaf = 0.025
LeafAbs = 0.5
Vcmax=34
EaV=62307
EdVC=2e5
delsC=639
Jmax = 60
EaJ=33115
EdVJ=2e5
delsJ=635
Rd0 = 0.92
kmax_25 = 0.5 
net = TRUE
netOrig = TRUE

# Ps = -0.5 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = 1500, VPD = 1.5, Ps = 0.5)

# Generate predictions with both models
preds.ProfitMax = make_pred_fn2(model = "Sperry", 
                       df = Tair_sim.df,
                       Tcrit = Tcrit, T50 = T50,
                       P50 = P50, P88 = P88,Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
                       Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
                       Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ, Rd0 = Rd0,
                       kmax_25 = kmax_25)
preds.ProfitMax_constrCi_varkmax = make_pred_fn2(model = "Sperry + CGnet_constrCi", 
                       df = Tair_sim.df,
                       Tcrit = Tcrit, T50 = T50,
                       P50 = P50, P88 = P88,Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
                       Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
                       Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ, Rd0 = Rd0,
                       kmax_25 = kmax_25)
preds.ProfitMax_net = make_pred_fn2(model = "Sperry + CGnet", 
                                         df = Tair_sim.df,
                                         Tcrit = Tcrit, T50 = T50,
                                         P50 = P50, P88 = P88,Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
                                         Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
                                         Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ, Rd0 = Rd0,
                                         kmax_25 = kmax_25)
bind_rows(preds.ProfitMax, preds.ProfitMax_constrCi_varkmax, preds.ProfitMax_net) %>% 
 # filter(Tair > 40 & Tair < 45) %>% 
  ggplot(aes(x = Tair, y = gs, color = Model)) +
  geom_point() +
  theme_classic()

HC_constkmax = hydraulic_cost(P, b, c, kmax_25, Tair = 44, constant_kmax = TRUE)
CG_gross = C_gain_corr(P, b, c, kmax_25 = kmax_25,
                       Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
                       Tair = 44, PPFD = PPFD, 
                       VPD = VPD, Tcrit = Tcrit, T50 = T50, 
                       Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
                       Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
                       Rd0 = Rd0,
                       constant_kmax = TRUE, net = FALSE, new_JT = FALSE)
CG_constrCi = C_gain_alt(
  P, b, c, Amax = NULL, kmax_25 = kmax_25,
  Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
  Tair = 44, PPFD = PPFD, 
  VPD = VPD, Tcrit = Tcrit, T50 = T50, 
  Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
  Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
  constant_kmax = TRUE, new_JT = FALSE)

range(CG_constrCi)
plot(CG_constrCi)
points(CG_gross)