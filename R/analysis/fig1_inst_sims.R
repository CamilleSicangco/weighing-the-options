# Fig 1: Overview of instantaneous simulation models
# by Camille Sicangco
# Created 23 May 2025

# Set variables/params ---------------------------------------------------------
# Environment
Tair = 46
Ps = 0.5
VPD = 1.5
PPFD = 1500

# Thermal damage
Tcrit = 43.4
T50 = 49.6

# Hydraulics
Weibull = fit_Weibull(P50 = 4.47, P88 = 5.50)
b = Weibull[1,1]
c = Weibull[1,2]
kmax_25 = 0.5
Pcrit = calc_Pcrit(b, c)
P = Ps_to_Pcrit(Ps, Pcrit)

# CGnet normalized with or without including the point where E = 0
CG_net = C_gain_corr(P, b, c, Amax = NULL, kmax_25, Tair = 48, VPD, PPFD, 
                     Wind = 8, Wleaf = 0.02, LeafAbs = 0.5,
                     Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
                     Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92,
                     constant_kmax = FALSE, net = TRUE, netOrig = TRUE)
CG_net_uncorr = C_gain(P, b, c, Amax = NULL, kmax_25, Tair = 48, VPD, PPFD, 
                       Wind = 8, Wleaf = 0.02, LeafAbs = 0.5, 
                       Vcmax=34,EaV=51780,EdVC=2e5,delsC=639, 
                       Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=635, Rd0 = 0.92,
                       constant_kmax = TRUE, net = TRUE, netOrig = TRUE)
plot(CG_net)
plot(CG_net_uncorr)

# Plots ------------------------------------------------------------------------

# Tair = 30
cost_gain30 = calc_costgain_netorig(P, b, c, kmax_25 = kmax_25, 
                                  Tair = 30, PPFD = PPFD, VPD = VPD,
                                  Tcrit = Tcrit, T50 = T50, #constant_kmax = TRUE,
                                  Wind = 8, Wleaf = 0.02, LeafAbs = 0.5,
                                  Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
                                  Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92)
composite_plot(cost_gain30)
ggsave(filename = "figs/composite_plot_30deg.pdf", width = 12.25, height = 6.5)

# Tair = 40
cost_gain41.3 = calc_costgain_netorig(P, b, c, kmax_25 = kmax_25, 
                                  Tair = 41.3, PPFD = PPFD, VPD = VPD,
                                  Tcrit = Tcrit, T50 = T50, #constant_kmax = TRUE,
                                  Wind = 8, Wleaf = 0.02, LeafAbs = 0.5,
                                  Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
                                  Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92)
composite_plot(cost_gain41.3)
ggsave(filename = "figs/composite_plot_41.3deg.pdf", width = 12.25, height = 6.5)

# Tair = 48
cost_gain48 = calc_costgain_netorig(P, b, c, kmax_25 = kmax_25, 
                                    Tair = 45, PPFD = PPFD, VPD = VPD,
                                    Tcrit = Tcrit, T50 = T50, #constant_kmax = FALSE,
                                    Wind = 8, Wleaf = 0.02, LeafAbs = 0.5,
                                    Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
                                    Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92
                                    )
composite_plot(cost_gain48)
ggsave(filename = "figs/composite_plot_48deg.pdf", width = 12.25, height = 6.5)

make_pred_fn(Tair = 60, Ps = 0.5, VPD = 1.5, PPFD = 1500, model = "Sperry", 
             Wind = 8, Wleaf = 0.02, LeafAbs = 0.5,
             Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
             Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92,
             Tcrit = 43.4, T50 = 44.4,
             kmax_25 = 2.5, constant_kmax = FALSE)
# Debugging --------------------------------------------------------------------
# Solve optimization
i = if (model == "Sperry") {
  which.max(cost_gain$CG_gross - cost_gain$HC)
} else if (model == "Sicangco") {
  which.max(cost_gain$CG_net - (cost_gain$HC + cost_gain$TC))
} else if (model == "Sperry + CGnet") {
  which.max(cost_gain$CG_net - cost_gain$HC)
} else if (model == "Sperry + TC") {
  which.max(cost_gain$CG_gross - (cost_gain$HC + cost_gain$TC))
} else {
  stop()
}


p1 = (cost_gain$CG_net - (cost_gain$HC + cost_gain$TC))
p2 = (CG_net_uncorr - (cost_gain$HC + cost_gain$TC))

plot(p1)
plot(p2)
P = P[i]
E = E_vec[i]

E = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax = TRUE)
Tleaf = calc_Tleaf(Tair = Tair, E = E, VPD = VPD, PPFD = PPFD, Wind = Wind, 
                   Wleaf = Wleaf, LeafAbs = LeafAbs)
Dleaf = VPDairToLeaf(Tleaf = Tleaf, Tair = Tair, VPD = VPD)
gs = calc_gw(E = E, Tleaf = Tleaf, Tair = Tair, VPD = VPD, 
             PPFD = PPFD, Wind = Wind, Wleaf = Wleaf)
A = calc_A(Tair = Tair, E = E, VPD = VPD, net = TRUE, netOrig = TRUE,
           PPFD = PPFD, Wind = 8, 
           Wleaf = 0.02, LeafAbs = 0.5,
           Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
           Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92)
plot(P,A)
plot(cost_gain$cost_gain[cost_gain$ID == "CG_net"])

Photosyn_custom(Tleaf = Tleaf[1], GS = gs[1], VPD = VPD, 
                PPFD = PPFD, 
                Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
                Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92)
