# Test of E calculation
# Created 12 Aug 2025
# by Camille Sicangco

# Environment
Tair_vec = seq(20,60, by = 1)
PPFD = 1500
VPD = 1.5

# Ps = -0.5 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 0.5)
#out_Ps0.5_constVPD = get_predictions(Tair_sim.df)

# Hydraulics
P50 = 4.07
P88 = 5.50

Weibull = fit_Weibull(P50, P88)
b = Weibull[1,1]
c = Weibull[1,2]
Pcrit = calc_Pcrit(b, c)
P = Ps_to_Pcrit(Ps, Pcrit, pts = 600)


# Calculate costs and gains
cost_gain = calc_costgain(P, b, c, kmax_25 = kmax_25, Tair = Tair, VPD = VPD, 
                          PPFD = PPFD, Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs, 
                          Tcrit = Tcrit, T50 = T50, 
                          Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
                          Jmax = Jmax,EaJ=EaJ,EdVJ=EdVC,delsJ=delsJ,
                          Rd0 = Rd0, constr_Ci = constr_Ci)

costgain_df = bind_rows(lapply(1:nrow(Tair_sim.df), 
                               function(i) calc_costgain(
                                 P, b, c,
                                 Tair = Tair_sim.df$Tair[i], Ps = Tair_sim.df$Ps[i], 
                                 VPD = Tair_sim.df$VPD[i], PPFD = Tair_sim.df$PPFD[i])))
costgain_df = calc_costgain(P, b, c,
                            Tair = Tair_sim.df$Tair, Ps = Tair_sim.df$Ps, 
                            VPD = Tair_sim.df$VPD, PPFD = Tair_sim.df$PPFD,
                            Tcrit = 46.5, T50 = 50.4, #P50 = 4.07, P88 = 5.50,
                            Wind = 5, Wleaf = 0.025, LeafAbs = 0.5,
                            Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
                            Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92,
                            constr_Ci = FALSE)
i = which.max(costgain_df$CG_gross_constkmax - costgain_df$HC_constkmax)

Weibull = fit_Weibull(P50, P88)
b = Weibull[1,1]
c = Weibull[1,2]

E_vec = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax = TRUE)


P = P[i]
E = E_vec[i]
Tleaf = calc_Tleaf(Tair = Tair, E = E, VPD = VPD, PPFD = PPFD, Wind = Wind, 
                   Wleaf = Wleaf, LeafAbs = LeafAbs)
Dleaf = VPDairToLeaf(Tleaf = Tleaf, Tair = Tair, VPD = VPD)
gs = calc_gw(E = E, Tleaf = Tleaf, Tair = Tair, VPD = VPD, 
             PPFD = PPFD, Wind = Wind, Wleaf = Wleaf)

new_JT = if (model %in% c("Sperry + TC", "Sperry + CGnet + TC")) {
  TRUE
} else {
  FALSE
}

A = calc_A(Tleaf = Tleaf, g_w = gs, VPD = VPD, net = TRUE, netOrig = TRUE,
           PPFD = PPFD, Wind = Wind, 
           Wleaf = Wleaf, LeafAbs = LeafAbs,
           Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
           Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ, new_JT = new_JT,
           ...)

out = c(model, Tair, Ps, P, E, Tleaf, Dleaf, gs, A)