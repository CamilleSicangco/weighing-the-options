# Fig 1: Overview of instantaneous simulation models
# by Camille Sicangco
# Created 23 May 2025

# Set variables/params ---------------------------------------------------------

# Environment
Tair = 46
Ps = 0.5
VPD = RHtoVPD(RH = 60, Tair)
PPFD = 1500

# Thermal damage
Tcrit = 46.5
T50 = 50.4

# Hydraulics
Weibull = fit_Weibull(P50 = 4.07, P88 = 5.50)
b = Weibull[1,1]
c = Weibull[1,2]
kmax_25 = 0.5
Pcrit = calc_Pcrit(b, c)
P = Ps_to_Pcrit(Ps, Pcrit, pts = 600)

# Plots ------------------------------------------------------------------------

# Tair = 40
Tair = 40
VPD = RHtoVPD(RH = 60, Tair)
cost_gain40 = calc_costgain(P, b, c, kmax_25 = kmax_25, 
                            Tair = Tair, PPFD = PPFD, VPD = VPD,
                            Tcrit = Tcrit, T50 = T50,
                            Wind = 5, Wleaf = 0.025, LeafAbs = 0.5,
                            Vcmax=100.52,EaV=62307,EdVC=2e5,delsC=639,
                            Jmax = 165.53,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92)


Fig2 = composite_plot(cost_gain40)
ggsave(filename = "figs/Fig2_40deg_inst_sim.tiff", Fig2, 
       width = 12.25, height = 6.5, bg = "white")

# Tair = 48
Tair = 48
VPD = RHtoVPD(RH = 60, Tair)
cost_gain48 = calc_costgain(P, b, c, kmax_25 = kmax_25, 
                            Tair = Tair, PPFD = PPFD, VPD = VPD,
                            Tcrit = Tcrit, T50 = T50,
                            Wind = 5, Wleaf = 0.02, LeafAbs = 0.5,
                            Vcmax=100.52,EaV=62307,EdVC=2e5,delsC=639,
                            Jmax = 165.53,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92
)
Fig3 = composite_plot(cost_gain48)
ggsave(filename = "figs/Fig3_48deg_inst_sim.tiff", Fig3, 
       width = 12.25, height = 6.5, bg = "white")

# Tair = 30
Tair = 30
VPD = RHtoVPD(RH = 60, Tair)
cost_gain30 = calc_costgain(P, b, c, kmax_25 = kmax_25, 
                            Tair = Tair, PPFD = PPFD, VPD = VPD,
                            Tcrit = Tcrit, T50 = T50,
                            Wind = 5, Wleaf = 0.02, LeafAbs = 0.5,
                            Vcmax=100.52,EaV=62307,EdVC=2e5,delsC=639,
                            Jmax = 165.53,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92
)
FigS4 = composite_plot(cost_gain30)
ggsave(filename = "figs/FigS4_30deg_inst_sim.tiff", FigS4, 
       width = 12.25, height = 6.5, bg = "white")
