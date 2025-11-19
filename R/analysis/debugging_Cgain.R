# Debugging of C gain function over a temperature range
# by Camille Sicangco
# Created 19 Nov 2025

# INSTANTANEOUS

# Environment
Tair = 46
Ps = 0.5
VPD = RHtoVPD(RH = 60, Tair)
PPFD = 1500

# Thermal damage
Tcrit = 46.5
T50 = 50.4

# Hydraulics
Weibull = fit_Weibull(P50 = 4.47, P88 = 5.50)
b = Weibull[1,1]
c = Weibull[1,2]
kmax_25 = 0.5
Pcrit = calc_Pcrit(b, c)
P = Ps_to_Pcrit(Ps, Pcrit, pts = 600)


E = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax = TRUE)
Tleaf = calc_Tleaf(Tair = Tair, VPD = VPD, PPFD = PPFD, 
                   E = E)
g_w = calc_gw(E, Tleaf, Tair = Tair, VPD = VPD, PPFD = PPFD)
Dleaf = plantecophys::VPDairToLeaf(Tleaf = Tleaf, Tair = Tair, 
                                   VPD = VPD)
A_FvCB = calc_A(Tair, VPD, PPFD, E = E, net = TRUE, netOrig = TRUE)

A_FvCB = as.numeric(Photosyn_out[2, ])
Ci_FvCB = as.numeric(Photosyn_out[1, ])

Ci_Ficks = 420.0 - A_FvCB * 1E-4 / (g_w/1.57) * 100^2
A_Ficks = g_w/1.57*(420 - Ci_FvCB)

head(A_FvCB)
head(A_Ficks)

head(Ci_FvCB)
head(Ci_Ficks)

head(A_FvCB)
head(Ci_Ficks)

Ci_FvCB = as.numeric(mapply(plantecophys::Photosyn, VPD = Dleaf,
                            Ca = 420,
                      PPFD = PPFD, Tleaf = Tleaf,
                      GS = g_w, Jmax = 100, Vcmax = 50, g1 = 2.9, g0 = 0.003)[1, ])
A_Ficks = g_w/1.57*(420 - Ci_FvCB)
head(Ci_FvCB)
head(A_Ficks)

# RUN SIMULATIONS ##############################################################

# Environment
Tair_vec = seq(20, 60, by = 5)
PPFD = 1500
VPD = 1.5

### Constant VPD -------------------------

# Ps = -0.5 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 0.5)
out_Ps0.5_constVPD = get_predictions(Tair_sim.df)
out_Ps0.5_constVPD_constrCi = get_predictions(Tair_sim.df, constr_Ci = TRUE)


# Environment
Tair = 46
Ps = 0.5
VPD = RHtoVPD(RH = 60, Tair)
PPFD = 1500

# Thermal damage
Tcrit = 46.5
T50 = 50.4

# Hydraulics
Weibull = fit_Weibull(P50 = 4.47, P88 = 5.50)
b = Weibull[1,1]
c = Weibull[1,2]
kmax_25 = 0.5
Pcrit = calc_Pcrit(b, c)
P = Ps_to_Pcrit(Ps, Pcrit, pts = 500)
C_gain_alt_new(P)


Ci <- 400.0 - out_Ps0.5_constVPD$A * 1E-4 / (1.57 * out_Ps0.5_constVPD$gs) * 100^2
out_Ps0.5_constVPD$Ci_Fick = Ci
print(out_Ps0.5_constVPD)
write.csv(out_Ps0.5_constVPD, "Ps0.5_constVPD.csv")

# Ps = -2 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 2)
out_Ps2_constVPD = get_predictions(Tair_sim.df)
Ci <- 400.0 - out_Ps2_constVPD$A * 1E-4 / (1.57 * out_Ps2_constVPD$gs) * 100^2
out_Ps2_constVPD$Ci_Fick = Ci
#print(out_Ps2_constVPD)
write.csv(out_Ps2_constVPD, "Ps2_constVPD.csv")

# Ps = -4 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 4)
out_Ps4_constVPD = get_predictions(Tair_sim.df)
Ci <- 400.0 - out_Ps4_constVPD$A * 1E-4 / (1.57 * out_Ps4_constVPD$gs) * 100^2
out_Ps4_constVPD$Ci_Fick = Ci
#print(out_Ps4_constVPD)
write.csv(out_Ps4_constVPD, "Ps4_constVPD.csv")

save(out_Ps0.5_constVPD, out_Ps2_constVPD, out_Ps4_constVPD,
     file = "data/out/theoretical_sims_constVPD.Rdata")
load("data/out/theoretical_sims_constVPD.Rdata")

### Constant RH --------------------------

VPD = RHtoVPD(RH = 60, TdegC = Tair_vec)

# Ps = -0.5 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 0.5)
out_Ps0.5_constRH = get_predictions(Tair_sim.df)
Ci <- 400.0 - out_Ps0.5_constRH$A * 1E-4 / (1.57 * out_Ps0.5_constRH$gs) * 100^2
out_Ps0.5_constRH$Ci_Fick = Ci
#print(out_Ps0.5_constRH)
write.csv(out_Ps0.5_constRH, "Ps0.5_constRH.csv")

# Ps = -2 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 2)
out_Ps2_constRH = get_predictions(Tair_sim.df)
Ci <- 400.0 - out_Ps2_constRH$A * 1E-4 / (1.57 * out_Ps2_constRH$gs) * 100^2
out_Ps2_constRH$Ci_Fick = Ci
#print(out_Ps2_constRH)
write.csv(out_Ps2_constRH, "Ps2_constRH.csv")

# Ps = -4 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 4)
out_Ps4_constRH = get_predictions(Tair_sim.df)
Ci <- 400.0 - out_Ps4_constRH$A * 1E-4 / (1.57 * out_Ps4_constRH$gs) * 100^2
out_Ps4_constRH$Ci_Fick = Ci
#print(out_Ps4_constRH)
write.csv(out_Ps4_constRH, "Ps4_constRH.csv")

save(out_Ps0.5_constRH, out_Ps2_constRH, out_Ps4_constRH,
     file = "data/out/theoretical_sims_constRH.Rdata")
load("data/out/theoretical_sims_constRH.Rdata")

# PLOTTING #####################################################################

## Constant RH -----------------------

# Combine dataframes into one list
outputs_constRH = list(out_Ps0.5_constRH, out_Ps2_constRH, out_Ps4_constRH)

# Generate composite plots
gsAPleafvsTleaf_constRH.plt = plot_composite_Trange(outputs_constRH, vars = "gs, A, Pleaf")
EDvsTleaf_constRH.plt = plot_composite_Trange(outputs_constRH, vars = "E, Dleaf")

# Save plots
ggsave(gsAPleafvsTleaf_constRH.plt,
       filename = "figs/Fig1_TheoreticalSims_gsAPleafvsTleaf_constRH.tiff", width = 12.25, height = 10)
ggsave(EDvsTleaf_constRH.plt,
       filename = "figs/TheoreticalSims_EDvsTleaf_constRH.tiff", width = 12.25, height = 6.5)

## Constant VPD ----------------------

# Combine dataframes into one list
outputs_constVPD = list(out_Ps0.5_constVPD, out_Ps2_constVPD, out_Ps4_constVPD)

# Generate composite plots
gsAPleafvsTleaf_constVPD.plt = plot_composite_Trange(outputs_constVPD, vars = "gs, A, Pleaf")
EDvsTleaf_constVPD.plt = plot_composite_Trange(outputs_constVPD, vars = "E, Dleaf")

# Save plots
ggsave(gsAPleafvsTleaf_constVPD.plt,
       filename = "figs/FigS3_TheoreticalSims_gsAPleafvsTleaf_constVPD.tiff",
       width = 12.25, height = 10, bg = "white")
ggsave(EDvsTleaf_constVPD.plt,
       filename = "figs/TheoreticalSims_EDvsTleaf_constVPD.tiff",
       width = 12.25, height = 6.5, bg = "white")

## THIS BIT FROM HERE DOES NOT WORK ----------------------

# Plot Tleaf vs Tair
Tleaf_plts.l = lapply(outputs_constRH, function(df) {
  Tleaf.plt = plot_physio_vs_Tleaf(df, all = FALSE, yvar = "Tleaf", xvar = "Tair") %>%
    add_Tthreshold_lines()
  return(Tleaf.plt)
})
Tleaf_comb_plt = plot_grid(Tleaf_plts.l[[1]] +
                             ggtitle(expression(bold(psi[s]*" = -0.5 MPa"))) +
                             theme(legend.position = "none",
                                   axis.title.x = element_blank(),
                                   plot.title = element_text(size = 12, face = "bold")),
                           Tleaf_plts.l[[2]]  +
                             ggtitle(expression(bold(psi[s]*" = -2 MPa"))) +
                             theme(legend.position = "none",
                                   axis.title.x = element_blank(),
                                   axis.title.y = element_blank(),
                                   axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   plot.title = element_text(size = 12, face = "bold")),
                           Tleaf_plts.l[[3]] +
                             ggtitle(expression(bold(psi[s]*" = -4 MPa"))) +
                             theme(legend.position = "none",
                                   axis.title.x = element_blank(),
                                   axis.title.y = element_blank(),
                                   axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   plot.title = element_text(size = 12, face = "bold")),
                           nrow = 1, rel_widths = c(1, 0.85, 0.85))
TleafvsTair.plt = plot_grid(legend,
                            annotate_figure(Tleaf_comb_plt, bottom = text_grob(expression("T"[air]*" (\u00B0C)"), size = 16)),
                            nrow = 2, rel_heights = c(.1,1)) +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
ggsave(TleafvsTair.plt,
       filename = "figs/TheoreticalSims_TleafvsTair.tiff", width = 12.25, height = 5.5)
