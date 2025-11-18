# Model testing over a temperature range
# by Camille Sicangco

# RUN SIMULATIONS ##############################################################

### Constant VPD -------------------------

# Environment
Tair_vec = seq(20,60, by = 1)
PPFD = 1500
VPD = 1.5

# Medlyn model: fit b
b_USO = log(((0.05 - 1e-5)/1.6*400/5.89 - 1) * sqrt(1.5) / 2.9)/(-0.5 + 0.2)
g1_alt = ((0.05 - 1e-5)/1.6*400/5.89 - 1)*sqrt(1.5) / exp(0.55*(-0.5 + 0.3))

# Ps = -0.5 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 0.5)
out_Ps0.5_constVPD = get_predictions(Tair_sim.df, g1 = g1_alt)
test = get_predictions(Tair_sim.df, kmax_25 = 2)
test = make_pred(kmax_25 = 4, Ps = 0.5,
                 Tair = 30, PPFD = PPFD, VPD = VPD,
                 Tcrit = Tcrit, T50 = T50,
                 Wind = 8, Wleaf = 0.025, LeafAbs = 0.5,
                 Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
                 Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92)
test$E*(pi * (3.25*3/8)^2)/LeafArea_df$LeafArea
  
# Ps = -2 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 2)
out_Ps2_constVPD = get_predictions(Tair_sim.df, g1 = g1_alt)

# Ps = -4 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 4)
out_Ps4_constVPD = get_predictions(Tair_sim.df, g1 = g1_alt)

save(out_Ps0.5_constVPD, out_Ps2_constVPD, out_Ps4_constVPD, 
     file = "data/out/theoretical_sims_constVPD.Rdata")
load("data/out/theoretical_sims_constVPD.Rdata")

### Constant RH --------------------------

VPD = RHtoVPD(RH = 60, TdegC = Tair_vec)

# Ps = -0.5 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 0.5)
out_Ps0.5_constRH = get_predictions(Tair_sim.df, g1 = g1_alt)

test= get_predictions(Tair_sim.df, kmax_25 = 4)
test$E0 = test$E
test$E = test$E0*(pi * (3.25*3/8)^2)/mean(LeafArea_df$LeafArea)
calc_kmax(0.5, 40) * mean(LeafArea_df$LeafArea)/(pi * (3.25*3/8)^2)
6*(pi * (3.25*3/8)^2)/mean(LeafArea_df$LeafArea)
plot_physio_vs_Tleaf(test, all = FALSE, yvar = "E") # Leaf area basis
plot_physio_vs_Tleaf(filter(test, Model == "Sperry"), all = FALSE, yvar = "E") # Canopy area basis

test2= get_predictions(Tair_sim.df, kmax_25 = 0.5)
plot_physio_vs_Tleaf(out_Ps0.5_constRH_newb, all = FALSE, yvar = "E")

# Ps = -2 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 2)
out_Ps2_constRH = get_predictions(Tair_sim.df, g1 = g1_alt)

# Ps = -4 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 4)
out_Ps4_constRH = get_predictions(Tair_sim.df, g1 = g1_alt)

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
