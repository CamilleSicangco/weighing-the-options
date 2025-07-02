# Model testing over a temperature range
# by Camille Sicangco

# RUN SIMULATIONS ##############################################################

### Constant VPD -------------------------

# Environment
Tair_vec = seq(20,60, by = 1)
PPFD = 1500
VPD = 1.5

# Ps = -0.5 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 0.5)
out_Ps0.5_constVPD = get_predictions(Tair_sim.df)

# Ps = -2 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 2)
out_Ps2_constVPD = get_predictions(Tair_sim.df)

# Ps = -4 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 4)
out_Ps4_constVPD = get_predictions(Tair_sim.df)

save(out_Ps0.5_constVPD, out_Ps2_constVPD, out_Ps4_constVPD, 
     file = "data/out/theoretical_sims_constVPD.Rdata")
load("data/out/theoretical_sims_constVPD.Rdata")

### Constant RH --------------------------

VPD = RHtoVPD(RH = 60, TdegC = Tair_vec)

# Ps = -0.5 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 0.5)
out_Ps0.5_constRH = get_predictions(Tair_sim.df)

# Ps = -2 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 2)
out_Ps2_constRH = get_predictions(Tair_sim.df)

# Ps = -4 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = PPFD, VPD = VPD, Ps = 4)
out_Ps4_constRH = get_predictions(Tair_sim.df)

save(out_Ps0.5_constRH, out_Ps2_constRH, out_Ps4_constRH, 
     file = "data/out/theoretical_sims_constRH.Rdata")
load("data/out/theoretical_sims_constRH.Rdata")

# PLOTTING #####################################################################

## Constant VPD ----------------------

# Combine dataframes into one list
outputs_constVPD = list(out_Ps0.5_constVPD, out_Ps2_constVPD, out_Ps4_constVPD)

# Generate composite plots
gsAPleafvsTleaf_constVPD.plt = plot_composite_Trange(outputs_constVPD, vars = "gs, A, Pleaf")
EDvsTleaf_constVPD.plt = plot_composite_Trange(outputs_constVPD, vars = "E, Dleaf")

# Save plots
ggsave(gsAPleafvsTleaf_constVPD.plt,
       filename = "figs/FigS3_TheoreticalSims_gsAPleafvsTleaf_constVPD.tiff", width = 12.25, height = 10)
ggsave(EDvsTleaf_constVPD.plt,
       filename = "figs/TheoreticalSims_EDvsTleaf_constVPD.pdf", width = 12.25, height = 6.5)


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
       filename = "figs/TheoreticalSims_EDvsTleaf_constRH.pdf", width = 12.25, height = 6.5)

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
       filename = "figs/TheoreticalSims_TleafvsTair.pdf", width = 12.25, height = 5.5)
