# Write pdf file with all figures
# by Camille Sicangco
# Created 3 July 2025

# Main
main_plots = list(
  gsAPleafvsTleaf_constRH.plt, # Fig 1
  Fig2,
  Fig3,
  Fig4 = plts.l[[1]],
  AEvT.plt, # Fig 5
  Tleaf_pred_obs.plt, # Fig 6
  Fig7_PleafvT)

plots_uniform <- lapply(main_plots, function(p) {
  ggdraw(p)  # or add margins manually
})

pdf("figs/main_figures.pdf", width = 11, height = 7)
for (p in plots_uniform) {
  print(p)
}
dev.off()

# Supplementary
supp_plots = list(FigS1,
                  FigS2,
                  gsAPleafvsTleaf_constVPD.plt,
                  FigS4,
                  FigS5,
                  FigS6 = plts.l[[3]])
supp_plots_uniform <- lapply(supp_plots, function(p) {
  ggdraw(p)  # or add margins manually
})

pdf("figs/supp_figures.pdf", width = 11, height = 7)
for (p in supp_plots_uniform) {
  print(p)
}
dev.off()
