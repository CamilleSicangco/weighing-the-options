# Sensitivity analysis of TC model to Tcrit and T50
# Created 27 May 2025

# Set environmental variables, with and without constant VPD
Tair_vec = seq(30,60, by = 1)
Tair_sim_constVPD.df = data.frame(Tair = Tair_vec, PPFD = 1500, VPD = 1.5,
                         Ps = 0.5)

VPD = RHtoVPD(RH = 60, TdegC = Tair_vec)
Tair_sim_constRH.df = data.frame(Tair = Tair_vec, PPFD = 1500, VPD,
                                  Ps = 0.5)
# GENERATE PREDICTIONS #########################################################

# Hold T50 constant, vary Tcrit
preds_varTcrit_constVPD = make_pred_Tthresholds(Tair_sim_constVPD.df)
preds_varTcrit_constRH = make_pred_Tthresholds(Tair_sim_constRH.df)

# Hold Tcrit constant, vary T50
preds_varT50_constVPD = make_pred_Tthresholds(Tair_sim_constVPD.df, 
                                     hold_Tcrit = TRUE,
                                     Thold_val = 46.5,
                                     Tvar_vals = c(47.5, 48.5, 49.5, 50.4))
preds_varT50_constRH = make_pred_Tthresholds(Tair_sim_constRH.df, 
                                              hold_Tcrit = TRUE,
                                              Thold_val = 46.5,
                                              Tvar_vals = c(47.5, 48.5, 49.5, 50.4))

# List all predictions
preds.l = list(preds_varTcrit_constRH, preds_varTcrit_constVPD,
                        preds_varT50_constRH, preds_varT50_constVPD)
names(preds.l) = c("varTcrit_constRH", "varTcrit_constVPD",
                   "varT50_constRH", "varT50_constVPD")
# PLOTTING #####################################################################

# Generate plots
plts.l = lapply(1:length(preds.l), function(i) {
  
  palette = scales::seq_gradient_pal("slategray1","darkslateblue", "Lab")(seq(0,1,length.out=4))
  
  # Specify if Tcrit or T50 are varied
  if (isTRUE(grepl("Tcrit", names(preds.l)[i]))) {
    plt = preds.l[[i]] %>% 
      ggplot(aes(x = Tleaf, y = gs, linetype = Tcrit, color = Tcrit))
    legend_label = expression("T"[crit]*" (\u00B0C)")
  } else {
    plt = preds.l[[i]] %>% 
      ggplot(aes(x = Tleaf, y = gs, linetype = T50, color = T50))
    legend_label = expression("T"[50]*" (\u00B0C)")
  }
  # Create plots
  plt = plt + 
    geom_line(linewidth = 1) + 
    theme_classic() +
    ylab(expression("g"[s]*" (mol m"^-2*"s"^-1*")")) +
    xlab(expression("T"[leaf]*" (\u00B0C)")) +
    guides(linetype = guide_legend(title = legend_label),
           color = guide_legend(title = legend_label)) +
    scale_colour_manual(values = palette) +
    theme(axis.title = element_text(size = 14))
})
names(plts.l) = names(preds.l)

ggsave("figs/Fig4_SA_Tcrit_constRH.tiff", plts.l[[1]], height = 7, width = 11, bg = "white")
ggsave("figs/SA_Tcrit_constVPD.tiff", plts.l[[2]], height = 7, width = 11, bg = "white")
ggsave("figs/SA_T50_constRH.tiff", plts.l[[3]], height = 7, width = 11, bg = "white")
ggsave("figs/SA_T50_constVPD.tiff", plts.l[[4]], height = 7, width = 11, bg = "white")