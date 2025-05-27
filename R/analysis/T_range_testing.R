# Model testing over a temperature range
# by Camille Sicangco

# Environment
Tair_vec = seq(30,50, by = 1)
Ps_vec = seq(0, 4, by = 0.2)

get_preds_theoretical_sims = 
  function(
    df, Tcrit = 43.4, T50 = 49.6, P50 = 4.07, P88 = 5.50,
    Wind = 8, Wleaf = 0.02, LeafAbs = 0.5, 
    Vcmax=34,EaV=51780,EdVC=2e5,delsC=639, 
    Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=635, Rd0 = 0.92,
    constant_kmax = TRUE, net = TRUE, netOrig = TRUE,
    g1 = 2.9,g0=0.003
  ) {
    out = make_pred(df = df, models = "final",
                    Tcrit_hw = Tcrit, T50_hw = T50, P50 = P50, P88 = P88,
                    Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs, 
                    Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC, 
                    Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ, Rd0 = Rd0,
                    constant_kmax = constant_kmax, net = net, netOrig = netOrig) 
    
    # Make Medlyn predictions
    Medlyn_preds = plantecophys::PhotosynEB(Tair=df$Tair, VPD=df$VPD, PPFD = df$PPFD,
                                            Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs, 
                                            g1 = g1,g0=g0,
                                            Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC, 
                                            Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ, Rd0 = Rd0)
    Medlyn_preds$P = get_Pleaf_Medlyn(Medlyn_preds, P50 = P50, P88 = P88, 
                                      Ps = df$Ps[1], kmax_25 = kmax_25, constant_kmax = constant_kmax)
    
    Medlyn_preds_sim = Medlyn_preds %>%
      rename(E = ELEAF, Dleaf = VPDleaf, gs = GS, A = ALEAF) %>%
      mutate(Model = "Medlyn", .before = 1) %>%
      select(Model, Tair, E, Tleaf, Dleaf, gs, A, P)
    
    # Combine all predictions
    out_all = bind_rows(out, Medlyn_preds_sim)
    return(out_all)
  }


get_Pleaf_Medlyn = function(pred_df, P50, P88, Ps, kmax_25, constant_kmax) {
  
  # Hydraulics
  Weibull = fit_Weibull(P50, P88)
  b = Weibull[1,1]
  c = Weibull[1,2]
  Pcrit = calc_Pcrit(b, c)
  P = Ps_to_Pcrit(Ps, Pcrit)
  
  E = Medlyn_preds$ELEAF
  Tair = Medlyn_preds$Tair
  
  Pleaf_vec = sapply(1:length(E), function(i) {
    E_vec = trans_from_vc(P, kmax_25, Tair[i], b, c, constant_kmax)
    j = which.min(abs(E[i] - E_vec))
    Pleaf = P[j]
    return(Pleaf)
  })
  
  return(Pleaf_vec)
}

## Run simulations -------------------------------------------------------------

# Ps = -0.5 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = 1500, VPD = 1.5, kmax = 0.7,
                         Ps = 0.5, HWtrt = "HW")
out_Ps0.5 = get_preds_theoretical_sims(Tair_sim.df, constant_kmax = TRUE)#,
                                               #T50 = 45.4)

# Ps = -2 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = 1500, VPD = 1.5, kmax = 0.7,
                         Ps = 2, HWtrt = "HW")
out_Ps2 = get_preds_theoretical_sims(Tair_sim.df, constant_kmax = TRUE)#,
                                             #Tcrit = 48.6)

# Ps = -4 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = 1500, VPD = 1.5, kmax = 0.7,
                         Ps = 4, HWtrt = "HW")
out_Ps4 = get_preds_theoretical_sims(Tair_sim.df, constant_kmax = TRUE)#,
                                             #Tcrit = 48.6)

# Plotting #####################################################################

plot_physio_vs_Tleaf = function(df, all = TRUE, yvar = NULL, xvar = "Tleaf") {
  
  linetype = c(Sperry = "solid", "Sperry + CGnet" = "dashed", "Sperry + CGnet + TC" = "dotted",
               Medlyn = "dotdash")
  shapes = c(Sperry = 1, "Sperry + CGnet" = 2, "Sperry + CGnet + TC" = 3,
            Medlyn = 4)
  palette = c(Sperry = "#88CCEE", "Sperry + CGnet" = "#332288", "Sperry + CGnet + TC" = "#DDCC77",
              Medlyn = "#CC6677")
  
  if (isTRUE(all)) {
    p = df %>% pivot_longer(
      cols = P:A,
      names_to = "var",
      values_to = "pred"
    ) %>% 
      ggplot(aes(x = Tair, y = pred, color = Model, linetype = Model)) + 
      geom_point() + geom_line() +
      facet_wrap(vars(var), scales = "free") + 
      theme_classic()  +
      scale_color_manual(values = palette, 
                         labels = c("Medlyn",
                                    "Sperry",
                                    expression("Sperry + CG"[net]),
                                    expression("Sperry + CG"[net]*" + TC")
                         )) +
      scale_linetype_manual(values = linetype, 
                            labels = c("Medlyn",
                                       "Sperry",
                                       expression("Sperry + CG"[net]),
                                       expression("Sperry + CG"[net]*" + TC")
                            )) +
      ylab("Predictions") +
      xlab(expression("T"[air]*" (\u00B0C)"))
  } else {
    ylabel = 
      if (yvar == "gs") {
        expression("g"[s]*" (mol m"^-2*"s"^-1*")")
      } else if (yvar == "P") {
        expression(psi[leaf]*" (-MPa)")
      } else if (yvar == "Tleaf") {
        expression("T"[leaf]*" (\u00B0C)")
      } else if (yvar == "A") {
        expression("A (" * mu * "mol m"^-2*"s"^-1*")")
      } else if (yvar == "Dleaf") {
        expression("VPD"[leaf]*" (kPa)")
      } else if (yvar == "E") {
        expression("E"*" (mmol m"^-2*"s"^-1*")")
      }
    xlabel = 
      if (xvar == "Tleaf") {
        expression("T"[leaf]*" (\u00B0C)")
      } else if (xvar == "Tair") {
        expression("T"[air]*" (\u00B0C)")
      } else {
        "oops!"
      }
    
    p = df %>%
      #filter(Tleaf > 40) %>% 
      ggplot(aes(x = !!sym(xvar), y = !!sym(yvar), color = Model, linetype = Model)) +
      #geom_point(size = 3) + 
      geom_line(size = 1) +
      theme_classic() +
      ylab(ylabel) +
      xlab(xlabel) +
      scale_color_manual(values = palette, 
                         labels = c("Medlyn",
                                    "Sperry",
                                    expression("Sperry + CG"[net]),
                                    expression("Sperry + CG"[net]*" + TC")
                         )) +
      scale_linetype_manual(values = linetype, 
                            labels = c("Medlyn",
                                       "Sperry",
                                       expression("Sperry + CG"[net]),
                                       expression("Sperry + CG"[net]*" + TC")
                            )) +
      theme(text = element_text(size = 14),
            axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
            axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
  }
  return(p)
}

add_Tthreshold_lines = function(plt, Tcrit = 43.4, T50 = 49.6) {
  #r = 2/(T50 - Tcrit)
  #T95 = log(1/0.95-1)/-r + T50
  
  p = plt +
    geom_hline(yintercept = Tcrit, linetype = "dashed") +
    annotate("text", x=40, y=Tcrit, label="Tcrit", vjust = -0.5) +
    geom_hline(yintercept = T50, linetype = "dotted", col = "darkorange2") +
    annotate("text", x=40, y=T50, label="T50", vjust = -0.5, col = "darkorange2") +
    #geom_hline(yintercept = T95, linetype = "dotdash", col = "red") +
    #annotate("text", x=40, y=T95, label="T95", vjust = -0.5, col = "red") +
    geom_abline(slope = 1, intercept = 0, col = "grey") +
    annotate("text", x = 40, y = 40, label = "1:1 line", col = "grey", vjust = 3)
  return(p)
}

physio_vars = c("gs", "P", "A", "E", "Dleaf")
outputs.l = list(out_Ps0.5, out_Ps2, out_Ps4)

# Create list of plots of each physio variable, for Ps = -0.5, -2, and -4 MPa
plts.l = lapply(physio_vars, function(var) {
  lapply(outputs.l, function(df) {
    plt = plot_physio_vs_Tleaf(df, all = FALSE, var)
  })
})
names(plts.l) = physio_vars

# Extract legend for plots
legend = get_legend(plts.l[[1]][[1]] +
                      guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
                      theme(legend.position = "bottom",
                            legend.text = element_text(size = 12),
                            legend.title = element_text(size = 14)))

# Combine all plots for each physio variable
comb_plts.l = lapply(seq_along(plts.l), function(i) {
  
  widths = if (names(plts.l)[i] %in% c("gs", "P")) {
    c(1,0.83,0.83)
  } else {
    c(1,0.83,0.83)
  }
  
  if (names(plts.l)[i] %in% c("gs", "A")) {
    plts.l[[i]][[1]] = plts.l[[i]][[1]] +
      ggtitle(expression(bold(psi[s]*" = -0.5 MPa"))) +
      theme(plot.title = element_text(size = 12, face = "bold"))
    
    plts.l[[i]][[2]] = plts.l[[i]][[2]] +
      ggtitle(expression(bold(psi[s]*" = -2 MPa"))) +
      theme(plot.title = element_text(size = 12, face = "bold"))
    
    plts.l[[i]][[3]] = plts.l[[i]][[3]] +
      ggtitle(expression(bold(psi[s]*" = -4 MPa"))) +
      theme(plot.title = element_text(size = 12, face = "bold"))
  }
  
  comp_plt = plot_grid(plts.l[[i]][[1]] + theme(legend.position = "none",
                                       axis.title.x = element_blank()), 
                       plts.l[[i]][[2]] + theme(legend.position = "none",
                                       axis.title.x = element_blank(),
                                       axis.title.y = element_blank(),
                                       axis.text.y = element_blank(),
                                       axis.ticks.y = element_blank()), 
                       plts.l[[i]][[3]] + theme(legend.position = "none",
                                       axis.title.x = element_blank(),
                                       axis.title.y = element_blank(),
                                       axis.text.y = element_blank(),
                                       axis.ticks.y = element_blank()), 
                       nrow = 1, rel_widths = widths
  )
  return(comp_plt)
})

# Plot gs and Pleaf vs Tleaf
gsPleafvsTleaf.plt = 
  plot_grid(
    legend, 
    comb_plts.l[[1]],
    annotate_figure(comb_plts.l[[2]], bottom = text_grob(expression("T"[leaf]*" (\u00B0C)"), size = 16)), 
    nrow = 3, rel_heights = c(.1,1.05,1)
  ) +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
ggsave(gsPleafvsTleaf.plt,
       filename = "figs/TheoreticalSims_gsPleafvsTleaf.pdf", width = 12.25, height = 6.5)

# Plot A, E, and Dleaf vs Tleaf
AEDvsTleaf.plt = 
  plot_grid(
    legend, 
    comb_plts.l[[3]], comb_plts.l[[4]],
    annotate_figure(comb_plts.l[[5]], bottom = text_grob(expression("T"[leaf]*" (\u00B0C)"), size = 16)), 
    nrow = 4, rel_heights = c(.1,1.05,1,1)) +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
ggsave(AEDvsTleaf.plt,
       filename = "figs/TheoreticalSims_AEDvsTleaf.pdf", width = 12.25, height = 10)

# Plot Tleaf vs Tair
Tleaf_plts.l = lapply(outputs.l, function(df) {
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
