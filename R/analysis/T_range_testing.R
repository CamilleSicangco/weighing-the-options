# Model testing over a temperature range
# by Camille Sicangco

# Environment
Tleaf_seq = seq(20,70, by = 1)

Vcs = TVcmax(Tleaf_seq,EaV=62307,EdVC=2e5,delsC=639)
Js = TVcmax(Tleaf_seq, EaV=33115,EdVC=2e5,delsC=635)

Vcs_2 = TVcmax_updated(Tleaf_seq,EaV=62307,EdVC=2e5,delsC=639)
Js_2 = TVcmax_updated(Tleaf_seq, EaV=33115,EdVC=2e5,delsC=635)

plot(Tleaf_seq, Vcs)
plot(Tleaf_seq, Js)

plot(Tleaf_seq, Vcs_2)
plot(Tleaf_seq, Js_2)



## Run simulations -------------------------------------------------------------

# Environment
Tair_vec = seq(20,60, by = 1)
Ps_vec = seq(0, 4, by = 0.2)

# Ps = -0.5 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = 1500, VPD = 1.5, Ps = 0.5)
out_Ps0.5 = get_preds_theoretical_sims(Tair_sim.df, kmax_25 = 0.5)

# Ps = -2 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = 1500, VPD = 1.5, Ps = 2)
out_Ps2 = get_preds_theoretical_sims(Tair_sim.df, kmax_25 = 0.5)

# Ps = -4 MPa
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = 1500, VPD = 1.5, Ps = 4)
out_Ps4 = get_preds_theoretical_sims(Tair_sim.df, kmax_25 = 0.5)

# Plotting #####################################################################

# Test for one data frame
out_Ps0.5 %>%
  #filter(Model != "Sperry + CGnet") %>% 
  ggplot(aes(x = Tleaf, y = gs, color = Model)) +
  geom_point(size = 3, shape = 1) + 
  theme_classic()
plot_physio_vs_Tleaf(out_Ps0.5, all = FALSE, yvar = "gs")

# Composite plots
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
legend = ggpubr::get_legend(plts.l[[1]][[1]] +
                      guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
                      theme(legend.position = "bottom",
                            legend.text = element_text(size = 12),
                            legend.title = element_text(size = 14)))

# Combine all plots for each physio variable
comb_plts.l = lapply(seq_along(plts.l), function(i) {
  
  widths = if (names(plts.l)[i] %in% c("gs")) {
    c(1,0.83,0.83)
  } else {
    c(1,0.83,0.83)
  }
  
  labs = 
    if (names(plts.l)[i] %in% c("gs", "E")) {
      c("(a)", "(b)", "(c)")
    } else if  (names(plts.l)[i] %in% c("A", "Dleaf")) {
      c("(d)", "(e)", "(f)")
    } else if  (names(plts.l)[i]  == "P") {
      c("(g)", "(h)", "(i)")
    }
  
  plts.l[[i]][[1]] = plts.l[[i]][[1]] +
    ggtitle(labs[1]) +
    theme(plot.title = element_text(size = 12, face = "bold"))
  
  plts.l[[i]][[2]] = plts.l[[i]][[2]] +
    ggtitle(labs[2]) +
    theme(plot.title = element_text(size = 12, face = "bold"))
  
  plts.l[[i]][[3]] = plts.l[[i]][[3]] +
    ggtitle(labs[3]) +
    theme(plot.title = element_text(size = 12, face = "bold"))
  
  if (names(plts.l)[i] %in% c("gs", "E")) {
    plts.l[[i]][[1]] = plts.l[[i]][[1]] +
      labs(title = expression(psi[s]*" = -0.5 MPa"),
           subtitle = expression(bold("(a)"))) +
      theme(plot.subtitle = element_text(size = 12),
            plot.title = element_text(size = 12))
    
    plts.l[[i]][[2]] = plts.l[[i]][[2]]+
      labs(title = expression(psi[s]*" = -2 MPa"),
           subtitle = expression(bold("(b)"))) +
      theme(plot.subtitle = element_text(size = 12),
            plot.title = element_text(size = 12))
    
    plts.l[[i]][[3]] = plts.l[[i]][[3]] +
      labs(title = expression(psi[s]*" = -4 MPa"),
           subtitle = expression(bold("(c)"))) +
      theme(plot.subtitle = element_text(size = 12),
            plot.title = element_text(size = 12))
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
gsAPleafvsTleaf.plt = 
  plot_grid(
    legend, 
    comb_plts.l[[1]],
    comb_plts.l[[3]],
    annotate_figure(comb_plts.l[[2]], bottom = text_grob(expression("T"[leaf]*" (\u00B0C)"), size = 16)), 
    nrow = 4, rel_heights = c(.2,1.05,1,1)
  ) +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
ggsave(gsAPleafvsTleaf.plt,
       filename = "figs/TheoreticalSims_gsAPleafvsTleaf_nopts.pdf", width = 12.25, height = 10)

# Plot A, E, and Dleaf vs Tleaf
EDvsTleaf.plt = 
  plot_grid(
    legend, 
    comb_plts.l[[4]],
    annotate_figure(comb_plts.l[[5]], bottom = text_grob(expression("T"[leaf]*" (\u00B0C)"), size = 16)), 
    nrow = 3, rel_heights = c(.1,1.05,1)) +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
ggsave(EDvsTleaf.plt,
       filename = "figs/TheoreticalSims_EDvsTleaf.pdf", width = 12.25, height = 6.5)

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
