# Fig 1: Overview of instantaneous simulation models
# by Camille Sicangco
# Created 23 May 2025

# Plotting functions -----------------------------------------------------------

# Create composite plot of costs/gains and profit
composite_plot = function(df) {
  
  ymin = min(df$cost_gain[!is.na(df$cost_gain)])
  ymax = max(df$cost_gain[!is.na(df$cost_gain)])
  
  min_P = min(df$P)
  df = df %>%  
    filter(P != min_P)
  
  # Sperry et al. 2017
  p_A = df %>% 
    filter(ID %in% c("HC", "CG_gross")) %>% 
    plot_costgain() + 
    ggtitle("Sperry") +
    theme(plot.title = element_text(size = 11, face = "bold")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3), limits = c(ymin, ymax)) 
  
  # CGnet with unconstrained Ci
  p_C = df %>% 
    filter(ID %in% c("HC", "CG_net_uncorr")) %>% 
    plot_costgain() + 
    ggtitle(expression(bold("Sperry + CG"[net]))) +
    theme(plot.title = element_text(size = 11, face = "bold")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3), limits = c(ymin, ymax)) 
  
  # CGnet with unconstrained Ci plus TC
  p_D = df %>% 
    filter(ID %in% c("HC", "CG_net_uncorr", "TC")) %>%     
    plot_costgain() + 
    ggtitle(expression(bold("Sperry + CG"[net]*" + TC"))) +
    theme(plot.title = element_text(size = 11, face = "bold")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3), limits = c(ymin, ymax)) 
  
  # CGnet with unconstrained Ci plus TC and RC
  p_E = df %>% 
    filter(ID %in% c("HC", "CG_net_uncorr", "TC", "RC")) %>%     # remove point with E = 0
    plot_costgain() + 
    ggtitle(expression(bold("Sperry + CG"[net]*" + TC + RC"))) +
    theme(plot.title = element_text(size = 11, face = "bold")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3), limits = c(ymin, ymax)) 
  
  # Extract legend from last plot
  legend = get_legend(p_E + 
                        guides(color = guide_legend(nrow = 1)) +
                        theme(legend.position = "bottom",
                              legend.text = element_text(size = 10),
                              legend.title = element_text(size = 12)))
  
  # Combine cost/gain plots
  comb_plt = ggarrange(p_A + theme(legend.position = "none",
                                   axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank()), 
                       p_C + theme(legend.position = "none",
                                   axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank()), 
                       p_D + theme(legend.position = "none"), 
                       p_E + theme(legend.position = "none",
                                   axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank()), 
                       nrow = 2, ncol = 2)
  comb_plt = annotate_figure(comb_plt,
                             bottom = "Leaf water potential (-MPa)", 
                             left = "Costs/gains")
  comb_plt = plot_grid(legend, comb_plt, nrow = 2, rel_heights = c(.05,1))
  
  # Plot profit
  profit.plt = plot_profit(df)
  
  # Combine cost/gain and profit plots
  p = ggarrange(comb_plt + ggtitle("A)") +
                  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0)), 
                profit.plt + ggtitle("B)") +
                  theme(plot.title = element_text(size = 15, face = "bold"), hjust = 0), 
                ncol = 2) +
    theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
  return(p)
}

# Plot costs and gains
plot_costgain = function(df, separate = TRUE)
{
  palette = c(HC = "#1E88E5", TC = "#FFC107", RC = "grey", CG_net = "#D81B60", CG_net_uncorr = "#D81B60", CG_gross = "#D81B60")
  linetype = c(HC = "solid", TC = "solid", RC = "solid", CG_net = "solid", CG_net_uncorr = "dotted", CG_gross = "dashed")
  
  if (isFALSE(separate)) {
    p = ggplot(df, aes(x = P, y = cost_gain)) + 
      geom_line(aes(color = ID, linetype = ID), linewidth = 1) + 
      theme_classic() + 
      scale_color_manual(values = palette, name = "Cost/gain") +
      xlab("Leaf water potential (-MPa)") +
      ylab("HC, TC, RC, CG") +
      expand_limits(x = 0) +
      scale_linetype_manual(values = linetype, "Cost/gain") +
      theme(text = element_text(size = 20)) +
      geom_hline(yintercept = 0, linetype = "dashed")
  } else {
    p = df %>% 
      ggplot(aes(x = P, y = cost_gain)) + 
      geom_line(aes(color = ID), linewidth = 1) + 
      theme_classic() + 
      scale_color_manual(values = palette, name = "Cost/gain",
                         labels = c("CG", "HC", "RC", "TC")) +
      expand_limits(x = 0) +
      theme(text = element_text(size = 20)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) + 
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank())
  }
  
  return(p)
}

# Plot net profit
plot_profit = function(df) {
  linetype = c(Sperry = "solid", CGnet = "dashed", CGnet_TC = "dotted",
               CGnet_TC_RC = "dotdash")
  
  min_P = min(df$P)
  df = df %>%  
    filter(P != min_P) %>% 
    pivot_wider(names_from = "ID", values_from = "cost_gain") %>% 
    mutate(Sperry = CG_gross - HC,
           CGnet = CG_net_uncorr - HC,
           CGnet_TC = CG_net_uncorr - (HC +TC),
           CGnet_TC_RC = CG_net_uncorr - (HC + TC + RC)) %>%
    select(P, Sperry:CGnet_TC_RC) %>% 
    pivot_longer(cols = Sperry:CGnet_TC_RC, names_to = "Model", values_to = "Profit") 
  
  
  df$Model = factor(df$Model, 
                    levels = c("Sperry", "CGnet", "CGnet_TC", "CGnet_TC_RC")
                    )
  
  p = df %>% 
    ggplot(aes(x = P, y = Profit, linetype = Model)) + 
    scale_linetype_manual(values = linetype, 
                          labels = c("Sperry",
                                     expression("Sperry + CG"[net]),
                                     expression("Sperry + CG"[net]*" + TC"),
                                     expression("Sperry + CG"[net]*" + TC + RC"))) +
    geom_line(linewidth = 1) + 
    theme_classic() +
    xlab("Leaf water potential (-MPa)")
  return(p)
}

# Set variables/params ---------------------------------------------------------
# Environment
Tair = 25
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
kmax_25 = 0.7
Pcrit = calc_Pcrit(b, c)
P = Ps_to_Pcrit(Ps, Pcrit)

# CGnet normalized with or without including the point where E = 0
CG_net = C_gain_corr(P, b, c, Amax = NULL, kmax_25, Tair, VPD, PPFD, 
                     Wind = 8, Wleaf = 0.02, LeafAbs = 0.5,
                     Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
                     Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92,
                     constant_kmax = TRUE, net = TRUE, netOrig = TRUE)
CG_net_uncorr = C_gain(P, b, c, Amax = NULL, kmax_25, Tair, VPD, PPFD, 
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
                                  Tcrit = Tcrit, T50 = T50, constant_kmax = TRUE,
                                  Wind = 8, Wleaf = 0.02, LeafAbs = 0.5,
                                  Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
                                  Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92)
composite_plot(cost_gain30)
ggsave(filename = "figs/composite_plot_30deg.pdf", width = 12.25, height = 6.5)

# Tair = 40
cost_gain40 = calc_costgain_netorig(P, b, c, kmax_25 = kmax_25, 
                                  Tair = 40, PPFD = PPFD, VPD = VPD,
                                  Tcrit = Tcrit, T50 = T50, constant_kmax = TRUE,
                                  Wind = 8, Wleaf = 0.02, LeafAbs = 0.5,
                                  Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
                                  Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92)
composite_plot(cost_gain40)
ggsave(filename = "figs/composite_plot_40deg.pdf", width = 12.25, height = 6.5)

# Tair = 48
cost_gain48 = calc_costgain_netorig(P, b, c, kmax_25 = kmax_25, 
                                  Tair = 48, PPFD = PPFD, VPD = VPD,
                                  Tcrit = Tcrit, T50 = T50, constant_kmax = TRUE,
                                  Wind = 8, Wleaf = 0.02, LeafAbs = 0.5,
                                  Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
                                  Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92)
composite_plot(cost_gain48)
ggsave(filename = "figs/composite_plot_48deg.pdf", width = 12.25, height = 6.5)

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
