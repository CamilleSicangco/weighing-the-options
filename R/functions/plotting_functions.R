# Plotting functions

# Instantaneous simulations ####################################################
# Create composite plot of costs/gains and profit
composite_plot = function(df) {
  df = df %>% 
    pivot_longer(cols = HC_constkmax:CG_net_newJT, 
                 values_to = "cost_gain",
                 names_to = "ID")
  ymin = min(df$cost_gain[!is.na(df$cost_gain)])
  ymax = max(df$cost_gain[!is.na(df$cost_gain)])
  
  min_P = min(df$P)
  df = df %>%  
    filter(P != min_P)
  
  # Sperry et al. 2017
  p_A = df %>% 
    filter(ID %in% c("HC_constkmax", "CG_gross_constkmax")) %>% 
    plot_costgain() + 
    ggtitle(expression(bold("(a)")*" ProfitMax")) +
    theme(plot.title = element_text(size = 12, face = "bold")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3), limits = c(ymin, ymax)) 

  # Sperry + variable kmax
  p_B = df %>% 
    filter(ID %in% c("HC_varkmax", "CG_gross_varkmax")) %>% 
    plot_costgain() + 
    ggtitle(expression(bold("(b)")*" ProfitMax"[k[max](T)])) +
    theme(plot.title = element_text(size = 12, face = "bold")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3), limits = c(ymin, ymax)) 
  
  
  # CGnet with unconstrained Ci
  p_C = df %>% 
    filter(ID %in% c("HC_varkmax", "CG_net")) %>% 
    plot_costgain() + 
    ggtitle(expression(bold("(c)")*" ProfitMax"[net])) +
    theme(plot.title = element_text(size = 12, face = "bold")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3), limits = c(ymin, ymax)) 
  
  # CGnet with unconstrained Ci plus TC
  p_D = df %>% 
    filter(ID %in% c("HC_varkmax", "CG_net_newJT", "TC")) %>%     
    plot_costgain() + 
    ggtitle(expression(bold("(d)")*" ProfitMax"[TC])) +
    theme(plot.title = element_text(size = 12, face = "bold")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3), limits = c(ymin, ymax)) 
  
  # CGnet with unconstrained Ci plus TC and RC
  #p_E = df %>% 
  #  filter(ID %in% c("HC", "CG_net", "TC", "RC")) %>%     # remove point with E = 0
  #  plot_costgain() + 
  #  ggtitle(expression(bold("Sperry + CG"[net]*" + TC + RC"))) +
  #  theme(plot.title = element_text(size = 11, face = "bold")) +
  #  scale_y_continuous(breaks = scales::pretty_breaks(n = 3), limits = c(ymin, ymax)) 
  
  # Extract legend from last plot
  legend = ggpubr::get_legend(p_D + 
                                guides(color = guide_legend(nrow = 1)) +
                                theme(legend.position = "bottom",
                                      legend.text = element_text(size = 10),
                                      legend.title = element_text(size = 12)))
  
  # Combine cost/gain plots
  comb_plt = ggarrange(p_A + theme(legend.position = "none",
                                   axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank()), 
                       p_B + theme(legend.position = "none",
                                   axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.text.x = element_blank(),
                                   axis.ticks.x = element_blank()), 
                       p_C + theme(legend.position = "none"), 
                       p_D + theme(legend.position = "none",
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
  p = ggarrange(comb_plt, 
                profit.plt + ggtitle("(e)") +
                  theme(plot.title = element_text(size = 12, face = "bold"), hjust = 0), 
                ncol = 2#, width = 12.25, height = 5.5
                ) +
    theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
  return(p)
}

## Helpers ---------------------------------------------------------------------

# Plot costs and gains
plot_costgain = function(df, separate = TRUE)
{
  palette = c(HC_constkmax = "#88CCEE", HC_varkmax = "#88CCEE", 
              TC = "#DDCC77", RC = "#117733",
              CG_net = "#CC6677", CG_net_newJT = "#CC6677", CG_gross_constkmax = "#CC6677", CG_gross_varkmax = "#CC6677")
  # palette = c(HC = "#1E88E5", TC = "#FFC107", RC = "grey", CG_net = "#D81B60", CG_net_uncorr = "#D81B60", CG_gross = "#D81B60")
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
                         labels = c("CG", "HC", "TC")) +
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
  palette = c(Sperry = "#1E88E5", Sperry_varkmax = "#332288", CGnet = "#e7d87d", CGnet_TC = "orange")
  #palette = c(Sperry = "#1E88E5", "Sperry + CGnet" = "#332288", "Sperry + CGnet + TC" = "#FFC107",
  #            Medlyn = "#D81B60")
  linetype = c(Sperry = "solid", Sperry_varkmax = "dotdash", CGnet = "dashed", CGnet_TC = "dotted")
  
  min_P = min(df$P)
  df = df %>%  
    filter(P != min_P) %>% 
    pivot_wider(names_from = "ID", values_from = "cost_gain") %>% 
    mutate(Sperry = CG_gross_constkmax - HC_constkmax,
           Sperry_varkmax = CG_gross_varkmax - HC_varkmax,
           CGnet = CG_net - HC_varkmax,
           CGnet_TC = CG_net_newJT - (HC_varkmax + TC)) %>%
    select(P, Sperry:CGnet_TC) %>% 
    pivot_longer(cols = Sperry:CGnet_TC, names_to = "Model", values_to = "Profit") 
  
  
  df$Model = factor(df$Model, 
                    levels = c("Sperry", "Sperry_varkmax", "CGnet", "CGnet_TC")
  )
  
  p = df %>% 
    ggplot(aes(x = P, y = Profit, linetype = Model, color = Model)) + 
    geom_line(linewidth = 1) + 
    scale_color_manual(values = palette, 
                       labels = c("ProfitMax",
                                  expression("ProfitMax"[k[max](T)]),
                                  expression("ProfitMax"[net]),
                                  expression("ProfitMax"[TC])
                       )) +
    scale_linetype_manual(values = linetype, 
                          labels = c("ProfitMax",
                                     expression("ProfitMax"[k[max](T)]),
                                     expression("ProfitMax"[net]),
                                     expression("ProfitMax"[TC])
                          )) +
    theme_classic() +
    xlab("Leaf water potential (-MPa)")
  
  return(p)
}
# T-range testing ##############################################################
plot_composite_Trange = function(outputs, vars = c("gs, A, Pleaf", "E, Dleaf")) {
  
  physio_vars = c("gs", "P", "A", "E", "Dleaf")
  
  # Create individual plots
  plts.l = lapply(physio_vars, function(var) {
    lapply(outputs, function(df) {
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
  
  # Format individual plots
  plts.l = lapply(seq_along(plts.l), function(i) {
    labs = 
      if (names(plts.l)[i] %in% c("gs", "E")) {
        c("(a)", "(b)", "(c)")
      } else if  (names(plts.l)[i] %in% c("A", "Dleaf")) {
        c("(d)", "(e)", "(f)")
      } else if  (names(plts.l)[i]  == "P") {
        c("(g)", "(h)", "(i)")
      }
    ylimits = c(layer_scales(plts.l[[i]][[1]])$y$get_limits(),
                layer_scales(plts.l[[i]][[2]])$y$get_limits(),
                layer_scales(plts.l[[i]][[3]])$y$get_limits())
    ymin = min(ylimits)
    ymax = max(ylimits)
    
    plts.l[[i]][[1]] = plts.l[[i]][[1]] + 
      ylim(c(ymin, ymax))  +
      ggtitle(labs[1]) +
      theme(plot.title = element_text(size = 12, face = "bold"))
    
    plts.l[[i]][[2]] = plts.l[[i]][[2]] + 
      ylim(c(ymin, ymax))  +
      ggtitle(labs[2]) +
      theme(plot.title = element_text(size = 12, face = "bold"))
    
    plts.l[[i]][[3]] = plts.l[[i]][[3]] + 
      ylim(c(ymin, ymax))  +
      ggtitle(labs[3]) +
      theme(plot.title = element_text(size = 12, face = "bold"))
    
    if (names(plts.l)[i] %in% c("gs", "E")) {
      plts.l[[i]][[1]] = plts.l[[i]][[1]]+
        labs(title = expression(psi[s]*" = -0.5 MPa"),
             subtitle = expression(bold("(a)"))) +
        theme(plot.subtitle = element_text(size = 12),
              plot.title = element_text(size = 12))
      
      plts.l[[i]][[2]] = plts.l[[i]][[2]] +
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
    
    return(plts.l[[i]])
  })
  # Align the leftmost panels
  aligned_plots = cowplot::align_plots(plts.l[[1]][[1]] +
                                         theme(legend.position = "none",
                                               axis.title.x = element_blank()), 
                                       plts.l[[2]][[1]] +
                                         theme(legend.position = "none",
                                               axis.title.x = element_blank()),
                                       plts.l[[3]][[1]] +
                                         theme(legend.position = "none",
                                               axis.title.x = element_blank()), 
                                       plts.l[[4]][[1]] +
                                         theme(legend.position = "none",
                                               axis.title.x = element_blank()), 
                                       plts.l[[5]][[1]] +
                                         theme(legend.position = "none",
                                               axis.title.x = element_blank()), 
                                       align = 'v', axis = 'l')
  
  comb_plt = 
    if (vars == "gs, A, Pleaf") {
      comb_plt = cowplot::plot_grid(
        #legend, 
        aligned_plots[[1]], 
        plts.l[[1]][[2]] +  theme(legend.position = "none",
                                  axis.title.x = element_blank(),
                                  axis.title.y = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank()), 
        plts.l[[1]][[3]] +  theme(legend.position = "none",
                                  axis.title.x = element_blank(),
                                  axis.title.y = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank()),
        aligned_plots[[3]], 
        plts.l[[3]][[2]] +  theme(legend.position = "none",
                                  axis.title.x = element_blank(),
                                  axis.title.y = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank()), 
        plts.l[[3]][[3]] +  theme(legend.position = "none",
                                  axis.title.x = element_blank(),
                                  axis.title.y = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank()),
        aligned_plots[[2]], 
        plts.l[[2]][[2]] +  theme(legend.position = "none",
                                  axis.title.x = element_blank(),
                                  axis.title.y = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank()), 
        plts.l[[2]][[3]] +  theme(legend.position = "none",
                                  axis.title.x = element_blank(),
                                  axis.title.y = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank()), 
        nrow = 3, rel_widths = c(1, 0.75, 0.75)
      ) 
    } else {
      cowplot::plot_grid(
        #legend, 
        aligned_plots[[4]], 
        plts.l[[4]][[2]] +  theme(legend.position = "none",
                                  axis.title.x = element_blank(),
                                  axis.title.y = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank()), 
        plts.l[[4]][[3]] +  theme(legend.position = "none",
                                  axis.title.x = element_blank(),
                                  axis.title.y = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank()),
        aligned_plots[[5]], 
        plts.l[[5]][[2]] +  theme(legend.position = "none",
                                  axis.title.x = element_blank(),
                                  axis.title.y = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank()), 
        plts.l[[5]][[3]] +  theme(legend.position = "none",
                                  axis.title.x = element_blank(),
                                  axis.title.y = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank()),
        nrow = 2, rel_widths = c(1, 0.75, 0.75)
      ) 
    }
  
  comb_plt = plot_grid(legend, comb_plt, nrow = 2, rel_heights = c(1,25)) +
    theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
  
  return(comb_plt)
}


plot_physio_vs_Tleaf = function(df, all = TRUE, yvar = NULL, xvar = "Tleaf") {
  
  df$Model = factor(df$Model, levels = c("Medlyn", "Sperry", "Sperry + varkmax", "Sperry + CGnet",  "Sperry + CGnet + TC"))
  
  linetype = c(Sperry = "solid", "Sperry + varkmax" = "longdash", "Sperry + CGnet" = "dashed", "Sperry + CGnet + TC" = "dotted",
               Medlyn = "dotdash")
  shapes = c(Sperry = 1, "Sperry + varkmax" = 5, "Sperry + CGnet" = 2, "Sperry + CGnet + TC" = 3,
             Medlyn = 4)
  palette = c(Sperry = "#1E88E5", "Sperry + varkmax" = "#332288", "Sperry + CGnet" = "#e7d87d", "Sperry + CGnet + TC" = "orange",
              Medlyn = "#D81B60")
  
  
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
                         labels = c("USO",
                                    "ProfitMax",
                                    expression("ProfitMax"[k[max](T)]),
                                    expression("ProfitMax"[net]),
                                    expression("ProfitMax"[TC])
                         )) +
      scale_linetype_manual(values = linetype, 
                            labels = c("USO",
                                       "ProfitMax",
                                       expression("ProfitMax"[k[max](T)]),
                                       expression("ProfitMax"[net]),
                                       expression("ProfitMax"[TC])
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
                         labels = c("USO",
                                    "ProfitMax",
                                    expression("ProfitMax"[k[max](T)]),
                                    expression("ProfitMax"[net]),
                                    expression("ProfitMax"[TC])
                         )) +
      scale_linetype_manual(values = linetype, 
                            labels = c("USO",
                                       "ProfitMax",
                                       expression("ProfitMax"[k[max](T)]),
                                       expression("ProfitMax"[net]),
                                       expression("ProfitMax"[TC])
                            )) +
      theme(text = element_text(size = 14),
            axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
            axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
  }
  return(p)
}

## WTC simulations #############################################################

# Plot A or E vs Tleaf
plot_AEvT_WTC = function(
    GAM = NULL,
    df,
    yvar = c("E", "A", "gs", "Dleaf", "Pleaf")
) {
  
  palette = c(Sperry = "#1E88E5", "Sperry + varkmax" = "#332288", "Sperry + CGnet" = "#e7d87d", "Sperry + CGnet + TC" = "orange", 
              Medlyn = "#D81B60", observed = "grey50")
  
  ylabel = 
    if (yvar == "E") {
      expression("E"*" (mmol m"^-2*"s"^-1*")")
    } else if (yvar == "A") {
      expression("A (" * mu * "mol m"^-2*"s"^-1*")")
    } else if (yvar == "gs") {
      expression("g"[s]*" (mol m"^-2*"s"^-1*")")
    } else if (yvar == "Pleaf") {
      expression(psi[leaf]*" (-MPa)")
    } else if (yvar == "Tleaf") {
      expression("T"[leaf]*" (\u00B0C)")
    } else if (yvar == "Dleaf") {
      expression("VPD"[leaf]*" (kPa)")
    } 
  
  plt = if (is.null(GAM)) {
    NULL
  } else {
    plotGAM(GAM, smooth.c = "Tcan") +
      geom_point(data = df, aes(x = Tleaf, y = !!sym(yvar), color = Model), shape = 1, size = .5) +
      theme_classic() +
      scale_color_manual(
        values = palette,
        labels = c("Observations",
                   "USO",
                   "ProfitMax",
                   expression("ProfitMax"[k[max](T)]),
                   expression("ProfitMax"[net]),
                   expression("ProfitMax"[TC])
        )) +
      xlab(expression("T"[leaf]*" (\u00B0C)")) +
      ylab(ylabel) + 
      guides(color = guide_legend(override.aes = list(shape = 19, size = 2)),
             linetype = "none") +
      theme(plot.title = element_blank())
  }
  
  return(plt)
}

# Plot evaluation metrics
plot_evals = function(eval_df,
                      trt = c("control", "heatwave"),
                      metric = c("R2", "MAE"),
                      binned = FALSE) {
  palette = c(Sperry = "#1E88E5", "Sperry + varkmax" = "#332288", 
              "Sperry + CGnet" = "#e7d87d", "Sperry + CGnet + TC" = "orange", 
              Medlyn = "#D81B60")
  
  plt = eval_df %>% 
    filter(treatment == trt) %>% 
    ggplot(aes(x = Model, y = !!sym(metric), fill = Model)) +
    geom_col() +
    scale_fill_manual(values = palette,
                      labels = c("Observations",
                                 "USO",
                                 "ProfitMax",
                                 expression("ProfitMax"[k[max](T)]),
                                 expression("ProfitMax"[net]),
                                 expression("ProfitMax"[TC])
                      )) +
    theme_classic() +
    ylab(metric) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          strip.background = element_blank(), 
          strip.text = element_text(face = "bold"))
  
  plt = if (isTRUE(binned)) {
    plt + facet_wrap(fit_group ~ variable, scales = "free", ncol = 5)
  } else {
    plt + facet_wrap(vars(variable), scales = "free", ncol = 5) 
  }
  
  return(plt)
  
}

## Helpers ---------------------------------------------------------------------
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
