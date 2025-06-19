# Sensitivity analysis of TC model to Tcrit and T50
# Created 27 May 2025

# Make predictions with different Tcrit/T50 values
make_pred_Tthresholds = function(df, Wind = 8, Wleaf = 0.01, LeafAbs=0.5,
                     hold_Tcrit = FALSE,
                     Thold_val = 49.6,
                     Tvar_vals = c(43.4, 45.5, 47.5, 48.5),
                     ...) {
  if(isTRUE(hold_Tcrit)) {
    Tcrit = Thold_val
    T50s = Tvar_vals
    
    sims.l = lapply(T50s, function(x,...) make_pred_fn2(model = "Sperry + CGnet + TC", df = df,
                               Tcrit = Thold_val, T50 = x,
                               ...))
    names(sims.l) = paste0("T50_", Tvar_vals)
    preds = bind_rows(sims.l, .id = "ID") %>% 
      mutate(Tcrit = Thold_val,
             T50 = str_sub(ID, 5)) %>% 
      select(!ID)
  } else {
    T50 = Thold_val
    Tcrits = Tvar_vals
    
    sims.l = lapply(Tcrits, function(x,...) make_pred_fn2(model = "Sperry + CGnet + TC", df = df, 
                                                              Tcrit = x, T50 = T50,
                                                              ...))
    names(sims.l) = paste0("T50_", Tvar_vals)
    preds = bind_rows(sims.l, .id = "ID") %>% 
      mutate(T50 = Thold_val,
             Tcrit = str_sub(ID, 5)) %>% 
      select(!ID)
  }
  
  return(preds)
}

Tair_vec = seq(30,60, by = 1)
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = 1500, VPD = 1.5,
                         Ps = 0.5)
#EaV=51780,EdVC=2e5,delsC=640,
#Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633, Rd = 0)
test1 = make_pred_fn2(model = "Sperry + CGnet + TC", df = Tair_sim.df, 
                     Wind = 8, Wleaf = 0.02, LeafAbs = 0.5, 
                     Tcrit = 43.6, T50 = 49.6,
                     kmax_25 = 0.5, constant_kmax = FALSE)#,Jmax = 100, Vcmax = 50,
                     #EaV=51780, EaJ=21640, delsC=640, delsJ=633)
test2 = make_pred_fn2(model = "Sperry + CGnet + TC", df = Tair_sim.df, 
                     Wind = 8, Wleaf = 0.02, LeafAbs = 0.5, 
                     Tcrit = 48.6, T50 = 49.6,
                     kmax_25 = 0.5, constant_kmax = FALSE#,Jmax = 100, Vcmax = 50
                     )
#test1 %>% 
bind_rows(test1, test2, .id = "test") %>% 
  ggplot(aes(x = Tair, y = gs, color = test
             )) + 
  geom_point(shape = 1) + 
  theme_classic() +
  ylab(expression("g"[s]*" (mol m"^-2*"s"^-1*")")) +
  xlab(expression("T"[leaf]*" (\u00B0C)"))


preds_varTcrit = make_pred_Tthresholds(Tair_sim.df, 
                              Wind = 8, Wleaf = 0.02, LeafAbs = 0.5, 
                              kmax_25 = 0.5, constant_kmax = FALSE)

palette = scales::seq_gradient_pal("slategray1","darkslateblue", "Lab")(seq(0,1,length.out=4))
Tcrit.plt = preds_varTcrit %>% 
  ggplot(aes(x = Tleaf, y = gs, linetype = Tcrit, color = Tcrit)) + 
  geom_line(linewidth = 1) + 
  theme_classic() +
  ylab(expression("g"[s]*" (mol m"^-2*"s"^-1*")")) +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  guides(linetype = guide_legend(title = expression("T"[crit]*" (\u00B0C)")),
         color = guide_legend(title = expression("T"[crit]*" (\u00B0C)"))) +
  scale_colour_manual(values = palette) +
  theme(axis.title = element_text(size = 14))
Tcrit.plt
ggsave("figs/SA_Tcrit.pdf", Tcrit.plt, height = 7, width = 11)

