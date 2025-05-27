# Sensitivity analysis of TC model to Tcrit and T50
# Created 27 May 2025

# Make predictions with different Tcrit/T50 values
make_pred_Tthresholds = function(df, Wind = 8, Wleaf = 0.01, LeafAbs=0.5,
                     hold_Tcrit = FALSE,
                     Thold_val = 49.6,
                     Tvar_vals = c(43.4, 47.4, 48, 48.6),
                     ...) {
  if(isTRUE(hold_Tcrit)) {
    Tcrit = Thold_val
    T50s = Tvar_vals
    
    sims.l = lapply(T50s, function(T50,...) make_pred_fn2(model = "Sperry + CGnet + TC", df = df, 
                               Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs, 
                               Tcrit_hw = Thold_val, T50_hw = T50,
                               HWtrt = "HW",...))
    names(sims.l) = paste0("T50_", Tvar_vals)
    preds = bind_rows(sims.l, .id = "ID") %>% 
      mutate(Tcrit = Thold_val,
             T50 = str_sub(ID, 5)) %>% 
      select(!ID)
  } else {
    T50 = Thold_val
    Tcrits = Tvar_vals
    
    sims.l = lapply(Tcrits, function(x,...) make_pred_fn2(model = "Sperry + CGnet + TC", df = df, 
                                                              Tcrit_hw = x, T50_hw = T50,
                                                              HWtrt = "HW",...))
    names(sims.l) = paste0("T50_", Tvar_vals)
    preds = bind_rows(sims.l, .id = "ID") %>% 
      mutate(T50 = Thold_val,
             Tcrit = str_sub(ID, 5)) %>% 
      select(!ID)
  }
  
  return(preds)
}

test = make_pred_fn2(model = "Sperry + CGnet + TC", df = Tair_sim.df, 
                     Wind = 8, Wleaf = 0.02, LeafAbs = 0.5, 
                     Tcrit_hw = 48.6, T50_hw = 49.6,
                     HWtrt = "HW", kmax_25 = 1.5, constant_kmax = TRUE)
test %>% 
  ggplot(aes(x = Tleaf, y = gs)) + 
  geom_point() + 
  theme_classic() +
  ylab(expression("g"[s]*" (mol m"^-2*"s"^-1*")")) +
  xlab(expression("T"[leaf]*" (\u00B0C)"))

Tair_vec = seq(30,50, by = 1)
Tair_sim.df = data.frame(Tair = Tair_vec, PPFD = 1500, VPD = 1.5,
                         Ps = 0.5, HWtrt = "HW")
preds = make_pred_Tthresholds(Tair_sim.df, 
                              Wind = 8, Wleaf = 0.02, LeafAbs = 0.5, 
                              kmax_25 = 1.5, constant_kmax = TRUE)

preds %>% 
  ggplot(aes(x = Tair, y = gs, linetype = Tcrit)) + 
  geom_line(linewidth = 1) + 
  theme_classic() +
  ylab(expression("g"[s]*" (mol m"^-2*"s"^-1*")")) +
  xlab(expression("T"[air]*" (\u00B0C)"))
