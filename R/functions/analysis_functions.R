# Functions for analysis
# by Camille Sicangco
# Created 30 Sept 2024

# Make predictions with Sperry and Sicangco models for multiple observations
make_pred = function(df, Wind = 8, Wleaf = 0.01, LeafAbs=0.86,...) {
  # Generate Sicangco model predictions
  preds_Sicangco = mapply(make_pred_fn, Tair = df$Tair_al, Ps = df$Ps, VPD = df$VPD,
                          PPFD = df$PAR, kmax_25 = df$kmax, Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
                          model = "Sicangco",...)
  preds_Sicangco = data.frame(t(preds_Sicangco))
  names(preds_Sicangco) = c("Model", "Tair", "P", "E", "Tleaf", "Dleaf", "gs", "A")
  preds_Sicangco = preds_Sicangco %>% mutate(across(Tair:A, as.numeric))
  
  # Generate Sperry model predictions
  preds_Sperry = mapply(make_pred_fn, Tair = df$Tair_al, Ps = df$Ps, VPD = df$VPD,
                        PPFD = df$PAR,kmax_25 = df$kmax,  Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
                        model = "Sperry",...)
  preds_Sperry = data.frame(t(preds_Sperry))
  names(preds_Sperry) = c("Model", "Tair", "P", "E", "Tleaf", "Dleaf", "gs", "A")
  preds_Sperry = preds_Sperry %>% mutate(across(Tair:A, as.numeric))
  
  preds = bind_rows(preds_Sicangco, preds_Sperry)
  return(preds)
}

# Make predictions with Sperry or Sicangco models for a single observation
make_pred_fn = function(Tair, 
                     Ps, 
                     VPD, 
                     PPFD, 
                     Wind, 
                     LeafAbs = 0.5,
                     Wleaf = 0.025,
                     Tcrit = 48.574,
                     T50 = 50.17,
                     P50 = 4.31,
                     P88 = 6.64,
                     kmax_25 = 10.116744,
                     Vcmax=34,EaV=51780,EdVC=2e5,delsC=640,
                     Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633,
                     model) {
  

  # Hydraulics
  Weibull = fit_Weibull(P50, P88)
  b = Weibull[1,1]
  c = Weibull[1,2]
  Pcrit = calc_Pcrit(b, c)
  P = Ps_to_Pcrit(Ps, Pcrit)
  E_vec = trans_from_vc(P, kmax_25, Tair, b, c)
  

  # Compute costs and gains
  cost_gain = calc_costgain(P, b, c, kmax_25 = kmax_25, Wind = Wind, 
                            Wleaf = Wleaf, LeafAbs = LeafAbs,
                            Tair = Tair, PPFD = PPFD, 
                            VPD = VPD, Tcrit = Tcrit, T50 = T50, 
                            Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
                            Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
                            constant_kmax = TRUE)
  
  # Compute marginal costs and gains
  marg_df = marginal_gaincost(cost_gain)
  
  # Solve optimization
  i = if (model == "Sperry") {
    which.min(abs(marg_df$CG_gross - marg_df$HC))
  } else if (model == "Sicangco") {
    which.min(abs(marg_df$CG_net - marg_df$CC))
  } else {
    stop()
  }
  
  P = marg_df$P[i]
  E = E_vec[i]
  Tleaf = calc_Tleaf(Tair = Tair, E = E, VPD = VPD, PPFD = PPFD, Wind = Wind, 
                     Wleaf = Wleaf, LeafAbs = LeafAbs)
  Dleaf = VPDairToLeaf(Tleaf = Tleaf, Tair = Tair, VPD = VPD)
  gs = calc_gw(E = E, D_leaf = Dleaf)
  net = if (model == "Sperry") {
    FALSE
  } else if (model == "Sicangco") {
    TRUE
  } 
  A = calc_A(Tair = Tair, E = E, VPD = VPD, net = net, PPFD = PPFD, Wind = Wind, 
             Wleaf = Wleaf, LeafAbs = LeafAbs,
             Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
             Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ)
  
  out = c(model, Tair, P, E, Tleaf, Dleaf, gs, A)
  names(out) = c("model", "Tair", "P", "E", "Tleaf", "Dleaf", "gs", "A")
  return(out)
}


# Calculate all marginal gains and costs
marginal_gaincost = function(df)
{
  df = df %>% pivot_wider(names_from = ID, values_from = cost_gain) %>% 
    mutate(CC = HC + TC)
  HC = marginal_GainCost(df$P, df$HC)
  CC = marginal_GainCost(df$P, df$CC)
  CG_gross = marginal_GainCost(df$P, df$CG_gross)
  CG_net = marginal_GainCost(df$P, df$CG_net)
  
  df_out = data.frame(P = df$P[-1], CG_gross, CG_net, HC, CC)
  return(df_out)
}

# Calculate marginal cost or gain 
marginal_GainCost = function(x = P, y)
{
  dy_dx = diff(y) / diff(x)
  return(dy_dx)
}