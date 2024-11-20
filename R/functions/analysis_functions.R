# Functions for analysis
# by Camille Sicangco
# Created 30 Sept 2024

# Make predictions with prescribed values of E
prescribedE_pred = function(df, Wind = 8, Wleaf = 0.01, LeafAbs=0.86, Rd0 = 0.92, TrefR = 25,...)
{
  Tleaf = try(calc_Tleaf(Tair = df$Tair, E = df$E, VPD = df$VPD, PPFD = df$PPFD, Wind = Wind, 
                     Wleaf = Wleaf, LeafAbs = LeafAbs))
  Tleaf[grep("Error", Tleaf)] = NA
  Tleaf = as.numeric(Tleaf)
  Dleaf = try(VPDairToLeaf(Tleaf = Tleaf, Tair = df$Tair, VPD = df$VPD))
  gs = try(calc_gw(E = df$Trans, D_leaf = Dleaf))
  Agr =  mapply(plantecophys::Photosyn, VPD = df$VPD, 
              Ca = 420, PPFD = df$PPFD, Tleaf = df$Tcan,  
              GS = gs,
              Vcmax=34,EaV=51780,EdVC=2e5,delsC=640,
              Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633, Rd0 = 0)
  Rd = Rd0 * exp(0.1012 * (Tleaf - TrefR) - 5e-04 * (Tleaf^2 - 
                                                       TrefR^2))
  A = as.numeric(Agr[2, ]) - Rd
  out = data.frame(Tleaf.pred = Tleaf, Dleaf.pred = Dleaf, gs.pred = gs, A.pred = A,
                   Tleaf.obs = df$Tcan, Dleaf.obs = df$Dleaf, gs.obs = df$gs, A.obs = df$A, E.obs = df$E)
  out = out %>% 
    mutate(Tleaf.diff = Tleaf.pred - Tleaf.obs, 
           Dleaf.diff = Dleaf.pred - Dleaf.obs,
           gs.diff = gs.pred - gs.obs,
           A.diff = A.pred - A.obs)
  return(out)
}

# Make predictions with Sperry and Sicangco models (or intermediate model versions) 
# for multiple observations
make_pred = function(models = c("final", "intermediate"),
                     df, Wind = 8, Wleaf = 0.01, LeafAbs=0.86,
                     Tcrit_c = 42.6, Tcrit_hw = 44.4,
                     T50_c = 48.6, T50_hw = 50.4,...) {
  if (models == "final") {
    model1 = "Sicangco"
    model2 = "Sperry"
  } else if (models == "intermediate") {
    model1 = "Sperry + CGnet"
    model2 = "Sperry + TC"
  }
  # Generate Sicangco (or Sperry + CGnet) model predictions
  preds1 = make_pred_fn2(model = model1, df = df, Wind = Wind, Wleaf = Wleaf,
                         LeafAbs = LeafAbs, Tcrit_c = Tcrit_c, Tcrit_hw = Tcrit_hw,
                         T50_c = T50_c, T50_hw = T50_hw,...)
  
  # Generate Sperry (or Sperry + TC) model predictions
  preds2 = make_pred_fn2(model = model2, df = df, Wind = Wind, Wleaf = Wleaf,
                         LeafAbs = LeafAbs, Tcrit_c = Tcrit_c, Tcrit_hw = Tcrit_hw,
                         T50_c = T50_c, T50_hw = T50_hw,...)
  
  preds = bind_rows(preds1, preds2)
  return(preds)
}

# Generate predictions for multiple observations for one model
make_pred_fn2 = function(
    model = c("Sicangco", "Sperry", "Sperry + CGnet", "Sperry + TC"),
    df, 
    Wind = 8,
    Wleaf = 0.01,
    LeafAbs=0.86,
    Tcrit_c = 42.6,
    Tcrit_hw = 44.4,
    T50_c = 48.6,
    T50_hw = 50.4,
    ...) {
  # Generate model predictions
  preds = mapply(make_pred_fn, Tair = df$Tair, Ps = df$Ps, VPD = df$VPD,
                 PPFD = df$PPFD, kmax_25 = df$kmax, Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
                 model = model, 
                 Tcrit = ifelse(df$HWtrt == "HW", Tcrit_hw, Tcrit_c), 
                 T50 = ifelse(df$HWtrt == "HW", T50_hw, T50_c),
                 ...)
  preds = data.frame(t(preds))
  names(preds) = c("Model", "Tair", "P", "E", "Tleaf", "Dleaf", "gs", "A")
  preds = preds %>% mutate(across(Tair:A, as.numeric))
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
                        kmax_25 = 2,
                        Vcmax=34,EaV=51780,EdVC=2e5,delsC=640,
                        Jmax = 60,EaJ=21640,EdVJ=2e5,delsJ=633,
                        model) {
  

  # Hydraulics
  Weibull = fit_Weibull(P50, P88)
  b = Weibull[1,1]
  c = Weibull[1,2]
  Pcrit = calc_Pcrit(b, c)
  P = Ps_to_Pcrit(Ps, Pcrit)
  E_vec = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax = TRUE)
  

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
  } else if (model == "Sperry + CGnet") {
    which.min(abs(marg_df$CG_net - marg_df$HC))
  } else if (model == "Sperry + TC") {
    which.min(abs(marg_df$CG_gross - marg_df$CC))
  } else {
    stop()
  }
  
  P = marg_df$P[i]
  E = E_vec[i]
  Tleaf = calc_Tleaf(Tair = Tair, E = E, VPD = VPD, PPFD = PPFD, Wind = Wind, 
                     Wleaf = Wleaf, LeafAbs = LeafAbs)
  Dleaf = VPDairToLeaf(Tleaf = Tleaf, Tair = Tair, VPD = VPD)
  gs = calc_gw(E = E, D_leaf = Dleaf)
  #net = if (model == "Sperry") {
  #  FALSE
  #} else if (model == "Sicangco") {
  #  TRUE
  #} 
  A = calc_A(Tleaf = Tleaf, g_w = gs, VPD = VPD, net = TRUE, PPFD = PPFD, Wind = Wind, 
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