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
                     df, Wind = 8, Wleaf = 0.01, LeafAbs=0.5,
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

# Custom version of Photosyn function using Heskel et al. 2017 equation for R(T)
Photosyn_custom <- function(VPD=1.5, 
                     Ca=400, 
                     PPFD=1500,
                     Tleaf=25,
                     Patm=100,
                     RH=NULL,
                     
                     gsmodel=c("BBOpti","BBLeuning","BallBerry","BBdefine"),
                     g1=4,
                     g0=0, 
                     gk=0.5,
                     vpdmin=0.5,
                     D0=5,
                     GS=NULL,
                     BBmult=NULL,
                     
                     alpha=0.24, 
                     theta=0.85, 
                     Jmax=100, 
                     Vcmax=50, 
                     gmeso=NULL,
                     TPU=1000,
                     alphag=0,
                     
                     Rd0 = 1.115,
                     Q10 = 1.92,
                     Rd=NULL,
                     TrefR = 25,
                     Rdayfrac = 1.0,
                     
                     EaV = 58550,
                     EdVC = 200000,
                     delsC = 629.26,
                     
                     EaJ = 29680,
                     EdVJ = 200000,
                     delsJ = 631.88,
                     
                     GammaStar = NULL,
                     Km = NULL,
                     
                     Ci = NULL,
                     Tcorrect=TRUE,  
                     returnParsOnly=FALSE,
                     whichA=c("Ah","Amin","Ac","Aj")){
  
  
  whichA <- match.arg(whichA)
  gsmodel <- match.arg(gsmodel)
  if(gsmodel == "BBdefine" && is.null(BBmult)){
    Stop("When defining your own BB multiplier, set BBmult.")
  }
  inputCi <- !is.null(Ci)
  inputGS <- !is.null(GS)
  
  if(inputCi & inputGS)Stop("Cannot provide both Ci and GS.")
  
  if(is.null(TPU))TPU <- 1000
  
  if(is.null(VPD) && !is.null(RH)){
    VPD <- RHtoVPD(RH, Tleaf)
  } 
  if(is.null(VPD) && is.null(RH)){
    Stop("Need one of VPD, RH.")
  }
  
  #---- Constants; hard-wired parameters.
  Rgas <- .Rgas()
  GCtoGW <- 1.57     # conversion from conductance to CO2 to H2O
  
  
  #---- Do all calculations that can be vectorized
  
  # g1 and g0 are input ALWAYS IN UNITS OF H20
  # G0 must be converted to CO2 (but not G1, see below)
  g0 <- g0/GCtoGW
  
  # Leaf respiration
  if(is.null(Rd)){
    Rd = Rd0 * exp(0.1178 * (Tleaf - TrefR) - 7.017e-4 * (Tleaf^2 - TrefR^2))
  }
  
  # CO2 compensation point in absence of photorespiration
  if(is.null(GammaStar))GammaStar <- TGammaStar(Tleaf,Patm)
  
  # Michaelis-Menten coefficient
  if(is.null(Km))Km <- TKm(Tleaf,Patm)
  
  #-- Vcmax, Jmax T responses
  if(Tcorrect){
    Vcmax <- Vcmax * TVcmax(Tleaf,EaV, delsC, EdVC)
    Jmax <- Jmax * TJmax(Tleaf, EaJ, delsJ, EdVJ)
  }
  
  # Electron transport rate
  J <- Jfun(PPFD, alpha, Jmax, theta)
  VJ <- J/4
  
  #--- Stop here if only the parameters are required
  if(returnParsOnly){
    return(list(Vcmax=Vcmax, Jmax=Jmax, Km=Km, GammaStar=GammaStar, VJ=VJ))
  }
  
  # Medlyn et al. 2011 model gs/A. NOTE: 1.6 not here because we need GCO2!
  if(gsmodel == "BBOpti"){
    vpduse <- VPD
    vpduse[vpduse < vpdmin] <- vpdmin
    GSDIVA <- (1 + g1/(vpduse^(1-gk)))/Ca
  }
  
  # Leuning 1995 model, without gamma (CO2 compensation point)
  if(gsmodel == "BBLeuning"){
    GSDIVA <- g1 / Ca / (1 + VPD/D0)
    GSDIVA <- GSDIVA / GCtoGW   # convert to conductance to CO2
  }
  
  # Original Ball&Berry 1987 model.
  if(gsmodel == "BallBerry"){
    if(is.null(RH))RH <- VPDtoRH(VPD, Tleaf)
    RH <- RH / 100
    GSDIVA <- g1 * RH / Ca
    GSDIVA <- GSDIVA / GCtoGW   # convert to conductance to CO2
  } 
  
  # Multiplier is user-defined.
  if(gsmodel == "BBdefine"){
    GSDIVA <- BBmult / GCtoGW
  }
  
  if(inputGS){
    
    GC <- GS / GCtoGW
    
    if(GS > 0){
      
      # Solution when Rubisco activity is limiting
      A <- 1./GC
      B <- (Rd - Vcmax)/GC - Ca - Km
      C <- Vcmax * (Ca - GammaStar) - Rd * (Ca + Km)
      Ac <- QUADM(A,B,C)
      
      # Photosynthesis when electron transport is limiting
      B <- (Rd - VJ)/GC - Ca - 2*GammaStar
      C <- VJ * (Ca - GammaStar) - Rd * (Ca + 2*GammaStar)
      Aj <- QUADM(A,B,C)
      
      # NOTE: the solution above gives net photosynthesis, add Rd
      # to get gross rates (to be consistent with other solutions).
      Ac <- Ac + Rd
      Aj <- Aj + Rd
    } else {
      Ac <- Aj <- 0
    }
    
  } else {
    
    # If CI not provided, calculate from intersection between supply and demand
    if(!inputCi){
      
      #--- non-vectorized workhorse
      getCI <- function(VJ,GSDIVA,PPFD,VPD,Ca,Tleaf,vpdmin,g0,Rd,
                        Vcmax,Jmax,Km,GammaStar){
        
        if(identical(PPFD, 0) | identical(VJ, 0)){
          vec <- c(Ca,Ca)
          return(vec)
        }
        
        # Taken from MAESTRA.
        # Following calculations are used for both BB & BBL models.
        # Solution when Rubisco activity is limiting
        A <- g0 + GSDIVA * (Vcmax - Rd)
        B <- (1. - Ca*GSDIVA) * (Vcmax - Rd) + g0 * 
          (Km - Ca)- GSDIVA * (Vcmax*GammaStar + Km*Rd)
        C <- -(1. - Ca*GSDIVA) * (Vcmax*GammaStar + Km*Rd) - g0*Km*Ca
        
        CIC <- QUADP(A,B,C)
        
        # Solution when electron transport rate is limiting
        A <- g0 + GSDIVA * (VJ - Rd)
        B <- (1 - Ca*GSDIVA) * (VJ - Rd) + g0 * (2.*GammaStar - Ca)- 
          GSDIVA * (VJ*GammaStar + 2.*GammaStar*Rd)
        C <- -(1 - Ca*GSDIVA) * GammaStar * (VJ + 2*Rd) - 
          g0*2*GammaStar*Ca
        
        CIJ <- QUADP(A,B,C)
        return(c(CIJ,CIC))
      }
      
      # get Ci
      x <- mapply(getCI, 
                  VJ=VJ,
                  GSDIVA = GSDIVA,
                  PPFD=PPFD,
                  VPD=VPD,
                  Ca=Ca,
                  Tleaf=Tleaf,
                  vpdmin=vpdmin,
                  g0=g0,
                  Rd=Rd,
                  Vcmax=Vcmax,
                  Jmax=Jmax,
                  Km=Km,
                  GammaStar=GammaStar)
      
      CIJ <- x[1,]
      CIC <- x[2,]
    } else {
      
      # Rare case where one Ci is provided, and multiple Tleaf (Jena bug).
      if(length(Ci) == 1){
        Ci <- rep(Ci, length(Km))
      }
      
      # Ci provided (A-Ci function mode)
      CIJ <- Ci
      
      if(length(GammaStar) > 1){
        CIJ[CIJ <= GammaStar] <- GammaStar[CIJ < GammaStar]
      } else {
        CIJ[CIJ <= GammaStar] <- GammaStar
      }
      
      CIC <- Ci
      
    }
    
    # Photosynthetic rates, without or with mesophyll limitation
    if(is.null(gmeso) || gmeso < 0){
      # Get photosynthetic rate  
      Ac <- Vcmax*(CIC - GammaStar)/(CIC + Km)
      Aj <- VJ * (CIJ - GammaStar)/(CIJ + 2*GammaStar)
      
    } else {
      # Ethier and Livingston (2004) (Equation 10).
      A <- -1/gmeso
      BC <- (Vcmax - Rd)/gmeso + CIC + Km
      CC <- Rd*(CIC+Km)-Vcmax*(CIC-GammaStar)
      Ac <- mapply(QUADP, A=A,B=BC,C=CC)
      
      BJ <- (VJ - Rd)/gmeso + CIC + 2.0*GammaStar
      CJ <- Rd*(CIC+2.0*GammaStar) - VJ*(CIC - GammaStar)
      Aj <- mapply(QUADP, A=A,B=BJ,C=CJ)
      
      Ac <- Ac + Rd
      Aj <- Aj + Rd
      
    }
    
    
    # When below light-compensation points, assume Ci=Ca.
    if(!inputCi){
      lesslcp <- vector("logical", length(Aj))
      lesslcp <- Aj <= Rd + 1E-09
      
      if(length(Ca) == 1)Ca <- rep(Ca, length(CIJ))
      if(length(GammaStar) == 1)GammaStar <- rep(GammaStar, length(CIJ))
      if(length(VJ) == 1)VJ <- rep(VJ, length(CIJ))
      
      CIJ[lesslcp] <- Ca[lesslcp]
      Aj[lesslcp] <- VJ[lesslcp] * (CIJ[lesslcp] - GammaStar[lesslcp]) / 
        (CIJ[lesslcp] + 2*GammaStar[lesslcp])
      
      Ci <- ifelse(Aj < Ac, CIJ, CIC)
      
    }
  }
  
  # Limitation by triose-phosphate utilization
  if(!is.null(Ci)){
    Ap <- 3 * TPU * (Ci - GammaStar)/(Ci - (1 + 3*alphag)*GammaStar)
    Ap[Ci < 400] <- 1000  # avoid nonsense
  } else {
    Ap <- 1000  # This is when inputGS = TRUE; 
  }  
  
  
  # Hyperbolic minimum.
  Am <- -mapply(QUADP, A = 1 - 1E-04, B = Ac+Aj, C = Ac*Aj)
  
  # Another hyperbolic minimum with the transition to TPU
  tpulim <- any(Ap < Am)
  if(!is.na(tpulim) && tpulim){
    Am <- -mapply(QUADP, A = 1 - 1E-07, B = Am+Ap, C = Am*Ap)
  }
  
  # Net photosynthesis
  Am <- Am - Rd
  
  # Calculate conductance to CO2
  if(!inputCi && !inputGS){
    if(whichA == "Ah")GS <- g0 + GSDIVA*Am
    if(whichA == "Aj")GS <- g0 + GSDIVA*(Aj-Rd)
    if(whichA == "Ac")GS <- g0 + GSDIVA*(Ac-Rd)
  } 
  if(inputCi) {
    if(whichA == "Ah")GS <- Am/(Ca - Ci)
    if(whichA == "Aj")GS <- (Aj-Rd)/(Ca - Ci)
    if(whichA == "Ac")GS <- (Ac-Rd)/(Ca - Ci)
  }
  
  # Extra step here; GS can be negative
  GS[GS < g0] <- g0
  
  # Output conductance to H2O
  if(!inputGS){
    GS <- GS*GCtoGW
  }
  
  # Calculate Ci if GS was provided as input.
  if(inputGS){
    Ci <- Ca - Am/GC
    
    # For zero GC:
    Ci[!is.finite(Ci)] <- Ca
    # Stomata fully shut; Ci is not really Ca but that's 
    # how we like to think about it.
  }
  
  # Chloroplastic CO2 concentration
  if(!is.null(gmeso)){
    Cc <- Ci - Am/gmeso
  } else {
    Cc <- Ci
  }
  
  # Transpiration rate assuming perfect coupling.
  # Output units are mmol m-2 s-1
  E <- 1000*GS*VPD/Patm
  
  df <- data.frame( Ci=Ci,
                    ALEAF=Am,
                    GS=GS,
                    ELEAF=E,
                    Ac=Ac,
                    Aj=Aj,
                    Ap=Ap,
                    Rd=Rd,
                    VPD=VPD,
                    Tleaf=Tleaf,
                    Ca=Ca,
                    Cc=Cc,
                    PPFD=PPFD,
                    Patm=Patm)
  
  return(df)
}

# Custom version of Photosyn function using Heskel et al. 2017 equation for R(T)
# and with alternative implementation of Ac/Aj calculation 
Photosyn_custom2 <- function(VPD=1.5, 
                            Ca=400, 
                            PPFD=1500,
                            Tleaf=25,
                            Patm=100,
                            RH=NULL,
                            
                            gsmodel=c("BBOpti","BBLeuning","BallBerry","BBdefine"),
                            g1=4,
                            g0=0, 
                            gk=0.5,
                            vpdmin=0.5,
                            D0=5,
                            GS=NULL,
                            BBmult=NULL,
                            
                            alpha=0.24, 
                            theta=0.85, 
                            Jmax=100, 
                            Vcmax=50, 
                            gmeso=NULL,
                            TPU=1000,
                            alphag=0,
                            
                            Rd0 = 1.115,
                            Q10 = 1.92,
                            Rd=NULL,
                            TrefR = 25,
                            Rdayfrac = 1.0,
                            
                            EaV = 58550,
                            EdVC = 200000,
                            delsC = 629.26,
                            
                            EaJ = 29680,
                            EdVJ = 200000,
                            delsJ = 631.88,
                            
                            GammaStar = NULL,
                            Km = NULL,
                            
                            Ci = NULL,
                            Tcorrect=TRUE,  
                            returnParsOnly=FALSE,
                            whichA=c("Ah","Amin","Ac","Aj")){
  
  
  whichA <- match.arg(whichA)
  gsmodel <- match.arg(gsmodel)
  if(gsmodel == "BBdefine" && is.null(BBmult)){
    Stop("When defining your own BB multiplier, set BBmult.")
  }
  inputCi <- !is.null(Ci)
  inputGS <- !is.null(GS)
  
  if(inputCi & inputGS)Stop("Cannot provide both Ci and GS.")
  
  if(is.null(TPU))TPU <- 1000
  
  if(is.null(VPD) && !is.null(RH)){
    VPD <- RHtoVPD(RH, Tleaf)
  } 
  if(is.null(VPD) && is.null(RH)){
    Stop("Need one of VPD, RH.")
  }
  
  #---- Constants; hard-wired parameters.
  Rgas <- .Rgas()
  GCtoGW <- 1.57     # conversion from conductance to CO2 to H2O
  
  
  #---- Do all calculations that can be vectorized
  
  # g1 and g0 are input ALWAYS IN UNITS OF H20
  # G0 must be converted to CO2 (but not G1, see below)
  g0 <- g0/GCtoGW
  
  # Leaf respiration
  if(is.null(Rd)){
    Rd = Rd0 * exp(0.1178 * (Tleaf - TrefR) - 7.017e-4 * (Tleaf^2 - TrefR^2))
  }
  
  # CO2 compensation point in absence of photorespiration
  if(is.null(GammaStar))GammaStar <- TGammaStar(Tleaf,Patm)
  
  # Michaelis-Menten coefficient
  if(is.null(Km))Km <- TKm(Tleaf,Patm)
  
  #-- Vcmax, Jmax T responses
  if(Tcorrect){
    Vcmax <- Vcmax * TVcmax(Tleaf,EaV, delsC, EdVC)
    Jmax <- Jmax * TJmax(Tleaf, EaJ, delsJ, EdVJ)
  }
  
  # Electron transport rate
  J <- Jfun(PPFD, alpha, Jmax, theta)
  VJ <- J/4
  
  #--- Stop here if only the parameters are required
  if(returnParsOnly){
    return(list(Vcmax=Vcmax, Jmax=Jmax, Km=Km, GammaStar=GammaStar, VJ=VJ))
  }
  
  # Medlyn et al. 2011 model gs/A. NOTE: 1.6 not here because we need GCO2!
  if(gsmodel == "BBOpti"){
    vpduse <- VPD
    vpduse[vpduse < vpdmin] <- vpdmin
    GSDIVA <- (1 + g1/(vpduse^(1-gk)))/Ca
  }
  
  # Leuning 1995 model, without gamma (CO2 compensation point)
  if(gsmodel == "BBLeuning"){
    GSDIVA <- g1 / Ca / (1 + VPD/D0)
    GSDIVA <- GSDIVA / GCtoGW   # convert to conductance to CO2
  }
  
  # Original Ball&Berry 1987 model.
  if(gsmodel == "BallBerry"){
    if(is.null(RH))RH <- VPDtoRH(VPD, Tleaf)
    RH <- RH / 100
    GSDIVA <- g1 * RH / Ca
    GSDIVA <- GSDIVA / GCtoGW   # convert to conductance to CO2
  } 
  
  # Multiplier is user-defined.
  if(gsmodel == "BBdefine"){
    GSDIVA <- BBmult / GCtoGW
  }
  
  if(inputGS){
    
    GC <- GS / GCtoGW
    
    if(GS > 0){
      # If CI not provided, calculate from intersection between supply and demand
      if(!inputCi){
        
        #--- non-vectorized workhorse
        getCI <- function(VJ,GSDIVA,PPFD,VPD,Ca,Tleaf,vpdmin,g0,Rd,
                          Vcmax,Jmax,Km,GammaStar){
          
          if(identical(PPFD, 0) | identical(VJ, 0)){
            vec <- c(Ca,Ca)
            return(vec)
          }
          
          
          # Taken from MAESTRA.
          # Following calculations are used for both BB & BBL models.
          # Solution when Rubisco activity is limiting
          A <- g0 + GSDIVA * (Vcmax - Rd)
          B <- (1. - Ca*GSDIVA) * (Vcmax - Rd) + g0 * 
            (Km - Ca)- GSDIVA * (Vcmax*GammaStar + Km*Rd)
          C <- -(1. - Ca*GSDIVA) * (Vcmax*GammaStar + Km*Rd) - g0*Km*Ca
          
          CIC <- QUADP(A,B,C)
          
          # Solution when electron transport rate is limiting
          A <- g0 + GSDIVA * (VJ - Rd)
          B <- (1 - Ca*GSDIVA) * (VJ - Rd) + g0 * (2.*GammaStar - Ca)- 
            GSDIVA * (VJ*GammaStar + 2.*GammaStar*Rd)
          C <- -(1 - Ca*GSDIVA) * GammaStar * (VJ + 2*Rd) - 
            g0*2*GammaStar*Ca
          
          CIJ <- QUADP(A,B,C)

          return(c(CIJ,CIC))
        }
        
        # get Ci
        x <- mapply(getCI, 
                    VJ=VJ,
                    GSDIVA = GSDIVA,
                    PPFD=PPFD,
                    VPD=VPD,
                    Ca=Ca,
                    Tleaf=Tleaf,
                    vpdmin=vpdmin,
                    g0=g0,
                    Rd=Rd,
                    Vcmax=Vcmax,
                    Jmax=Jmax,
                    Km=Km,
                    GammaStar=GammaStar)
        
        CIJ <- x[1,]
        CIC <- x[2,]
      } else {
        
        # Rare case where one Ci is provided, and multiple Tleaf (Jena bug).
        if(length(Ci) == 1){
          Ci <- rep(Ci, length(Km))
        }
        
        # Ci provided (A-Ci function mode)
        CIJ <- Ci
        
        if(length(GammaStar) > 1){
          CIJ[CIJ <= GammaStar] <- GammaStar[CIJ < GammaStar]
        } else {
          CIJ[CIJ <= GammaStar] <- GammaStar
        }
        
        CIC <- Ci
        
      }
      
      # Solution when Rubisco activity is limiting
      A <- 1.
      B <- Rd - Vcmax - (CIC + Km) * GC / (Patm * 1.e3)
      C <- (Vcmax * (CIC - GammaStar) - (CIC + Km) * Rd) * GC / (Patm * 1.e3)
      Ac <- QUADM(A,B,C)
      
      # Photosynthesis when electron transport is limiting
      B <- Rd - VJ - (CIC + 2. * GammaStar) * GC / (Patm * 1.e3)
      C <- (VJ * (CIC - GammaStar) - (CIC + 2. * GammaStar) * Rd) * GC / (Patm * 1.e3)
      Aj <- QUADM(A,B,C)
      
      # Solution when Rubisco activity is limiting
      #A <- 1./GC
      #B <- (Rd - Vcmax)/GC - Ca - Km
      #C <- Vcmax * (Ca - GammaStar) - Rd * (Ca + Km)
      #Ac <- QUADM(A,B,C)
      
      # Photosynthesis when electron transport is limiting
      #B <- (Rd - VJ)/GC - Ca - 2*GammaStar
      #C <- VJ * (Ca - GammaStar) - Rd * (Ca + 2*GammaStar)
      #Aj <- QUADM(A,B,C)
      
      # NOTE: the solution above gives net photosynthesis, add Rd
      # to get gross rates (to be consistent with other solutions).
      Ac <- Ac + Rd
      Aj <- Aj + Rd
    } else {
      Ac <- Aj <- 0
      CIC = NA
      CIJ = NA
    }
    
  } else {
    
    # If CI not provided, calculate from intersection between supply and demand
    if(!inputCi){
      
      #--- non-vectorized workhorse
      getCI <- function(VJ,GSDIVA,PPFD,VPD,Ca,Tleaf,vpdmin,g0,Rd,
                        Vcmax,Jmax,Km,GammaStar){
        
        if(identical(PPFD, 0) | identical(VJ, 0)){
          vec <- c(Ca,Ca)
          return(vec)
        }
        
        
        # Modification using Medlyn model
        # Solution when Rubisco activity is limiting
        A <- (Rd - Vcmax)*(GSDIVA - g0) - g0
        B <- (Ca - Km)*(g0 - Rd * (GSDIVA - g0)) + 
          Vcmax * (GSDIVA - g0) * (Ca + GammaStar) - Vcmax
        C <- Km * Ca * (g0 - Rd * (GSDIVA - g0)) - 
          Vcmax * Ca * GammaStar * (GSDIVA - g0) + Vcmax * GammaStar
        
        CIC <- QUADP(A,B,C)
        
        # Solution when electron transport rate is limiting
        A <- (Rd - VJ)*(GSDIVA - g0) - g0
        B <- (Ca - 2. * GammaStar)*(g0 - Rd * (GSDIVA - g0)) + 
          VJ * (GSDIVA - g0) * (Ca + GammaStar) - VJ
        C <- 2. * GammaStar * Ca * (g0 - Rd * (GSDIVA - g0)) - 
          VJ * Ca * GammaStar * (GSDIVA - g0) + VJ * GammaStar
        
        CIJ <- QUADP(A,B,C)
        
        # Taken from MAESTRA.
        # Following calculations are used for both BB & BBL models.
        # Solution when Rubisco activity is limiting
        #A <- g0 + GSDIVA * (Vcmax - Rd)
        #B <- (1. - Ca*GSDIVA) * (Vcmax - Rd) + g0 * 
        #  (Km - Ca)- GSDIVA * (Vcmax*GammaStar + Km*Rd)
        #C <- -(1. - Ca*GSDIVA) * (Vcmax*GammaStar + Km*Rd) - g0*Km*Ca
        
        #CIC <- QUADP(A,B,C)
        
        # Solution when electron transport rate is limiting
        #A <- g0 + GSDIVA * (VJ - Rd)
        #B <- (1 - Ca*GSDIVA) * (VJ - Rd) + g0 * (2.*GammaStar - Ca)- 
        #  GSDIVA * (VJ*GammaStar + 2.*GammaStar*Rd)
        #C <- -(1 - Ca*GSDIVA) * GammaStar * (VJ + 2*Rd) - 
        #  g0*2*GammaStar*Ca
        
        #CIJ <- QUADP(A,B,C)
        return(c(CIJ,CIC))
      }
      
      # get Ci
      x <- mapply(getCI, 
                  VJ=VJ,
                  GSDIVA = GSDIVA,
                  PPFD=PPFD,
                  VPD=VPD,
                  Ca=Ca,
                  Tleaf=Tleaf,
                  vpdmin=vpdmin,
                  g0=g0,
                  Rd=Rd,
                  Vcmax=Vcmax,
                  Jmax=Jmax,
                  Km=Km,
                  GammaStar=GammaStar)
      
      CIJ <- x[1,]
      CIC <- x[2,]
    } else {
      
      # Rare case where one Ci is provided, and multiple Tleaf (Jena bug).
      if(length(Ci) == 1){
        Ci <- rep(Ci, length(Km))
      }
      
      # Ci provided (A-Ci function mode)
      CIJ <- Ci
      
      if(length(GammaStar) > 1){
        CIJ[CIJ <= GammaStar] <- GammaStar[CIJ < GammaStar]
      } else {
        CIJ[CIJ <= GammaStar] <- GammaStar
      }
      
      CIC <- Ci
      
    }
    
    # Photosynthetic rates, without or with mesophyll limitation
    if(is.null(gmeso) || gmeso < 0){
      # Get photosynthetic rate  
      Ac <- Vcmax*(CIC - GammaStar)/(CIC + Km)
      Aj <- VJ * (CIJ - GammaStar)/(CIJ + 2*GammaStar)
      #A <- 1
      #BC <- Rd - Vcmax - (CIC + Km) * GSDIVA / (Patm * 1e3)
      #CC <- (Vcmax * (CIC - GammaStar) - (CIC + Km) * Rd) * GSDIVA / (Patm * 1e3)
      #Ac <- mapply(QUADP, A=A,B=BC,C=CC)
      
      #BJ <- Rd - Jmax - (CIJ + 2. * GammaStar) * GSDIVA / (Patm * 1e3)
      #CJ <- (Jmax * (CIJ - GammaStar) - (CIJ + 2. * GammaStar) *
      #         Rd) * GSDIVA / (Patm * 1e3)
      #Aj <- mapply(QUADP, A=A,B=BJ,C=CJ)
      
      #Ac <- Ac + Rd
      #Aj <- Aj + Rd
      
    } else {
      # Ethier and Livingston (2004) (Equation 10).
      A <- -1/gmeso
      BC <- (Vcmax - Rd)/gmeso + CIC + Km
      CC <- Rd*(CIC+Km)-Vcmax*(CIC-GammaStar)
      Ac <- mapply(QUADP, A=A,B=BC,C=CC)
      
      BJ <- (VJ - Rd)/gmeso + CIC + 2.0*GammaStar
      CJ <- Rd*(CIC+2.0*GammaStar) - VJ*(CIC - GammaStar)
      Aj <- mapply(QUADP, A=A,B=BJ,C=CJ)
      
      Ac <- Ac + Rd
      Aj <- Aj + Rd
      
    }
    
    
    # When below light-compensation points, assume Ci=Ca.
    if(!inputCi){
      lesslcp <- vector("logical", length(Aj))
      lesslcp <- Aj <= Rd + 1E-09
      
      if(length(Ca) == 1)Ca <- rep(Ca, length(CIJ))
      if(length(GammaStar) == 1)GammaStar <- rep(GammaStar, length(CIJ))
      if(length(VJ) == 1)VJ <- rep(VJ, length(CIJ))
      
      CIJ[lesslcp] <- Ca[lesslcp]
      Aj[lesslcp] <- VJ[lesslcp] * (CIJ[lesslcp] - GammaStar[lesslcp]) / 
        (CIJ[lesslcp] + 2*GammaStar[lesslcp])
      
      Ci <- ifelse(Aj < Ac, CIJ, CIC)
      
    }
  }
  
  # Limitation by triose-phosphate utilization
  if(!is.null(Ci)){
    Ap <- 3 * TPU * (Ci - GammaStar)/(Ci - (1 + 3*alphag)*GammaStar)
    Ap[Ci < 400] <- 1000  # avoid nonsense
  } else {
    Ap <- 1000  # This is when inputGS = TRUE; 
  }  
  
  
  # Hyperbolic minimum.
  Am <- -mapply(QUADP, A = 1 - 1E-04, B = Ac+Aj, C = Ac*Aj)
  
  # Another hyperbolic minimum with the transition to TPU
  tpulim <- any(Ap < Am)
  if(!is.na(tpulim) && tpulim){
    Am <- -mapply(QUADP, A = 1 - 1E-07, B = Am+Ap, C = Am*Ap)
  }
  
  # Net photosynthesis
  Am <- Am - Rd
  
  # Calculate conductance to CO2
  if(!inputCi && !inputGS){
    if(whichA == "Ah")GS <- g0 + GSDIVA*Am
    if(whichA == "Aj")GS <- g0 + GSDIVA*(Aj-Rd)
    if(whichA == "Ac")GS <- g0 + GSDIVA*(Ac-Rd)
  } 
  if(inputCi) {
    if(whichA == "Ah")GS <- Am/(Ca - Ci)
    if(whichA == "Aj")GS <- (Aj-Rd)/(Ca - Ci)
    if(whichA == "Ac")GS <- (Ac-Rd)/(Ca - Ci)
  }
  
  # Extra step here; GS can be negative
  GS[GS < g0] <- g0
  
  # Output conductance to H2O
  if(!inputGS){
    GS <- GS*GCtoGW
  }
  
  # Calculate Ci if GS was provided as input.
  if(inputGS){
    Ci <- Ca - Am/GC
    
    # For zero GC:
    Ci[!is.finite(Ci)] <- Ca
    # Stomata fully shut; Ci is not really Ca but that's 
    # how we like to think about it.
  }
  
  # Chloroplastic CO2 concentration
  if(!is.null(gmeso)){
    Cc <- Ci - Am/gmeso
  } else {
    Cc <- Ci
  }
  
  # Transpiration rate assuming perfect coupling.
  # Output units are mmol m-2 s-1
  E <- 1000*GS*VPD/Patm
  
  df <- data.frame( Ci=Ci,
                    ALEAF=Am,
                    GS=GS,
                    ELEAF=E,
                    Ac=Ac,
                    Aj=Aj,
                    Ap=Ap,
                    Rd=Rd,
                    VPD=VPD,
                    Tleaf=Tleaf,
                    Ca=Ca,
                    Cc=Cc,
                    PPFD=PPFD,
                    Patm=Patm,
                    CIC=CIC,
                    CIJ=CIJ
                    )
  
  return(df)
}

# Modified cost gain function
calc_costgain_netorig = function (P = NULL, b = -2.5, c = 2, Amax_gross = NULL, Amax_net = NULL, 
          kmax_25 = 4, Tair = 25, VPD = 1.5, ratiocrit = 0.05, PPFD = 1000, 
          Patm = 101.325, Wind = 2, Wleaf = 0.01, LeafAbs = 0.5, Tcrit = 50, 
          T50 = 51, Ca = 420, Jmax = 100, Vcmax = 50, constant_kmax = FALSE, 
          Rd0 = 0.92, TrefR = 25, ...) 
{
  if(is.null(P)) {
    return(NULL)
  }
  HC = hydraulic_cost(P, b, c, kmax_25, Tair, ratiocrit, constant_kmax)
  TC = thermal_cost(P, b, c, kmax_25, Tair, VPD, PPFD, Patm, 
                    Wind, Wleaf, LeafAbs, Tcrit, T50, constant_kmax)
  CG_net = C_gain(P, b, c, Amax_net, kmax_25, Tair, VPD, PPFD, 
                  Patm, Wind, Wleaf, LeafAbs, Ca, Jmax, Vcmax, constant_kmax, 
                  net = TRUE, Rd0, TrefR, netOrig = TRUE, ...)
  CG_gross = C_gain(P, b, c, Amax_gross, kmax_25, Tair, VPD, 
                    PPFD, Patm, Wind, Wleaf, LeafAbs, Ca, Jmax, Vcmax, constant_kmax, 
                    net = FALSE, Rd0, TrefR, ...)
  cost_gain = c(HC, TC, CG_net, CG_gross)
  ID = c(rep("HC", length(HC)), rep("TC", length(TC)), rep("CG_net", 
                                                           length(CG_net)), rep("CG_gross", length(CG_gross)))
  df = data.frame(P, ID, cost_gain)
  return(df)
}
