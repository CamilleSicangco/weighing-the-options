# Functions for analysis
# by Camille Sicangco
# Created 30 Sept 2024

# Make predictions with prescribed values of E
prescribedE_pred = function(df, Wind = 5, Wleaf = 0.025, LeafAbs=0.5, 
                            Rd0 = 0.92, TrefR = 25,...)
{
  Tleaf = try(calc_Tleaf(Tair = df$Tair, E = df$E, VPD = df$VPD, PPFD = df$PPFD, Wind = Wind, 
                     Wleaf = Wleaf, LeafAbs = LeafAbs))
  Tleaf[grep("Error", Tleaf)] = NA
  Tleaf = as.numeric(Tleaf)
  Dleaf = try(VPDairToLeaf(Tleaf = Tleaf, Tair = df$Tair, VPD = df$VPD))
  gs = try(calc_gw(E = df$Trans, D_leaf = Dleaf))
  #gs = try(calc_gw_diff(E = df$E, Tleaf = Tleaf, Tair = df$Tair, VPD = df$VPD, 
  #                 PPFD = df$PPFD, Wind = Wind, Wleaf = Wleaf))
  Agr =  mapply(plantecophys::Photosyn, VPD = df$VPD, 
              Ca = 400, PPFD = df$PPFD, Tleaf = df$Tcan, Tair = df$Tair,
              GS = gs,
              Vcmax=100.52,EaV=51780,EdVC=2e5,delsC=640,
              Jmax = 165.53,EaJ=21640,EdVJ=2e5,delsJ=633, Rd0 = 0)
  Rd = Rd0 * exp(0.1178 * (Tleaf - TrefR) - 7.017e-4 * (Tleaf^2 - TrefR^2)) 
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


# Generate predictions for all models
make_pred = function(
    Tair, Ps, VPD, PPFD, 
    P50 = 4.07, P88 = 5.5,
    Wind = 5,
    Wleaf = 0.025,
    LeafAbs=0.5,
    Tcrit = 46.5,
    T50 = 50.4,
    kmax_25 = 1.5,
    Vcmax=100.52,EaV=62307,EdVC=2e5,delsC=639,
    Jmax = 165.53,EaJ=33115,EdVJ=2e5,delsJ=635,
    Rd0 = 0.92,
    constr_Ci = FALSE,
    ...) {
  
  # Hydraulics
  Weibull = fit_Weibull(P50, P88)
  b = Weibull[1,1]
  c = Weibull[1,2]
  Pcrit = calc_Pcrit(b, c)
  P = Ps_to_Pcrit(Ps, Pcrit, pts = 500)
  
  
  # Calculate costs and gains
  cost_gain = calc_costgain(P, b, c, kmax_25 = kmax_25, Tair = Tair, VPD = VPD, 
                            PPFD = PPFD, Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs, 
                            Tcrit = Tcrit, T50 = T50, 
                            Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
                            Jmax = Jmax,EaJ=EaJ,EdVJ=EdVC,delsJ=delsJ,
                            Rd0 = Rd0, constr_Ci = constr_Ci, ...)

  models = c("Sperry", "Sperry + varkmax", "Sperry + CGnet", "Sperry + CGnet + TC")

  
  # Generate model predictions
  preds = sapply(models, 
                 function(x) make_pred_fn(model = x, 
                                          cost_gain = cost_gain, 
                                          P = P, b = b, c = c,
                                          Tair = Tair, Ps = Ps, VPD = VPD, 
                                          PPFD = PPFD, kmax_25 = kmax_25, 
                                          Wind = Wind, Wleaf = Wleaf, 
                                          LeafAbs = LeafAbs,
                                          Tcrit = Tcrit, T50 = T50,
                                          Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
                                          Jmax = Jmax,EaJ=EaJ,EdVJ=EdVC,delsJ=delsJ,
                                          Rd0 = Rd0,
                                          ...)
  )
  preds = data.frame(t(preds))
  rownames(preds) = NULL
  preds = preds %>% mutate(across(Tair:idx, as.numeric))
  return(preds)
}

# Make predictions with Sperry or Sicangco models for a single observation
make_pred_fn = function(P, b, c,
                         cost_gain,
                         Tair, 
                         Ps, 
                         VPD, 
                         PPFD, 
                         Wind, 
                         LeafAbs = 0.5,
                         Wleaf = 0.025,
                         Tcrit = 46.5,
                         T50 = 50.4,
                         P50 = 4.07,
                         P88 = 5.50,
                         kmax_25 = 1.5,
                         Vcmax=100.52,EaV=62307,EdVC=2e5,delsC=639,
                         Jmax = 165.53,EaJ=33115,EdVJ=2e5,delsJ=635,
                         Rd0 = 0.92,
                         model,
                         constr_Ci = FALSE,
                        Patm = 101.325,
                        matrix_gross = FALSE,
                        gs_feedback = TRUE,
                         #constant_kmax = TRUE,
                         ...) {
  
  
  # Solve optimization
  if (model == "Sperry") {
    i = which.max(cost_gain$CG_gross_constkmax - cost_gain$HC_constkmax)
  } else if (model == "Sperry + CGnet") {
    i = if(isTRUE(matrix_gross)) {
      which.max(cost_gain$CG_net_constkmax - cost_gain$HC_constkmax)
    } else {
      which.max(cost_gain$CG_net - cost_gain$HC_varkmax)
      }
  } else if (model == "Sperry + varkmax") {
    i = if(isTRUE(matrix_gross)) {
      which.max(cost_gain$CG_net - cost_gain$HC_varkmax)
    } else {
      which.max(cost_gain$CG_gross_varkmax - cost_gain$HC_varkmax)
    }
  } else if (model == "Sperry + CGnet + TC") {
    i = which.max(cost_gain$CG_net_newJT - (cost_gain$HC_varkmax + cost_gain$TC))
  } else if (model == "Sperry + TC") {
    i = which.max(cost_gain$CG_gross_varkmax - (cost_gain$HC_varkmax + cost_gain$TC))
  } else {
    stop()
  }
  
  # Hydraulics
  Weibull = fit_Weibull(P50, P88)
  b = Weibull[1,1]
  c = Weibull[1,2]
    
    E_vec = if (model %in% c("Sperry", "Sperry + CGnet_constrCi")) {
      trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax = TRUE)
    } else {
      trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax = FALSE)
    }
    
    P = P[i]
    E = E_vec[i]
    Tleaf = calc_Tleaf(Tair = Tair, E = E, VPD = VPD, PPFD = PPFD, Wind = Wind, 
                       Wleaf = Wleaf, LeafAbs = LeafAbs, gs_feedback = gs_feedback)
    Dleaf = VPDairToLeaf(Tleaf = Tleaf, Tair = Tair, VPD = VPD)
    gs = calc_gw(E = E, Tleaf = Tleaf, Tair = Tair, VPD = VPD, 
                 PPFD = PPFD, Wind = Wind, Wleaf = Wleaf)
    
    # Recalculate E with Penman-Monteith
    #E = gs*1000*Dleaf/101.325
    #H2OLV0 <- 2501000
    #H2OMW <- 0.018
    #CPAIR <- 1010
    #AIRMA <- 0.029
    #UMOLPERJ <- 4.57
    #DHEAT <- 2.15e-05
    #Tair_k <- Tair + 273.15
    #LHV <- (H2OLV0 - 2365 * Tair) * H2OMW
    #GAMMA <- CPAIR * AIRMA * Patm * 1000/LHV
    #SLOPE <- (plantecophys::esat(Tair + 0.1, Pa = Patm) - plantecophys::esat(Tair, 
    #                                                                         Pa = Patm))/0.1
    #CMOLAR <- Patm * 1000/(8.314 * Tair_k)
    #Gbhforced <- 0.003 * sqrt(Wind/Wleaf) * CMOLAR
    #GRASHOF <- 1.6e+08 * abs(Tleaf - Tair) * (Wleaf^3)
    #Gbhfree <- 0.5 * DHEAT * (GRASHOF^0.25)/Wleaf * CMOLAR
    #Gbh <- 2 * (Gbhfree + Gbhforced)
    #Rsol <- 2 * PPFD/UMOLPERJ
    #E = (SLOPE * Rsol + Dleaf * 1000 * Gbh * 
    #    CPAIR * AIRMA)/(GAMMA * Gbh / gs + SLOPE)/LHV*1000
    
    # Keep this re-updated Pleaf 
    #Pleaf_idx = which.min(abs(E - E_vec))
    #Pleaf_alt = P[Pleaf_idx]
    new_JT = if (model %in% c("Sperry + TC", "Sperry + CGnet + TC")) {
      TRUE
    } else {
      FALSE
    }
    
    A = calc_A(Tleaf = Tleaf, g_w = gs, VPD = VPD, net = TRUE, netOrig = TRUE,
               PPFD = PPFD, Wind = Wind, 
               Wleaf = Wleaf, LeafAbs = LeafAbs,
               Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
               Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ, new_JT = new_JT,
               ...)
    Ci = calc_Ci(Tleaf = Tleaf, g_w = gs, VPD = VPD, net = TRUE, netOrig = TRUE,
               PPFD = PPFD, Wind = Wind, 
               Wleaf = Wleaf, LeafAbs = LeafAbs,
               Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
               Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ, new_JT = new_JT,
               ...)
    opt = c(model, Tair, Ps, P, E, Tleaf, Dleaf, gs, A, Ci, i/length(E_vec)#, Ci, 
            #i, cost_gain$CG_gross_constkmax[i], cost_gain$HC_constkmax[i]
            )
    
    names(opt) = c("Model", "Tair", "Ps", "P", "E", "Tleaf", "Dleaf", "gs", "A","Ci", "idx"#, 
                   #"Ci", "idx", "CG", "HC"
                   )
    
    #out = list(opt, P, E = E_vec, 
    #           CG = cost_gain$CG_gross_constkmax, HC = cost_gain$HC_constkmax)
    
    return(opt)
}

# Calculate Pleaf for the Medlyn model
calc_Pleaf0 = function(E, Ps, P50 = 4.07, P88 = 5.50) {
  # Hydraulics
  Weibull = fit_Weibull(P50, P88)
  b = Weibull[1,1]
  c = Weibull[1,2]
  Pcrit = calc_Pcrit(b, c)
  P = Ps_to_Pcrit(Ps, Pcrit, pts = 600)
  E_vec = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax = TRUE)
  
  # Get corresponding Pleaf
  i = which.min(abs(E - E_vec))
  Pleaf = P[i]
  
  return(Pleaf)
}

# Custom version of Photosyn function using Heskel et al. 2017 equation for R(T)
Photosyn_custom <- function(VPD=1.5, 
                     Ca=400, 
                     PPFD=1500,
                     Tleaf=25,
                     Patm=101.325,
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
                     whichA=c("Ah","Amin","Ac","Aj"),
                     
                     Tcrit = 46.5,
                     T50 = 50.4,
                     
                     new_JT = TRUE,
                     
                     b_USO = 0.55, # sensitivity of g1 to SWP
                     Ps = 0.5, # soil water potential, -MPa
                     
                     Tair = 25, # deg C
                     Wind = 5, # Wind speed, m/s
                     Wleaf = 0.025 # Leaf width, m
                     ){
  
  
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
    Jmax <- 
      if(isTRUE(new_JT)) {
        Jmax * TJmax_updated(Tleaf, EaJ, delsJ, EdVJ, Tcrit, T50)
      } else {
        Jmax * TJmax(Tleaf, EaJ, delsJ, EdVJ)  
        }
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
    
    g1 = g1 * exp(b_USO * (-Ps + 0.2))
    
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
    
    # Convert stomatal conductance to water vapor -> carbon
    Gsc <- GS / GCtoGW
    
    # Calculate boundary layer conductance
    DHEAT <- 2.15e-05
    Tair_k <- Tair + 273.15
    CMOLAR <- Patm * 1000/(8.314 * Tair_k)
    Gbhforced <- 0.003 * sqrt(Wind/Wleaf) * CMOLAR
    GRASHOF <- 1.6e+08 * abs(Tleaf - Tair) * (Wleaf^3)
    Gbhfree <- 0.5 * DHEAT * (GRASHOF^0.25)/Wleaf * CMOLAR
    Gbh <- 2 * (Gbhfree + Gbhforced)
    Gbw = Gbh / 0.93
    Gbc <- Gbw / 1.37 # Aphalo and Jarvis 1993
    #Gbc = 1E9
    GC = (Gsc * Gbc) / (Gsc + Gbc)
    
    #GC = Gsc
    
    if(GS > 0){
      
      # Solution when Rubisco activity is limiting
      #A <- 1./Gsc
      #B <- (Rd - Vcmax)/Gsc - Ca - Km
      #C <- Vcmax * (Ca - GammaStar) - Rd * (Ca + Km)
      #Ac <- QUADM(A,B,C)
      
      # Photosynthesis when electron transport is limiting
      #B <- (Rd - VJ)/Gsc - Ca - 2*GammaStar
      #C <- VJ * (Ca - GammaStar) - Rd * (Ca + 2*GammaStar)
      #Aj <- QUADM(A,B,C)
      
      # NOTE: the solution above gives net photosynthesis, add Rd
      # to get gross rates (to be consistent with other solutions).
      #Ac <- Ac + Rd
      #Aj <- Aj + Rd
      
      # Solution when Rubisco activity is limiting
      A <- 1./GC
      B <- -((Vcmax + Rd)/GC + Ca + Km)
      C <- Vcmax * (Ca - GammaStar + Rd/GC)
      Ac <- QUADM(A,B,C)
      
      # Photosynthesis when electron transport is limiting
      B <- -((VJ + Rd)/GC + Ca + 2*GammaStar)
      C <- VJ * (Ca - GammaStar + Rd/GC)
      Aj <- QUADM(A,B,C)
      
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
      lesslcp <- Aj <= Rd + 1E-1
      
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


# Function to fit temperature response of Vcmax 
# by Dushan Kumarathunge
# inputs; dat is a dataframe with Vcmax, Tleaf in Kelvin (named as "TsK")
# return; "Arr" for fitted Arrhenius model parameters, "Peak" for fitted Peak Arrhenius parameters
fitpeaked<-function(dat,return=c("Arr","Peak"),start=list(k25=100, Ea=60, delS = 0.64)){
  
  return <- match.arg(return)
  
  try(Vc<-nls(fVc, start=start, data=dat))
  Vc2<-summary(Vc)
  res1<-data.frame(cbind(Vc2$coefficients[1],Vc2$coefficients[2],Vc2$coefficients[3],Vc2$coefficients[4],Vc2$coefficients[5],
                         Vc2$coefficients[6]))
  
  names(res1)[1:6]<-c("Vcmax25","EaV","delsV","Vcmax25.se","EaV.se","delsV.se")
  
  try(Vr<-nls(fVc.a, start = list(k25=100, Ea=40), data = dat))
  Vr2<-summary(Vr)
  res<-data.frame(cbind(Vr2$coefficients[1],Vr2$coefficients[2],Vr2$coefficients[3],Vr2$coefficients[4]))
  
  
  names(res)[1:4]<-c("Vcmax25","EaV","Vcmax25.se","EaV.se")
  
  
  an<-anova(Vr,Vc)
  AIC_Arr<-AIC(Vr)
  AIC_Peak<-AIC(Vc)
  prob_V<-an[[6]][[2]]
  
  
  r1<-cor(fitted(Vc),dat$Vcmax)
  R2_Peak<-r1*r1
  
  r2<-cor(fitted(Vr),dat$Vcmax)
  R2_Arr<-r2*r2
  
  #test for normality of residuals
  rest<-residuals(Vc)
  norm<-shapiro.test(rest)
  s<-norm$statistic
  pvalue<-norm$p.value
  
  topt<-200/(res1[[3]]-0.008314*log(res1[[2]]/(200-res1[[2]])))
  Topt<-topt-273.15
  
  Topt.se <- topt * 
    sqrt(res1[[6]]**2 + 
           (0.008314 * ((res1[[5]] / res1[[2]])**2 + (res1[[5]] / (200 - res1[[2]]))**2))**2) /
    (res1[[3]]-0.008314*log(res1[[2]]/(200-res1[[2]])))
  
  
  param_peak<-cbind(res1,Topt,Topt.se,R2_Arr,R2_Peak,AIC_Arr,AIC_Peak,prob_V)
  names(param_peak)[1:13]<-c("Vcmax25","EaV","delsV","Vcmax25.se","EaV.se","delsV.se","ToptV",
                             "ToptV.se","R2_Arr","R2_Peak","AIC_Arr","AIC_Peak","prob_V")
  
  param_arr<-cbind(res,R2_Arr,R2_Peak,AIC_Arr,AIC_Peak,prob_V)
  names(param_arr)[1:9]<-c("Vcmax25","EaV","Vcmax25.se","EaV.se","R2_Arr",
                           "R2_Peak","AIC_Arr","AIC_Peak","prob_V")
  
  
  if(return == "Peak")return(param_peak)
  if(return == "Arr")return(param_arr)
  
}

# Function to fit temperature response of Jmax 
# by Dushan Kumarathunge
# inputs; dat is a dataframe with Vcmax, Tleaf in Kelvin (named as "TsK")
# return; "Arr" for fitted Arrhenius model parameters, "Peak" for fitted Peak Arrhenius parameters

fitpeakedJ<-function(dat,return=c("Arr","Peak"),start=list(k25=100, Ea=60, delS = 0.64)){
  return <- match.arg(return)
  dat$TsK <- dat$Tleaf + 273.15
  
  try(Vj<-nls(fVJ, start = start, data = dat))
  Vj2<-summary(Vj)
  res1<-Vj2$coefficients[1:6]
  names(res1)[1:6]<-c("Jmax25","Ea","delsJ","Jmax25.se","EaJ.se","delsJ.se")
  
  try(Vr<-nls(fVJ.a, start = list(k25=100, Ea=40), data = dat))
  Vr2<-summary(Vr)
  res<-Vr2$coefficients[1:4]
  names(res)[1:4]<-c("Vcmax25","EaV","Vcmax25.se","EaV.se")
  
  
  an<-anova(Vr,Vj)
  AIC_Arr<-AIC(Vr)
  AIC_Peak<-AIC(Vj)
  prob_J<-an[[6]][[2]]
  
  
  r1<-cor(fitted(Vj),dat$Jmax)
  R2_Peak<-r1*r1
  
  r2<-cor(fitted(Vr),dat$Jmax)
  R2_Arr<-r2*r2
  
  
  
  #test for normality of residuals
  rest<-residuals(Vj)
  norm<-shapiro.test(rest)
  s<-norm$statistic
  pvalue<-norm$p.value
  
  topt<-200/(res1[[3]]-0.008314*log(res1[[2]]/(200-res1[[2]])))
  Topt<-topt-273.15
  Topt.se <- topt * 
    sqrt(res1[[6]]**2 + 
           (0.008314 * ((res1[[5]] / res1[[2]])**2 + (res1[[5]] / (200 - res1[[2]]))**2))**2) /
    (res1[[3]]-0.008314*log(res1[[2]]/(200-res1[[2]])))
  
  param_peak<-c(res1,Topt,Topt.se,R2_Arr,R2_Peak,AIC_Arr,AIC_Peak,prob_J)
  names(param_peak)[1:13]<-c("Jmax25","EaJ","delsJ","Jmax25.se","EaJ.se","delsJ.se","ToptJ","ToptJ.se","R2_Arr","R2_Peak",
                             "AIC_Arr","AIC_Peak","prob_J")
  
  param_arr<-c(res,R2_Arr,R2_Peak,AIC_Arr,AIC_Peak,prob_J)
  names(param_arr)[1:9]<-c("Jmax25","EaJ","Jmax25.se","EaJ.se","R2_Arr",
                           "R2_Peak","AIC_Arr","AIC_Peak","prob_J")
  
  
  if(return == "Peak")return(param_peak)
  if(return == "Arr")return(param_arr)
  
  #return(c(res,r2,s,pvalue))
  #return(list(res))
}

# define peak Arrhenius model for Vcmax
fVc <- as.formula(Vcmax ~ k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK)) * 
                    (1+exp((298.15*delS - 200)/(298.15*0.008314))) / 
                    (1+exp((TsK*delS-200)/(TsK*0.008314))))

# define peak Arrhenius model for Jmax
fVJ <- as.formula(Jmax ~ k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK)) * 
                    (1+exp((298.15*delS - 200)/(298.15*0.008314))) / 
                    (1+exp((TsK*delS-200)/(TsK*0.008314))))

# define standard Arrhenius model for Vcmax

fVc.a<-as.formula(Vcmax ~ k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK)))


# define standard Arrhenius model for Vcmax

fVJ.a<-as.formula(Jmax ~ k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK)))

# Generate predictions with all models
get_predictions = 
  function(
    df, Tcrit = 46.5, T50 = 50.4, P50 = 4.07, P88 = 5.50,
    Wind = 5, Wleaf = 0.025, LeafAbs = 0.5,
    Vcmax=100.52,EaV=62307,EdVC=2e5,delsC=639,
    Jmax = 165.53,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92,
    kmax_25 = 1.5, #net = TRUE, netOrig = TRUE,
    g1 = 2.9, g0=1.e-5, b_USO = 0.55,
    constr_Ci = FALSE, scaling_factor = 2.9,
    ...
  ) {
    out = bind_rows(lapply(1:nrow(df), 
                 function(i) make_pred(
                   Tair = df$Tair[i], Ps = df$Ps[i], VPD = df$VPD[i], PPFD = df$PPFD[i],
                   Tcrit = Tcrit, T50 = T50, P50 = P50, P88 = P88,
                   Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs, 
                   Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC, 
                   Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ, 
                   Rd0 = Rd0, kmax_25 = kmax_25, constr_Ci = constr_Ci,
                   ...)))
    
    # Make Medlyn predictions
    Medlyn_preds = plantecophys::PhotosynEB(Tair=df$Tair, VPD=df$VPD, PPFD = df$PPFD,
                                            Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs, 
                                            g1 = g1,g0=g0,
                                            Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC, 
                                            Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ, Rd0 = Rd0,
                                            b_USO = b_USO, Ps = df$Ps
                                            )
    get_Pleaf_Medlyn_V = Vectorize(get_Pleaf_Medlyn, c("Ps", "Tair", "E"))
    Medlyn_preds$P = get_Pleaf_Medlyn_V(
      Ps = df$Ps, Tair = df$Tair, E = Medlyn_preds$ELEAF, 
      P50 = P50, P88 = P88, kmax_25 = kmax_25, constant_kmax = FALSE)
    
    Medlyn_preds_sim = Medlyn_preds %>%
      rename(E = ELEAF, Dleaf = VPDleaf, gs = GS, A = ALEAF) %>%
      mutate(Model = "Medlyn", .before = 1) %>%
      #mutate(P = NA) %>% 
      select(Model, Tair, E, Tleaf, Dleaf, gs, A, P, Ci
             )
    
    # Combine all predictions
    #out$Ci = as.numeric(out$Ci)
    #Medlyn_preds_sim$Ci = as.numeric(Medlyn_preds_sim$Ci)
    out_all = bind_rows(out, Medlyn_preds_sim)
    
    out_all$E = out_all$E/scaling_factor
    out_all$A = out_all$A/scaling_factor
    out_all$gs = out_all$gs/scaling_factor
    
    return(out_all)
  }

# Estimate Pleaf for the USO model
get_Pleaf_Medlyn = function(Ps, E, Tair, P50, P88, kmax_25, constant_kmax) {
  
  # Hydraulics
  Weibull = fit_Weibull(P50, P88)
  b = Weibull[1,1]
  c = Weibull[1,2]
  Pcrit = calc_Pcrit(b, c)
  
  P = Ps_to_Pcrit(Ps, Pcrit, pts = 500)
  E_vec = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax)
  j = which.min(abs(E - E_vec))
  Pleaf = P[j]
  
  return(Pleaf)
}


# Make predictions with different Tcrit/T50 values
make_pred_Tthresholds = function(df, Wind = 5, Wleaf = 0.025, LeafAbs=0.5,
                                 hold_Tcrit = FALSE,
                                 Thold_val = 50.4,
                                 Tvar_vals = c(46.5, 47.5, 48.5, 49.5),
                                 constr_Ci = FALSE,
                                 ...) {
  if(isTRUE(hold_Tcrit)) {
    Tcrit = Thold_val
    T50s = Tvar_vals
    
    sims.l = lapply(T50s, function(T50,...) {
      preds = get_predictions(df, Tcrit = Tcrit, T50 = T50) %>% 
        filter(Model == "Sperry + CGnet + TC")
    }
    )
    names(sims.l) = paste0("T50_", Tvar_vals)
    preds = bind_rows(sims.l, .id = "ID") %>% 
      mutate(Tcrit = Thold_val,
             T50 = str_sub(ID, 5)) %>% 
      select(!ID)
  } else {
    T50 = Thold_val
    Tcrits = Tvar_vals
    
    sims.l = lapply(Tcrits, function(Tcrit,...) {
      preds = get_predictions(df, Tcrit = Tcrit, T50 = T50) %>% 
        filter(Model == "Sperry + CGnet + TC")
    }
    )
    names(sims.l) = paste0("T50_", Tvar_vals)
    preds = bind_rows(sims.l, .id = "ID") %>% 
      mutate(T50 = Thold_val,
             Tcrit = str_sub(ID, 5)) %>% 
      select(!ID)
  }
  preds = preds %>% mutate(across(Tair:A, as.numeric),
                           across(T50:Tcrit, as.factor))
  return(preds)
}

# From plantecophys
# Jmax temperature response (Arrhenius)
TJmax <- function(Tleaf, EaJ, delsJ, EdVJ){
  J1 <- 1+exp((298.15*delsJ-EdVJ)/8.314/298.15)
  J2 <- 1+exp(((Tleaf + 273.15)*delsJ-EdVJ)/8.314/(Tleaf + 273.15))
  exp(EaJ/8.314*(1/298.15 - 1/(Tleaf + 273.15)))*J1/J2
}


# Calculate Ci
calc_Ci = function (Tair = 25, VPD = 1.5, PPFD = 1000, Patm = 101.325, 
                    E = 2, Wind = 5, Wleaf = 0.025, LeafAbs = 0.5, Ca = 400, 
                    Jmax = 100, Vcmax = 50, net = FALSE, Rd0 = 0.92, TrefR = 25, 
                    netOrig = TRUE, g1 = 2.9, g0 = 0.003, g_w = NULL, Tleaf = NULL, 
                    ...) 
{
  if (is.null(g_w) & is.null(Tleaf)) {
    Tleaf = calc_Tleaf(Tair = Tair, VPD = VPD, PPFD = PPFD, 
                       E = E, Wind = Wind, Patm = Patm, Wleaf = Wleaf, LeafAbs = LeafAbs)
    g_w = calc_gw(E, Tleaf, Patm, Tair, VPD, PPFD, Wind, 
                  Wleaf)
  }
  Dleaf = plantecophys::VPDairToLeaf(Tleaf = Tleaf, Tair = Tair, 
                                     VPD = VPD, Pa = Patm)
  if (net == FALSE) {
    Photosyn_out = mapply(plantecophys::Photosyn, VPD = Dleaf, 
                          Ca = Ca, PPFD = PPFD, Tleaf = Tleaf, Patm = Patm, 
                          GS = g_w, Rd = 0, Jmax = Jmax, Vcmax = Vcmax, g1 = g1, 
                          g0 = g0, ...)
  }
  else {
    if (isTRUE(netOrig)) {
      Photosyn_out = mapply(plantecophys::Photosyn, VPD = Dleaf, 
                            Ca = Ca, PPFD = PPFD, Tleaf = Tleaf, Patm = Patm, 
                            GS = g_w, Jmax = Jmax, Vcmax = Vcmax, g1 = g1, 
                            g0 = g0, ...)
    }
    else {
      Rd = Rd0 * exp(0.1178 * (Tleaf - TrefR) - 7.017e-4 * (Tleaf^2 - TrefR^2)) 
      Photosyn_out = mapply(plantecophys::Photosyn, VPD = Dleaf, 
                            Ca = Ca, PPFD = PPFD, Tleaf = Tleaf, Patm = Patm, 
                            GS = g_w, Rd = 0, Jmax = Jmax, Vcmax = Vcmax, 
                            g1 = g1, g0 = g0, ...)
    }
  }
  Ci = as.numeric(Photosyn_out[1, ])
  return(Ci)
}

# Generate predictions for all models
make_pred_compEB = function(
    Tair, Ps, VPD, PPFD, 
    P50 = 4.07, P88 = 5.5,
    Wind = 5,
    Wleaf = 0.025,
    LeafAbs=0.5,
    Tcrit = 46.5,
    T50 = 50.4,
    kmax_25 = 1.5,
    Vcmax=100.52,EaV=62307,EdVC=2e5,delsC=639,
    Jmax = 165.53,EaJ=33115,EdVJ=2e5,delsJ=635,
    Rd0 = 0.92,
    constr_Ci = FALSE,
    ...) {
  
  # Hydraulics
  Weibull = fit_Weibull(P50, P88)
  b = Weibull[1,1]
  c = Weibull[1,2]
  Pcrit = calc_Pcrit(b, c)
  P = Ps_to_Pcrit(Ps, Pcrit, pts = 500)
  
  
  # Calculate costs and gains
  cost_gain = calc_costgain(P, b, c, kmax_25 = kmax_25, Tair = Tair, VPD = VPD, 
                            PPFD = PPFD, Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs, 
                            Tcrit = Tcrit, T50 = T50, 
                            Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
                            Jmax = Jmax,EaJ=EaJ,EdVJ=EdVC,delsJ=delsJ,
                            Rd0 = Rd0, constr_Ci = constr_Ci, ...)
  
  # Calculate costs and gains
  cost_gain_nogsfeedback = calc_costgain(P, b, c, kmax_25 = kmax_25, Tair = Tair, VPD = VPD, 
                            PPFD = PPFD, Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs, 
                            Tcrit = Tcrit, T50 = T50, 
                            Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
                            Jmax = Jmax,EaJ=EaJ,EdVJ=EdVC,delsJ=delsJ,
                            Rd0 = Rd0, constr_Ci = constr_Ci, 
                            gs_feedback = FALSE, ...)
  
  CG_dfs = list(cost_gain, cost_gain_nogsfeedback)
  
  # Generate model predictions
  preds = sapply(CG_dfs, 
                 function(CG_df) make_pred_fn(model = "Sperry + CGnet", 
                                          cost_gain = CG_df, 
                                          P = P, b = b, c = c,
                                          Tair = Tair, Ps = Ps, VPD = VPD, 
                                          PPFD = PPFD, kmax_25 = kmax_25, 
                                          Wind = Wind, Wleaf = Wleaf, 
                                          LeafAbs = LeafAbs,
                                          Tcrit = Tcrit, T50 = T50,
                                          Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
                                          Jmax = Jmax,EaJ=EaJ,EdVJ=EdVC,delsJ=delsJ,
                                          Rd0 = Rd0,
                                          ...)
  )
  preds = data.frame(t(preds))
  rownames(preds) = NULL
  preds = preds %>% mutate(across(Tair:idx, as.numeric))
  return(preds)
}

# Corrected Tleaf equation for plantecophys LeafEnergyBalance
LeafEnergyBalance_custom <- function(Tleaf = 21.5, Tair = 20, 
                              gs = 0.15,
                              PPFD = 1500, VPD = 2, Patm = 101,
                              Wind = 2, Wleaf = 0.02, 
                              StomatalRatio = 1,   # 2 for amphistomatous
                              LeafAbs = 0.5,   # in shortwave range, much less than PAR
                              returnwhat=c("balance","fluxes")){   
  
  
  returnwhat <- match.arg(returnwhat)
  
  # Constants
  Boltz <- 5.67 * 10^-8     # w M-2 K-4
  Emissivity <- 0.95        # -
  LatEvap <- 2.54           # MJ kg-1
  #CPAIR <- 1010.0           # J kg-1 K-1
  CPAIR <- 29.29            # J mol-1 K-1
  
  H2OLV0 <- 2.501e6         # J kg-1
  H2OMW <- 18e-3            # J kg-1
  AIRMA <- 29.e-3           # mol mass air (kg/mol)
  AIRDENS <- 1.204          # kg m-3
  UMOLPERJ <- 4.57
  DHEAT <- 21.5e-6          # molecular diffusivity for heat
  
  
  
  # Density of dry air
  AIRDENS <- Patm*1000/(287.058 * Tk(Tair))
  
  # Latent heat of water vapour at air temperature (J mol-1)
  LHV <- (H2OLV0 - 2.365E3 * Tair) * H2OMW
  
  # Const s in Penman-Monteith equation  (Pa K-1)
  SLOPE <- (esat(Tair + 0.1) - esat(Tair)) / 0.1
  
  # Radiation conductance (mol m-2 s-1)
  Gradiation <- 4.*Boltz*Tk(Tair)^3 * Emissivity / (CPAIR * AIRMA)
  
  # See Leuning et al (1995) PC&E 18:1183-1200 Appendix E
  # Boundary layer conductance for heat - single sided, forced convection
  CMOLAR <- Patm*1000 / (8.314 * Tk(Tair))   # .Rgas() in package...
  Gbhforced <- 0.003 * sqrt(Wind/Wleaf) * CMOLAR
  
  # Free convection
  GRASHOF <- 1.6E8 * abs(Tleaf-Tair) * (Wleaf^3) # Grashof number
  Gbhfree <- 0.5 * DHEAT * (GRASHOF^0.25) / Wleaf * CMOLAR
  
  # Total conductance to heat (both leaf sides)
  Gbh <- 2*(Gbhfree + Gbhforced)
  
  # Heat and radiative conductance
  Gbhr <- Gbh + 2*Gradiation
  
  # Boundary layer conductance for water (mol m-2 s-1)
  Gbw <- StomatalRatio * 1.075 * Gbh  # Leuning 1995
  gw <- gs*Gbw/(gs + Gbw)
  
  # Longwave radiation
  # (positive flux is heat loss from leaf)
  Rlongup <- Emissivity*Boltz*Tk(Tleaf)^4
  
  # Rnet
  Rsol <- 2*PPFD/UMOLPERJ   # W m-2
  Rnet <- LeafAbs*Rsol - Rlongup   # full
  
  # Isothermal net radiation (Leuning et al. 1995, Appendix)
  ea <- esat(Tair) - 1000*VPD
  ema <- 0.642*(ea/Tk(Tair))^(1/7)
  Rnetiso <- LeafAbs*Rsol - (1 - ema)*Boltz*Tk(Tair)^4 # isothermal net radiation
  
  # Isothermal version of the Penmon-Monteith equation
  GAMMA <- CPAIR*AIRMA*Patm*1000/LHV
  ET <- (1/LHV) * (SLOPE * Rnetiso + 1000*VPD * Gbh * CPAIR * AIRMA) / (SLOPE + GAMMA * Gbhr/gw)
  
  # Latent heat loss
  lambdaET <- LHV * ET
  
  # Heat flux calculated using Gradiation (Leuning 1995, Eq. 11)
  Y <- 1/(1 + Gradiation/Gbh)
  H2 <- Y*(Rnetiso - lambdaET)
  
  # Heat flux calculated from leaf-air T difference.
  # (positive flux is heat loss from leaf)
  #H <- -CPAIR * AIRDENS * (Gbh/CMOLAR) * (Tair - Tleaf)
  H <- (CPAIR *Gbh) * (Tleaf - Tair)
  
  # Leaf-air temperature difference recalculated from energy balance.
  # (same equation as above!)
  #Tleaf2 <- Tair + H2/(CPAIR * AIRDENS * (Gbh/CMOLAR))
  Tleaf2 <- Tair + H2/(CPAIR * Gbh)
  
  # Difference between input Tleaf and calculated, this will be minimized.
  EnergyBal <- Tleaf - Tleaf2
  
  if(returnwhat == "balance")return(EnergyBal)
  
  if(returnwhat == "fluxes"){
    
    l <- data.frame(ELEAFeb=1000*ET, Gradiation=Gradiation, Rsol=Rsol, Rnetiso=Rnetiso, Rlongup=Rlongup, H=H, lambdaET=lambdaET, gw=gw, Gbh=Gbh, H2=H2, Tleaf2=Tleaf2)
    return(l)
  }
}