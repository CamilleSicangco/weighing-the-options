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
Photosyn_custom_oldTresp <- function(VPD=1.5, 
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
                                     whichA=c("Ah","Amin","Ac","Aj"),
                                     
                                     Tcrit = 43.4,
                                     T50 = 48.6){
  
  
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
    Vcmax <- #ifelse(Rd != 0,
      #      Vcmax * TVcmax_updated(Tleaf,EaV, delsC, EdVC, Tcrit),
      Vcmax * TVcmax(Tleaf,EaV, delsC, EdVC)
    #)
    Jmax <- #ifelse(Rd != 0,
      #      Jmax * TJmax_updated(Tleaf, EaJ, delsJ, EdVJ, Tcrit),
      Jmax * TJmax(Tleaf, EaJ, delsJ, EdVJ)
    #)
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
    #if(!inputCi){
    # lesslcp <- vector("logical", length(Aj))
    #lesslcp <- Aj <= Rd + 1E-09
    
    #if(length(Ca) == 1)Ca <- rep(Ca, length(CIJ))
    #if(length(GammaStar) == 1)GammaStar <- rep(GammaStar, length(CIJ))
    #if(length(VJ) == 1)VJ <- rep(VJ, length(CIJ))
    
    #CIJ[lesslcp] <- Ca[lesslcp]
    #Aj[lesslcp] <- VJ[lesslcp] * (CIJ[lesslcp] - GammaStar[lesslcp]) / 
    # (CIJ[lesslcp] + 2*GammaStar[lesslcp])
    
    #Ci <- ifelse(Aj < Ac, CIJ, CIC)
    
    #  }
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



# Modified cost gain function
calc_costgain_corr = function (P = NULL, b = -2.5, c = 2, Amax_gross = NULL, Amax_net = NULL, 
                               kmax_25 = 4, Tair = 25, VPD = 1.5, ratiocrit = 0.05, PPFD = 1000, 
                               Patm = 101.325, Wind = 2, Wleaf = 0.025, LeafAbs = 0.5, Tcrit = 43.4, 
                               T50 = 49.6, Ca = 400, Rd0 = 0.92, TrefR = 25, 
                               Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
                               Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635,
                               constr_Ci = FALSE, ...) 
{
  if(is.null(P)) {
    return(NULL)
  }
  
  # Hydraulic costs
  HC_constkmax = hydraulic_cost(P, b, c, kmax_25, Tair, constant_kmax = TRUE)
  HC_varkmax = hydraulic_cost(P, b, c, kmax_25, Tair, constant_kmax = FALSE)
  
  # Thermal cost
  TC = thermal_cost(P, b, c, kmax_25, Tair, VPD, PPFD, Patm, 
                    Wind, Wleaf, LeafAbs, Tcrit, T50, constant_kmax = FALSE)
  
  # Carbon gains
  if(isFALSE(constr_Ci)) {
    CG_gross_constkmax = C_gain_corr(
      P, b, c, kmax_25 = kmax_25,
      Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
      Tair = Tair, PPFD = PPFD, 
      VPD = VPD, Tcrit = Tcrit, T50 = T50, 
      Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
      Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
      constant_kmax = TRUE, net = FALSE, new_JT = FALSE, 
      ...)
    CG_gross_varkmax = C_gain_corr(
      P, b, c, kmax_25 = kmax_25,
      Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
      Tair = Tair, PPFD = PPFD, 
      VPD = VPD, Tcrit = Tcrit, T50 = T50, 
      Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
      Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
      constant_kmax = FALSE, net = FALSE, new_JT = FALSE, 
      ...)
    CG_net = C_gain_corr(
      P, b, c, kmax_25 = kmax_25,
      Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
      Tair = Tair, PPFD = PPFD, 
      VPD = VPD, Tcrit = Tcrit, T50 = T50, 
      Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
      Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
      constant_kmax = FALSE, net = TRUE, netOrig = TRUE, new_JT = FALSE,  
      ...)
    CG_net_newJT = C_gain_corr(
      P, b, c, kmax_25 = kmax_25,
      Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
      Tair = Tair, PPFD = PPFD, 
      VPD = VPD, Tcrit = Tcrit, T50 = T50, 
      Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
      Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
      constant_kmax = FALSE, net = TRUE, netOrig = TRUE, new_JT = TRUE,  
      ...)
  } else {
    CG_gross_constkmax = C_gain_alt(
      P, b, c, Amax = NULL, kmax_25 = kmax_25,
      Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
      Tair = Tair, PPFD = PPFD, 
      VPD = VPD, Tcrit = Tcrit, T50 = T50, 
      Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
      Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
      constant_kmax = TRUE, new_JT = FALSE, net = FALSE,
      ...)
    CG_gross_varkmax = C_gain_alt(
      P, b, c, Amax = NULL, kmax_25 = kmax_25,
      Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
      Tair = Tair, PPFD = PPFD, 
      VPD = VPD, Tcrit = Tcrit, T50 = T50, 
      Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
      Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
      constant_kmax = FALSE, new_JT = FALSE, net = FALSE,
      ...)
    CG_net = C_gain_alt(
      P, b, c, Amax = NULL, kmax_25 = kmax_25,
      Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
      Tair = Tair, PPFD = PPFD, 
      VPD = VPD, Tcrit = Tcrit, T50 = T50, 
      Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
      Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
      constant_kmax = FALSE, new_JT = FALSE, net = TRUE,
      ...)
    CG_net_newJT = C_gain_alt(
      P, b, c, Amax = NULL, kmax_25 = kmax_25,
      Wind = Wind, Wleaf = Wleaf, LeafAbs = LeafAbs,
      Tair = Tair, PPFD = PPFD, 
      VPD = VPD, Tcrit = Tcrit, T50 = T50, 
      Vcmax=Vcmax,EaV=EaV,EdVC=EdVC,delsC=delsC,
      Jmax = Jmax,EaJ=EaJ,EdVJ=EdVJ,delsJ=delsJ,
      constant_kmax = FALSE, new_JT = TRUE, net = TRUE,
      ...)
  }
  
  df = data.frame(P, 
                  HC_constkmax, HC_varkmax, 
                  TC,
                  CG_gross_constkmax, CG_gross_varkmax,
                  CG_net, CG_net_newJT)
  return(df)
}

# Updated V, J T-responses
TVcmax_updated <- function(Tleaf, EaV, delsC, EdVC, Tcrit = 43.4){
  
  if(EdVC > 0){
    V1 <- 1+exp((delsC*(25 + 273.15)-EdVC)/(.Rgas()*(25 + 273.15)))
    V2 <- 1+exp((delsC*(Tleaf+273.15)-EdVC)/(.Rgas()*(Tk(Tleaf))))
    f <- V1/V2
  } else f <- 1
  
  V = ifelse(Tleaf < Tcrit, 
             exp((Tleaf-25)*EaV/(.Rgas()*Tk(Tleaf)*Tk(25))) * f,
             0)
  return(V)
}

TJmax_updated <- function(Tleaf, EaJ, delsJ, EdVJ, Tcrit = 43.4, T50 = 49.6){
  J1 <- 1+exp((298.15*delsJ-EdVJ)/.Rgas()/298.15)
  J2 <- 1+exp((Tk(Tleaf)*delsJ-EdVJ)/.Rgas()/Tk(Tleaf))
  
  r = 2/(T50 - Tcrit)
  functionality = 1 - 1/(1 + exp(-r * (Tleaf - T50)))
  
  J = functionality * (exp(EaJ/.Rgas()*(1/298.15 - 1/Tk(Tleaf)))*J1/J2)
  
  return(J)
}

# Original V, J T-responses
# Vcmax temperature response (Arrhenius)
TVcmax <- function(Tleaf, EaV, delsC, EdVC){
  if(EdVC > 0){
    V1 <- 1+exp((delsC*(25 + 273.15)-EdVC)/(.Rgas()*(25 + 273.15)))
    V2 <- 1+exp((delsC*(Tleaf+273.15)-EdVC)/(.Rgas()*(Tk(Tleaf))))
    f <- V1/V2
  } else f <- 1
  
  exp((Tleaf-25)*EaV/(.Rgas()*Tk(Tleaf)*Tk(25))) * f
}

# Jmax temperature response (Arrhenius)
TJmax <- function(Tleaf, EaJ, delsJ, EdVJ){
  J1 <- 1+exp((298.15*delsJ-EdVJ)/.Rgas()/298.15)
  J2 <- 1+exp((Tk(Tleaf)*delsJ-EdVJ)/.Rgas()/Tk(Tleaf))
  exp(EaJ/.Rgas()*(1/298.15 - 1/Tk(Tleaf)))*J1/J2
}

# Correction to C_gain
C_gain_corr = function (P, b = -2.5, c = 2, Amax = NULL, kmax_25 = 4, Tair = 25, 
                        VPD = 1.5, PPFD = 1000, Patm = 101.325, Wind = 2, Wleaf = 0.01, 
                        LeafAbs = 0.5, Ca = 400, Jmax = 100, Vcmax = 50, constant_kmax = FALSE, 
                        net = FALSE, Rd0 = 0.92, TrefR = 25, netOrig = TRUE, ...) 
{
  E = trans_from_vc(P, kmax_25, Tair, b, c, constant_kmax)
  A = calc_A_corr(Tair, VPD, PPFD, Patm, E, Wind, Wleaf, LeafAbs, 
                  Ca, Jmax, Vcmax, net, Rd0, TrefR, netOrig, ...)
  E = ifelse(E == 0, NA, E)
  A = ifelse(E == 0, NA, A)
  Amax = if (is.null(Amax)) {
    max(abs(A[!is.na(A)]))
  } else {
    Amax
  }
  
  gain = 
    if (!is.na(Amax) & Amax != 0) {
      A/Amax
    } else {
      rep(0, length = length(A))
    }
  return(gain)
}


calc_A_corr = function (Tair = 25, VPD = 1.5, PPFD = 1000, Patm = 101.325, 
                        E = 2, Wind = 2, Wleaf = 0.01, LeafAbs = 0.5, Ca = 400, 
                        Jmax = 100, Vcmax = 50, net = FALSE, Rd0 = 0.92, TrefR = 25, 
                        netOrig = TRUE, g1 = 2.9, g0 = 0.003, g_w = NULL, Tleaf = NULL, 
                        new_JT = FALSE,
                        ...) 
{
  if (is.null(g_w) & is.null(Tleaf)) {
    Tleaf = calc_Tleaf(Tair = Tair, VPD = VPD, PPFD = PPFD, 
                       E = E, Wind = Wind, Patm = Patm, Wleaf = Wleaf, 
                       LeafAbs = LeafAbs)
    g_w = calc_gw(E, Tleaf, Patm, Tair, VPD, PPFD, Wind, 
                  Wleaf)
  }
  if (net == FALSE) {
    Photosyn_out = mapply(plantecophys::Photosyn, VPD = VPD, 
                          Ca = Ca, PPFD = PPFD, Tleaf = Tleaf, Patm = Patm, 
                          GS = g_w, Rd = 0, Jmax = Jmax, Vcmax = Vcmax, g1 = g1, 
                          g0 = g0, new_JT = new_JT, ...)
    Anet = as.numeric(Photosyn_out[2, ])
    Rd = as.numeric(Photosyn_out[8, ])
    A = Anet + Rd
  }
  else {
    if (isTRUE(netOrig)) {
      Photosyn_out = mapply(plantecophys::Photosyn, VPD = VPD, 
                            Ca = Ca, PPFD = PPFD, Tleaf = Tleaf, Patm = Patm, 
                            GS = g_w, Jmax = Jmax, Vcmax = Vcmax, g1 = g1, 
                            g0 = g0, new_JT = new_JT, ...)
      A = as.numeric(Photosyn_out[2, ])
    }
    else {
      Rd = Rd0 * exp(0.1012 * (Tleaf - TrefR) - 5e-04 * 
                       (Tleaf^2 - TrefR^2))
      Photosyn_out = mapply(plantecophys::Photosyn, VPD = VPD, 
                            Ca = Ca, PPFD = PPFD, Tleaf = Tleaf, Patm = Patm, 
                            GS = g_w, Rd = 0, Jmax = Jmax, Vcmax = Vcmax, 
                            g1 = g1, g0 = g0, new_JT = new_JT, ...)
      A = as.numeric(Photosyn_out[2, ]) - Rd
    }
  }
  return(A)
}