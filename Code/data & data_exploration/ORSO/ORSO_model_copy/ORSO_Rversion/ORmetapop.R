# Oregon sea otter meta-population model
# EXPLANATION
# Model simulates growth of recently introduced sea otter population(s) in Oregon.
# Dynamics modeled for P coastal blocks (sub-populations) linked by dispersal. 
# Stage-structured matrix model includes density-dependence, environmental stochasticity,
# demographic stochasticity, immigration and dispersal, range expansion via diffusion,
# and habitat-based variation in local equilibrium densities.
# User provides model parameters, which are informed by analyses of data from
# sea otter populations elsewhere (mostly California)
# Spatial scenarios (for each, assume 180 otters total):
# ) Rogue Reef/Crook Point, 2) Humbug/Cape Blanco, 3) Cape Arago/Coos Bay, 4) Yaquina Bay/Cape Foulweather.
#  =      S3 & S4,              S5 & S6                 S9 & S10 & SE2        C7 & C8 & CE4
rm(list = ls())
# Set User Parameters  ---------------------------------------------------------
reps = 10000         # Number replications for population sims (should use at least 500)
StartYear = 2026    # Initial year
Nyrs = 100           # Number of years to project population dynamics
Mltplscn = 0        # set to 0 for single scenario, 1 for range of initN, 2 for range of areas
Initseg = c("N3","C7","S6")  # List of coastal segments for introduction: e.g. c(1) = Block 1 only
                             # Pacific city "N3", New Port "C7", Port Ordford "S6"
Initpop = c(30,30,30)  # Number of animals introduced to each block)
AddOt = c(0,0,0)      # Additional otters per year introduced after initial intro: 0 = none
AddYrs = 0         # Number of years for supplementary otter additions (e.g. re-habs) added  
PpnFem = .65        # Proportion of females in introduced pop
PpnAd = .25         # Proportion of adults in introduced pop
Ephase = 10         # Avg Years for population to become "established" 
AddMortE = .14      # Excess annual mortality during establishment phase (0.15 for San Nic)
ProbMv = .6         # Probability of adult otters dispersing from initial translocation site 
EstryMvAdj = .5     # Proportional decrease in dispersal prob for estuaries relative to outer coast
AgeMvAdj = .5       # Proportional decrease in dispersal probability for subadults relative to adults 
Mvmort = .75        # proportion of initially dispersing animals that die/disappear from Oregon
sig = 0.1           # Environmental stochasticity (std. deviation in log-lambda)
rmax = 0.18         # Maximum rate of growth: default = log(1.2), or 20% per year
theta = .9          # theta parameter for theta-logistic (1 = Ricker model, >1 = delayed DD) 
V_sp = 2.5          # Population front asymptotic wavespeed, km/yr, minimum  
# ~~~~~~END User parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
plotmap = 1
plotrngexp = 1
savetable = 1
Historical_comp = 0
Scenario_name = c("Norm_Succes") # set to c("") for no scenario
# Scenario_name = c("SC3_") # set to c("") for no scenario
#
# Load necessary libraries -------------------------------------------------
# NOTE: Ensure following packages are installed
library(gtools)
library(mvtnorm)
library(boot)
library(ggplot2)
# library(rgdal) # deprecated
library(terra)
library(dplyr) 
library(ggrepel)
# library(ggsn) # deprecated (but not needed?  for map scale bar)
library(reshape2)
require(parallel)
require(doParallel)
require(foreach)
require(parallelly)
library(abind)
require(openxlsx)
forloop_packages<-c('gtools','mvtnorm')
#
#
# Load files ----------------------------------------------------------------
data = read.csv("./data/OR_pops.csv", header = TRUE)  # Data for coastal blocks
Demdat = read.csv("./data/RandDem.csv", header = TRUE)    # Stochastic vital rates
Distmat =  as.matrix(read.csv("./data/Distmat.csv", header = FALSE)) # Inter-blk LCP distances
# Probabilities of dispersal from each block for each age/sex class
DispP = read.csv("./data/ORDispProb.csv", header = TRUE) 
# Inter-pop movement matrices: pairwise prob of dispersal based on 
#  pairwise LCP distances and dispersal kernels for each age/sex class
destJF = read.csv("./data/DispMatJF.csv", header = FALSE) 
destAF = read.csv("./data/DispMatAF.csv", header = FALSE);
destJM = read.csv("./data/DispMatJM.csv", header = FALSE);
destAM = read.csv("./data/DispMatAM.csv", header = FALSE);
# Load GIS data for plotting results (NAD_1983_Albers)
load("./data/GISdata.rdata")
#
# Process data ---------------------------------------------------------------
#
if(Scenario_name == "" | is.na(Scenario_name)){Scenario_name = "Sim"}
P = nrow(Distmat)  # number blocks (or sub-populations)
datasort = order(as.numeric(data$Ycoord))
data$Sortord = rep(0,P)
data$Sortord[datasort] = seq(1,P)
Nfin = numeric()
Nfin_lo = numeric()
Nfin_hi = numeric()
if(Mltplscn==0){
  Nscn = 1
}else if (Mltplscn==1){
  Ipopsz = seq(50,600,by = 25)
  Nscn = length(Ipopsz)
}else if (Mltplscn==2){
  Initsgs = data$Segment[datasort][1:40]
  Nscn = length(Initsgs)
}
NcsSc = array(0,dim = c(P,Nyrs,Nscn))
DcsSc = array(0,dim = c(P,Nyrs,Nscn))
NcsSc_lo = array(0,dim = c(P,Nyrs,Nscn))
NcsSc_hi = array(0,dim = c(P,Nyrs,Nscn))
#
if(Nyrs>100){
  rmax = 0.2
  theta = 1.2
  Kbase = 0.9*data$K_value
}else{
  Kbase = data$K_value
}

# Set up for parallel processing
acomb <- function(...) abind(..., along=3)
# ncores = min(40,detectCores()-4)
# cl <- makeCluster(ncores)
# registerDoParallel(cl)
ncores = min(40,availableCores()-4)
cl = parallelly::makeClusterPSOCK(ncores, autoStop = F)
registerDoParallel(cl,cores=ncores)
# *** Outer loop *** 
for (z in 1:Nscn){
  if (Mltplscn==1){
    Initpop = rep( round(Ipopsz[z]/length(Initseg)), length(Initseg))
  }
  if (Mltplscn==2){  
    Initseg = Initsgs[z]   
  }  
  Yr1 = StartYear # as.numeric(format(Sys.Date(), "%Y"))+1
  alpha = 1/V_sp
  rmax = min(log(1.2),rmax)
  Years = c(Yr1:(Yr1+Nyrs-1))  
  Yrs = seq(1:Nyrs)
  Years = Yrs-1+Yr1
  # Set age/sex distribution of initial translocation
  stgdist = c(PpnFem*(1-PpnAd),PpnFem*PpnAd,
              (1-PpnFem)*(1-PpnAd),(1-PpnFem)*PpnAd)
  # set sex distribution of additional otters added (assume juvenile)
  juvdist = c(PpnFem,0,(1-PpnFem),0)
  #
  # Initialize population vector
  Initblk = numeric(length = length(Initseg))
  initblks = rep(0,P)
  Nadd = matrix(0,nrow = Nyrs,ncol = P)
  N0 = rep(0,P)
  for (i in 1:length(Initseg)){
    Initblk[i] = which(data$Segment==Initseg[i])
    initblks[Initblk[i]] = 1
    N0[Initblk[i]] = Initpop[i]
    if(AddOt[i]>0){
      for(y in 2:(AddYrs+1)){
        Nadd[y,Initblk[i]] = AddOt[i]
      }
    }
  }
  # Mean K for each Block (use average value, possibly add stochasticity?): 
  K = Kbase
  Ktot = sum(data$K_value)
  Areahab = data$Area_km2
  Kdns = K/Areahab
  Etag = data$Estry   # Tag the segments that are estuary habitat
  # initial mortality beta params
  mm = max(.001,AddMortE)
  # Assume CV of .3 for annual survival
  vv = (0.3*mm)^2
  aa = -1*(mm*(vv+mm^2-mm))/vv
  bb = (vv+mm^2-mm)*(mm-1)/vv
  # 
  # Environmental stochasticity: calc SIGMA for correlated random effects
  SIGMA = as.matrix((sig^2)*exp(-.005*Distmat))
  MU_Z = rep(0,P)
  zvec = matrix(data = 0,nrow = 4, ncol = 1)
  #
  # Dispersal probabilities
  disp = matrix(data = NA,nrow = 4, ncol = P)
  disp[1,] = DispP$Jf
  disp[2,] = DispP$Af
  disp[3,] = DispP$Jm
  disp[4,] = DispP$Am
  #
  # *Create a default matrix with lambda ~1.01 (slight positive growth)
  br = 1; wr = .6
  fsj = 0.78; fsa = 0.88; 
  msj = 0.76; msa = 0.86;
  gf = .315; gm = .302;
  FF=(br/2)*wr*fsa
  FM=FF         
  # Construct matrix               
  A = matrix(c(
    fsj*(1-gf),  FF,      0,           0,    
    fsj*gf,      fsa,     0,           0,
    0,           FM,      msj*(1-gm),  0,
    0,           0,       msj*gm,      msa),nrow=4,ncol=4,byrow = T)
  #
  W=eigen(A)$vector;          # W=matrix of right eigenvectors 
  lambdas=eigen(A)$values;    # lambdas=vector of eigenvalues
  lambda=max(Re(lambdas))		# lambda1=dominant eigenvalue, real part only
  w=abs(W[,1])					      # w=stable distribution, unscaled
  sad = w/sum(w)                # w=stable distribution, scaled
  #
  #~~~~~~~~~~~~~~ End data processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  # Run Simulations -----------------------------------------------------------
  # Initialize arrays to store results:
  # Cycle through reps
  N_sim <- foreach(r=1:reps, .combine='acomb', .multicombine=TRUE,
                   .packages=forloop_packages) %dopar% {
                     N = matrix(0,nrow = P,ncol = Nyrs)
                     # Number of years of population establishment (before range expansion begins)
                     E = pmax(5,round(rnorm(1,Ephase,1)))
                     # Determine wavespeed
                     # V = runif(1,V_mn,V_mx)
                     # Reinitialize Block occupation matrix (keeps track of which blocks occupied each year)
                     BlokOcc = matrix(data = 0,nrow = P, ncol = Nyrs)
                     # Re-Initialize population vector at year 1 (using N0) for each rep
                     n = array(data = 0, c(4,P,Nyrs))
                     nt = matrix(data = 0,nrow = 4, ncol = 1) 
                     nd = matrix(data = 0,nrow = 4, ncol = 1)
                     N0r = N0
                     # Allow for initial mortality and dispersal from reintroduction site 
                     #  (initial translocated otters that "jump" to new location)
                     pdsp = ProbMv * (AgeMvAdj * (1-PpnAd) + PpnAd)
                     V = .03
                     a =-1*(pdsp*(V+pdsp^2-pdsp))/V
                     b = (V+pdsp^2-pdsp)*(pdsp-1)/V
                     probdisp = rbeta(1,a,b)
                     # hist(probdisp)
                     for (i in 1:P){
                       if(N0[i]>0){
                         # Additional initial mortality 
                         initmort =max(.05,min(0.9, (rbeta(1,aa,bb)) ))
                         N0r[i] = round(N0[i]*(1-initmort))
                         # Determine number of animals that move as stochastic binomial variable
                         mv = rbinom(1,N0r[i], probdisp * (EstryMvAdj * (Etag[i]) + (1-Etag[i])))
                         # Ensure at least some animals left in areas selected for suppl. additions
                         if (Nadd[2,i]>0){
                           mv = min(mv,round(0.9*N0r[i]))
                         }    
                         N0r[i] = N0r[i] - mv
                         # select recipient location (outer coast)
                         noc = which(N0r==0 & Etag==0)
                         imv = noc[which(rmultinom(1,1,Kdns[noc]/sum(Kdns[noc]))==1)]
                         # Adjust for proportion of movers that die or move out of Oregon
                         N0r[imv] = round(mv*(1-Mvmort))
                       }
                     }
                     # N0r
                     for (i in 1:P){
                       if(N0r[i]>0){
                         n[1:4, i, 1] = rmultinom(1, N0r[i], stgdist)
                         BlokOcc[i,1] = 1
                       }else{
                         BlokOcc[i,1] = 0
                         n[1:4, i, 1] = c(0,0,0,0)
                       }
                     }
                     # Total otters by block in year 1 
                     N[,1] = colSums(n[, , 1])
                     # env. stoch., correlated random annual deviations 
                     eps = rmvnorm(Nyrs, MU_Z, SIGMA)
                     for (y in 2:Nyrs) {
                       #  First, determine if any new blocks occupied (allows for range expansion)
                       oc = which(BlokOcc[,y-1]==1)
                       BlokOcc[oc,y] = 1      
                       # Find all unoccupied cells
                       noc = which(BlokOcc[,y-1]==0)
                       # S_noc = data$Block[noc]
                       # Loop through unoccupied cells, see if they could be colonized
                       # by "assymptotic" range expansion after colonization phase
                       if(sum(noc) > 0 & y > E){
                         for(k in 1:length(noc)){
                           j = noc[k]
                           iib = which(Distmat[,j]<25 & Distmat[,j]>0 & BlokOcc[,(y-1)]==1)
                           nAd = length(iib)
                           # for each neighbouring block, see if its duration of occupation
                           # is greater than alpha*Distance between centroids (ie would expect 
                           #  range expansion into this new habitat)
                           if(nAd==0){
                             BlokOcc[j,y] = 0
                           }else{
                             for (a in 1:nAd){
                               Adblk = iib[a]
                               dst = Distmat[Adblk,j]  
                               # Probability of colonization: logit fxn of years occupied relative to 
                               #   movement of population front given assymptotic wave speed
                               # NOTE: estuaries take longer to become colonized from outer coast            
                               if(Etag[j]==1){
                                 probexp = inv.logit(0.5*(sum(BlokOcc[Adblk,1:(y-1)])-alpha*dst - 5))
                               }else{
                                 probexp = inv.logit(1.5*(sum(BlokOcc[Adblk,1:(y-1)])-alpha*dst))
                               }  
                               BlokOcc[j,y] = max(BlokOcc[j,y],rbinom(1,1,probexp))
                             }
                           }
                         }
                       }
                       # Next, step through blocks and compute dynamics for all occupied blocks
                       for (i in 1:P){
                         if (BlokOcc[i, y] == 1){
                           if ( sum(n[1:4,i,y-1])>0 ){
                             Nt = N[i,y-1]
                             # NOTE: account for minimal growth in "population establishment" phase 
                             # and slower growth at very sall pop sizes (allee effect):
                             if (y <= E){
                               # Note: use demographic schedule for 1.18 to get appropriate relative survival
                               # rates by age/sex, then correct for growth below (age-indep. adj. x 1/1.18 )
                               lamstoch = max(.94,min(1.20, round(1.18 + 0.5*eps[y,i],2)))
                             }else if(y > E && Nt < 5){ # Allee effect, reduced growth at very low pop sizes
                               lamstoch = max(.94,min(1.20, round(1.1 + eps[y,i],2)))
                             }else{
                               # Calculate D-D lambda (with stochasticity) and select appropriate vital rates
                               lamstoch = max(.94,min(1.20, round(exp(rmax*(1-(Nt/(K[i]))^theta)+eps[y,i]),2)))
                             } 
                             idxs = which(Demdat$Lam==lamstoch)
                             j = sample(idxs,1)
                             br = Demdat$br2[j]; wr = Demdat$wr2[j];
                             fsj = Demdat$fs1[j]; fsa = Demdat$fs2[j]; 
                             msj = Demdat$ms1[j]; msa = Demdat$ms2[j];
                             gf = Demdat$Gf[j]; gm = Demdat$Gm[j];
                             nt[1:4,1] = n[1:4,i,y-1]  
                             # Add birth rate mod if not enough males relative to females
                             br = br * min(1, 2 * (nt[4]/(nt[2]+.1))^.5)
                             # Add growth mod if age ratio too skewed to adults 
                             gf = min(1, gf / min(1, ((2.5*nt[1])/(nt[2]+.1))^.5))
                             gm = min(1, gm / min(1, ((2.5*nt[3])/(nt[4]+.1))^.5))
                             FF=(br/2)*wr*fsa
                             FM=FF         
                             # Construct matrix               
                             AP = matrix(c(
                               fsj*(1-gf),  FF,      0,           0,    
                               fsj*gf,      fsa,     0,           0,
                               0,           FM,      msj*(1-gm),  0,
                               0,           0,       msj*gm,      msa),nrow=4,ncol=4,byrow = T)
                             #
                             # Next lines do matrix multiplication (within-block demog transitions)
                             nt1 = AP%*%nt
                             # Apply demographic stochastic and 
                             # additional mortality during establishment phase
                             # and random extinctions for small pops (<5)
                             if (y <= E){
                               Emort = rbeta(1,aa,bb)
                               nt1 = nt1*.848*(1-Emort)
                               if(sum(nt1) < 5){
                                 nt1 = rbinom(1,1,.5) * nt1 
                               }           
                               nt1 = floor(nt1) + rbinom(4, 1, asin(sqrt(nt1 - floor(nt1)))/asin(1) )
                             }else if(y > E && sum(nt1) < 2){
                               nt1 = rbinom(1,1,.9) * nt1
                               nt1 = floor(nt1) + rbinom(4, 1, asin(sqrt(nt1 - floor(nt1)))/asin(1) )
                             }else{
                               # This equation applies a demographic stochasticity adjustment,
                               # randomly rounding up or down to nearest integer, with relative prob of
                               # rounding up or down determined by decimal value
                               nt1 = floor(nt1) + rbinom(4, 1, asin(sqrt(nt1 - floor(nt1)))/asin(1) )
                             } 
                           }else{
                             nt1 = zvec
                           }
                           # NEXT LINES ACCOUNT FOR ADDITIONAL RE-INTRODUCTIONS 
                           #  (randomly assign sex but assume that additional otters are juvs)
                           if (Nadd[y,i]>0) {
                             ni = c(1,0,0,0) + rmultinom(1, Nadd[y,i]-1, juvdist)              
                           } else {
                             ni = zvec
                           }
                           # Next, Calculate number of dispersers (with stochasticity)
                           # ****NOTE: no dispersal happens until pop established, 
                           # must be other occupied blocks nearby,
                           #  some individuals of each age class must remain resident
                           recip_poss = length(which(BlokOcc[,y]*Distmat[,i]>0 & BlokOcc[,y]*Distmat[,i]<150))
                           if (y > E && recip_poss>0 && sum(nt1)>10){
                             nd[1] = min(floor(.5*nt1[1]),rpois(1,nt1[1]*disp[1,i]))
                             nd[2] = min(floor(.33*nt1[2]),rpois(1,nt1[2]*disp[2,i]))
                             nd[3] = min(floor(.9*nt1[3]),rpois(1,nt1[3]*disp[3,i]))
                             nd[4] = min(floor(.5*nt1[4]),rpois(1,nt1[4]*disp[4,i]))
                           } else {
                             nd[1:4] = zvec
                           }
                           n[1:4,i,y] = pmax(zvec, n[1:4,i,y] + nt1 + ni - nd)
                           #
                           # Now distribute dispersers randomly among "currently occupied" blocks 
                           # with probabilities determined appropriately for each age/sex class
                           if (sum(nd) > 0){  
                             # JF dispersal
                             JFprob = BlokOcc[,y]*destJF[,i] 
                             JFprob = JFprob/sum(JFprob)
                             ndJF = rmultinom(1, nd[1], JFprob)
                             n[1,,y] = n[1,,y] + t(ndJF)
                             # AF dispersal
                             AFprob = BlokOcc[,y]*destAF[,i] 
                             AFprob = AFprob/sum(AFprob)
                             ndAF = rmultinom(1, nd[2], AFprob)
                             n[2,,y] = n[2,,y] + t(ndAF)     
                             # JM dispersal
                             JMprob = BlokOcc[,y]*destJM[,i] 
                             JMprob = JMprob/sum(JMprob)
                             ndJM = rmultinom(1, nd[3], JMprob)
                             n[3,,y] = n[3,,y] + t(ndJM)
                             # AM dispersal
                             AMprob = BlokOcc[,y]*destAM[,i] 
                             AMprob = AMprob/sum(AMprob)
                             ndAM = rmultinom(1, nd[4], AMprob)
                             n[4,,y] = n[4,,y] + t(ndAM)            
                           }
                         }
                       }
                       # Tabulate sub-population abundance in each block 
                       N[,y] = colSums(n[,,y])
                       # Nmn[,y] = Nmn[,y] + N[,y,r]/reps
                     }
                     return(N)
                   }
  N = N_sim
  D = N_sim
  for (y in 1:Nyrs){
    for (i in 1:reps){
      D[, y, i] = N_sim[, y, i] / Areahab
    }}
  Nmn = apply(N_sim,c(1,2),mean)
  Nlo = apply(N_sim,c(1,2),quantile,prob=0.1)
  Nhi = apply(N_sim,c(1,2),quantile,prob=0.9)
  Dmn = Nmn
  for (y in 1:Nyrs){
    Dmn[,y] = Nmn[,y]/Areahab
  }
  NcsSc[,,z] = Nmn
  DcsSc[,,z] = Dmn
  NcsSc_lo[,,z] = Nlo
  NcsSc_hi[,,z] = Nhi
  Nfin[z] = sum(Nmn[,Nyrs])
  Nfin_lo[z] = sum(Nlo[,Nyrs])
  Nfin_hi[z] = sum(Nhi[,Nyrs])
  
  for(p in 1:P){
    write.csv(data.frame(Years, D[p, ,]), 
              paste0('./results/replicates/Scenario-', 
                     Scenario_name, 'Segment-', data$Segment[p],
                     '.csv'),
              row.names = FALSE)
  }
}
# save(Dmn, file = './results/DensityProjections.Rdata')
# stopCluster(cl)


# if(Mltplscn == 1){
#   df_target = data.frame(Abund_25 = c(100,200,300),
#                          Init_req = rep(0,3))
#   df_tmp = data.frame(N_introduced = Ipopsz,
#                       N25_mn = Nfin,N25_lo = Nfin_lo,N25_hi = Nfin_hi)
#   ggplot(df_tmp,aes(x=N_introduced,y=N25_mn)) +
#     geom_ribbon(aes(ymin=N25_lo,ymax=N25_hi),alpha=0.3) +
#     geom_line() +
#     labs(x="Number introduced",y="Number after 25 years",
#          title="Population size vs. number introduced") +
#     theme_classic()
#   # Target 100
#   ii = which(abs(df_tmp$N25_mn-100)==min(abs(df_tmp$N25_mn-100)))
#   adj = df_tmp$N25_mn[ii]/100
#   df_target$Init_req[1] = ceiling(df_tmp$N_introduced[ii]/adj)
#   df_sec_Tseries_100 = round(NcsSc[,,ii]/adj,digits = 1); 
#   #adj = sum(df_sec_Tseries_100[,Nyrs])/100
#   df_sec_Tseries_100 = cbind(data$Segment,as.data.frame(df_sec_Tseries_100))
#   colnames(df_sec_Tseries_100) = c("Section",as.character(Years))
#   # Target 200
#   ii = which(abs(df_tmp$N25_mn-200)==min(abs(df_tmp$N25_mn-200)))
#   adj = df_tmp$N25_mn[ii]/200
#   df_target$Init_req[2] = ceiling(df_tmp$N_introduced[ii]/adj)
#   df_sec_Tseries_200 = round(NcsSc[,,ii]/adj,digits = 1)
#   df_sec_Tseries_200 = cbind(data$Segment,as.data.frame(df_sec_Tseries_200))
#   colnames(df_sec_Tseries_200) = c("Section",as.character(Years))
#   # Target 300
#   ii = which(abs(df_tmp$N25_mn-300)==min(abs(df_tmp$N25_mn-300)))
#   adj = df_tmp$N25_mn[ii]/300
#   df_target$Init_req[3] = ceiling(df_tmp$N_introduced[ii]/adj)
#   df_sec_Tseries_300 = round(NcsSc[,,ii]/adj,digits = 1)
#   df_sec_Tseries_300 = cbind(data$Segment,as.data.frame(df_sec_Tseries_300))
#   colnames(df_sec_Tseries_300) = c("Section",as.character(Years))
#   
#   wb = createWorkbook()
#   addWorksheet(wb, "Target_100") 
#   writeData(wb, sheet = "Target_100", x = df_sec_Tseries_100, startCol = 1)
#   addWorksheet(wb, "Target_200") 
#   writeData(wb, sheet = "Target_200", x = df_sec_Tseries_200, startCol = 1)
#   addWorksheet(wb, "Target_300") 
#   writeData(wb, sheet = "Target_300", x = df_sec_Tseries_300, startCol = 1)
#   saveWorkbook(wb, file=paste0("./results/SctnTrends_IntroSct_",paste(Initseg,collapse="_"),
#                                ".xlsx"),overwrite = TRUE)
# }
# if(Mltplscn == 2){
#   df_tmp = data.frame(InitSeg = factor(Initsgs,levels = Initsgs),
#                       N25_mn = Nfin,N25_lo = Nfin_lo,N25_hi = Nfin_hi)
#   df_tmp = df_tmp[-40,]
#   ggplot(df_tmp,aes(x=InitSeg,y=N25_mn)) +
#     geom_errorbar(aes(ymin=N25_lo,ymax=N25_hi),alpha=0.3) +
#     geom_point() +
#     labs(x="Coastal Segment of Intriduction",y="Number after 25 years",
#          title="Population size vs. number introduced") +
#     theme_classic()
# }
# 
# tryCatch(stopCluster(cl), error = function(e1) e1 = print("Cluster terminted"))
# tryCatch(rm(cl), error = function(e1) e1 = print("Cluster removed"),
#          warning = function(e1) e1 = print("Cluster removed"))
# 
# if(Nyrs >100){
#   Ntot = t(apply(N_sim,c(2,3),sum))
#   ii = which(Ntot[,Nyrs]>0); Ntot = Ntot[ii,]
#   Nmn_tot = apply(Ntot,2,mean)
#   Nlo_tot = apply(Ntot,2,quantile,prob=0.025)
#   Nhi_tot = apply(Ntot,2,quantile,prob=0.975)
#   Nmn_tot = as.numeric(smooth.spline(Nmn_tot,spar=.3)$y)
#   Nlo_tot = as.numeric(smooth.spline(Nlo_tot,spar=.3)$y)
#   Nhi_tot = as.numeric(smooth.spline(Nhi_tot ,spar=.3)$y)
#   Nyrs_to_K = which(abs(Nhi_tot - 1*Ktot)==min(abs(Nhi_tot - 1*Ktot)))
#   # Nyrs_to_K = which(abs(Nmn_tot - .5*Ktot)==min(abs(Nmn_tot - .5*Ktot)))
#   df_Y2K = data.frame(Year = seq(1,Nyrs),
#                       Ott_est_mn = Nmn_tot,
#                       Ott_est_lo = Nlo_tot,
#                       Ott_est_hi = Nhi_tot)
#   plt_Y2K = ggplot(df_Y2K,aes(x=Year,y=Ott_est_mn)) +
#     geom_ribbon(aes(ymin=Ott_est_lo,ymax=Ott_est_hi),alpha = 0.2) +
#     geom_line(size=1.1) + 
#     geom_hline(yintercept = Ktot, color="red", linetype = "dashed") +
#     labs(x="Year of Projection",y="Estimated abundance") +
#     ggtitle("Long-term population projection",
#             subtitle = paste(c(paste0("Initial translocation of ", 
#                                       sum(Initpop), " introduced to sections "),Initseg),collapse=" ")) +
#     theme_classic()
#   print(paste0("Number years to reach K statewide = ",Nyrs_to_K))
#   print(plt_Y2K)
# }
# 
# #  Do some plots ---------------------------------------------------------
# # Heatmap of Density vs Coastal Block
# if(Mltplscn == 0 & Nyrs <= 100){
#   
#   if(length(Initblk)==1){
#     plotlab = paste(c(paste0("Initial translocation of ",Initpop, " located in section"),Initseg),
#                     collapse=" ")
#     
#     if (sum(AddOt) > 0 ){
#       plotlab = paste(c(plotlab, paste0("supplemented by ",AddOt," subadults per year for ",AddYrs,"yrs")), 
#                       collapse=" ")
#     }
#   }else{
#     plotlab = paste(c(paste0("Initial translocation of ", sum(Initpop), " divided among sections"),
#                       Initseg),
#                     collapse=", ")
#     if (sum(AddOt) > 0 ){
#       plotlab = paste(c(plotlab, paste0("supplemented by ",sum(AddOt)," subadults per year for ",AddYrs,"yrs")), 
#                       collapse=" ")
#     }
#   }
#   
#   dfDens = data.frame(Blocks = data$Segment,Sort = data$Sortord,
#                       Area = data$Area_km2, Dmn)
#   colnames(dfDens) = c('Block','Order','Area',as.character(Years))
#   df_Dens <- melt(dfDens, id.vars = c("Block","Order","Area"))
#   names(df_Dens)[4:5] <- c("Year", "Density")
#   df_Dens$Number = df_Dens$Density * df_Dens$Area
#   dfDens = df_Dens[with(df_Dens,order(Order,Year)),]; rm(df_Dens)
#   dfDens$Block = factor(dfDens$Block,levels = data$Segment[order(data$Sortord)])
#   maxD <- ceiling(100*max(dfDens$Number))/100
#   plt1 = ggplot(dfDens, aes(Year, Block)) +
#     geom_tile(aes(fill = Number), color = "white") +
#     scale_fill_gradient(low = "white", high = "steelblue",limits=c(0, maxD)) +
#     xlab("Year in Future") +
#     ylab("Coastal Section #") +
#     theme(legend.title = element_text(size = 12),
#           legend.text = element_text(size = 12),
#           plot.title = element_text(size=14,face="bold"),
#           axis.title=element_text(size=12),
#           axis.text.y = element_text(size=8),
#           axis.text.x = element_text(angle = 90, hjust = 1)) +
#     labs(fill = "Mean Expected Number",
#          title=paste0("Projected Number by Section after ", Nyrs," Years"),
#          subtitle=plotlab) +
#     ggtitle(paste0("Projected Number by Section after ", Nyrs," Years"),
#             subtitle=plotlab)
#   print(plt1)
#   
#   # Trend plot of abundance over time
#   # Calculate mean trend and CI (use bootstrap CI, 1000 reps)
#   # First, entire Metapopulation (all of SEAK):
#   mean.fun <- function(dat, idx) mean(dat[idx], na.rm = TRUE)
#   CIL.fun <- function(dat, idx) quantile(dat[idx], 0.025, na.rm = TRUE)
#   CIH.fun <- function(dat, idx) quantile(dat[idx], 0.975, na.rm = TRUE)
#   Extnct.fun <- function(dat, idx) 100*(length(dat[idx[dat[idx]<1]])/length(dat[idx]))
#   means = numeric(length=Nyrs)
#   SEmean = numeric(length=Nyrs)
#   Lo = numeric(length=Nyrs)
#   Hi = numeric(length=Nyrs)
#   ExtnctP = numeric(length=Nyrs)
#   CImn = matrix(0,nrow = Nyrs,ncol=2)
#   Nsum = colSums(N[,1,])
#   means[1] = mean(Nsum)
#   bootobj = boot(Nsum, CIL.fun, R=1000, sim="ordinary")
#   Lo[1] = mean(bootobj$t)
#   bootobj = boot(Nsum, CIH.fun, R=1000, sim="ordinary")
#   Hi[1] = mean(bootobj$t)
#   CImn[1,1:2] = sum(Nmn[,1])
#   for(y in 2:Nyrs){
#     Nsum = colSums(N[,y,])
#     bootobj = boot(Nsum, mean.fun, R=1000, sim="ordinary")
#     means[y] = median(bootobj$t)
#     SEmean[y] = sd(bootobj$t)
#     CImn[y,1] = means[y] - 1.96*SEmean[y]
#     CImn[y,2] = means[y] + 1.96*SEmean[y]
#     # tmp = boot.ci(bootobj, type="bca", conf = 0.90); CImn[y,] = tmp$bca[2:3]
#     bootobj = boot(Nsum, CIL.fun, R=1000, sim="ordinary")
#     Lo[y] = mean(bootobj$t)
#     bootobj = boot(Nsum, CIH.fun, R=1000, sim="ordinary")
#     Hi[y] = mean(bootobj$t)  
#     bootobj = boot(Nsum, Extnct.fun, R=1000, sim="ordinary")
#     ExtnctP[y] = mean(bootobj$t)  
#   }
#   Pop_Overall <- data.frame(Year=Years,Mean=means,lower=Lo,upper=Hi,
#                             SEmean=SEmean,CImeanLo=CImn[,1],CImeanHi=CImn[,2],
#                             ExtnctP = ExtnctP)
#   if (Historical_comp==1){
#     titletxt = paste0("Projected Sea Otter Population: Avg. after ", Nyrs," years = ",
#                       round(Pop_Overall$Mean[Nyrs]),
#                       ", Probability of Extinction = ",round(Pop_Overall$ExtnctP[Nyrs]),"%")
#     
#   }else{
#     titletxt = paste0("Projected Sea Otter Population, ", Nyrs," Years: mean = ",
#                       round(Pop_Overall$Mean[Nyrs]),
#                       "(95% CI ",round(Pop_Overall$CImeanLo[Nyrs]),
#                       " - ", round(Pop_Overall$CImeanHi[Nyrs]),")")
#   }
#   maxN = ceiling( max(Pop_Overall$upper)/100)*100
#   plt2 = (ggplot(Pop_Overall, aes(Year, Mean))+
#             geom_line(data=Pop_Overall)+
#             geom_ribbon(data=Pop_Overall,aes(ymin=lower,ymax=upper),alpha=0.2)+
#             geom_ribbon(data=Pop_Overall,aes(ymin=CImeanLo,ymax=CImeanHi),alpha=0.3)+
#             ylim(0,maxN) +  
#             xlab("Year") +
#             ylab("Expected Abundance") +
#             ggtitle(titletxt, subtitle=plotlab)) + theme_classic(base_size = 12)
#   
#   if (Historical_comp==1){
#     minyr = min(Pop_Overall$Year)
#     tmp = data.frame(Year = minyr + c(2,3,4,5,6,7,8,11),
#                      Count = c(21,23,21,13,12,4,4,1))
#     plt2 = plt2 + geom_point(data=tmp,aes(x=Year,y=Count),color="red",size=2)
#   }
#   
#   print(plt2)
#   
#   if (plotrngexp==1){
#     # Compute observed rate of range spread:
#     # NOTES: 
#     # - if starting from one focal area (e.g. south end of island), then there are 
#     #  two "mostly independent" coastlines (east and west coast) that range front
#     #  is moving along, so observed rate of range spread should be ~ 2x V_sp
#     # - block considered "occupied" when >2 otters present, on average
#     # - range extent corrected for coastline complexity, to approximate 1-D coast
#     #
#     Range = numeric(length = Nyrs) 
#     for (i in 1:Nyrs){
#       ii = which(Nmn[,i]>2.5 & Etag==0)
#       tmp = numeric()
#       for(j in 1:length(ii)){
#         tmp[j] = mean(sort(Distmat[,ii[j]],decreasing=F)[2:3])
#       }
#       ii = which(Nmn[,i]>2.5 & Etag==1)
#       tmp2 = numeric()
#       for(j in 1:length(ii)){
#         tmp2[j] = .5
#       }
#       Range[i] = sum(tmp) # + sum(tmp2) #*CoastCrct
#     }
#     df_Rngspr = data.frame(Year = Years[(Ephase+1):min(80,Nyrs)],
#                            Range_ext = Range[(Ephase+1):min(80,Nyrs)])
#     ftV = lm(Range_ext ~ Year, data=df_Rngspr)
#     summary(ftV)
#     WScrct = 1/(length(Initpop)*2) 
#     plt3 = ggplot(data = df_Rngspr, aes(x=Year,y=Range_ext)) +
#       geom_point() +
#       stat_smooth(method = "lm", col = "red") + 
#       labs(x="Year",y="Coastal range extent (km)",
#            title = paste("Observed range spread over", Nyrs,"Years,", "V_sp = ",V_sp, "km/yr"),
#            subtitle = paste("Adj R2 = ",signif(summary(ftV)$adj.r.squared, 3),
#                             ", Slope =",signif(ftV$coef[[2]], 3),
#                             ", Approx Observed Wave Spd (correcting for # initial pops) =", 
#                             signif(ftV$coef[[2]]*WScrct,2))) +
#       theme_classic()
#     print(plt3)
#   }
#   #
#   # Output summary tables ---------------------------------------------------
#   # Simulation summary (Overall):
#   
#   if (savetable==1 | plotmap == 1){
#     # Calculate summary For each block 
#     randsmp = sample(seq(1,reps),1000,replace = TRUE)
#     meansALL = numeric(length=P*Nyrs)
#     LoALL = numeric(length=P*Nyrs)
#     HiALL = numeric(length=P*Nyrs)
#     YearsALL = numeric(length=P*Nyrs)
#     BlockIDs = character(length=P*Nyrs)
#     BlockArea = numeric(length=P*Nyrs)
#     for (p in 1:P){
#       tmp = matrix(nrow=1000,ncol=Nyrs)
#       for(r in 1:1000){
#         tmp[r,] <- (N[p,,randsmp[r]])
#       }
#       Lo = numeric(length=Nyrs)
#       Hi = numeric(length=Nyrs)
#       means = numeric(length=Nyrs)
#       means[1] = mean(tmp[,1])
#       Lo[1] = mean(tmp[,1])
#       Hi[1] = mean(tmp[,1])
#       for(y in 2:Nyrs){
#         means[y] = mean(tmp[,y])
#         # Optional: load library BMS, and fit density curve
#         # dens = as.numeric(density(tmp[,y],adjust = 3))
#         # CI =  quantile.density(dens, probs=c(.05, .95))
#         # Lo[y] <- max(0,as.numeric(CI[1]))
#         # Lo[y] <- max(0,as.numeric(CI[1]))
#         Lo[y] <- quantile(tmp[,y], probs=c(.05))
#         Hi[y] <- quantile(tmp[,y], probs=c(.95))
#       }
#       meansALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = means
#       LoALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = Lo
#       HiALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = Hi
#       YearsALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = Years
#       BlockIDs[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = as.character(rep(data$Segment[p],Nyrs))
#       BlockArea[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = rep(Areahab[p],Nyrs)
#     }
#     Hab_Blocks <- data.frame(Block = BlockIDs, Area = BlockArea, 
#                              Year=YearsALL, Mean=meansALL, 
#                              lower=LoALL,upper=HiALL,Density=dfDens$Density,
#                              DensityLO = LoALL/BlockArea, DensityHI = HiALL/BlockArea)
#     Hab_Blocks_Fin = Hab_Blocks[which(Hab_Blocks$Year==max(Hab_Blocks$Year)),]
#     # Output block summaries of simulation results 
#     if (savetable==1){
#       savename = paste0("_Sct_",paste0(Initseg,collapse = "_"),"_Nott_",
#                         paste(Initpop,collapse = "_"))
#       write.csv(Hab_Blocks,paste0('./results/',Scenario_name,
#                                   'Results_ORblks',savename,'.csv'),row.names = FALSE)
#       write.csv(Pop_Overall,paste0('./results/',Scenario_name,
#                                    'Results_tot_',savename,'.csv'),row.names = FALSE)
#     }
#   }
#   #
#   if (plotmap == 1){
#     Cdata <- ORgrd %>% 
#       group_by(Segment,id) %>%
#       summarize(Kval_mn = mean(Kval_mn, na.rm = TRUE))
#     CellID = character()
#     BlockID = character()
#     Celldens = numeric()
#     for (i in 1:P){
#       ii = which(Cdata$Segment==data$Segment[i])
#       CellID = c(CellID, Cdata$id[ii] )
#       BlockID = c(BlockID, Cdata$Segment[ii])
#       cellscale = Cdata$Kval_mn[ii] / sum(Cdata$Kval_mn[ii])
#       cellscale = cellscale / max(cellscale)
#       Celldens = c(Celldens, Hab_Blocks_Fin$Mean[i]*cellscale)
#     }
#     CellDensProject = data.frame(CellID =CellID,
#                                  BlockID=BlockID,
#                                  Celldens = Celldens)
#     write.csv(CellDensProject,
#               paste0("./results/",Scenario_name,"Results_CellDensProject.csv"),
#               row.names = FALSE)
#     #
#     # MAP OUTPUT:  ------------------------------------------------------
#     endvals = Hab_Blocks_Fin[,c(1,4,7)]
#     ORblk = merge(ORblk, endvals, by.x='Segment', by.y = 'Block')
#     plt4 = ggplot() + 
#       geom_sf(data=ORblk, aes(fill = Mean, color=Mean),
#                    alpha = 1,size = 1) +
#       scale_fill_continuous(low = "#fff7ec", high = "#7F0000") + 
#       scale_color_continuous(guide = "none", low = "#fff7ec", high = "#7F0000") + 
#       scale_x_continuous(name = "Longitude", breaks = seq(-125,-123)) + # , breaks = NULL, labels = NULL
#       scale_y_continuous(name = "Latitude") + # , breaks = NULL, labels = NULL
#       geom_sf(data = ORlnd, aes(fill=piece),
#                    color="wheat4", fill="cornsilk1",size = 0.1) +   
#       # north(HGlnd,location = "topright") +
#       # scalebar(HGlnd, dist = 50, dist_unit = "km", st.size = 3.5, 
#       #         transform = FALSE, location = "bottomleft") +
#       ggtitle("Projected Distribution") +
#       # coord_equal(ratio=1) + 
#       coord_sf() +
#       theme_minimal()
#     print(plt4)
#   }
# }  

