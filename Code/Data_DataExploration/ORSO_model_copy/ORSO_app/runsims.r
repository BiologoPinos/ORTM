runsims <- function(){
  # Load files ---------------------------------------------------------
  data = read.csv("./data/OR_pops.csv", header = TRUE)  # Data for coastal blocks
  Demdat = read.csv("./data/RandDem.csv", header = TRUE)    # Stochastic vital rates
  Distmat =  as.matrix(read.csv("./data/Distmat.csv", header = FALSE)) # Inter-blk LCP distances
  # Probabilities of dispersal from each block for each age/sex class
  DispP = read.csv("./data/ORDispProb.csv", header = TRUE) 
  # Inter-pop movement matrices: pairwise prob of dispersal based on 
  #  pairwise LCP distances and dispersal kernels for each age/sex class
  destJF = read.csv("./data/DispMatJF.csv", header = FALSE) 
  destAF = read.csv("./data/DispMatAF.csv", header = FALSE)
  destJM = read.csv("./data/DispMatJM.csv", header = FALSE)
  destAM = read.csv("./data/DispMatAM.csv", header = FALSE)
  # Load GIS data for plotting results (NAD_1983_Albers)
  load("./data/GISdata.rdata")
  # User Parameters  ---------------------------------------------------------
  # (set via UI inputs)
  PpnFem = input$PpnFem    # Proportion of females in introduced pop
  PpnAd = input$PpnAd      # Proportion of adults in introduced pop
  reps = input$Reps        # Number replications for population sims (should use at least 100)
  Nyrs = input$Nyrs        # Number of years to project population dynamics
  AddYrs = input$AddYrs    # Number of years that "supplemental" otters added to introduced pops
  Ephase = input$Ephase      # Maximum years before pop "established" (before range expansion begins)
  ProbMv = input$ProbMv        # Probability of substantial % of otters moving in establishment phase
  EstryMvAdj = input$EstryMvAdj   # Proportional dispersal prob for estuaries relative to outer coast
  AgeMvAdj = input$AgeMvAdj       # Proportional dispersal probability for subadults relative to adults 
  Mvmort = input$Mvmort           # proportion of initially dispersing animals that die/disappear from Oregon
  AddMortE = input$AddMortE      # Excess annual mortality during establishment phase (0.13 for San Nic)  
  V_sp = input$V_sp        # Population front asymptotic wavespeed, km/yr, minimum  
  sig = input$Estoch       # Environmental stochasticity (std. deviation in log-lambda)
  rmax = input$Rmax        # Maximum rate of growth: default = log(1.22), or 22% per year
  theta = input$Theta      # theta parameter for theta-logistic (1 = Ricker model, >1 = delayed DD) 
  Init_df = values$Intro_df # Data frame of introduction params (where and how many)
  if(nrow(Init_df)==1){
    Initseg = as.character(Init_df[1]) # List of coastal sections where otters introduced (Initial translocation)
    Initpop = as.integer(Init_df[2])    # Number of animals in initial introduction to each section  
    AddOt = as.integer(Init_df[3])      # Number of animals per year in supplemental additions to each section 
  }else{
    Initseg = as.character(Init_df[,1]) # List of coastal sections where otters introduced (Initial translocation)
    Initpop = as.integer(Init_df[,2])    # Number of animals in initial introduction to each section  
    AddOt = as.integer(Init_df[,3])      # Number of animals per year in supplemental additions to each section 
  }
  #
  # Process data and set up sims --------------------------------------
  #
  Yr1 = as.numeric(format(Sys.Date(), "%Y"))+1
  alpha = 1/V_sp
  rmax = min(log(1.2),rmax)
  Years = c(Yr1:(Yr1+Nyrs-1))  
  Yrs = seq(1:Nyrs)
  Years = Yrs-1+Yr1
  P = nrow(Distmat)  # number blocks (or sub-populations)
  datasort = order(as.numeric(data$Ycoord))
  data$Sortord = rep(0,P)
  data$Sortord[datasort] = seq(1,P)
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
  K = data$K_value
  Areahab = data$Area_km2
  Kdns = K/Areahab
  Etag = data$Estry   # Tag the segments that are estuary habitat
  #
  # Assume CV of .25 for annual survival
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
  N = array(data = 0, c(P,Nyrs,reps))
  # Nmn = matrix(0,nrow=P,ncol=Nyrs)
  # Nmn[,1] = N0
  zvec = matrix(data = 0,nrow = 4, ncol = 1)
  #
  withProgress(message = 'Running simulations, please be patient...', 
               session=session, 
  { # for loop for simulations
    for (r in 1:reps){
      incProgress(0.8*(1/reps),session=session)
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
      # Allow for initial high mortality and the possibility that a high proportion 
      #  of initial translocated otters "jump" to new location  
      pdsp = ProbMv * (AgeMvAdj * (1-PpnAd) + PpnAd)
      V = .03
      a =-1*(pdsp*(V+pdsp^2-pdsp))/V
      b = (V+pdsp^2-pdsp)*(pdsp-1)/V
      probdisp = rbeta(1,a,b)
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
      N[,1,r] = colSums(n[, , 1])
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
        if(sum(noc) > 0 && y > E){
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
              Nt = N[i,y-1,r]
              # NOTE: account for minimal growth in "population establishment" phase 
              # and slower growth at very sall pop sizes (allee effect):
              if (y <= E){
                lamstoch = max(.94,min(1.20, round(1.18 + 0.5*eps[y,i],2)))
              }else if(y > E & Nt <= 5){
                lamstoch = max(.94,min(1.20, round(1.05 + eps[y,i],2)))
              }else{
                # Calculate D-D lambda (with stochasticity) and select appropriate vital rates
                lamstoch = max(.94,min(1.20, round(exp(rmax*(1-(Nt/K[i])^theta)+eps[y,i]),2)))
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
              # Add growth mod if age ratio too skewed (too few juvs)
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
              # Apply demographic stichasticity and 
              # additional mortalty during establishment phae
              # and random extinctions for small pops (<5)
              if (y <= E){
                Emort = rbeta(1,aa,bb)
                nt1 = nt1*.848*(1-Emort)
                if(sum(nt1) < 5){
                  nt1 = rbinom(1,1,.5) * nt1 
                }           
                nt1 = floor(nt1) + rbinom(4, 1, asin(sqrt(nt1 - floor(nt1)))/asin(1) )
              }else if(y > E && sum(nt1) < 5){
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
        N[,y,r] = colSums(n[,,y])
        # Nmn[,y] = Nmn[,y] + N[,y,r]/reps
      }
    }
    # Process results --------------------------------------------------
    # 
    # Mean abundance over time
    Nmn = apply(N,c(1,2),mean)
    # Mean density over time
    Dmn = Nmn
    for (y in 1:Nyrs){
      Dmn[,y] = Nmn[,y]/Areahab
    }
    #
    dfDens = data.frame(Blocks = data$Segment,Sort = data$Sortord,
                        Area = data$Area_km2, Dmn)
    colnames(dfDens) = c('Block','Order','Area',as.character(Years))
    df_Dens <- melt(dfDens, id.vars = c("Block","Order","Area"))
    names(df_Dens)[4:5] <- c("Year", "Density")
    df_Dens$Number = df_Dens$Density * df_Dens$Area
    dfDens = df_Dens[with(df_Dens,order(Order,Year)),]; rm(df_Dens)
    dfDens$Block = factor(dfDens$Block,levels =  data$Segment[order(data$Sortord)])
    # Trend plot of abundance over time
    # Calculate mean trend and CI (use bootstrap CI, 1000 reps)
    mean.fun <- function(dat, idx) mean(dat[idx], na.rm = TRUE)
    CIL.fun <- function(dat, idx) quantile(dat[idx], 0.025, na.rm = TRUE)
    CIH.fun <- function(dat, idx) quantile(dat[idx], 0.975, na.rm = TRUE)
    means = numeric(length=Nyrs)
    SEmean = numeric(length=Nyrs)
    Lo = numeric(length=Nyrs)
    Hi = numeric(length=Nyrs)
    CImn = matrix(0,nrow = Nyrs,ncol=2)
    Nsum = colSums(N[,1,])
    means[1] = mean(Nsum)
    bootobj = boot(Nsum, CIL.fun, R=100, sim="ordinary")
    Lo[1] = mean(bootobj$t)
    bootobj = boot(Nsum, CIH.fun, R=100, sim="ordinary")
    Hi[1] = mean(bootobj$t)
    CImn[1,1:2] = sum(Nmn[,1])
    for(y in 2:Nyrs){
      Nsum = colSums(N[,y,])
      bootobj = boot(Nsum, mean.fun, R=100, sim="ordinary")
      means[y] = median(bootobj$t)
      SEmean[y] = sd(bootobj$t)
      CImn[y,1] = means[y] - 1.96*SEmean[y]
      CImn[y,2] = means[y] + 1.96*SEmean[y]
      # tmp = boot.ci(bootobj, type="bca", conf = 0.90); CImn[y,] = tmp$bca[2:3]
      bootobj = boot(Nsum, CIL.fun, R=100, sim="ordinary")
      Lo[y] = mean(bootobj$t)
      bootobj = boot(Nsum, CIH.fun, R=100, sim="ordinary")
      Hi[y] = mean(bootobj$t)  
    }
    Pop_Overall <- data.frame(Year=Years,Mean=means,lower=Lo,upper=Hi,
                              SEmean=SEmean,CImeanLo=CImn[,1],CImeanHi=CImn[,2])
    #
    # Calculate summary of abundance and density by block/year, with CI
    randsmp = sample(seq(1,reps),1000,replace = TRUE)
    meansALL = numeric(length=P*Nyrs)
    LoALL = numeric(length=P*Nyrs)
    HiALL = numeric(length=P*Nyrs)
    YearsALL = numeric(length=P*Nyrs)
    BlockIDs = character(length=P*Nyrs)
    BlockArea = numeric(length=P*Nyrs)
    for (p in 1:P){
      tmp = matrix(nrow=1000,ncol=Nyrs)
      for(r in 1:1000){
        tmp[r,] <- (N[p,,randsmp[r]])
      }
      Lo = numeric(length=Nyrs)
      Hi = numeric(length=Nyrs)
      means = numeric(length=Nyrs)
      means[1] = mean(tmp[,1])
      Lo[1] = mean(tmp[,1])
      Hi[1] = mean(tmp[,1])
      for(y in 2:Nyrs){
        means[y] = mean(tmp[,y])
        # Optional: load library BMS, and fit density curve
        # dens = as.numeric(density(tmp[,y],adjust = 3))
        # CI =  quantile.density(dens, probs=c(.05, .95))
        # Lo[y] <- max(0,as.numeric(CI[1]))
        # Lo[y] <- max(0,as.numeric(CI[1]))
        Lo[y] <- quantile(tmp[,y], probs=c(.05))
        Hi[y] <- quantile(tmp[,y], probs=c(.95))
      }
      meansALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = means
      LoALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = Lo
      HiALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = Hi
      YearsALL[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = Years
      BlockIDs[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = as.character(rep(data$Segment[p],Nyrs))
      BlockArea[((p-1)*Nyrs+1):((p-1)*Nyrs+Nyrs)] = rep(Areahab[p],Nyrs)
    }
    Hab_Blocks <- data.frame(Block = BlockIDs, Area = BlockArea, 
                             Year=YearsALL, Mean=meansALL, 
                             lower=LoALL,upper=HiALL,Density=dfDens$Density,
                             DensityLO = LoALL/BlockArea, DensityHI = HiALL/BlockArea)
    Hab_Blocks_Fin = Hab_Blocks[which(Hab_Blocks$Year==max(Hab_Blocks$Year)),]
    Hab_Blocks_Fin = Hab_Blocks_Fin[data$Sortord,]
    #
    datasim = Hab_Blocks_Fin[,c(1,4,7)]
    ORblk = merge(ORblk, datasim, by.x='id', by.y = 'Block')
    titletxt <- paste0("Projected Sea Otter Population, ", Nyrs," Years")
    mapplot2 = ggplot() + 
      geom_polygon(data=ORblk, aes(x = long, y = lat, fill = Mean, color=Mean, group = group),
                   alpha = 1,size = 1) +
      scale_fill_continuous(low = "#fff7ec", high = "#7F0000") + 
      scale_color_continuous(guide = FALSE, low = "#fff7ec", high = "#7F0000") + 
      scale_x_continuous(name = "Longitude", breaks = seq(-125,-123)) + # , breaks = NULL, labels = NULL
      scale_y_continuous(name = "Latitude") + # , breaks = NULL, labels = NULL
      geom_polygon(data = ORlnd, aes(x=long,y=lat,fill=piece,group=group),
                   color="wheat4", fill="cornsilk1",size = 0.1) +   
      # north(HGlnd,location = "topright") +
      # scalebar(HGlnd, dist = 50, dist_unit = "km", st.size = 3.5, 
      #         transform = FALSE, location = "bottomleft") +
      ggtitle(titletxt) +
      # coord_equal(ratio=1) + 
      coord_map("conic", lat0 = 18, xlim = c(-125, -123), ylim=c(41.8, 46.3)) +
      theme_minimal()
  })
  Tab1 <- Pop_Overall
  colnames(Tab1) <- c("Year",
                      "Average Number",
                      "Lower Estimate (CI)",
                      "Upper Estimate (CI)",
                      "Estimation Uncertainty (SE)",
                      "Lower 95% CI for the Mean",
                      "Upper 95% CI for the Mean")
  colnames(Hab_Blocks_Fin) <- c("Coastal Section",
                      "Area (km2)",
                      "Year",
                      "Average Number",
                      "Lower Estimate (CI)",
                      "Upper Estimate (CI)",
                      "Density (#/km2)",
                      "Lower Density Est.(#/km2)",
                      "Upper Density Est.(#/km2)")
  # Update "values" structure
  values$Pop_Overall <- Pop_Overall
  values$Tab1 <- Tab1
  values$dfDens <- dfDens
  values$Hab_Blocks <- Hab_Blocks
  values$Hab_Blocks_Fin <- Hab_Blocks_Fin
  # 
  ggsave("./www/mapdens.png",plot=mapplot2,width = 6,height = 8,dpi=300)
  return(datasim)
}
