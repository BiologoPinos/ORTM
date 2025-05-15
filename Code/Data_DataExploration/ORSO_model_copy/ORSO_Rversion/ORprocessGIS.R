# Script to process spatial data for Gwaii Haanas sea otter populaiton model
#
# Load necessary libraries
library(gdistance)
library(marmap)
library(raster)
library(terra) 
library(sp) 
library(sf)
library(fields)
library(fitdistrplus)
#
# Load Data ----------------------------------------------------------------------
# Data for coastl blocks and habitat cells around Haida Gwaii
Bdata = read.csv("./data/OR_pops.csv", header = TRUE)
BLK = as.matrix(cbind(Bdata$Xcoord,Bdata$Ycoord))
pts = st_as_sf(Bdata, coords = c("Xcoord","Ycoord"))
ORblkpts = vect("./GIS/Sea_otter_hab_sct_pts.shp", layer="Sea_otter_hab_sct_pts")
# Read in polygon layer
ORlnd_sp = vect("./GIS/OR_polygon.shp", layer="OR_polygon")
# Create raster layer and rasterize polygon
rs <- rast(ncol=2626, nrow=9550) # NOTE: gives grid size of 50m
ext(rs) <- ext(ORlnd_sp)
rp <- rasterize(ORlnd_sp, rs, 'OBJECTID')
values(rp)[values(rp) > 0] <- 100
values(rp)[is.na(values(rp))] <- 1
# rp@data@values[rp@data@values==100] = NA
# Cdata ... # NOTE: cell data is in GISdata.rdata  
# Cdata$PU_ID = Cdata$Cell_ID
# matrix of x-y coordinates of coastl block centroids
# matrix of x-y coordinates of grid cell centroids
# CELL = as.matrix(cbind(Cdata$Xcoord,Cdata$Ycoord))
# Load Raster Map of Haida Gwaii and process for LCP distances ------------------
# ORlnd <- raster("./GIS/OR_poly.tif") # Raster of Oregon, UTM Zone 10 NAD_1983
ORlnd = rp
ORlnd <- raster(ORlnd)
plot(ORlnd) #SEAK Raster
points(BLK) #plots Block Centroid points
# Process data for estimating least cost path (LCP) distances between points
trOR <- transition(ORlnd,mean, directions = 8) ##8 directions 
#correct distances (map distortion) for diagonal movements and large area of study site
# trOR=geoCorrection(trOR) 
# Compute euclidean distance and LCP distance matrices 
EUC=rdist(BLK) #euclidian distances
LCD=costDistance(trOR,BLK,BLK) #calculates LCD matrix, pairwise block-block distances 
Distmat = (LCD*52)/1000   # Convert pairwise LCP distances to km
# Save matrix as csv file to load in metapopulaiton model
write.table(Distmat, file = "./data/Distmat2.csv",row.names=FALSE, 
            na="",col.names=FALSE, sep=",")

# Compute Dispersal parameters--------------------------------------------------
# Import movement data from radio tagged sea otters in southern SE Alaska
# to use for fitting sea otter dispersal kernels:
# SEmoves = read.csv("./data/SE_Ottermoves.csv", header = TRUE)
CAmoves = read.csv("./data/ottmovesBSMB.csv", header = TRUE)
ottmoves = data.frame(dist = CAmoves$move,Sex=CAmoves$Sex+1, 
                      Age=CAmoves$Ageclass+1)
#
Distmat =  as.matrix(read.csv("./data/Distmat.csv", header = FALSE)) # Inter-blk LCP distances
# Fit sea otter dispersal kernels for 4 age/sex classes (M/F juv and adult)
#  and use these to compute a) probability of emmigration from each block,
#  and b) dispersal probability matrix, cell i,j = prob that otter emmigrating from 
#  block i ends up in block j (0's along diagonal)
#  ****CODE THIS AND SAVE RESULTS
Wpar = matrix(0,nrow = 4,ncol = 2)
Xpar = matrix(0,nrow = 4,ncol = 1)
cntr = 0
for (i in 1:2){
  for (j in 1:2){
    cntr = cntr+1
    # ii = which(SEmoves$SexN==i & SEmoves$Ageglass==j)
    # xi = pmax(0.1,SEmoves$LCD_km[ii])
    ii = which(ottmoves$Sex==i &ottmoves$Age==j)
    xi = pmax(0.1,ottmoves$dist[ii])
    fitd = fitdist(xi,"weibull")
    fitx = fitdist(xi,"exp")
    plot(fitd)
    # cdfcomp(list(fitd,fitx), addlegend=TRUE)
    Wpar[cntr,1] = fitd$estimate[1]
    Wpar[cntr,2] = fitd$estimate[2]  
    Xpar[cntr,1] = fitx$estimate  
  }
}
# Save DispP to GHDispProb.csv
NBlk = length(Bdata$Segment)
DispP = data.frame(Block = Bdata$Segment)
DispP$Jf = numeric(length =NBlk)
DispP$Af = numeric(length =NBlk) 
DispP$Jm = numeric(length =NBlk)
DispP$Am = numeric(length =NBlk)
DispPx = DispP
for (i in 1:NBlk){
  mndst = 0.75*mean(sort(Distmat[,i],decreasing=F)[2:3])
  DispP$Jf[i] = 1-pweibull(mndst,Wpar[1,1],Wpar[1,2])
  DispP$Af[i] = 1-pweibull(mndst,Wpar[2,1],Wpar[2,2])
  DispP$Jm[i] = 1-pweibull(mndst,Wpar[3,1],Wpar[3,2])
  DispP$Am[i] = 1-pweibull(mndst,Wpar[4,1],Wpar[4,2])
  DispPx$Jf[i] = 1-pexp(mndst,Xpar[1,1])
  DispPx$Af[i] = 1-pexp(mndst,Xpar[2,1])
  DispPx$Jm[i] = 1-pexp(mndst,Xpar[3,1])
  DispPx$Am[i] = 1-pexp(mndst,Xpar[4,1])
}
Disp = DispPx
write.csv(Disp,"./data/ORDispProbLR.csv",row.names = FALSE)
# Inter-pop dispersal matrices for each age/sex class
#  (save each as matrix using write.table)

ORDispMatJF = 1-pexp(Distmat,Xpar[1,1])
diag(ORDispMatJF) = 0
ORDispMatAF = 1-pexp(Distmat,Xpar[2,1])
diag(ORDispMatAF) = 0
ORDispMatJM = 1-pexp(Distmat,Xpar[3,1])
diag(ORDispMatJM) = 0
ORDispMatAM = 1-pexp(Distmat,Xpar[4,1])
diag(ORDispMatAM) = 0

SumJF = colSums(ORDispMatJF)
SumAF = colSums(ORDispMatAF)
SumJM = colSums(ORDispMatJM)
SumAM = colSums(ORDispMatAM)
for (i in 1:NBlk){
  ORDispMatJF[,i] = ORDispMatJF[,i]/SumJF[i]
  ORDispMatAF[,i] = ORDispMatAF[,i]/SumAF[i]
  ORDispMatJM[,i] = ORDispMatJM[,i]/SumJM[i]
  ORDispMatAM[,i] = ORDispMatAM[,i]/SumAM[i]  
}
write.table(ORDispMatJF, file = "./data/DispMatJF.csv",row.names=FALSE, 
            na="",col.names=FALSE, sep=",")
write.table(ORDispMatAF, file = "./data/DispMatAF.csv",row.names=FALSE, 
            na="",col.names=FALSE, sep=",")
write.table(ORDispMatJM, file = "./data/DispMatJM.csv",row.names=FALSE, 
            na="",col.names=FALSE, sep=",")
write.table(ORDispMatAM, file = "./data/DispMatAM.csv",row.names=FALSE, 
            na="",col.names=FALSE, sep=",")
# Calculate pairwise LCP distances from each hab cell to each block centroid -----
#  (for use in interpolating densities by weighted averaging) 
# NOTE: NO LONGER USING THIS
# LCD=costDistance(trHAIDA,CELL,BLK) #calculates LCD matrix, distances cell to block centroid 
# DistmatHab = as.data.frame(LCD/1000)
# #
# # Create "HabAvg" matrix: uses inverse distance^2 weighting,
# # for row i, col j, value represents proportional contribution for block j 
# # on habitat call i
# tmp = DistmatHab^2
# HabAvg = 1/(tmp+0.000001)
# Divisor = rowSums(HabAvg)
# for (i in 1:NBlk){
#   HabAvg[,i] = HabAvg[,i]/Divisor
# }
# HabAvg = data.frame(cbind(Cdata$PU_ID,HabAvg))
# colnames(HabAvg) = c('PU_ID',as.character(Bdata$BlockID))  
# write.csv(HabAvg,"HabAvg.csv",row.names = FALSE)  
detach(package:gdistance)
detach(package:raster)