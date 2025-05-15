require(readxl)
require(openxlsx)
require(ggplot2)
require(gridExtra)
tdr_all_feeding_bouts <- read_excel("~/Otterdat/tdr_all_feeding_bouts.xlsx")
tdr_all_feeding_bouts$dpthmean = pmin(100,tdr_all_feeding_bouts$dpthmean)
datF = tdr_all_feeding_bouts[which(tdr_all_feeding_bouts$sex=="F"),]
datM = tdr_all_feeding_bouts[which(tdr_all_feeding_bouts$sex=="M"),]
datF_ca = datF[which(datF$area=="ca2" | datF$area=="mba" | datF$area=="sni"),]
datM_ca = datM[which(datM$area=="ca2" | datM$area=="mba" | datM$area=="sni"),]
dat_ca = rbind(datF_ca,datM_ca)
#
brks = c(0,5,10,15,20,25,30,35,40,60,100)
F_dep_hist = hist(datF_ca$dpthmean, brks)
M_dep_hist = hist(datM_ca$dpthmean, brks)
All_dep_hist = hist(dat_ca$dpthmean, brks)
#
df_Divedepsum = data.frame(Mn_Depth = All_dep_hist$mids,
                           Depbin_lo = All_dep_hist$breaks[1:length(All_dep_hist$mids)],
                           Depbin_hi = All_dep_hist$breaks[2:length(All_dep_hist$breaks)],
                           All_otts = round(All_dep_hist$density,digits = 5),
                           Males = round(M_dep_hist$density,digits = 5),
                           Females= round(F_dep_hist$density,digits = 5))
wb = createWorkbook()
addWorksheet(wb, "Depth_Use_Sum") 
writeData(wb, sheet = "Depth_Use_Sum", x = df_Divedepsum, startCol = 1)
saveWorkbook(wb, file=("./results/CA_foragedepth_sum.xlsx"),overwrite = TRUE)

dat_ca$sex = factor(dat_ca$sex,levels = "F","M")

plt_F = ggplot(datF_ca,aes(x=dpthmean)) +
  geom_histogram(breaks = seq(0,100,by=2)) +
  labs(x = "Depth bin (m)", y = "Frequency of feeding dives", 
       title = "A) Females") +
  xlim(0,60) +
  theme_classic()

plt_M = ggplot(datM_ca,aes(x=dpthmean)) +
  geom_histogram(breaks = seq(0,100,by=2)) +
  labs(x = "Depth bin (m)", y = "Frequency of feeding dives", 
       title = "B) Males") +
  xlim(0,60) +
  theme_classic()

plt_A = ggplot(dat_ca,aes(x=dpthmean)) +
  geom_histogram(breaks = seq(0,100,by=2)) +
  labs(x = "Depth bin (m)", y = "Frequency of feeding dives", 
       title = "C) All otters") +
  xlim(0,60) +
  theme_classic()

grid.arrange(plt_F,plt_M,plt_A)

