% PredUrchinKelp_ImplicitCC_MngtScen_Plotting.m
% May 2025

% Authorship: 
% Andrés Pinos-Sánchez | andres.pinos.sanchez@gmail.com
% Co-authors: Jess Hopf, Leif Rasmuson, Mark Novak, Will White


% %% persistance and extension of kelp across simulations
% % (last 4*20 yrs avg)
% % multipled by 0.006 to get per transect values
% KelpRR = median(mean(sum(kts(1:2,(end-80):end,1,1,kelp_avg>0)),2)) *0.006;
% DriftRR = median(mean(sum(kts(3,(end-80):end,1,1,kelp_avg>0)),2)) *0.006;
% UrchinsRR = median(mean(sum(uts(2:3,(end-80):end,1,1,kelp_avg>0)),2)) *0.006;
% 
% how many persist (don't go extinct) at end?
sum(kelp_avg>0)
% 
% % which ones dont/do go extinct 
% persti = find(kelp_avg>0);
% exti = find(kelp_avg==0);

% %% plot persistence over time (for multiple reps, single scenario)
% 
% figure
% hold on
% plot((1:T2+1),sum(kt2(2,:,:)>0,3),'k','LineWidth',1)
% % xline(dist.yrs,'--r')
% % xline(dist.yrs(end)+8*4+1,'--k')
% % xline(dist.yrs(1)+mngt.time+(1:mngt.length),':k')
% xlabel('Time (seasons)')
% ylabel('Number of sims persisting')
% ylim([0,RR])


%% plot distributions of mean biomasses for all replicates

% over last 20 years

figure
% subplot(2,1,1)
% histogram(mean(PBE((end-80):end,1,1,kelp_avg>0))*0.006,20)
% xlabel('Pred biomass')
subplot(2,1,1)
histogram(mean(sum(kts(1:2,(end-80):end,1,1,kelp_avg>0)),2)*0.006,20)
xlabel('Kelp biomass kg')
ylabel('Number of sims')
subplot(2,1,2)
histogram(mean(sum(uts(2:3,(end-80):end,1,1,kelp_avg>0)),2)*0.006,20)
xlabel('Urchin biomass kg')
ylabel('Number of sims')


% %% plot kelp biomass v urchin biomass
% % 
% % kts_all = reshape(sum(kts(1:2,:,1,1,:)),1,[])*0.006;
% % uts_all = reshape(sum(uts(1:2,:,1,1,:)),1,[])*0.006;
% kts_all = reshape(mean(sum(kts(1:2,(end-4):end,1,1,:)),2),1,[])*0.006;
% uts_all = reshape(mean(sum(uts(2:3,(end-4):end,1,1,:)),2),1,[])*0.006;
% 
% % kts_all_k = reshape(mean(sum(kts(1:2,(end-4):end,1,1,kelp_avg>0)),2),1,[])*0.006/60;
% % uts_all_k = reshape(mean(sum(uts(2:3,(end-4):end,1,1,kelp_avg>0)),2),1,[])*0.006/60;
% 
% figure
% hold on
% scatter(uts_all, kts_all, 'k', 'filled', 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
% % scatter(uts_all_k, kts_all_k, 'k', 'filled', 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
% xline(1.3*60, '--r') % approx cali kg threshold for shift (Ling et al 2015)
% xline(2.5*60, '--r') % approx global kg threshold for shift (Ling et al 2015)
% xlabel('Urchin biomass (kg.60m^2)')
% ylabel('Kelp biomass (kg.60m^2)')
% xlim([0,600])
% ylim([0,80*60])


%% single rep run: urchin-kelp figure and outputs

% which rep? 
R = 2;

    % calculate if persisting or not (mean over last 1 year)
    % kelp_state = double(kelp_avg(1,1,1,R) > 1)
    

% (4*20 yrs avg)
% preds = mean(PBE);
Kelp = [max(sum(kts(1:2,(end-80):end,1,1,R))), min(sum(kts(1:2,(end-80):end,1,1,R))), mean(sum(kts(1:2,(end-80):end,1,1,R)))];
Drift = [max(kts(3,(end-80):end,1,1,R)), min(kts(3,(end-40):end,1,1,R)), mean(kts(3,(end-80):end,1,1,R))];
Urchins = [max(sum(uts(2:3,(end-80):end,1,1,R))), min(sum(uts(2:3,(end-80):end,1,1,R))), mean(sum(uts(2:3,(end-80):end,1,1,R)))];
 
% preds_transect = preds*0.006
Kelp_transect = Kelp*0.006;
Drift_transect = Drift*0.006;
Urchins_transect = Urchins*0.006;


figure
subplot(2,1,1)
hold on
plot(repmat((1:(T2+1))'./4,1,3), kts(:,:,1,R)')
yline(8.3*10^4,'--k')
xline(dist.yrs./4,'--r')
grid minor
xlim([10,(T2+1)/4])
ylabel('Kelp density (kg.ha)')
legend('Juv Kelp','Adult kelp','Drift kelp')

subplot(2,1,2)
hold on
plot(repmat((1:(T2+1))'./4,1,2), uts(2:3,:,1,R)')
xline(dist.yrs./4,'--r')
grid minor
xlim([10,(T2+1)/4])
legend('Hiding adults','Exposed adult')
ylabel('Urchin density (kg.ha)')
xlabel('Time (seasons)')

text(0.01,0.3,...
    "Predators: avg(bioeat) = " + mean(PBE(:,:,:,:,R)) + newline +...
    "Kelp: RK = " + kelp.RK + ",  mu = " + kelp.mu + ",  rS = " + kelp.rS + ", g = " + kelp.g + ",  c = " + kelp.c + ",  d = " + kelp.d + newline +...
    "Urchins: RU = " + urchin.RU + ",  MJ = " + urchin.MJ + ",  MH = " + urchin.MH + ",  ME = " + urchin.ME + newline +...
    "Swtiching: w1 = " + urchin.w1 + ",  w2 = " + urchin.w2 + newline +...
    "Step func: kmin = " + urchin.kmin + newline +...
    "Init: kt = [" + num2str(kts(:,1)') + "],     ut = [" + num2str(uts(:,1)') +"]" + newline +...
    "Final numbers: Kelp avg = " + Kelp(3) + ", Drift avg = " + Drift(3) + newline +...
    "               Urchins (adults) = " + Urchins(3) + newline +...
    "Scenarios: Disturbace = " + (dist.yrs(1)>0)*1 + ", F = " + pred.F + ", Temp F = " + pred.fish + ", Culling =" + urchin.culling + ", Restoration =" + kelp.restore)
% % 
axis off