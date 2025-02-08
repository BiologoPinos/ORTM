% PredUrchinKelp_ImplicitCC_v6_MngtScen_v1_Plotting2.m
% May 2023
% 
% Jess Hopf
% jess.k.hopf@gmail.com 

% updated plotting code from PredUrchinKelp_ImplicitCC_v6_MngtScen_v1_Plotting.m
% goes with PredUrchinKelp_ImplicitCC_v6_MngtScen_v1.m



%% medians over multiple reps
% (last 4*20 yrs avg)
% multipled by 0.006 to get per transect values
predsRR = median(mean(PBE((end-80):end,1,1,kelp_avg>0))) *0.006
KelpRR = median(mean(sum(kts(1:2,(end-80):end,1,1,kelp_avg>0)),2)) *0.006
DriftRR = median(mean(sum(kts(3,(end-80):end,1,1,kelp_avg>0)),2)) *0.006
UrchinsRR = median(mean(sum(uts(2:3,(end-80):end,1,1,kelp_avg>0)),2)) *0.006

% how many persist (don't go extinct) at end?
sum(kelp_avg>0)

% which ones dont/do go extinct 
persti = find(kelp_avg>0);
exti = find(kelp_avg==0);

%% plot persistence over time (for multiple reps, single scenario)

figure(234524)
hold on
plot((1:T2+1),sum(kt2(2,:,:)>0,3),'k','LineWidth',1)
% xline(dist.yrs,'--r')
% xline(dist.yrs(end)+8*4+1,'--k')
% xline(dist.yrs(1)+mngt.time+(1:mngt.length),':k')
xlabel('Time (seasons)')
ylabel('Number of sims persisting')
ylim([0,10000])


%% plot distributions of mean biomasses for all replicates

% over last 20 years

figure
subplot(3,1,1)
histogram(mean(PBE((end-80):end,1,1,kelp_avg>0))*0.006,20)
xlabel('Pred biomass')
subplot(3,1,2)
histogram(mean(sum(kts(1:2,(end-80):end,1,1,kelp_avg>0)),2)*0.006,20)
xlabel('Kelp biomass')
subplot(3,1,3)
histogram(mean(sum(uts(2:3,(end-80):end,1,1,kelp_avg>0)),2)*0.006,20)
xlabel('Urchin biomass')


%% plot kelp biomass v urchin biomass
% 
% kts_all = reshape(sum(kts(1:2,:,1,1,:)),1,[])*0.006;
% uts_all = reshape(sum(uts(1:2,:,1,1,:)),1,[])*0.006;


kts_all = reshape(mean(sum(kts(1:2,(end-4):end,1,1,:)),2),1,[])*0.006;
uts_all = reshape(mean(sum(uts(2:3,(end-4):end,1,1,:)),2),1,[])*0.006;

% kts_all_k = reshape(mean(sum(kts(1:2,(end-4):end,1,1,kelp_avg>0)),2),1,[])*0.006/60;
% uts_all_k = reshape(mean(sum(uts(2:3,(end-4):end,1,1,kelp_avg>0)),2),1,[])*0.006/60;

figure
hold on
scatter(uts_all, kts_all, 'k', 'filled', 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
% scatter(uts_all_k, kts_all_k, 'k', 'filled', 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
xline(1.3*60, '--r') % approx cali kg threshold for shift (Ling et al 2015)
xline(2.5*60, '--r') % approx global kg threshold for shift (Ling et al 2015)
xlabel('Urchin biomass (kg.60m^2)')
ylabel('Kelp biomass (kg.60m^2)')
xlim([0,600])
ylim([0,80*60])




%% Single rep run: urchin-kelp figure and outputs

% which rep? 
R = 1;

    % calculate if persisting or not (mean over last 1 year)
    % kelp_state = double(kelp_avg(1,1,1,R) > 1)
    

% (4*20 yrs avg)
% preds = mean(PBE);
Kelp = [max(sum(kts(1:2,(end-80):end,1,1,R))), min(sum(kts(1:2,(end-80):end,1,1,R))), mean(sum(kts(1:2,(end-80):end,1,1,R)))];
Drift = [max(kts(3,(end-80):end,1,1,R)), min(kts(3,(end-40):end,1,1,R)), mean(kts(3,(end-80):end,1,1,R))];
Urchins = [max(sum(uts(2:3,(end-80):end,1,1,R))), min(sum(uts(2:3,(end-80):end,1,1,R))), mean(sum(uts(2:3,(end-80):end,1,1,R)))];
 
% preds_transect = preds*0.006
Kelp_transect = Kelp*0.006
Drift_transect = Drift*0.006
Urchins_transect = Urchins*0.006


figure
subplot(3,1,1)
hold on
plot(repmat((1:(T2+1))'./4,1,3), kts(:,:,1,R)')
yline(8.3*10^4,'--k')
xline(dist.yrs./4,'--r')
grid minor
xlim([10,(T2+1)/4])
ylabel('Kelp density (kg.ha)')
legend('Juv Kelp','Adult kelp','Drift kelp')

subplot(3,1,2)
hold on
plot(repmat((1:(T2+1))'./4,1,2), uts(2:3,:,1,R)')
xline(dist.yrs./4,'--r')
grid minor
xlim([10,(T2+1)/4])
legend('Hiding adults','Exposed adult')
ylabel('Urchin density (kg.ha)')
xlabel('Time (seasons)')

subplot(3,1,3)
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
    
axis off

%% Single run: predator figure

figure
hold on
plot((1:T2),sum(nb(:,(end-T2+1):end),1),'k')
plot((1:T2),PBE,'r')

%% heatmap scenarios for rep runs (single action)

% sample year 
samp_t = dist.yrs(end)+8*4+1; % 26;

% baseline dist (from scenario runs)
bd = repmat(squeeze(kelp_pers_t(:,samp_t,:,1,1)),1,length(mngt.length));

% baseline undist
load('Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_nodist_20241011.mat');
bud = baseline_persit_overtime(samp_t);

kelp_samp_rec2undist = round((squeeze(kelp_pers_t(:,samp_t,:,:))-bd)./(bud-bd).*100,0);

edges = [0 25 50 75 100 200];
kelp_samp_rec2undist_dist = discretize(kelp_samp_rec2undist, edges, edges(1:(end-1)));

figure
subplot(2,1,1)
heatmap(mngt.time./4, flip(mngt.length(2:6)./4), flip(kelp_samp_rec2undist(:,2:6)'),...
       'MissingDataColor',[0.9,0.9,0.9])
xlabel('Start of management action (yr relative to start of heatwave)')
ylabel('Length of management action (yrs)')
title("%Increase in prob of persistence: " + mngt_scen)
% colormap summer
colormap(cbrewer('seq','Greens',100,'linear'))   % Blues  Greens Purples
clim([0 100]) 

subplot(2,1,2)
heatmap(mngt.time./4, flip(mngt.length(2:6)./4), flip(kelp_samp_rec2undist_dist(:,2:6)'),...
       'MissingDataColor',[0.9,0.9,0.9])
xlabel('Start of management action (yr relative to start of heatwave)')
ylabel('Length of management action (yrs)')
title("%Increase in prob of persistence: " + mngt_scen)
colormap(cbrewer('seq','Greens',100,'linear'))   % Blues  Greens Purples


%% heatmap scenarios for rep runs (single action, var sample year)

% sample year 
% samp_t = dist.yrs(end)+8*4+1; % 26;

% load baseline undist vals
load('Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_nodist_20241011.mat');

% for loop to select year x after end of mngt action

for i = 1:length(mngt.time)
    for j = 1:length(mngt.length)

        % selection time (end of disturbance or end of action + 1yr, which is last)
        samp_t = max(dist.yrs(1) + mngt.time(i) + mngt.length(j), dist.yrs(end)) +4;

        % baseline dist (from scenario runs)
        bd = squeeze(kelp_pers_t(:,samp_t,i,1,1));
        
        % baseline undist
        bud = baseline_persit_overtime(samp_t);

        kelp_samp_rec2undist(i,j) = round((kelp_pers_t(:,samp_t,i,j)-bd)./(bud-bd).*100,0);
    end
end


edges = [0 25 50 75 100 200];
kelp_samp_rec2undist_dist = discretize(kelp_samp_rec2undist, edges, edges(1:(end-1)));

figure
% subplot(3,1,2);
heatmap(mngt.time./4, flip(mngt.length(2:6)./4), flip(kelp_samp_rec2undist(:,2:6)'),...
       'MissingDataColor',[0.9,0.9,0.9])
xlabel('Start of management action (yr relative to start of heatwave)')
ylabel('Length of management action (yrs)')
title("%Increase in prob of persistence: " + mngt_scen)
% colormap summer
colormap(cbrewer('seq','Greens',100,'linear'))   % Blues  Greens Purples
clim([0 100]) 

% subplot(2,1,2)
% heatmap(mngt.time./4, flip(mngt.length(2:6)./4), flip(kelp_samp_rec2undist_dist(:,2:6)'),...
%        'MissingDataColor',[0.9,0.9,0.9])
% xlabel('Start of management action (yr relative to start of heatwave)')
% ylabel('Length of management action (yrs)')
% title("%Increase in prob of persistence: " + mngt_scen)
% colormap(cbrewer('seq','Greens',100,'linear'))   % Blues  Greens Purples

%% Lineplot: Effects of length and degree of mngt action (single action)

% sample year  (dist.yrs(end))
samp_t = 88+8*4+1; % 26;

% time raltive to heatwave 
trh = 4; % -4

% baseline undist
load('Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_nodist_20241011.mat')
bud = baseline_persit_overtime(samp_t);

% load culling values
load("Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_culling_20241011_degreesLRG.mat")
% time raltive to heatwave 
timeC = find(mngt.time == trh); 
% baseline dist
bd_c = repmat(squeeze(kelp_pers_t(:,samp_t,timeC,1,1)),length(mngt.length),length(mngt.degreeC));
proppersist_cull = round((squeeze(kelp_pers_t(:,samp_t,timeC,:,:,:))-bd_c)./(bud-bd_c).*100,0);
degreeC = mngt.degreeC;

% load resoration values
load("Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_restoration_20241011_degreesLRG.mat")
% time raltive to heatwave 
timeR = find(mngt.time == trh);
% baseline dist
bd_r = repmat(squeeze(kelp_pers_t(:,samp_t,timeR,1,1)),length(mngt.length),length(mngt.degreeR));
proppersist_rest = round((squeeze(kelp_pers_t(:,samp_t,timeR,:,:,:))-bd_r)./(bud-bd_r).*100,0);
degreeR = mngt.degreeR;

% fishery closure benefits
load("Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_fishing_20241008.mat")
% time raltive to heatwave 
timeF = find(mngt.time == trh); 
% baseline dist (from scenario runs)
bd = squeeze(kelp_pers_t(:,samp_t,timeF,1,1));
proppersist_fish = round((squeeze(kelp_pers_t(:,samp_t,timeF,2:end))-bd)./(bud-bd).*100,0);

% load pre-heatwave values
load("Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_none_20241011.mat",...
    "kelpJ_avg_pre", "kelpA_avg_pre", "urchinA_avg_pre")
kelpJA_avg_pre = kelpJ_avg_pre + kelpA_avg_pre;

% increase v degree w/length
figure(11)
hold on
subplot(3,1,3)
hold on
plot((degreeC./urchinA_avg_pre).*100, proppersist_cull(2:end,:), 'LineWidth', 2)
% plot(degreeCF./100, proppersist_cullfish(2:end,:),'--', 'LineWidth', 2)
% plot(degreeCF./100, proppersist_cullfish(2:end,:)-proppersist_cull(2:end,2:end), 'LineWidth', 2)
xlabel('% Urchins removed')
ylabel('% Increase in probability of persistence')
mycolors =viridis(5);
ax = gca; 
ax.ColorOrder = mycolors;
% lgd = legend(num2str(mngt.length(2:end)'./4));
% title(lgd,'Length of action (yrs)')
% title("culling")
yline(proppersist_fish,':b')

figure(12)
hold on
subplot(3,1,3)
hold on
plot(degreeR./kelpJA_avg_pre.*100, proppersist_rest(2:end,:), 'LineWidth', 2)
% plot(degreeRF./(1.3*10^5).*100, proppersist_restfish(2:end,:),'--', 'LineWidth', 2)
% plot(degreeRF./100, proppersist_cullfish(2:end,:)-proppersist_rest(2:end,2:end), 'LineWidth', 2)
xlabel('% Kelp biomass added')
ylabel('% Increase in probability of persistence')
mycolors =viridis(5);
ax = gca; 
ax.ColorOrder = mycolors;
% lgd = legend(num2str(mngt.length(2:end)'./4));
% title(lgd,'Length of action (yrs)')
% title("restoration")
% ylim([0,120])
yline(proppersist_fish,':b')


%% Heatmap: combined 

% join kelp & urchin datasets
clear

% heatmap
% sample year 
samp_t = 88+8*4+1; % 26;

% baseline undist
load('Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_nodist_20241011.mat');
bud = baseline_persit_overtime(samp_t);

% baseline disturbed
load("Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_none_20241011.mat")
bd = baseline_persit_overtime(samp_t);

% load first dataset
load('Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_cull&rest_20241014_cull01.mat')
kelp_pers_all = kelp_pers_t;

% time realtive to heatwave [-4 0 4]
trh = 0;
trhi = find(mngt.time == trh);

% mngt action length [12]
mngtl = 12;
mngtli = find(mngt.length == mngtl);

proppersist = round((squeeze(kelp_pers_t(:,samp_t,trhi,mngtli,:))-bd)./(bud-bd).*100,0)';

% load rest and add on
load('Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_cull&rest_20241014_cull05.mat')
proppersist(2,:) = round((squeeze(kelp_pers_t(:,samp_t,trhi,mngtli,:))-bd)./(bud-bd).*100,0)';
load('Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_cull&rest_20241014_cull10.mat')
proppersist(3,:) = round((squeeze(kelp_pers_t(:,samp_t,trhi,mngtli,:))-bd)./(bud-bd).*100,0)';
load('Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_cull&rest_20241014_cull25.mat')
proppersist(4,:) = round((squeeze(kelp_pers_t(:,samp_t,trhi,mngtli,:))-bd)./(bud-bd).*100,0)';
load('Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_cull&rest_20241014_cull50.mat')
proppersist(5,:) = round((squeeze(kelp_pers_t(:,samp_t,trhi,mngtli,:))-bd)./(bud-bd).*100,0)';
load('Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_cull&rest_20241014_cull75.mat')
proppersist(6,:) = round((squeeze(kelp_pers_t(:,samp_t,trhi,mngtli,:))-bd)./(bud-bd).*100,0)';
load('Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_cull&rest_20241014_cull100.mat')
proppersist(7,:) = round((squeeze(kelp_pers_t(:,samp_t,trhi,mngtli,:))-bd)./(bud-bd).*100,0)';

% individual kelp (to add)
    % load resoration values
    load("Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_restoration_20241011_degreesLRG.mat")
    
    % mngt action length [12]
    mngtli = find(mngt.length == 12);
    
    % time raltive to heatwave [-4 0 4]
    trhi = find(mngt.time == trh);
    
    % baseline dist
    bd_r = repmat(squeeze(kelp_pers_t(:,samp_t,trhi,mngtli,1)),1,length(mngt.degreeR));
    proppersist_rest = round((squeeze(kelp_pers_t(:,samp_t,trhi,mngtli,:))'-bd_r)./(bud-bd_r).*100,0);
    degreeR = mngt.degreeR;

% flip combinedand add single mngt action values
proppersistAdd = [flip(proppersist); proppersist_rest];
% 
% figure(34)
% % subplot(3,1,2)
% heatmap([0,1,5,10,25,50,75,100], flip([0,1,5,10,25,50,75,100]), proppersistAdd,...
%        'MissingDataColor',[0.9,0.9,0.9])
% xlabel('Degree of Kelp Seeding')
% ylabel('Degree of Urchin Removal')
% title("% Increase prob persistence. Time:" + trh + ", length: " + mngtl/4)
% % colormap summer
% % colormap(cbrewer('seq','Purples',100,'linear'))   %  Blues Greens
% clim([0 150]) 
% 
% % sum of individ actions 
% sumproppersist = repmat(proppersistAdd(:,1),1,length(degreeR)) + repmat(proppersist_rest,length(degreeR),1);
% diffproppersist = proppersistAdd./(sumproppersist);
% 
% figure(23)
% % subplot(3,1,1)
% heatmap([0,1,5,10,25,50,75,100], flip([0,1,5,10,25,50,75,100]), diffproppersist,...
%        'MissingDataColor',[0.9,0.9,0.9])
% xlabel('Degree of Kelp Seeding')
% ylabel('Degree of Urchin Removal')
% title("Differnce: rest&cull. Time: " + trh + ", length: " + mngtl/4)
% colormap(cbrewer('seq','Greys',100,'linear'))
% clim([0 1])

% combined fishing + cull/rest ------
% Note, uses stuff from previous

% sample year 
samp_t = dist.yrs(end)+8*4+1; % 26;

% time raltive to heatwave 
trh = 0; % -4

% time raltive to heatwave 
time = find(mngt.time == trh);


% Fishing + cull
load("Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_fish&cull_20241014.mat")
fishcull_rec2undist = round((squeeze(kelp_pers_t(:,samp_t,time,:,:))-bd)./(bud-bd).*100,0);

figure(240)
hold on
subplot(2,1,1)
heatmap([0,1,5,10,25,50,75,100], 1, fishcull_rec2undist', ...
       'MissingDataColor',[0.9,0.9,0.9])
xlabel('Degree of management action (%) (in addtion to closure)')
ylabel(mngt_scen)
title("%Increase in prob persistence. Time:" + trh)
clim([0 150]) 

% sum of individ actions 
sumproppersist = flip(proppersistAdd(:,1)) + fishcull_rec2undist(1);
diffproppersist = fishcull_rec2undist./(sumproppersist);

figure(241)
hold on
subplot(2,1,1)
heatmap([0,1,5,10,25,50,75,100], 1, diffproppersist', ...
       'MissingDataColor',[0.9,0.9,0.9])
xlabel('Degree of management action (%) (in addtion to closure)')
ylabel(mngt_scen)
title("%Increase in prob persistence. Time:" + trh)
colormap(cbrewer('seq','Greys',100,'linear'))
clim([0 1])


% Fishing + rest 
load("Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_fish&rest_20241014.mat")
fishrest_rec2undist = round((squeeze(kelp_pers_t(:,samp_t,time,:,:))-bd)./(bud-bd).*100,0);

figure(240)
subplot(2,1,2)
heatmap([0,1,5,10,25,50,75,100], 1, fishrest_rec2undist', ...
       'MissingDataColor',[0.9,0.9,0.9])
xlabel('Degree of management action (%) (in addtion to closure)')
ylabel(mngt_scen)
title("%Increase in prob persistence. Time:" + trh)
clim([0 150]) 

% sum of individ actions 
sumproppersist = proppersist_rest + fishrest_rec2undist(1);
diffproppersist = fishrest_rec2undist'./(sumproppersist);

figure(241)
subplot(2,1,2)
heatmap([0,1,5,10,25,50,75,100], 1, diffproppersist, ...
       'MissingDataColor',[0.9,0.9,0.9])
xlabel('Degree of management action (%) (in addtion to closure)')
ylabel(mngt_scen)
title("%Increase in prob persistence. Time:" + trh)
colormap(cbrewer('seq','Greys',100,'linear'))
clim([0 1])


%% Comparing combined actions (OLD)

% sample year 
samp_t = 30*4+1; % 26;

% time raltive to heatwave 
trh = -4; % -4

% baseline undist
load('Model outputs\MngtScenarios\Implicitv6_MngtScen_v1_nodist_20240806.mat')
bud = baseline_persit_overtime(samp_t);

% baseline disturbed (pre-disturbed values)
load("Model outputs\MngtScenarios\Implicitv6_MngtScen_v1_none_20240806.mat")
kelpJA_avg_pre = kelpJ_avg_pre + kelpA_avg_pre;

% load fishing
load("Model outputs\MngtScenarios\Implicitv6_MngtScen_v1_fishing_20240801.mat")
% time realtive to heatwave 
timer = find(mngt.time == trh); 
% baseline dist
bd_f = repmat(squeeze(kelp_pers_t(:,samp_t,timer,1,1)),length(mngt.length),1);
proppersist_fish = round((squeeze(kelp_pers_t(:,samp_t,timer,:,:,:))-bd_f)./(bud-bd_f).*100,0);

% load culling values
load("Model outputs\MngtScenarios\Implicitv6_MngtScen_v1_culling_25pc_20240816.mat")
% baseline dist
bd_c = repmat(squeeze(kelp_pers_t(:,samp_t,timer,1,1)),length(mngt.length),length(mngt.degreeC));
proppersist_cull = round((squeeze(kelp_pers_t(:,samp_t,timer,:,:,:))-bd_c)./(bud-bd_c).*100,0);

% load culling + fishing
load("Model outputs\MngtScenarios\Implicitv6_MngtScen_v1_fish&cull_20240823_25pc.mat")
% baseline dist
bd_cf = repmat(squeeze(kelp_pers_t(:,samp_t,timer,1,1)),length(mngt.length),length(mngt.degreeC));
proppersist_cullfish = round((squeeze(kelp_pers_t(:,samp_t,timer,:,:,:))-bd_cf)./(bud-bd_cf).*100,0);

% load resoration values
load("Model outputs\MngtScenarios\Implicitv6_MngtScen_v1_restoration_20240823_25pc.mat")
% baseline dist
bd_r = repmat(squeeze(kelp_pers_t(:,samp_t,timer,1,1)),length(mngt.length),length(mngt.degreeR));
proppersist_rest = round((squeeze(kelp_pers_t(:,samp_t,timer,:,:,:))-bd_r)./(bud-bd_r).*100,0);

% load restoation + fishing
load("Model outputs\MngtScenarios\Implicitv6_MngtScen_v1_fish&rest_20240823_25pc.mat")
% baseline dist
bd_rf = repmat(squeeze(kelp_pers_t(:,samp_t,timer,1,1)),length(mngt.length),length(mngt.degreeR));
proppersist_restfish = round((squeeze(kelp_pers_t(:,samp_t,timer,:,:,:))-bd_rf)./(bud-bd_rf).*100,0);

% load culling + restoration
load("Model outputs\MngtScenarios\Implicitv6_MngtScen_v1_cull&rest_20240823_25pc.mat")
% baseline dist
bd_cr = repmat(squeeze(kelp_pers_t(:,samp_t,timer,1,1)),length(mngt.length),length(mngt.degreeR));
proppersist_cullrest = round((squeeze(kelp_pers_t(:,samp_t,timer,:,:,:))-bd_cr)./(bud-bd_cr).*100,0);

% load all actions
load("Model outputs\MngtScenarios\Implicitv6_MngtScen_v1_fishrestcull_20240823_25pc.mat")
% baseline dist
bd_all = repmat(squeeze(kelp_pers_t(:,samp_t,timer,1,1)),length(mngt.length),length(mngt.degreeR));
proppersist_all = round((squeeze(kelp_pers_t(:,samp_t,timer,:,:,:))-bd_all)./(bud-bd_all).*100,0);

% combine
proppersist_actions = [proppersist_fish, proppersist_cull, proppersist_cullfish,...
    proppersist_rest, proppersist_restfish, proppersist_cullrest, proppersist_all];

% names
proppersist_names = ["Fishery closure", "Urchin removal", "Closure & removal",...
    "Kelp restoration", "Closure & restoration", "Culling & restoration", "All actions"];


%%
figure(16757)
hold on
subplot(3,1,2)
hold on
plot(1:5, proppersist_actions(2:end,:), 'LineWidth', 2)
% colororder('blue', 'purple','k', 'green')
ylabel('% Increase in probability of persistence')
xlabel('Length of action (yrs)')
xlim([0.5,5.5])
% legend(proppersist_names)


%% differences

d_cf = proppersist_cullfish - (proppersist_fish + proppersist_cull);
d_rf = proppersist_restfish - (proppersist_fish + proppersist_rest);
d_cr = proppersist_cullrest - (proppersist_cull + proppersist_rest);
d_all = proppersist_all - (proppersist_fish + proppersist_cull + proppersist_rest);

diffs = [d_cf, d_rf, d_cr, d_all];
diffs = diffs(2:end,:);

diff_names = categorical(["Closure & removal", "Closure & restoration",...
    "Culling & restoration", "All actions"], "Ordinal", true);

figure(244)
hold on
scatter(diffs', diff_names, repmat(1:5,4,1)*10, 'filled', 'r')
xline(0, '--k')

scatter(-proprest)



