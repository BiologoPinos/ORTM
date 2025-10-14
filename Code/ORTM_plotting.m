% ORTM_plotting.m
    % June 2025

% Authorship: 
    % Andrés Pinos-Sánchez | andres.pinos.sanchez@gmail.com
    % Co-authors: Jess Hopf, Leif Rasmuson, Mark Novak, Will White

% Coding notes:
    % multiply by 0.006 to get per transect values


%% In how/which many simulations (SIMS) does kelp persist over time?

% Persistence and extension of kelp across simulations (last twenty years "4*20" avg)
    KelpRR = median(mean(sum(kts(1:2,(end-80):end,1,1,kelp_avg>0)),2)) *0.006;
    DriftRR = median(mean(sum(kts(3,(end-80):end,1,1,kelp_avg>0)),2)) *0.006;
    UrchinsRR = median(mean(sum(uts(2:3,(end-80):end,1,1,kelp_avg>0)),2)) *0.006;

% How many persist (don't go extinct) at end?
    sum(kelp_avg>0)

% Which ones don't/do go extinct 
    find(kelp_avg>0) % persti = find(kelp_avg>0);
    find(kelp_avg==0) % exti = find(kelp_avg==0);

%%

% Plot SIMS persistence (Kelp) over time (for multiple reps, single scenario)
    figure(444)
    hold on
    plot((1:T2+1),sum(kt2(2,:,:)>0,3),'r','LineWidth',1)
    % xline(dist.yrs,'--r')
    % xline(dist.yrs(end)+8*4+1,'--k')
    % xline(dist.yrs(1)+mngt.time+(1:mngt.length),':k')
    xlabel('Timesteps (seasons)')
    ylabel('Proportion of simulations with persisting kelp forests')
    ylim([0,RR])


%% Plot distributions of mean biomass for (across) all replicates

figure()
    % 1) Kelp 🌿 (juvenile + adult only)
        subplot(3,1,1)
        histogram(mean(sum(kts(1:2,(end-80):end,1,1,kelp_avg>0)),2)*0.006,20)
        xlabel('Kelp biomass (kg)')
    % 2) Urchins 🟣 (hiding + exposed only)
        subplot(3,1,2)
        histogram(mean(sum(uts(2:3,(end-80):end,1,1,kelp_avg>0)),2)*0.006,20)
        xlabel('Urchin biomass (kg)')
    % 3) Crabs total 🦀 (males + females, all age classes)
        subplot(3,1,3)
        histogram(mean(sum(cfts(:,(end-80):end,1,1,kelp_avg>0)) + ...
                       sum(cmts(:,(end-80):end,1,1,kelp_avg>0)),2)*0.006,20)
        xlabel('Total crab biomass (kg)')

% figure % (crabs females & males separated)
%     % 1) Crab females 🦀♀ (all age classes)
%         subplot(2,1,1)
%         histogram(mean(sum(cfts(:,(end-80):end,1,1,kelp_avg>0)),2)*0.006,20)
%         xlabel('Female crab biomass (kg)')
% 
%     % 2) Crab males 🦀♂ (all age classes)
%         subplot(2,1,2)
%         histogram(mean(sum(cmts(:,(end-80):end,1,1,kelp_avg>0)),2)*0.006,20)
%         xlabel('Male crab biomass (kg)')


%% Plot kelp biomass vs urchin biomass

% edit this, only sample summer timestep, becuase comparable pisco data for
% kelp is only of summer (august)

% sample over 10 yrs (after 10 years)
kts_all = reshape(mean(sum(kts(1:2,41:80,1,1,:)),2),1,[])*0.006;
uts_all = reshape(mean(sum(uts(2:3,41:80,1,1,:)),2),1,[])*0.006;
% plot
figure
hold on
scatter(uts_all, kts_all, 'k', 'filled', 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
% limits
xlim([0 150])
ylim([0 800])
% axis
xlabel('Urchin biomass (kg.60m^2)')
ylabel('Kelp biomass (kg.60m^2)')


% % sample over last 20 yrs, only summer (season 3)
% summer_idx = 120:3:180; % summer seasons in the last 20 years
% kts_all = reshape(mean(kts(1:2, summer_idx, 1, 1, :), 2), 1, []) * 0.006;
% uts_all = reshape(mean(uts(2:3, summer_idx, 1, 1, :), 2), 1, []) * 0.006;
% % plot
% figure
% hold on
% scatter(uts_all, kts_all, 'k', 'filled', 'MarkerFaceAlpha', .2, 'MarkerEdgeAlpha', .2)
% % limits
% xlim([0 150])
% ylim([0 500])
% % axis
% xlabel('Urchin biomass (kg.60m^2)')
% ylabel('Kelp biomass (kg.60m^2)')


%% Single replicate run: urchin-kelp figure and outputs

% Which replicate? 
R = 2;

% Summary stats over last 20 years (~80 seasons)
Kelp = [max(sum(kts(1:2,(end-80):end,1,1,R))), min(sum(kts(1:2,(end-80):end,1,1,R))), mean(sum(kts(1:2,(end-80):end,1,1,R)))];
Drift = [max(kts(3,(end-80):end,1,1,R)), min(kts(3,(end-80):end,1,1,R)), mean(kts(3,(end-80):end,1,1,R))];
Urchins = [max(sum(uts(2:3,(end-80):end,1,1,R))), min(sum(uts(2:3,(end-80):end,1,1,R))), mean(sum(uts(2:3,(end-80):end,1,1,R)))];

% Summary per transect
Kelp_transect = Kelp*0.006;
Drift_transect = Drift*0.006;
Urchins_transect = Urchins*0.006;

% Plot dynamics
figure

% Kelp
subplot(2,1,1)
hold on
plot(repmat((1:(T2+1))'./4,1,3), kts(:,:,1,R)')
yline(8.3e4,'--k')                  % reference line
xline(dist.yrs./4,'--r')           % disturbance timing
grid minor
xlim([10,(T2+1)/4])
ylabel('Kelp density (kg/ha)')
legend('Juvenile','Adult','Drift')

% Urchins
subplot(2,1,2)
hold on
plot(repmat((1:(T2+1))'./4,1,2), uts(2:3,:,1,R)')
xline(dist.yrs./4,'--r')
grid minor
xlim([10,(T2+1)/4])
legend('Hiding adults','Exposed adults')
ylabel('Urchin density (kg/ha)')
xlabel('Time (seasons)')

% Add annotation
% text(0.01,0.3,...
%     "Kelp: RK = " + kelp.RK + ", mu = " + kelp.mu + ", rS = " + kelp.rS + ", g = " + kelp.g + ", c = " + kelp.c + ", d = " + kelp.d + newline +...
%     "Urchins: RU = " + urchin.RU + ", MJ = " + urchin.MJ + ", MH = " + urchin.MH + ", ME = " + urchin.ME + newline +...
%     "Switching: w1 = " + urchin.w1 + ", w2 = " + urchin.w2 + ", kmin = " + urchin.kmin + newline +...
%     "Initial: kt = [" + num2str(kts(:,1,1,1,R)') + "], ut = [" + num2str(uts(:,1,1,1,R)') + "]" + newline +...
%     "Final avg: Kelp = " + Kelp(3) + ", Drift = " + Drift(3) + ", Urchins (adults) = " + Urchins(3) + newline +...
%     "Scenario: Disturbance = " + (~isnan(dist.yrs(1))) + ", Fishing F = " + pred.F + ", Culling = " + urchin.culling + ", Restoration = " + kelp.restore)

%% Single replicate run: urchin-kelp-crab figure and outputs

% R = 2;
% 
% % Crabs
% CrabF_total = sum(cfts(:,:,1,R),1); % females summed over age
% CrabM_total = sum(cmts(:,:,1,R),1); % males summed over age
% Crab_total = CrabF_total + CrabM_total;
% 
% % Plot dynamics
% figure
% 
% % 3) Crabs total
% subplot(3,1,1)
% hold on
% plot((1:(T2+1))'./4, Crab_total')
% xline(dist.yrs./4,'--r')
% grid minor
% xlim([10,(T2+1)/4])
% ylabel('Crabs total (kg/ha)')
% 
% % 4) Crabs females by age
% subplot(3,1,2)
% hold on
% plot(repmat((1:(T2+1))'./4,1,5), CrabF_total')
% xline(dist.yrs./4,'--r')
% grid minor
% xlim([10,(T2+1)/4])
% ylabel('Crabs F (age)')
% legend(arrayfun(@(x) sprintf('F%d',x-1), 1:5, 'UniformOutput', false))
% 
% % 5) Crabs males by age
% subplot(3,1,3)
% hold on
% plot(repmat((1:(T2+1))'./4,1,5), CrabM_total')
% xline(dist.yrs./4,'--r')
% grid minor
% xlim([10,(T2+1)/4])
% ylabel('Crabs M (age)')
% xlabel('Time (seasons)')
% legend(arrayfun(@(x) sprintf('M%d',x-1), 1:5, 'UniformOutput', false))
