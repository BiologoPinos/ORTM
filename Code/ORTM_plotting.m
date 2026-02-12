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

%% Plot SIMS persistence (Kelp) over time (for multiple reps, single scenario)

    figure(4)
    hold on
    % plot((1:T2+1),sum(kt2(2,:,:)>0,3)./RR,'k','LineWidth',1) % over seasons
    plot(((0:T2)/4),sum(kt2(2,:,:)>0,3)./RR,'Color','b','LineWidth',2) % over years
    % plot(((0:T2)/4),sum(kt2(2,:,:)>0,3)./RR,'k','LineWidth',1) % over years
    
    % xline(dist.yrs(1),'--','k') % start of disturbance and/or management (seasons)
    
    % xline(dist.yrs(1)/4,'--','Start Mng','Color','k','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','FontSize',11) % start of disturbance and/or management (years)
    % xline(mngt.length(end),'--','End Mng','Color','k','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','FontSize',11) % start of disturbance and/or management (years)

    % xline(buffer/4,'--','Sea Otter Reintroduction','Color','k','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','FontSize',11) % Sea otter reintroduction (years)

    ylabel('Proportion of Kelp Forests Persisting','FontSize',16)
    % ylim([0,RR])
    ylim([0,1])
    
    % xlabel('Timesteps (seasons)')
    xlabel('Years','FontSize',16)
    xlim([0,40])    



%% Plot distributions of mean biomass for (across) all replicates (where kelp persists)

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
        histogram(mean(sum(cfts(2:11,(end-80):end,1,1,kelp_avg>0)) + ...
                       sum(cmts(2:11,(end-80):end,1,1,kelp_avg>0)),2)*0.006,20)
        xlabel('Crab biomass (kg)')




%% Plot kelp biomass vs urchin biomass

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



%% Single replicate run: urchin-kelp figure and outputs

% Which replicate? 
R = 9973;

% Plot dynamics
figure()

% Kelp
subplot(2,1,1)
hold on
plot(repmat((1:(T2+1))'./4,1,3), kts(:,:,1,R)')
% yline(8.3e4,'--k')                  % reference line
xline(dist.yrs./4,'--r')           % disturbance timing
% grid minor
xlim([5,(T2+1)/4])
ylim([0 12*10^4])
ylabel('Kelp density (kg/ha)')
legend('Juvenile','Adult','Drift')

% Urchins
subplot(2,1,2)
hold on
plot(repmat((1:(T2+1))'./4,1,2), uts(2:3,:,1,R)')
xline(dist.yrs./4,'--r')
% grid minor
xlim([5,(T2+1)/4])
ylim([0 10000])
legend('Hiding adults','Exposed adults')
ylabel('Urchin density (kg/ha)')
xlabel('Years')



%% Crab biomass distribution across replicates

% --- sample only winter seasons and bound results to 0 ---
T = size(uts,2);
winterIdx = 2:4:T;                      % winters
x = (winterIdx-1)/4;                    % years corresponding to winter seasons

% Prepare full matrices then subset to winter rows
C = reshape(squeeze(sum(cmts(4,:,:,:,:,:),1)), T, []);  % T x Nsim (crabs; use 5:11 for age 4+)

C = C(winterIdx,:);    % only winters

% Bound to zero
C(C < 0) = 0;

% --- Figure 1: mean + ribbon (2 rows x 1 col) ---
figure(1111111)

% Crabs (winter only)
hold on
m = mean(C,2,'omitnan'); s = std(C,0,2,'omitnan');
upper = m + s;
lower = max(m - s, 0);                     % ensure lower bound >= 0
patch([x fliplr(x)], [upper.' fliplr(lower.')], [0.8 1 0.8], 'EdgeColor','none', 'FaceAlpha',0.35);
plot(x, m, 'k-', 'LineWidth', 1.5);
xlabel('Years'); 
xlim([5 49]);
xline(buffer/4,'--') % Sea otter reintroduction (years)
ylabel('Pre-fishery male Dungeness crab male biomass distribution (Kg/ha)');



%% Crab biomass for each age class by sex across replicates

figure()

meanF = squeeze(mean(mean(cft2(:, 80:end, :), 2), 3));  % females
meanM = squeeze(mean(mean(cmt2(:, 80:end, :), 2), 3));  % males

stdF = squeeze(std(mean(cft2(:, 80:end, :), 2), 0, 3));
stdM = squeeze(std(mean(cmt2(:, 80:end, :), 2), 0, 3));


bar_data = [meanF(:) meanM(:)];
b = bar(1:11, bar_data, 'grouped');
b(1).FaceColor = [0.85 0.2 0.2]; % red for females
b(2).FaceColor = [0.2 0.4 0.85]; % blue for males
hold on

% Add error bars
ngroups = size(bar_data, 1);
nbars = size(bar_data, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    if i == 1
        errorbar(x, meanF, stdF, 'k.', 'LineWidth', 1)
    else
        errorbar(x, meanM, stdM, 'k.', 'LineWidth', 1)
    end
end
hold off

xlabel('Crab Age Class')
ylabel('Biomass (Kg/ha)')
title('Mean Biomass per Age Class (Last 80 Time Steps) - With Sea otters')
legend({'Females','Males'}, 'Location','best')
box on



%% Urchin and crab biomass distribution across replicates

% --- sample only winter seasons and bound results to 0 ---
T = size(uts,2);
winterIdx = 2:4:T;                      % winters
x = (winterIdx-1)/4;                    % years corresponding to winter seasons

% Prepare full matrices then subset to winter rows
U = reshape(squeeze(sum(uts(2:3,:,:,:,:,:),1)), T, []);    % T x Nsim (adult urchins)
C = reshape(squeeze(sum(cmts(4,:,:,:,:,:),1)), T, []);  % T x Nsim (crabs; use 5:11 for age 4+)

U = U(winterIdx,:);    % only winters
C = C(winterIdx,:);    % only winters

% Bound to zero
U(U < 0) = 0;
C(C < 0) = 0;

% --- Figure 1: mean + ribbon (2 rows x 1 col) ---
figure()

% Urchins (winter only)
subplot(2,1,1); hold on
m = mean(U, 2,'omitnan'); s = std(U,0,2,'omitnan');
upper = m + s;
lower = max(m - s, 0);                     % ensure lower bound >= 0
patch([x fliplr(x)], [upper.' fliplr(lower.')], [0.8 0.8 1], 'EdgeColor','none', 'FaceAlpha',0.35);
plot(x, m, 'r-', 'LineWidth', 1.5);
xlabel('Years'); ylabel('Adult urchin biomass');
title('Winter: adult urchin biomass \pm 1 SD'); hold off

% Crabs (winter only)
subplot(2,1,2); hold on
m = mean(C,2,'omitnan'); s = std(C,0,2,'omitnan');
upper = m + s;
lower = max(m - s, 0);                     % ensure lower bound >= 0
patch([x fliplr(x)], [upper.' fliplr(lower.')], [0.8 1 0.8], 'EdgeColor','none', 'FaceAlpha',0.35);
plot(x, m, 'r-', 'LineWidth', 1.5);
xlabel('Years'); ylabel('Male Dungeness biomass (age 4+)');
title('Winter: male Dungeness biomass \pm 1 SD'); 
