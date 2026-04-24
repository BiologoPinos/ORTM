% ORTM_plotting.m
    % June 2025

% Authorship: 
    % Andrés Pinos-Sánchez | andres.pinos.sanchez@gmail.com
    % Co-authors: Jess Hopf, Leif Rasmuson, Mark Novak, Will White

% Coding notes:
    % multiply by 0.006 to get per transect values


%% Single replicate run: urchin-kelp (validation)
    % 
    % % Which replicate? 
    %     R = 2;
    % 
    % % Plot dynamics
    %     figure()
    % 
    %     % Kelp 🌿 
    %         subplot(2,1,1)
    %         hold on
    %         plot(repmat((1:(T2+1))'./4,1,3), kts(:,:,1,R)')
    % 
    %         % Disturbance timing
    %         xline(dist.yrs./4,'--r')           
    % 
    %         % Limits
    %         xlim([5,(T2+1)/4])
    %         ylim([0 12*10^4])
    % 
    %         % Axis
    %         ylabel('Kelp density (kg/ha)')
    %         legend('Juvenile','Adult','Drift')
    % 
    %     % Urchins 🟣 
    %         subplot(2,1,2)
    %         hold on
    %         plot(repmat((1:(T2+1))'./4,1,2), uts(2:3,:,1,R)')
    % 
    %         % Disturbance timing
    %         xline(dist.yrs./4,'--r')
    % 
    %         % Limits
    %         xlim([5,(T2+1)/4])
    %         ylim([0 10000])
    % 
    %         % Axis
    %         legend('Hiding adults','Exposed adults')
    %         ylabel('Urchin density (kg/ha)')
    %         xlabel('Years')


%% Explore initial conditions for urchins 🟣 and crabs 🦀
    % 
    % % Mean values at a given point across replicates 
    % 
    %     % Urchins 🟣 
    %         disp(mean(ut2(1, 41, :), 3)); % recruits
    %         disp(mean(ut2(2, 41, :), 3)); % Hiding urchins
    %         disp(mean(ut2(3, 41, :), 3)); % Exposed urchins
    % 
    %     % Crabs females 🦀
    %         % disp(mean(cft2(1:44,161,:), 3)); % all age classes x season
    %         disp(mean(sum(reshape(cft2(1:44, 161, :),4,11,[]),1),3)); % grouped by every 4 rows
    % 
    %     % Crabs males 🦀
    %         % disp(mean(cmt2(1:44,161,:), 3));
    %         disp(mean(sum(reshape(cmt2(1:44, 161, :),4,11,[]),1),3)); % grouped by every 4 rows
    %
    % % Plot mean biomass over time across all replicates
    %     figure()
    % 
    %     % Urchins 🟣
    %         subplot(3,1,1)
    %         hold on    
    %         plot((0:T2), squeeze(mean(ut2(2,:,:), 3))) % Hiding urchins    
    %         plot((0:T2), squeeze(mean(ut2(3,:,:), 3))) % Exposed urchins
    % 
    %      % Crabs females 🦀
    %         subplot(3,1,2)
    %         hold on    
    %         C_female = squeeze(sum(reshape(cft2(1:44,:,:),4,11,size(cft2,2),[]),1));
    %         for a = 1:11
    %             plot(0:T2, squeeze(C_female(a,:,:)))
    %         end
    % 
    %     % Crabs males 🦀
    %         subplot(3,1,3)
    %         hold on    
    %         C_male = squeeze(sum(reshape(cmt2(1:44,:,:),4,11,size(cmt2,2),[]),1));
    %         for a = 1:11
    %             plot(0:T2, squeeze(C_male(a,:,:)))
    %         end


%% Plot kelp 🌿 biomass vs urchin 🟣 biomass (validation)
    % 
    % % Plot vs/vs biomasses
    %     figure()
    % 
    %     % sample summer/fall (start of fall to start of winter) over 10 yrs (after) - transect evaluation (*0.006)
    %         kts_all = reshape(mean(sum(kts(1:2,sort([(77-41):4:77,(78-41):4:78]),1,1,:)),2),1,[])*0.006;
    %         uts_all = reshape(mean(sum(uts(2:3,sort([(77-41):4:77,(78-41):4:78]),1,1,:)),2),1,[])*0.006;
    % 
    %     % Scatter plot        
    %         hold on
    %         scatter(uts_all, kts_all, 'k', 'filled', 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
    % 
    %     % Limits
    %     % xlim([0 200])
    %     % ylim([0 800])
    % 
    %     % Axis
    %     xlabel('Urchin biomass (kg.60m^2)')
    %     ylabel('Kelp biomass (kg.60m^2)')


%% Plot distributions of mean biomass across all replicates (where kelp persists)

    % Plot
        figure()

        % 1) Kelp 🌿 (juvenile + adult only)
            subplot(3,1,1)
            histogram(mean(sum(kts(1:2,(end-80):end,1,1,kelp_avg>0)),2)*0.006,20)
            xlabel('Kelp biomass (kg/ha)')

        % 2) Urchins 🟣 (hiding + exposed only)
            subplot(3,1,2)
            histogram(mean(sum(uts(2:3,(end-80):end,1,1,kelp_avg>0)),2)*0.006,20)
            xlabel('Urchin biomass (kg/ha)')

        % 3) Crabs total 🦀 (males + females, adult classes only)
            subplot(3,1,3)
            histogram(mean(sum(cfts(9:44,(end-80):end,1,1,kelp_avg>0)) + ...
                           sum(cmts(9:44,(end-80):end,1,1,kelp_avg>0)),2)*0.006,20)
            xlabel('Crab density (N/ha)')


%% In how/which many simulations (SIMS) does kelp 🌿 persist over time?
    % 
    % % Persistence and extension of kelp across simulations (last twenty years "4*20" avg)
    %     KelpRR = median(mean(sum(kts(1:2, (end-80):end, 1, 1, kelp_avg>0)), 2)) *0.006;
    %     DriftRR = median(mean(sum(kts(3, (end-80):end, 1, 1, kelp_avg>0)), 2)) *0.006;
    %     UrchinsRR = median(mean(sum(uts(2:3, (end-80):end, 1, 1, kelp_avg>0)), 2)) *0.006;
    % 
    % % How many persist (don't go extinct) at end?
    %     sum(kelp_avg>0)
    % 
    % % Which ones don't/do go extinct 
    %     find(kelp_avg>0) % persti = find(kelp_avg>0);
    %     find(kelp_avg==0) % exti = find(kelp_avg==0);


%% Plot SIMS persistence (Kelp 🌿) over time (for multiple reps, single scenario)
    % 
    % % Plot
    %     figure()
    %     hold on
    % 
    %     % Over seasons? (note: baseline always at the end)
    %         % plot((1:T2+1), sum(kt2(2, :, :)>0, 3)./RR, 'k', 'LineWidth', 2) % baseline
    %         % plot((1:T2+1), sum(kt2(2, :, :)>0, 3)./RR, 'b', 'LineWidth', 2) % management (note: choose color)
    % 
    %     % Over years
    %         plot(((0:T2)/4), sum(kt2(2, :, :)>0, 3)./RR, 'k', 'LineWidth',1) % baseline
    %         % plot(((0:T2)/4), sum(kt2(2, :, :)>0, 3)./RR, 'b', 'LineWidth',1) % management (note: choose color)
    % 
    %     % When does management starts?
    %         % xline(dist.yrs(1)/4,'--','Start Mng','Color','k','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','FontSize',11) % start of management (years)
    %         % xline(dist.yrs(1),'--','Start Mng','Color','k','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','FontSize',11) % start of management (seasons)
    %         % xline(buffer/4,'--','Sea Otter Reintroduction','Color','k','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','FontSize',11) % Sea otter reintroduction (years)
    % 
    %     % Limits
    %         xlim([0,40])  
    %         ylim([0,1])
    %         % ylim([0,RR])
    % 
    %     % Axis
    %         ylabel('Proportion of Kelp Forests Persisting','FontSize',16)
    %         xlabel('Years','FontSize',16)
    %         % xlabel('Timesteps (seasons)')          


%% URCHIN 🟣 & CRAB 🦀 BIOMASS DISTRIBUTION ACROSS REPLICATES (seasonal scale - loop)
    % 
    % figure()
    % 
    % t = 0:size(ut2,2)-1;   % seasonal time steps
    % 
    % % Panel definitions
    % titles = { ...
    %     'Exposed urchin densities Mean ± 5–95%', ...
    %     'Fishable crab densities Mean ± 5–95%', ...
    %     'Hiding urchin densities Mean ± 5–95%', ...
    %     'Adult crab densities Mean ± 5–95%'};
    % 
    % ylabels = { ...
    %     'Adult urchin biomass (kg/ha)', ...
    %     'Fishable crab biomass (N/ha)', ...
    %     'Adult urchin biomass (kg/ha)', ...
    %     'Adult crab biomass (N/ha)'};
    % 
    % colors = { ...
    %     [0.6 0.4 0.8], ... % urchins
    %     [1 0.6 0.2], ...   % crabs
    %     [0.6 0.4 0.8], ...
    %     [1 0.6 0.2]};
    % 
    % tiledlayout(2,2)
    % 
    % for k = 1:4
    % 
    %     nexttile; hold on
    % 
    %     % Extract biomass (R x T)
    %     switch k
    %         case 1 % Exposed urchins
    %             B = squeeze(sum(ut2(3,:,:),1))';
    % 
    %         case 2 % Fishable crabs (male only, ages 17–24)
    %             B = squeeze(sum(cmt2(17:24,:,:),1))';
    % 
    %         case 3 % Hiding urchins
    %             B = squeeze(sum(ut2(2,:,:),1))';
    % 
    %         case 4 % All adult crabs (female + male)
    %             B = squeeze(sum(cmt2(9:end,:,:) + cft2(9:end,:,:),1))';
    %     end
    % 
    %     % Stats (no smoothing)
    %     m   = mean(B,1);
    %     p95 = prctile(B,95,1);
    %     p05 = prctile(B,5,1);
    % 
    %     % Plot
    %     fill([t fliplr(t)], [p95 fliplr(p05)], ...
    %          colors{k}, 'EdgeColor','none','FaceAlpha',0.4);
    % 
    %     plot(t, m, 'k-','LineWidth',1.5);
    % 
    %     % Labels
    %     xlim([20 160])
    %     ylabel(ylabels{k})
    %     title(titles{k})
    % 
    %     if k > 2
    %         xlabel('Seasonal time step')
    %     end
    % 
    %     grid on
    % end


%% URCHIN 🟣 & CRAB 🦀 BIOMASS DISTRIBUTION ACROSS REPLICATES (seasonal scale)
    % 
    % figure()
    % 
    %     t = 0:size(ut2,2)-1;   % seasonal time steps
    % 
    %         % Exposed Urchin 🟣 plot
    %             subplot(2,2,1); hold on
    % 
    %             U = squeeze(sum(ut2(3,:,:),1))';    % R x T
    % 
    %             m   = mean(U,1);
    %             p95 = prctile(U,95,1);
    %             p05 = prctile(U,5,1);
    % 
    %             fill([t fliplr(t)], [p95 fliplr(p05)], ...
    %                  [0.6 0.4 0.8], 'EdgeColor','none','FaceAlpha',0.4);
    % 
    %             plot(t, m, 'k-','LineWidth',1.5);
    % 
    %             xlim([20 160])
    %             ylabel('Adult urchin biomass (kg/ha)')
    %             title('Exposed urchin densities Mean ± 5–95%')
    % 
    %         % Fishable Crab 🦀 plot
    %             subplot(2,2,2); hold on
    % 
    %             C = squeeze(sum(cmt2(17:24,:,:),1))';   % R x T
    % 
    %             m   = mean(C,1);
    %             p95 = prctile(C,95,1);
    %             p05 = prctile(C,5,1);
    % 
    %             fill([t fliplr(t)], [p95 fliplr(p05)], ...
    %                  [1 0.6 0.2], 'EdgeColor','none','FaceAlpha',0.4);
    % 
    %             plot(t, m, 'k-','LineWidth',1.5);
    % 
    %             xlim([20 160])
    %             ylabel('Fishable crab biomass (kg/ha)')
    %             xlabel('Seasonal time step')
    %             title('Fishable crab densities Mean ± 5–95%')
    % 
    %         % Hiding Urchin 🟣 plot
    %             subplot(2,2,3); hold on
    % 
    %             U = squeeze(sum(ut2(2,:,:),1))';    % R x T
    % 
    %             m   = mean(U,1);
    %             p95 = prctile(U,95,1);
    %             p05 = prctile(U,5,1);
    % 
    %             fill([t fliplr(t)], [p95 fliplr(p05)], ...
    %                  [0.6 0.4 0.8], 'EdgeColor','none','FaceAlpha',0.4);
    % 
    %             plot(t, m, 'k-','LineWidth',1.5);
    % 
    %             xlim([20 160])
    %             ylabel('Adult urchin biomass (kg/ha)')
    %             title('Hiding urchin densities Mean ± 5–95%')
    % 
    %          % All Adult Crab 🦀 plot
    %             subplot(2,2,4); hold on
    % 
    %             C = squeeze(sum(cmt2(9:end,:,:) + cft2(9:end,:,:),1))';  % R x T
    % 
    %             m   = mean(C,1);
    %             p95 = prctile(C,95,1);
    %             p05 = prctile(C,5,1);
    % 
    %             fill([t fliplr(t)], [p95 fliplr(p05)], ...
    %                  [1 0.6 0.2], 'EdgeColor','none','FaceAlpha',0.4);
    % 
    %             plot(t, m, 'k-','LineWidth',1.5);
    % 
    %             xlim([20 160])
    %             ylabel('Adult crab biomass (kg/ha)')
    %             xlabel('Seasonal time step')
    %             title('Adult crab densities Mean ± 5–95%')


%% URCHIN 🟣 & CRAB 🦀 BIOMASS DISTRIBUTION ACROSS REPLICATES (yearly smoothed - loop)

    figure()

    window = 4;                         % 1 year smoothing
    t = (0:size(ut2,2)-1)/4;            % time in years

    % Panel definitions
    titles = { ...
        'Exposed urchin densities Mean ± 5–95%', ...
        'Fishable crab densities Mean ± 5–95%', ...
        'Hiding urchin densities Mean ± 5–95%', ...
        'Adult crab densities Mean ± 5–95%'};

    ylabels = { ...
        'Adult urchin biomass (kg/ha)', ...
        'Fishable crab biomass (N/ha)', ...
        'Adult urchin biomass (kg/ha)', ...
        'Adult crab biomass (N/ha)'};

    colors = { ...
        [0.6 0.4 0.8], ... % urchins
        [1 0.6 0.2], ...   % crabs
        [0.6 0.4 0.8], ...
        [1 0.6 0.2]};

    % Plot
        tiledlayout(2,2)

        for k = 1:4

            nexttile; hold on

            % Stract biomasses
            switch k
                case 1 % Exposed urchins
                    B = squeeze(sum(ut2(3,:,:),1))';

                case 2 % Fishable crabs (male only, ages 17–24)
                    B = squeeze(sum(cmt2(17:24,:,:),1))';

                case 3 % Hiding urchins
                    B = squeeze(sum(ut2(2,:,:),1))';

                case 4 % All adult crabs (female + male)
                    B = squeeze(sum(cmt2(9:end,:,:) + cft2(9:end,:,:),1))';
            end

            % Stats
            m   = movmedian(mean(B,1), window);
            p95 = movmedian(prctile(B,95,1), window);
            p05 = movmedian(prctile(B,5,1), window);

            % Plot
            fill([t fliplr(t)], [p95 fliplr(p05)], ...
                 colors{k}, 'EdgeColor','none','FaceAlpha',0.4);

            plot(t, m, 'k--','LineWidth',1.5);

            % Labels
            xlim([20/4 160/4])
            ylabel(ylabels{k})
            title(titles{k})

            if k > 2
                xlabel('Years')
            end

            grid on
        end


%% URCHIN 🟣 & CRAB 🦀 BIOMASS DISTRIBUTION ACROSS REPLICATES (yearly smoothed)
    % 
    % figure()
    % 
    %     window = 4;                         % 1 year (4 seasonal steps)
    %     t = (0:size(ut2,2)-1)/4;            % time in years
    % 
    %         % Exposed Urchin 🟣 plot
    %             subplot(2,2,1); hold on
    % 
    %             U = squeeze(sum(ut2(3,:,:),1))';    % R x T
    % 
    %             m = movmedian(mean(U,1), window);
    %             p95 = movmedian(prctile(U,95,1), window);
    %             p05 = movmedian(prctile(U,5,1), window);
    % 
    %             fill([t fliplr(t)], [p95 fliplr(p05)], ...
    %                  [0.6 0.4 0.8], 'EdgeColor','none','FaceAlpha',0.4);
    % 
    %             plot(t, m, 'k--','LineWidth',1.5);
    % 
    %             xlim([20/4 160/4])
    %             ylabel('Adult urchin biomass (kg/ha)')
    %             title('Exposed urchin densities Mean ± 5–95%')
    % 
    %         % Fishable Crab 🦀 plot
    %             subplot(2,2,2); hold on
    % 
    %             C = squeeze(sum(cmt2(17:24,:,:),1))';  % R x T
    % 
    %             m = movmedian(mean(C,1), window);
    %             p95 = movmedian(prctile(C,95,1), window);
    %             p05 = movmedian(prctile(C,5,1), window);
    % 
    %             fill([t fliplr(t)], [p95 fliplr(p05)], ...
    %                  [1 0.6 0.2], 'EdgeColor','none','FaceAlpha',0.4);
    % 
    %             plot(t, m, 'k--','LineWidth',1.5);
    % 
    %             xlim([20/4 160/4])
    %             ylabel('Fishable crab biomass (kg/ha)')
    %             xlabel('Years')
    %             title('Fishable crab densities Mean ± 5–95%')
    % 
    %         % Hiding Urchin 🟣 plot
    %             subplot(2,2,3); hold on
    % 
    %             U = squeeze(sum(ut2(2,:,:),1))';    % R x T
    % 
    %             m = movmedian(mean(U,1), window);
    %             p95 = movmedian(prctile(U,95,1), window);
    %             p05 = movmedian(prctile(U,5,1), window);
    % 
    %             fill([t fliplr(t)], [p95 fliplr(p05)], ...
    %                  [0.6 0.4 0.8], 'EdgeColor','none','FaceAlpha',0.4);
    % 
    %             plot(t, m, 'k--','LineWidth',1.5);
    % 
    %             xlim([20/4 160/4])
    %             ylabel('Adult urchin biomass (kg/ha)')
    %             title('Hiding urchin densities Mean ± 5–95%')
    % 
    %         % All Adult Crab 🦀 plot
    %             subplot(2,2,4); hold on
    % 
    %             C = squeeze(sum(cmt2(9:end,:,:) + cft2(9:end,:,:),1))';  % R x T
    % 
    %             m = movmedian(mean(C,1), window);
    %             p95 = movmedian(prctile(C,95,1), window);
    %             p05 = movmedian(prctile(C,5,1), window);
    % 
    %             fill([t fliplr(t)], [p95 fliplr(p05)], ...
    %                  [1 0.6 0.2], 'EdgeColor','none','FaceAlpha',0.4);
    % 
    %             plot(t, m, 'k--','LineWidth',1.5);
    % 
    %             xlim([20/4 160/4])
    %             ylabel('Adult crab biomass (kg/ha)')
    %             xlabel('Years')
    %             title('Adult crab densities Mean ± 5–95%')


%% FUNCTIONAL RESPONSE COMPARISON --------------------------
    % 
    % % 1) Hiding urchins (state 2)
    %     prey_hide = squeeze(uts(2,:,:,:,:,:));
    %     prey_hide = prey_hide(:);
    %     prey_hide = prey_hide(~isnan(prey_hide));
    % 
    % % 2) Exposed urchins (state 3)
    %     prey_exposed = squeeze(uts(3,:,:,:,:,:));
    %     prey_exposed = prey_exposed(:);
    %     prey_exposed = prey_exposed(~isnan(prey_exposed));
    % 
    % % 3) Adult crab biomass (sum of adult classes)
    %     adult_idx = 9:44;   % <-- modify if needed
    %     crab_adult = squeeze(sum(cmt2(adult_idx,:,:),1));  
    %     crab_adult = crab_adult(:);
    %     crab_adult = crab_adult(~isnan(crab_adult));
    % 
    % % Prey ranges
    %     prey1 = linspace(0, max(prey_hide), 500);
    %     prey2 = linspace(0, max(prey_exposed), 500);
    %     prey3 = linspace(0, max(crab_adult), 500);
    % 
    % % Functional responses
    %     Mort1 = -Func_TypeII(urchin.aH, urchin.bH, prey1); % Hiding urchins
    %     Mort2 = -Func_TypeII(urchin.aE, urchin.bE, prey2); % Exposed urchins
    %     Mort3 = -Func_TypeII(crab.aC, crab.bC, prey3); % Adult crabs
    % 
    % % Plotting
    % 
    %     figure()
    % 
    %     % Hiding urchins
    %         subplot(1,3,1)
    %         plot(prey1, Mort1, 'LineWidth', 2)
    %         title('Hiding urchins')
    %         xlabel('Biomass')
    %         ylabel('Per capita mortality')
    %         grid on
    % 
    %     % Exposed urchins
    %         subplot(1,3,2)
    %         plot(prey2, Mort2, 'LineWidth', 2)
    %         title('Exposed urchins')
    %         xlabel('Biomass')
    %         grid on
    % 
    %     % Adult crabs
    %         subplot(1,3,3)
    %         plot(prey3, Mort3, 'LineWidth', 2)
    %         title('Adult crabs')
    %         xlabel('Biomass')
    %         grid on
    % 
    %     sgtitle('Type II Functional Response Across Prey Types')


%% PREY CONSUMED OVER TIME BY SEA OTTERS
    % 
    % % Time (years)
    % t = (1:T2)/4;
    % 
    % % Otter density (T x RR)
    % otter = squeeze(pred_forced(1:T2,:));
    % 
    % % Helper function
    % q = @(x,p) prctile(x,p,2);
    % 
    % % Otter stats (computed once)
    % mO  = mean(otter,2);
    % p5O = q(otter,5);
    % p95O= q(otter,95);
    % 
    % % Labels (order matters!)
    % titles = { ...
    %     'Hiding urchins', ...
    %     'Exposed urchins)', ...
    %     'Female crabs', ...
    %     'Male crabs'};
    % 
    % figure()
    % tiledlayout(2,2)
    % 
    % for k = 1:4
    % 
    %     nexttile; hold on
    % 
    %     % Extract prey consumed (RR x T)
    %     pc = squeeze(prey_consumed2(k,:,:))';
    % 
    %     % Stats
    %     m  = mean(pc,1);
    % 
    %     % Keep your custom quantiles (note: you used different ones in panel 1)
    %         p5  = prctile(pc,5,1);
    %         p95 = prctile(pc,95,1);
    % 
    %     % Ribbon
    %     fill([t fliplr(t)], [p5 fliplr(p95)], ...
    %          [0.7 0.7 1], 'EdgeColor','none','FaceAlpha',0.3)
    % 
    %     % Mean line
    %     plot(t, m, 'k', 'LineWidth',1.5)
    % 
    %     % Secondary axis (otters)
    %     yyaxis right
    %     plot(t, mO, 'r--', 'LineWidth',1.5)
    % 
    %     % Labels
    %     title(titles{k})
    %     xlabel('Time (years)')
    %     yyaxis left
    %     ylabel('Prey consumed (Kg/ha)')
    %     yyaxis right
    %     ylabel('Otter density (Kg/ha)')
    %     grid on
    % end


%% SCATTER: prey consumed vs otter density
    % 
    % % Extract otter (T x RR)
    % otter = squeeze(pred_forced(1:T2,:));
    % 
    % % Flatten otter
    % otter_vec = otter(:);
    % 
    % figure()
    % tiledlayout(2,2)
    % 
    % labels = {'Hiding urchins','Exposed urchins','Female crabs','Male crabs'};
    % 
    % for k = 1:4
    % 
    %     nexttile; hold on
    % 
    %     % Extract prey consumed (T x RR)
    %     pc = squeeze(prey_consumed2(k,:,:)); 
    % 
    %     % Flatten
    %     pc_vec = pc(:);
    % 
    %     % Scatter
    %     scatter(otter_vec, pc_vec, 10, 'filled', ...
    %         'MarkerFaceAlpha', 0.2)
    % 
    %     % Labels
    %     xlabel('Otter density (Kg/ha)')
    %     ylabel('Prey consumed (Kg/ha)')
    %     title(labels{k})
    % 
    %     grid on
    % end