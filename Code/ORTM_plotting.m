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
    %         disp(mean(cft2(1:44,161,:), 3)); % all age classes x season
    %         % disp(mean(sum(reshape(cft2(1:44, 161, :),4,11,[]),1),3)); % grouped by every 4 rows
    % 
    %     % Crabs males 🦀
    %         disp(mean(cmt2(1:44,161,:), 3));
    %         % disp(mean(sum(reshape(cmt2(1:44, 161, :),4,11,[]),1),3)); % grouped by every 4 rows
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
    % 
    % % Plot
    %     figure()
    % 
    %     % 1) Kelp 🌿 (juvenile + adult only)
    %         subplot(3,1,1)
    %         histogram(mean(sum(kts(1:2,(end-80):end,1,1,kelp_avg>0)),2)*0.006,20)
    %         xlabel('Kelp biomass (kg)')
    % 
    %     % 2) Urchins 🟣 (hiding + exposed only)
    %         subplot(3,1,2)
    %         histogram(mean(sum(uts(2:3,(end-80):end,1,1,kelp_avg>0)),2)*0.006,20)
    %         xlabel('Urchin biomass (kg)')
    % 
    %     % 3) Crabs total 🦀 (males + females, adult classes only)
    %         subplot(3,1,3)
    %         histogram(mean(sum(cfts(9:44,(end-80):end,1,1,kelp_avg>0)) + ...
    %                        sum(cmts(9:44,(end-80):end,1,1,kelp_avg>0)),2)*0.006,20)
    %         xlabel('Crab biomass (kg)')



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

    % Plot
        figure(1)
        hold on

        % Over seasons? (note: baseline always at the end)
            % plot((1:T2+1), sum(kt2(2, :, :)>0, 3)./RR, 'k', 'LineWidth', 2) % baseline
            % plot((1:T2+1), sum(kt2(2, :, :)>0, 3)./RR, 'b', 'LineWidth', 2) % management (note: choose color)

        % Over years
            plot(((0:T2)/4), sum(kt2(2, :, :)>0, 3)./RR, 'y', 'LineWidth',1) % baseline
            % plot(((0:T2)/4), sum(kt2(2, :, :)>0, 3)./RR, 'b', 'LineWidth',1) % management (note: choose color)

        % When does management starts?
            % xline(dist.yrs(1)/4,'--','Start Mng','Color','k','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','FontSize',11) % start of management (years)
            % xline(dist.yrs(1),'--','Start Mng','Color','k','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','FontSize',11) % start of management (seasons)
            % xline(buffer/4,'--','Sea Otter Reintroduction','Color','k','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','FontSize',11) % Sea otter reintroduction (years)

        % Limits
            xlim([0,40])  
            ylim([0,1])
            % ylim([0,RR])

        % Axis
            ylabel('Proportion of Kelp Forests Persisting','FontSize',16)
            xlabel('Years','FontSize',16)
            % xlabel('Timesteps (seasons)')          



%% Urchin 🟣 and crab 🦀 biomass distribution across replicates (raw - lots of peaks)
    % 
    % % Plot
    %     figure()
    %     hold on
    % 
    % % Time vector directly from ut2
    %     t = 0:size(ut2,2)-1;
    % 
    % % Urchin 🟣 plot
    % 
    %     % Urchin subplot
    %         subplot(2,1,1)
    %         hold on
    % 
    %     % Create a 5–95% envelope
    %         fill([t fliplr(t)], ...
    %              [squeeze(prctile(sum(ut2(2:3,:,:),1),95,3)) ...
    %               fliplr(squeeze(prctile(sum(ut2(2:3,:,:),1),5,3)))], ...
    %              [0.8 0.8 0.8], 'EdgeColor','none','FaceAlpha',0.4);
    % 
    %     % Plot the Mean
    %         plot(t, squeeze(mean(sum(ut2(2:3,:,:),1),3)), ...
    %              'k-','LineWidth',1.5);
    % 
    %     % Limits
    %         xlim([20,160])  
    % 
    %     % Axis
    %         xlabel('Seasonal time step');
    %         ylabel('Adult urchin biomass (kg/ha)');
    %         title('Adult urchin biomass: Mean ± 5–95% envelope');
    %         % box on; grid on;
    % 
    % % Crab 🦀 plot
    % 
    %     % crab subplot
    %         subplot(2,1,2)
    %         hold on
    % 
    %     % Create a 5–95% envelope
    %         fill([t fliplr(t)], ...
    %              [squeeze(prctile(sum(cmt2(13:24,:,:),1),95,3)) ...
    %               fliplr(squeeze(prctile(sum(cmt2(13:24,:,:),1),5,3)))], ...
    %              [0.8 0.8 0.8], 'EdgeColor','none','FaceAlpha',0.4);
    % 
    %     % Plot the Mean
    %         plot(t, squeeze(mean(sum(cmt2(13:24,:,:),1),3)), ...
    %              'k-','LineWidth',1.5);
    % 
    %     % Limits
    %         xlim([20,160])  
    % 
    %     % Axis
    %         xlabel('Seasonal time step');
    %         ylabel('Adult fishable crab biomass (kg/ha)');
    %         title('Adult fishable crab biomass: Mean ± 5–95% envelope');
    %         % box on; grid on;



%% Urchin 🟣 and crab 🦀 biomass distribution across replicates (yearly smoothed)

    % Plot
        figure ()
        hold on
        window = 4;   % 1 year (4 seasonal steps)
    
    % Time vector directly from ut2
            t = (0:size(ut2,2)-1)/4;

    % Exposed Urchin 🟣 plot

        % Urchin subplot
            subplot(2,2,1)
            hold on

        % Raw statistics
            meanUE = squeeze(mean(sum(ut2(3,:,:),1),3));
            p95U = squeeze(prctile(sum(ut2(3,:,:),1),95,3));
            p05U = squeeze(prctile(sum(ut2(3,:,:),1),5,3));

        % Yearly smoothing
            meanUE = movmedian(meanUE, window);
            p95U = movmedian(p95U, window);
            p05U = movmedian(p05U, window);

        % Force row vectors (fix for fill)
            meanUE = meanUE(:)'; 
            p95U = p95U(:)'; 
            p05U = p05U(:)';

        % Envelope
            fill([t fliplr(t)], ...
                 [p95U fliplr(p05U)], ...
                 [0.8 0.8 0.8], 'EdgeColor','none','FaceAlpha',0.4);

        % Mean
            plot(t, meanUE, 'k-','LineWidth',1.5);

        % Limits
            xlim([20/4 160/4])
            
        % Axis
            ylabel('Adult urchin biomass (kg/ha)');
            % xlabel('Seasonal time step');
            title('Exposed urchin densities Mean ± 5–95%');
            % grid on;

    % Hiding Urchin 🟣 plot

        % Urchin subplot
            subplot(2,2,3)
            hold on

        % Raw statistics
            meanUH = squeeze(mean(sum(ut2(2,:,:),1),3));
            p95U = squeeze(prctile(sum(ut2(2,:,:),1),95,3));
            p05U = squeeze(prctile(sum(ut2(2,:,:),1),5,3));

        % Yearly smoothing
            meanUH = movmedian(meanUH, window);
            p95U = movmedian(p95U, window);
            p05U = movmedian(p05U, window);

        % Force row vectors (fix for fill)
            meanUH = meanUH(:)'; 
            p95U = p95U(:)'; 
            p05U = p05U(:)';

        % Envelope
            fill([t fliplr(t)], ...
                 [p95U fliplr(p05U)], ...
                 [0.8 0.8 0.8], 'EdgeColor','none','FaceAlpha',0.4);

        % Mean
            plot(t, meanUH, 'k-','LineWidth',1.5);

        % Limits
            xlim([20/4 160/4])
            
        % Axis
            ylabel('Adult urchin biomass (kg/ha)');
            % xlabel('Seasonal time step');
            title('Hiding urchin densities Mean ± 5–95%');
            % grid on;


    % Fishable Crab 🦀 plot

        % crab subplot
            subplot(2,2,2)
            hold on

        % Raw statistics
            meanC = squeeze(mean(sum(cmt2(17:24,:,:),1),3));
            p95C = squeeze(prctile(sum(cmt2(17:24,:,:),1),95,3));
            p05C = squeeze(prctile(sum(cmt2(17:24,:,:),1),5,3));
        
        % Yearly smoothing
            meanC = movmedian(meanC, window);
            p95C = movmedian(p95C, window);
            p05C = movmedian(p05C, window);
        
        % Force row vectors
            meanC = meanC(:)';
            p95C = p95C(:)';
            p05C = p05C(:)';

        % Envelope
            fill([t fliplr(t)], ...
                 [p95C fliplr(p05C)], ...
                 [0.8 0.8 0.8], 'EdgeColor','none','FaceAlpha',0.4);
        
        % Mean
            plot(t, meanC, 'k-','LineWidth',1.5);
        
        % Limits
            xlim([20/4 160/4])

        % Axis
            ylabel('Fishable crab biomass (kg/ha)');
            xlabel('Years');
            title('Fishable crab densities Mean ± 5–95%');
            % grid on;

    % All adult Crab 🦀 plot

        % crab subplot
            subplot(2,2,4)
            hold on

        % Raw statistics
            meanC = squeeze(mean(sum(cmt2(9:end,:,:) + cft2(9:end,:,:),1),3));
            p95C = squeeze(prctile(sum(cmt2(9:end,:,:) + cft2(9:end,:,:),1),95,3));
            p05C = squeeze(prctile(sum(cmt2(9:end,:,:) + cft2(9:end,:,:),1),5,3));
        
        % Yearly smoothing
            meanC = movmedian(meanC, window);
            p95C = movmedian(p95C, window);
            p05C = movmedian(p05C, window);
        
        % Force row vectors
            meanC = meanC(:)';
            p95C = p95C(:)';
            p05C = p05C(:)';

        % Envelope
            fill([t fliplr(t)], ...
                 [p95C fliplr(p05C)], ...
                 [0.8 0.8 0.8], 'EdgeColor','none','FaceAlpha',0.4);
        
        % Mean
            plot(t, meanC, 'k-','LineWidth',1.5);
        
        % Limits
            xlim([20/4 160/4])

        % Axis
            ylabel('Adult crab biomass (kg/ha)');
            xlabel('Years');
            title('Adult crab densities Mean ± 5–95%');
            % grid on;


%% FUNCTIONAL RESPONSE COMPARISON --------------------------

% 1) Hiding urchins (state 2)
    prey_hide = squeeze(uts(2,:,:,:,:,:));
    prey_hide = prey_hide(:);
    prey_hide = prey_hide(~isnan(prey_hide));

% 2) Exposed urchins (state 3)
    prey_exposed = squeeze(uts(3,:,:,:,:,:));
    prey_exposed = prey_exposed(:);
    prey_exposed = prey_exposed(~isnan(prey_exposed));

% 3) Adult crab biomass (sum of adult classes)
    adult_idx = 9:44;   % <-- modify if needed
    crab_adult = squeeze(sum(cmt2(adult_idx,:,:),1));  
    crab_adult = crab_adult(:);
    crab_adult = crab_adult(~isnan(crab_adult));

% Prey ranges
    prey1 = linspace(0, max(prey_hide), 500);
    prey2 = linspace(0, max(prey_exposed), 500);
    prey3 = linspace(0, max(crab_adult), 500);

% Functional responses
    Mort1 = -Func_TypeII(urchin.aH, urchin.bH, prey1); % Hiding urchins
    Mort2 = -Func_TypeII(urchin.aE, urchin.bE, prey2); % Exposed urchins
    Mort3 = -Func_TypeII(crab.aC, crab.bC, prey3); % Adult crabs

% Plotting

    figure()

    % Hiding urchins
        subplot(1,3,1)
        plot(prey1, Mort1, 'LineWidth', 2)
        title('Hiding urchins')
        xlabel('Biomass')
        ylabel('Per capita mortality')
        grid on

    % Exposed urchins
        subplot(1,3,2)
        plot(prey2, Mort2, 'LineWidth', 2)
        title('Exposed urchins')
        xlabel('Biomass')
        grid on

    % Adult crabs
        subplot(1,3,3)
        plot(prey3, Mort3, 'LineWidth', 2)
        title('Adult crabs')
        xlabel('Biomass')
        grid on

    sgtitle('Type II Functional Response Across Prey Types')


%% 
% Example: choose one scenario index
h = 1; i = 1; j = 1;

t = 1:T2;

figure
tiledlayout(4,1,'TileSpacing','compact','Padding','compact')

labels = {'Hiding urchins', 'Exposed urchins', 'Female crabs', 'Male crabs'};

for k = 1:4
    nexttile
    
    % ----- LEFT AXIS (prey consumption) -----
    yyaxis left
    hold on
    
    fill([t fliplr(t)], ...
         [prctile(squeeze(prey_consumed(k,:,:,:,:,:)),25,2)' ...
          fliplr(prctile(squeeze(prey_consumed(k,:,:,:,:,:)),75,2)')], ...
         [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha',0.6)
    
    plot(t, mean(squeeze(prey_consumed(k,:,:,:,:,:)),2,'omitnan'), ...
         'k', 'LineWidth', 1.8)
    
    ylabel(labels{k})
    
    % ----- RIGHT AXIS (predators) -----
    yyaxis right
    hold on
    
    plot(t, mean(pred_forced(1:T2,:),2,'omitnan'), ...
         '--', 'LineWidth', 1.5)
    
    ylabel('Sea otters')
    
    % ----- Formatting -----
    if k == 1
        title('Realized prey consumption vs predator density')
    end
    if k == 4
        xlabel('Season')
    end
end


%%
pred = squeeze(pred_forced(1:T2,:));              % T x RR
pc   = squeeze(prey_consumed(:,:,h,i,j,:));       % 4 x T x RR

labels = {'Hiding urchins', 'Exposed urchins', 'Female crabs', 'Male crabs'};

figure
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

for k = 1:4
    nexttile
    
    y = squeeze(pc(k,:,:));   % T x RR
    x = pred;                 % T x RR
    
    % Time index for coloring
    c = repmat((1:T2)', 1, RR);   % T x RR
    
    % Flatten
    x = x(:);
    y = y(:);
    c = c(:);
    
    % Remove NaNs only (keep zeros — they are informative)
    idx = ~isnan(x) & ~isnan(y);
    x = x(idx);
    y = y(idx);
    c = c(idx);
    
    % Scatter with time coloring
    scatter(x, y, 10, c, 'filled', 'MarkerFaceAlpha', 0.25); hold on
    
    xlabel('Sea otter density')
    ylabel('Prey consumed')
    title(labels{k})
end

colorbar