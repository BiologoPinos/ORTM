function kelp = ParaKelp_Implicit(tmax, species)

% Description:
    % sets parameter (para) values for kelp portion of the model
    % relevant to ORTM_model_otter.m

% Validate species:
    species = validatestring(species, {'Giant_kelp'; 'Bull_kelp'; 'Bull_kelp_south'; 'Bull_kelp_north'}); 

% Row names:
    Species = {'Giant_kelp'; 'Bull_kelp'; 'Bull_kelp_south'; 'Bull_kelp_north'};

% Recruitment:

    % successful zoo-spore production
        RK = [4e4;             % Giant kelp
              6.21e3;          % Bull kelp
              6.21e3;          % Bull kelp south
              6.21e3];         % Bull kelp north

    % standard deviation of recruitment
        RKstdv = [0.389;           % Giant kelp
                  0.32;            % Bull kelp
                  0.30;            % Bull kelp south
                  0];              % Bull kelp north

    % strength of density dependence
        mu = [9e-5;            % Giant kelp
              2.5e4;           % Bull kelp
              2.5e4;           % Bull kelp south
              2.5e4];          % Bull kelp north

    % relative-per-capita effect on juvenile survival
        ddD = [0.01;            % Giant kelp
               1;               % Bull kelp
               1;               % Bull kelp south
               1];              % Bull kelp north

    % spatial variance in adult kelp densities
        muvar = [189090;          % Giant kelp
                 116535;          % Bull kelp
                 457553           % Bull kelp south
                 1899];           % Bull kelp north

    % recruitment timing function
        RTk = zeros(4, 4 * (tmax/4));
        RTk(1,:) = repmat([0.1 0.1 0.4 0.4], 1, tmax/4); % Giant kelp
        RTk(2,:) = repmat([0 1 0 0], 1, tmax/4);         % Bull kelp
        RTk(3,:) = repmat([0 1 0 0], 1, tmax/4);         % Bull kelp south
        RTk(4,:) = repmat([0 1 0 0], 1, tmax/4);         % Bull kelp north

    % recruitment lag
        lag = [0 0;     % Giant kelp
               2 3;     % Bull kelp
               2 3;     % Bull kelp south
               2 3];    % Bull kelp north

    % reproduction weightings
        reproWeight = [0.5 0.5;   % Giant kelp
                       0.9 0.1;   % Bull kelp
                       0.9 0.1;   % Bull kelp south
                       0.9 0.1];  % Bull kelp north

% Growth:

    % seasonal growth rate
        g = [6.825;     % Giant kelp
             29.65;     % Bull kelp
             29.65;     % Bull kelp south
             29.65];    % Bull kelp north

% Mortality/survival:

    % change in standing biomass
        lambda = zeros(4, 4 * (tmax/4));
        lambda(1,:) = repmat([1 1 1 1], 1, tmax/4);        % Giant kelp
        lambda(2,:) = repmat([0.1 0.8 1 0.9], 1, tmax/4);  % Bull kelp
        lambda(3,:) = repmat([0.1 0.8 1 0.9], 1, tmax/4);  % Bull kelp south
        lambda(4,:) = repmat([0.1 0.8 1 0.9], 1, tmax/4);  % Bull kelp north

    % kelp retention
        rS = [0.5688;    % Giant kelp
              0.6;       % Bull kelp
              0.6;       % Bull kelp south
              0.6];      % Bull kelp north

    % drift production
        c = [0.9;       % Giant kelp
             0.9;       % Bull kelp
             0.9;       % Bull kelp south
             0.9];      % Bull kelp north

    % drift retention
        rD = [0.7;       % Giant kelp
              0.7;       % Bull kelp
              0.7;       % Bull kelp south
              0.7];      % Bull kelp north

    % decomposition
        d = [0.1;       % Giant kelp
             0.1;       % Bull kelp
             0.1;       % Bull kelp south
             0.1];      % Bull kelp north

% Grazing:

    % attack rates (same for all urchins)
        aij_template = [0 0.5;
                        0 0.5;
                        1 0];
        aij = repmat({aij_template}, 4, 1);

    % max prey consumed (same for all urchins)
        bhij_template = [0 2.985;
                         0 2.985;
                         2.985 0];
        bhij = repmat({bhij_template}, 4, 1);

    % handling time with seasonal scaling (same for all urchins)
        scale = reshape([1 0.9 1.15 1.2], 1,1,4);
        hij = cell(4,1);
        for i = 1:4
            hij{i} = repmat(bhij{i} .* scale, 1, 1, 1, tmax/4);
        end

% Build parameter table:

    Paratable = table(RK, RKstdv, mu, ddD, muvar, RTk, lag, reproWeight, ...
                      g, lambda, rS, c, rD, d, aij, bhij, hij, ...
                      'RowNames', Species);

% Return only selected species
    kelp = Paratable(species, :);

    
end


% function kelp = ParaKelp_Implicit(tmax, species)
% 
% % Description:
%     % sets parameter (para) values for kelp portion of the model
%     % relevant to ORTM_model_otter.m
% 
% % Validate species:
%     species = validatestring(species, {'Giant_kelp'; 'Bull_kelp'});
% 
% % Row names:
%     Species = {'Giant_kelp'; 'Bull_kelp'};
% 
% % Recruitment (mean-successful settlers):
% 
%     % successful zoo-spore production, fertilization, and settlement
%         RK = [4*10^4; 6.21*10^3]; % 0;
% 
%     % temporal (norm) standard deviation (noise) of recruits
%         RKstdv = [0.389; 0.32]; % 0;
% 
%     % strength of density dependence 
%         mu = [9*10^-5; 2.5*10^4]; % (Tuning parameter)
% 
%         % relative-per-ca-pita effect on juvenile survival by adults and juveniles
%             ddD = [0.01; 1]; % 0 = inter-cohort DD | 1 = intracohort DD
% 
%         % spatial variance in adult kelp densities
%             muvar = [189090; 116535]; % 0;
% 
%     % recruitment timing function [winter, spring, summer, autumn]
%         RTk = repmat([0.1 0.1 0.4 0.4;...
%                       0 1 0 0],1,tmax/4); % 0;
% 
%         % recruitment lag
%             lag = [0 0; ... 
%                    2 3];    
% 
%         % reproduction weightings
%             reproWeight = [0.5 0.5; ... 
%                            0.9 0.1];    
% 
% % Growth/Maturation:
% 
%     % growth rate season to season (Max growth for OR)
%         g = [6.825; 29.65]; % 0; 
% 
% % Mortality/survival:
% 
%     % change in standing kelp biomass over the season
%         lambda = repmat([1 1 1 1;...
%                          0.1 0.8 1 0.9],1,tmax/4);
% 
%     % kelp retention (0 or 1)
%         rS = [0.5688; 0.6];       
% 
%     % drift production (0 or 1) | low = more kelp, high = more drift
%         c = [0.9; 0.9]; 
% 
%     % drift retention (0 or 1) | low = less retained, high = more retained
%         rD = [0.7; 0.7];
% 
%     % decomposition (0 or 1) | low = slower, high = faster
%         d = [0.1; 0.1];
% 
% % Grazing:
% 
%     % attack rates (urchin stage j on kelp stage i)
%         aij = [{[ 0  0.5;
%                   0  0.5;
%                   1  0 ]};...
%                {[ 0  0.5;
%                   0  0.5;
%                   1  0 ]}];     
% 
%     % max prey consumed (1/handling time) (of urchin stage j on kelp stage i)
%         bhij = [{[ 0  2.985;  
%                    0  2.985;
%                    2.985  0 ]};...
%                 {[ 0  2.985;  
%                    0  2.985;
%                    2.985  0 ]}];
% 
%         % including seasonal variation
%             scale = reshape([1 0.9 1.15 1.2],1,1,4);
%             hij = cell(2,1);
%             for i = 1:2
% 
%                 hij{i} = repmat(bhij{i} .* scale, 1,1,1,tmax/4);
% 
%             end
% 
%         % not including seasonal variation
%             % hij = repmat(hij.*reshape([1 1 1 1],1,1,4),1,1,1,tmax/4);
% 
% % Build Kelp species table:
% 
%     % join in table
%     Paratable = table(RK, RKstdv, mu, ddD, muvar, RTk, lag, reproWeight, ...
%                         g, lambda, rS, c, rD, d, aij, bhij, hij, ...
%                         'RowNames', Species);
% 
%     % select species
%     kelp = Paratable(species, :);
% 
% end