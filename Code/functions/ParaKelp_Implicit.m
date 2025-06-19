function kelp = ParaKelp_Implicit(tmax)

% Description:
% sets parameter (para) values for kelp portion of the model
% relevant to ORTM_model_otter.m

% Parameter values follow this order
    Species = {'Giant_kelp'; 'Bull_kelp'};

% Recruitment (successful settlers)
    % mean zoospore production, successful fertilization, and settlement of spore
    % (annual (mean) incoming settlers per kg adult standing kelp (in the absence of DD))
    RK = [4*10^4; 4.97*10^4]; % 0;

    % temporal (normalized) variation|standard deviation (noise) of recruits
    RKstdv = [0.389; 0.32]; % 3*10^-16;

    % strength of denisity dependence (shading by local adults & juveniles)
    mu = [9*10^-5; 1*10^10]; % Tuning parameter % 0;
        
        % relative-per-capita effect on juvenile survival by adults and juveniles
        % ddD = 0 for ricker (inter-cohort DD)
        % ddD = 1 for beverton-holt (intracohort DD)
        ddD = [0.01; 1]; % 0.5;

        % spatial variance in kelp densities (for use in scale transition to landscape)
        muvar = [189090; 156840]; % 0;

    % recruitment timing function
    % vector that dictates if, and how much, reproduction occurs in that time step
    % [winter, spring, summer, autumn]
    RTk = repmat([0.1 0.1 0.4 0.4;...
                  0 1 0 0],1,tmax/4); % 0;
        
        % recruitment lag
        lag = [0 0; ... % Giant kelp (no lag)
               2 3];    % Bull kelp (2 & 3 time-steps before spring)
        
        % reproduction weightings
        reproWeight = [0.5 0.5; ... % Giant kelp
                       0.9 0.1];    % Bull kelp

% Growth/Maturation 
    % growth rate season to season 
    % (for every 1kg of kelp this season, there will be 1kg*g next season)
    g = [6.825; 29.65]; % 0; 

% Mortality/survival
    % change in standing kelp biomass over the season
    % [winter, spring, summer, autumn]
    lambda = repmat([1 1 1 1;...
                     0.1 0.8 1 0.9],1,tmax/4);
    
    % kelp retention
    % (proportion of standing (juvenile + adult) kelp biomass that remain at the end of the season)
    rS = [0.5688; 0.5688]; % 0;       

    % drift production [0,1]
    % (proportion of standing kelp biomass converting to drift)
    c = [0.9; 0.9]; % low = more kelp, high = more drift

    % drift retention [0,1]
    % (proportion of drift biomass retained locally)
    rD = [0.7; 0.7]; %  low = less retained, high = more retained

    % decomposition [0,1]
    % (proportion of drift decomposing)
    d = [0.1; 0.1]; % low = slower, high = faster


% Grazing
    % attack rates (of urchin stage j on kelp stage i)
    % (p.c. grazing mortality in absence of conspecifics)
    aij = [{[ 0  0.5;
              0  0.5;
              1  0 ]};...
           {[ 0  0.5;
              0  0.5;
              1  0 ]}];     
                 
    % max prey consumed (1/handling time) (of urchin stage j on kelp stage i)
    bhij = [{[ 0  2.985;  
               0  2.985;
               2.985  0 ]};...
            {[ 0  2.985;  
               0  2.985;
               2.985  0 ]}];

        % including seasonal variation
        % [winter spring summer autum]
        scalingFactor = reshape([1 0.9 1.15 1.2], 1, 1, 4);
        
        hij = cell(size(bhij));  % preallocate cell array of same size
       
        for i = 1:numel(bhij)
            hij{i} = repmat(bhij{i} .* scalingFactor, 1, 1, 1, tmax/4);
        end
    
        % % not including seasonal variation
        % hij = repmat(hij.*reshape([1 1 1 1],1,1,4),1,1,1,tmax/4);
            
% Kelp species table and choice
    % join in table
    Paratable = table(RK, RKstdv, mu, ddD, muvar, RTk, lag, ...
                        reproWeight, g, lambda, rS, c, rD, ...
                        d, aij, bhij, hij,'RowNames', Species);

    % select species
    kelp = Paratable('Bull_kelp',:);

end