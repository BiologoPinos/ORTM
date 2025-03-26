function kelp = ParaKelp_Implicit(tmax)

% Description:
% sets para values for kelp portion of the model
% relevant to PredUrchinKelp_ImplicitCC.m

% para values follow this order
    Species = {'Giant_kelp'; 'Bull_kelp'};

% Recruitment
    % yearly incomming settlers per kg adult standing kelp
    % (zoospore production and successful fertilization and settlement of spore)
    % mean value
    %% Do the math again
    RK = [4*10^4; 4.97*10^4]; %4.97*10^4 

    % temporal stdev (noise)
    % normalised
    RKstdv = [0.389; 0.389]; % 

    % strength of denisity dependence (shading by local adults)
    % ricker -> smaller = weaker effect = higher numbers (more survival)
    % beverton-holt -> bigger = weaker effect = higher number (more survival)
    mu = [9*10^-5; 1*10^10]; %1*10^5  1*10^10
        
        % the relative per capita effect on juvenile survival by adults and juveniles
        % ddD = 0 for ricker (inter-cohort DD)
        % ddD = 1 for beverton-holt (intracohort DD)
        % we'll assume that one adult has the same effect as one juvenile on the survival for giant kelp
        ddD = [0.01; 1]; % 0.5;

        % spatial variance in kelp densitites 
        % (for use in scale transition to landscape)
        % value from PISCO_kelp-urchin-sheephead_data.Rmd
        muvar = [189090; 189090]; %

    % recruitment timing
    % vector that dictates if, and how much, reproduction occurs in that time step
    % [winter, spring, summer, autumn]
    RTk = repmat([0.1 0.1 0.4 0.4;...
                  0 1 0 0],1,tmax/4); %
        
        % recruitment lag
        lag = [0 0; ... % Giant kelp (no lag)
               2 3];    % Bull kelp (vector of lags, 2 and 3 time-steps before spring)
        
        % reproduction weightings
        reproWeight = [0.5 0.5;   % Giant kelp
                       0.9 0.1];  % Bull kelp
                        % 1 .11
            
% Growth/Maturation    
    % growth rate season to season
    %% Do the math again
    g = [6.825; 76.44]; %76.44

% Mortality/survival
    % change in standing kelp biomass; including seasonal variation
    % [winter, spring, summer, autumn]
    %% Change this
    lambda = repmat([1 1 1 1;...
                     1 1 1 1],1,tmax/4); %0.1 1 0.8 0.6; %0 1 0.8 0.6;
  
    % standing (juvenile + adult) kelp retention
    %% Change this
    rS = [0.5688; 0.5688]; % 0.83;       

    % drift production [0,1]
    % (standing kelp mortality)
    % low = more kelp, high = more drift
    c = [0.9; 0.9]; %

    % drift retention [0,1]
    % low = less retained, high = more retained
    rD = [0.7; 0.8]; % 

    % decomposition [0,1]
    % low = slower decomp, high = fast decomp
    d = [0.1; 0.1]; % 0.9;% 


% Grazing
    % attack rates (of urchin stage j on kelp stage i)
    % (p.c. grazing mortality in absence of conspecifics)
    aij = [{[ 0  0.5;
              0  0.5;
              1  0 ]};...
           {[ 0  0.5;
              0  0.5;
              1  0 ]}];     
                 
    % max prey consumed (1/handling time) 
    % (of urchin stage j on kelp stage i)
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
    
        % not seasonal
        % hij = repmat(hij.*reshape([1 1 1 1],1,1,4),1,1,1,tmax/4);
            

%% kelp species table and choice

% join in table
    Paratable = table(RK, RKstdv, mu, ddD, muvar, RTk, lag, reproWeight, g, lambda, rS, c, rD, d, aij, bhij, hij,'RowNames', Species);

% select species
    kelp = Paratable('Bull_kelp',:);

end