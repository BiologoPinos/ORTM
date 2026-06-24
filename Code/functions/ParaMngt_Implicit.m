function mngt = ParaMngt_Implicit(scenario)

% Description:
% sets para values for mngt scenarios portion of the model
% relevant to PredUrchinKelp_ImplicitCC.m


% set default scenarios off
mngt.fish = "N"; 
mngt.culling = "N"; 
mngt.restore= "N"; 
mngt.time = NaN;
mngt.length = NaN; 
mngt.degree = NaN;
mngt.season = NaN;

% Set mngt for continuous or periodic (mngt-time-gap)

    % Choose scenario (continuous = 0 vs mngt-time-gap = 1)
        mngt.strategy = 0; % Set to 0 to run without gaps
    
    % Define the gap in time 
        % mngt.time_vec = [0, 1, 2, 3, 8, 9, 10, 11]; % 1 year gap
        % mngt.time_vec = [0, 1, 2, 3, 12, 13, 14, 15]; % 2 year gap
        mngt.time_vec = [0, 1, 2, 3, 16, 17, 18, 19]; % 3 year gap


% Set time and length for any scenario (if one occurs)

if contains(scenario,{'cull', 'rest'}) 

    % Timing of when does the action start (timesteps). This is relative to start of disturbance/mng
        % before = -x, during = 0, after = x-1 (yr)
        mngt.time =  0;        

    % Length of mngt action (timesteps). For how long
        mngt.length = 5*4;    

end


% If the scenario contains any of the strings then the action will be implemented

if contains(scenario,'cull') 

    % Urchin Removal

    % What seasons does mngt happen in?
        % 1 = winter, 2 = spring, 3 = summer, 4 = fall
        mngt.season = 3; % 1:4; %

    % = reducing urchin biomass by removing a set amount 
        mngt.culling = "Y";  

    % Proportion of biomass removed in a season (kg.ha)
        mngt.degreeC = (1*4.6286e+03)/2; % Target is to remove 4.6286e+03 urchins total      

end


if contains(scenario,'rest') 

    % Kelp Restoration  

        % What seasons does mngt happen in?
            % 1 = winter, 2 = spring, 3 = summer, 4 = fall
            mngt.season = 2:3; % 1:4; %

        % = reseeding juveniles into the population
            mngt.restore = "Y"; 

        % Proportion of biomass of recruits added per season (kg.ha)
            mngt.degreeR =  (1*(1.0094e+04 * 4))/2; % Target is to add 4 times more the mean peak spore production (1.0094e+04 * 4)  

end


end 

