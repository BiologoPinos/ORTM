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

% set time and length for any scenario (if one occurs)

if contains(scenario,{'fish', 'cull', 'rest'}) 

    % timing of action start
    % relative to start of disturbance, in terms of seasons (timesteps)
    % before = -x, during = 0, after = x-1 (yr)
    mngt.time =  -4;%  [-4,0,4]; %(-8:2:20); %  

    % length of mngt action (timesteps)
    mngt.length = 12; % 0:4:20; %   [0,4,12,20]; %    

end


% if the scenario contains any of the strings then the action will be
% implimented

if contains(scenario,'fish') 
    % Temporarily reducing fishing pressure (or temp MPA)
    mngt.fish = "Y"; 

    % temporary fishing pressure (baseline for SH = 0.05 per season)
    mngt.degreeF = 0; % [0,0.05]; % 0:0.005:0.05; %  

end

if contains(scenario,'cull') 
    % Urchin Removal
    % = reducing urchin biomass by removing a set amount 
    mngt.culling = "Y";  

     % biomass removed in a season (kg.ha)
    % avg pre-disturbance urchin biomass, get value from:
        % model was run with disturbance & no mngt action, to get biomass
        % time-step before disturbance
        % load("Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_none_20241011.mat", "kelpJ_avg_pre", "kelpA_avg_pre", "urchinA_avg_pre")
        % urchinA_avg_pre = 9.5790*10^3
    mngt.degreeC = 1*9.5790*10^3; %  [0,0.01,0.05,0.1,0.25,0.5,0.75,1].*9.5790*10^3; %    [10,25,(50:50:1000)]; % linspace(1000,5000,5); % 2000; % [20000,2000,200]; % 50 1*10^4; %

end
    
if contains(scenario,'rest') 
    % Kelp restoration  
    % = reseeding juveniles into the population
    mngt.restore = "Y"; 

    
  % biomass of recruits added per season (kg.ha)
    % avg pre-disturbance juv + adult standing kelp biomass (fished state)
    % get value from:
        % model was run with disturbance & no mngt action, to get biomass
        % time-step before disturbance
        % load("Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_none_20241011.mat", "kelpJ_avg_pre", "kelpA_avg_pre", "urchinA_avg_pre")
        % calculate kelpJA_avg_pre = kelpJ_avg_pre + kelpA_avg_pre;
        % = 1.0848*10^5 kg
    mngt.degreeR = 1*1.0848*10^5; %  [0,0.01,0.05,0.1,0.25,0.5,0.75,1].*1.0848*10^5; %  [10,25,(50:50:1000)]; % logspace(1,3,20); % 10^4; % [10^6, 10^4, 10^2]; %  

end


end 