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

% set time and length for any scenario (if one occurs)

if contains(scenario,{'cull', 'rest'}) 

    % timing of action start
    % relative to start of disturbance or realtive to the start of the model when no disturbance is happening
    % in terms of seasons (timesteps)
    % before = -x, during = 0, after = x-1 (yr)
    mngt.time =  0; % [-4,0,4]; %  -8:4:20; %         

    % length of mngt action (timesteps)
    mngt.length = 5*4; %  0:4:20; %    [0,4,12,20]; %   


end


% if the scenario contains any of the strings then the action will be
% implimented

if contains(scenario,'cull') 
    % Urchin Removal
    % what seasons does mngt happen in?
    % 1 = winter, 2 = spring, 3 = summer, 4 = fall
    mngt.season = 3; % 1:4; %

    % = reducing urchin biomass by removing a set amount 
    mngt.culling = "Y";  

    % proportion of biomass removed in a season (kg.ha)
    mngt.degreeC = 1*4.2951e3; %  [0,0.01,0.05,0.1,0.25,0.5,0.75,1]; %    

end

if contains(scenario,'rest') 
    % Kelp restoration  
    % what seasons does mngt happen in?
    % 1 = winter, 2 = spring, 3 = summer, 4 = fall
    mngt.season = 2:3; % 2:3; %

    % = reseeding juveniles into the population
    mngt.restore = "Y"; 


    % proportion of biomass of recruits added per season (kg.ha)
    mngt.degreeR =  1*(5.7735e+04 * 0.1); %(4.1791e+03 * 0.01); %6.21*10^5; % (4.1791e+03 * 0.01); % 0.01(or other %)*adult kelp biomass  

end


end 




% 
% function mngt = ParaMngt_Implicit(scenario)
% 
% % Description:
% % sets para values for mngt scenarios portion of the model
% % relevant to PredUrchinKelp_ImplicitCC.m
% 
% 
% % set default scenarios off
% mngt.fish = "N"; 
% mngt.culling = "N"; 
% mngt.restore= "N"; 
% mngt.time = NaN;
% mngt.length = NaN; 
% mngt.degree = NaN;
% 
% % set time and length for any scenario (if one occurs)
% 
% if contains(scenario,{'fish', 'cull', 'rest'}) 
% 
%     % timing of action start
%     % relative to start of disturbance, in terms of seasons (timesteps)
%     % before = -x, during = 0, after = x-1 (yr)
%     mngt.time =  4; % -4 %  [-4,0,4]; %(-8:2:20); %  
% 
%     % length of mngt action (timesteps)
%     mngt.length = 100; % 12 % 0:4:20; %   [0,4,12,20]; %    
% 
% end
% 
% 
% % if the scenario contains any of the strings then the action will be
% % implimented
% 
% if contains(scenario,'fish') 
%     % Temporarily reducing fishing pressure (or temp MPA)
%     mngt.fish = "Y"; 
% 
%     % temporary fishing pressure (baseline for SH = 0.05 per season)
%     mngt.degreeF = 0; % [0,0.05]; % 0:0.005:0.05; %  
% 
% end
% 
% 
% 
% if contains(scenario,'cull') 
%     % Urchin Removal
%     % = reducing urchin biomass by removing a set amount 
%     mngt.culling = "Y";  
% 
%     % biomass removed in a season (kg.ha)
%     % avg biomass, get value from:
%         % model was run with no mngt action, to get biomass
%         % last time step
%     mngt.degreeC = 1*4.2952e3; % 9.5790*10^3; %
% end
% 
% if contains(scenario,'rest') 
%     % Kelp restoration  
%     % = reseeding juveniles into the population
%     mngt.restore = "Y"; 
% 
% 
%   % biomass of recruits added per season (kg.ha)
%     % avg pre-disturbance juv + adult standing kelp biomass (fished state)
%     % get value from:
%         % model was run with disturbance & no mngt action, to get biomass
%         % time-step before disturbance
%         % load("Model outputs\MngtScenarios\Implicitv6a_MngtScen_v1_none_20241011.mat", "kelpJ_avg_pre", "kelpA_avg_pre", "urchinA_avg_pre")
%         % calculate kelpJA_avg_pre = kelpJ_avg_pre + kelpA_avg_pre;
%         % = 1.0848*10^5 kg
%     mngt.degreeR = 1*1.0848*10^5; %  [0,0.01,0.05,0.1,0.25,0.5,0.75,1].*1.0848*10^5; %  [10,25,(50:50:1000)]; % logspace(1,3,20); % 10^4; % [10^6, 10^4, 10^2]; %  
% 
% end
% 
% 
% end 