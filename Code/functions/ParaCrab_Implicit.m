function crab = ParaCrab_Implicit(tmax, species)

% Description:
    % sets parameter (para) values for Dungeness crab portion of the model
    % relevant for ORTM_model_otter.m

% Validate species:
    species = validatestring(species, {'Dungeness_OR_Norm'; 'Dungeness_OR_LogNorm'});

% Row names:
    Species = {'Dungeness_OR_Norm'; 'Dungeness_OR_LogNorm'};
    
% Recruitment:

    % Recruitment (egg production) (Higgins et al., 1997)
        RC = [2000000; 2000000];  % eggs per female in fall

    % Temporal Standard deviation
        RCstdv = [1.26; 1.26]; %% NEED NUMBER

    % Recruitment distribution
        RCdist = {'normal'; 'log-normal'};

    % Recruitment timing function [winter, spring, summer, autumn] (Higgins et al., 1997)
        RTC = repmat([0 1 0 0; ...
                      0 1 0 0],1,tmax/4); % all age-0 recruits arrive in spring

% Mortality:

    % Natural instantaneous mortality rate (Botsford & Hobbs, 1995)
        MRC = [0.2; 0.2]; % new recruits age 0
        MJC = [0.2; 0.2]; % juveniles age 1
        MAC = [0.2; 0.2]; % adult age 2-11
        
    % Fishing rate (Botsford & Hobbs, 1995)
        FC = repmat([0.33 0.33 0.33 0; ...
                     0.33 0.33 0.33 0], 1, tmax/4); 

    % Predation (attack) instantaneous mortality rate of urchins (relevant for type II)
        PC = [0.3650; 0.3650]; %% NEED NUMBER

    % Handling time or max prey consumed "h = 1/handling time" (relevant for type II)
        HC = [1/0.001; 1/0.001]; %% NEED NUMBER

    % Cannibalism 
        % C = % Happens only on newly settle age 0 crabs

    % Planktonic larval duration (PLD)
        % PLD = [182; 182];
        % tauC = 1-PLD./91; 

% Build Dungeness populations table:

    % join in table
        Paratable = table(RC, RCstdv, RCdist, RTC, MRC, MJC, MAC, FC, PC, HC, ...
                  'RowNames', Species);
    
    % select species
        crab = Paratable(species,:);

end