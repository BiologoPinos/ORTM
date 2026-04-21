function urchin = ParaUrchin_Implicit(tmax, species)

% Description:
    % sets parameter (para) values for urchin portion of the model
    % relevant to ORTM_model_otter.m

% Validate species:
    species = validatestring(species, {'Urchins_CA'; 'Urchins_OR'});

% Row names:
    Species = {'Urchins_CA'; 'Urchins_OR'};

% Recruitment (mean-successful settlers):

    % Larval production, dispersal and settlement (assumes open population)
        RU = [3*10^5; 2*10^2]; % (Tuning parameter) %2*10^2
    
    % Temporal (norm) standard deviation (noise) of recruits
        RUstdv = [0.621; 1.89]; % 0;

    % Recruitment distribution
        RUdist = {'normal'; 'log-normal'};

    % Recruitment timing function [winter, spring, summer, autumn]
        RTU = repmat([0.05 0.54 0.36 0.05;...
                      0.2 0.2 0.5 0.1],1,tmax/4);

% Growth/Maturation:

    % Proportion maturing from juvenile to adults
        gJ = [1/(4*2); 1/(4*2)];

% Mortality rates (instantaneous):

    % Natural instantaneous mortality rate
        MJU = [0.1; 0.1]; % juveniles
        MHU = [0.1; 0.1]; % hiding adults
        MEU = [0.1; 0.1]; % exposed adults

    % Per kg mortality of urchins - search/attack rate (relevant for type II)
        aE = [0.615; 0.365]; % exposed 
        aH = aE.*0.5; % hiding

    % Max kg of urchins consumed by a kg of predator (1/handling time) (relevant for type II)
        bE = [1/0.0733; 1/0.0733]; % exposed
        bH = bE.*0.5; % hiding

    % Predator baseline preference (relevant for Yodzis functional response)
        wu = [0.5; 0.5];

    % Predator preference sensitivity to relative abundance of prey (relevant for Yodzis functional response)
        phi = [0; 0];

    % Fishing mortality rate (we are assuming no urchin fishery)
        FU = [0; 0]; %

    % Planktonic larval duration (PLD)
        PLD = [65; 91]; %
        tau = 1-PLD./91; % [1-PLD(1)/91; 1-PLD(2)/91]; % 0;

% Behavioral switching function (hiding <> exposed):

    % inflection point (drift density at which urchin hiding:exposed = 1:1)
        w1 = [1; 1]; % 
    % slope around inflection point
        w2 = [0.5; 0.5]; % 

% Minimum standing kelp (juvs + adults) biomass density threshold in which a barren state is declared:

    % Turns predation and fishing of exposed on/off
        kmin = [1170; 408.33]; %  408.33 % 0; 
        % 0 = events never off (predation on)
        % large number = events always off (predation off)

% Build urchin populations table:

    % join in table
        Paratable = table(RU, RUstdv, RUdist, RTU, gJ, MJU, MHU, MEU, aE, ...
                          aH, bE, bH, wu, phi, FU, PLD, tau, w1, w2, kmin, ...
                            'RowNames', Species);
    
    % select species
        urchin = Paratable(species,:);

end
