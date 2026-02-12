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
        RU = [3*10^5; 5*10^2]; % (Tuning parameter) %2*10^2
    
    % Temporal (norm) standard deviation (noise) of recruits
        RUstdv = [0.621; 1.26]; % 0;

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
            % strength of recruitment facilitation by adults (DD)            
            alpha = [1*10^-5; 1*10^-5]; % NOT BEING USED IN OR VERSION
            % variance in adult urchin densities
            alphavar = [5751518; 505]; % NOT BEING USED IN OR VERSION

    % Predation (attack) instantaneous mortality rate of urchins (relevant for type II)
        PE = [0.0065; 0.3650]; % exposed (lower = less mortality)
            % 0.3650 -> Burt et al., 2018
            % 0.0027 -> 27.38*0.21/23.02/91
        PH = PE.*0.5; % hiding  (lower = less mortality) % [PE(1)*0.5; PE(2)*0.5]; % 0;

    % Handling time or max prey consumed "h = 1/handling time" (relevant for type II)
        HE = [1/0.001; 1/0.001]; % exposed
            % 1/0.001
            % 1/1.25e-4 -> 30*0.21/23.02/24/91 
        HH = HE.*0.5; % hiding

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
        Paratable = table(RU, RUstdv, RUdist, RTU, gJ, MJU, alpha, alphavar, ...
                            MHU, MEU, PE, PH, HE, HH, FU, PLD, tau, w1, w2, kmin, ...
                            'RowNames', Species);
    
    % select species
        urchin = Paratable(species,:);

end
