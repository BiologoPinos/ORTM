function urchin = ParaUrchin_Implicit(tmax)

% Description:
% sets parameter (para) values for urchin portion of the model
% relevant to ORTM_model_otter.m

% Parameter values follow this order
    Species = {'Urchins_CA'; 'Urchins_OR'};

% Recruitment
    % Mean larval production (fecundity), dispersal and settlement for urchins (assumes open population)
    RU = [3*10^5; 1.5*10^5]; % Tuning parameter % 3*10^5; % 0;
    
    % Temporal (normalized) standard deviation (noise)
    RUstdv = [0.621; 1.89]; % 0;

    % Recruitment distribution
    RUdist = {'normal'; 'log-normal'};

    % Reproduction timing
    % vector that dictates if reproduction occurs in that time step
    % urchins reproduce 90% January-September (March-May = spring)
    RTu = repmat([0.05 0.54 0.36 0.05;...
                  0.2 0.2 0.5 0.1],1,tmax/4);

% Growth    
    % proportion maturing from juvenile to adults
    % 4 seasons in a year * 2 = 2 yrs maturation
    gJ = [1/(4*2); 1/(4*2)];


% Mortlaity rates (instantaneous)
    % natural mortality
        % juveniles
        MJ = [0.1; 0.1];
            % strength of recruitment facilitation by adults (DD) - NOT BEING USED IN OR VERSION            
            alpha = [0.00001; 0.00001]; % smaller = stronger DD effects
            % variance in adult urchin densitites - NOT BEING USED IN OR VERSION
            alphavar = [5751518; 505]; 

        % hiding adults
        MH = [0.1; 0.1];
        % exposed adults
        ME = [0.1; 0.1];

    % predation rate on urchins by predators (otters) (lower = less mortality)
        % exposed
        PE = [0.0065; 0.3650]; % 0;
        % hiding
        PH = PE.*0.5;
        % PH = [PE(1)*0.5; PE(2)*0.5]; % 0;

    % Fishing rate
        F = [0; 0];   

    % Discount rate for juv survival
        % PLD in days
        PLD = [65; 91];
        tau = 1-PLD./91;
        % tau = [1-PLD(1)/91; 1-PLD(2)/91];


% **behavioural switching function (hiding <> exposed)
    % inflection point (drift density at which urchin hiding:exposed = 1:1)
    w1 = [1; 1]; % 
    % slope around inflection point 
    % (larger values = steeper inflection = sharper transition)
    % (smaller values = flatter line = urchins out all the time)
    w2 = [0.5; 0.5]; % 

    
% **step function (turns predation/fishing of exposed on/off) 
% min density of all standing kelp (juvs + adults)
% 0 = events never off (predation on)
% large number = events always off (predation off)
    kmin = [1170; 1170]; % 0;% 10^10; %  


% urchin population (which US state) table and choice

% join in table
    Paratable = table(RU, RUstdv, RUdist, RTu, gJ, MJ, alpha, alphavar, ...
                        MH, ME, PE, PH, F, PLD, tau, w1, w2, kmin, ...
                        'RowNames', Species);

% select species
    urchin = Paratable('Urchins_OR',:);

end
