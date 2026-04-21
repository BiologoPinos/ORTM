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
        RC = [632; 632];

    % Temporal Standard deviation
        RCstdv = [1.26; 1.26]; %% NEED NUMBER

    % Recruitment distribution
        RCdist = {'normal'; 'log-normal'};

    % Recruitment timing function [winter, spring, summer, autumn] (Higgins et al., 1997)
        RTC = repmat([0 0.5 0.5 0; ...
                      0 0.5 0.5 0],1,tmax/4);

        beta = [1.58e-3; 1.58e-3];

% Mortality:

    % Natural instantaneous mortality rate (Botsford & Hobbs, 1995)
        MRC = [0.05; 0.05];
        MJC = [0.05; 0.05];
        MAC = [0.05; 0.05];
        
    % Fishing rate (Botsford & Hobbs, 1995)
        FC = repmat([0.33 0.33 0.33 0; ...
                     0.33 0.33 0.33 0], 1, tmax/4); 

    % Per kg mortality of crabs - search/attack rate (relevant for type II)
        aC = [0.615; 0.275]; 

    % Handling time or max prey consumed "h = 1/handling time" (relevant for type II)
        bC = [1/0.0545; 1/0.0545];

    % Predator baseline preference (relevant for Yodzis functional response)
        wc = [0.3; 0.5];

    % Normalized Cannibalism influence function
        NCIF = {[0.0028; 0.0028; 0.0028; 0.0028; 
                 0.0315; 0.0315; 0.0315; 0.0315;
                 0.0516; 0.0516; 0.0516; 0.0516; 
                 0.0772; 0.0772; 0.0772; 0.0772;
                 0.1085; 0.1085; 0.1085; 0.1085;
                 0.1457; 0.1457; 0.1457; 0.1457; 
                 0.1457; 0.1457; 0.1457; 0.1457; 
                 0.1457; 0.1457; 0.1457; 0.1457; 
                 0.1457; 0.1457; 0.1457; 0.1457; 
                 0.1457; 0.1457; 0.1457; 0.1457];
                [0.0028; 0.0028; 0.0028; 0.0028; 
                 0.0315; 0.0315; 0.0315; 0.0315;
                 0.0516; 0.0516; 0.0516; 0.0516; 
                 0.0772; 0.0772; 0.0772; 0.0772;
                 0.1085; 0.1085; 0.1085; 0.1085;
                 0.1457; 0.1457; 0.1457; 0.1457; 
                 0.1457; 0.1457; 0.1457; 0.1457; 
                 0.1457; 0.1457; 0.1457; 0.1457; 
                 0.1457; 0.1457; 0.1457; 0.1457; 
                 0.1457; 0.1457; 0.1457; 0.1457]};

% Build Dungeness populations table:

    % join in table
        Paratable = table(RC, RCstdv, RCdist, RTC, beta, MRC, MJC, MAC, FC, aC, bC, wc, NCIF,  ...
                  'RowNames', Species);
    
    % select species
        crab = Paratable(species,:);

end