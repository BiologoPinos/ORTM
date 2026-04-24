%% Life history parameters for Dungeness crab (Hobbs & Botsford 1989)
A = 1:10; % Age classes


%% Cannibalism influence function by age class (Hobbs & Botsford 1989)

    % carapace width (mm) of each age class 
        Wid = [40 120 150 180 210 240 240 240 240 240]; 
    
    % mass as a function of width
        Wgt = Wid.^2.76; 
    
    % Influence function of adult biomass on recruit cannibalism
        Infl_c = Wgt.^0.8; 
    
    % Normalizedwg
        Infl_c_n = Infl_c./sum(Infl_c); % SUM = 1


%% Recruitment at steady state

    % Mortality terms
        S = [0.2, 0.2, 0.2, 0.2, 0,2, 0.2, 0.2, 0.2, 0.2]; 
        % S = [0.0200, 0.0200, 0.0200, 0.0200, 0.0200, 0.0200, 0.0200, 0.0200, 0.0200, 0.0200];

        Sc = [1; cumprod(S(:))];

    % Set recruitment
        R = 1000; % 2000000; % number of recruits
        % R = R * 0.000316; % recruits mass in kg
        N = sum(R*Sc);

        R_Sc = R*Sc;

        % This is what we are aiming for:
            % Kg per age class: [0.1582 3.6307 45.1850 41.8250 23.0598 7.0576 6.1216 4.0811 2.0405 0 0]
            % Kg of the population: 133.15


%% Value of Beta for ricker type cannibalism function
    
    % Stable state for K
        k = -1;

    % Sum of the influence function
        phi = sum(Infl_c_n);

    % Kg of incoming recruits
        % R = 632; % Kg of incoming recruits
    
    % Solve for beta:
        b = -(k/(phi*R));


%% Initial conditions for crabs

    % Normal population in a hectare distributed by age 
        cft0 = [500, 250, 150, 75, 25, 5, 3, 2, 1, 0, 0]; % number
        cmt0 = [500, 250, 150, 75, 25, 0, 0, 0, 0, 0, 0]; % number
    
    % Carapace width (mm) per age class
        Wid = [10 40 120 150 180 210 240 240 240 240 240];

    % Individual mass (arbitrary mass units)
        Wgt = (Wid.^2.76)*0.00055; % unites of weight in Kg

    % Convert numbers to total biomass per hectare (Initial conditions)
        cfw0 = cft0 .* Wgt; % female biomass per hectare
        cmw0 = cmt0 .* Wgt; % male biomass per hectare

        c_w = cfw0+cmw0;
        c_wT = sum(c_w);    
