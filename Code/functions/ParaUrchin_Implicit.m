function urchin = ParaUrchin_Implicit(tmax)

% Description:
% sets para values for urchin portion of the model
% relevant to PredUrchinKelp_ImplicitCC.m


% Recruitment
    % **constant incoming settlers (mean value; inc fecundity & dispersal)
    % This value makes all the difference to what state the model is in! 
    % - set to 0 to turn off urchins in model
    % - threshold urchin biomass value = 2.64x10^4 kg.ha, above which urchins dominate) get this in non-MPA scenario with RU = 3.605*10^5  
    % - 1*10^6 gives urchin values seen at Anacapa East, with sheephead present (6.7x10^3 kg.ha), but urchins dominate when sheephead are
    %   under fishing pressure
    % Set so that kg.60m^2 is within the range seen in the Channel Islands
    % (~60 kg.60m2 outside MPAs, ~24 kg.60m2 in MPAs)
    urchin.RU = 1.5*10^5; %1.5*10^5; % 3*10^5; % 0;
    
    %% Tune this value (smaller?)

    % temporal stdv of recruitment (noise)
    urchin.RUstdv = 0.621; % 0;%  

    % Reproduction timing
    % vector that dictates if reproduction occurs in that time step
    % urchins reproduce 90% January-September (March-May = spring)
    urchin.RTu = repmat([0.2 0.2 0.5 0.1],1,tmax/4);


% Growth    
    % proportion maturing from juvenile to adults
    % 4 seasons in a year * 2 = 2 yrs maturation
    urchin.gJ = 1/(4*2);


% Mortlaity rates (instantaneous)
    % natural mortality
        % juveniles
        urchin.MJ = 0.1;
            % strength of DD 
            % smaller = stronger DD effects (survival increases slower with
            % inc adult biomass)
            urchin.alpha = 0.00001;
            % variance in adult urchin densitites 
            % (for use in scale transition to landscape)
            % value from PISCO_kelp-urchin-sheephead_data.Rmd
            urchin.alphavar = 16113; % 

        % hiding adults
        urchin.MH = 0.1;
        % exposed adults
        urchin.ME = 0.1;

    % predation rate on urchins by predators (otters) (lower = less mortality)
        % exposed
        urchin.PE = 0.3650; %(1.46/4) %0.013(Cali); %(1.46/4)
        % hiding   0.1825 %0.0065(Cali)
        urchin.PH = urchin.PE*0.5;
%% ASK this to Jess

    % fishing rate
        urchin.F = 0;   

    % discount rate for juv survival
    % PLD in days
    urchin.PLD = 91; %65(Cali);
    urchin.tau = 1-urchin.PLD/91;


% **behavioural switching function (hiding <> exposed)
    % inflection point (drift density at which urchin hiding:exposed = 1:1)
    urchin.w1 = 1; % 
    % slope around inflection point 
    % (larger values = steeper inflection = sharper transition)
    % (smaller values = flatter line = urchins out all the time)
    urchin.w2 = 0.5; % 

    
% **step function (turns predation/fishing of exposed on/off) 
% min density of all standing kelp (juvs + adults)
% 0 = events never off (predation on)
% large number = events always off (predation off)
    urchin.kmin = 1170;% 0;% 10^10; %  

end