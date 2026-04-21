function kelp = ParaKelp_Implicit(tmax, species)

% Description:
    % sets parameter (para) values for kelp portion of the model
    % relevant to ORTM_model_otter.m

% Validate species:
    species = validatestring(species, {'Giant_kelp'; 'Bull_kelp'; 'Bull_kelp_south'; 'Bull_kelp_north'}); 

% Row names:
    Species = {'Giant_kelp'; 'Bull_kelp'; 'Bull_kelp_south'; 'Bull_kelp_north'};

% Recruitment:

    % Successful zoo-spore production
        RK = [4*10^4; 6.21*10^3; 6.21*10^3; 6.21*10^3]; % [Giant; Bull; Bull_S; Bull_N]

    % Standard deviation of recruitment
        RKstdv = [0.389; 0.32; 0.30; 0]; % [Giant; Bull; Bull_S; Bull_N]

    % Strength of density dependence
        RKbeta = [9*10^-5; 3*10^-5; 3*10^-5; 3*10^-5]; % [Giant; Bull; Bull_S; Bull_N]
        % intra-cohort needs larger RKbeta, inter-cohort needs smaller RKbeta

    % Relative-per-capita effect on juvenile survival (Beverton-Holt=1 | Ricker=0)
        D = [0.01; 0.9; 0.9; 0.9]; % [Giant; Bull; Bull_S; Bull_N]
        % closer to 1 intra-cohort, closer to 0 inter-cohort

    % Spatial variance in adult kelp densities
        muvar = [189090; 116535; 457553; 1899]; % [Giant; Bull; Bull_S; Bull_N]

    % Recruitment timing function
        RTk = zeros(4, 4 * (tmax/4));
        RTk(1,:) = repmat([0.1 0.1 0.4 0.4], 1, tmax/4); % Giant kelp
        RTk(2,:) = repmat([0 1 0 0], 1, tmax/4);         % Bull kelp
        RTk(3,:) = repmat([0 1 0 0], 1, tmax/4);         % Bull kelp south
        RTk(4,:) = repmat([0 1 0 0], 1, tmax/4);         % Bull kelp north

    % Recruitment lag
        RKlag = [0  0;     % Giant kelp
                 2  3;     % Bull kelp
                 2  3;     % Bull kelp south
                 2  3];    % Bull kelp north

    % Reproduction weightings
        reproWeight = [0.5  0.5;   % Giant kelp
                       0.9  0.1;   % Bull kelp
                       0.9  0.1;   % Bull kelp south
                       0.9  0.1];  % Bull kelp north

% Growth:

    % Seasonal growth rate
        g = [6.825; 29.65; 29.65; 29.65]; % [Giant; Bull; Bull_S; Bull_N]

% Mortality/survival:

    % Change in standing biomass
        lambda = zeros(4, 4 * (tmax/4));
        lambda(1,:) = repmat([1 1 1 1], 1, tmax/4);        % Giant kelp
        lambda(2,:) = repmat([0.1 0.8 1 0.9], 1, tmax/4);  % Bull kelp
        lambda(3,:) = repmat([0.1 0.8 1 0.9], 1, tmax/4);  % Bull kelp south
        lambda(4,:) = repmat([0.1 0.8 1 0.9], 1, tmax/4);  % Bull kelp north

    % Kelp retention
        rS = [0.5688; 0.5688; 0.5688; 0.5688]; % [Giant; Bull; Bull_S; Bull_N]

    % Drift production
        c = [0.9; 0.9; 0.9; 0.9]; % [Giant; Bull; Bull_S; Bull_N]

    % Drift retention
        rD = [0.7; 0.7; 0.7; 0.7]; % [Giant; Bull; Bull_S; Bull_N]

    % Decomposition
        d = [0.1; 0.1; 0.1; 0.1]; % [Giant; Bull; Bull_S; Bull_N]

% Grazing:

    % Attack rates (same for all urchins)
        aij_template = [0   0.5;
                        0   0.5;
                        1   0];
        aij = repmat({aij_template}, 4, 1);

    % Max prey consumed (same for all urchins)
        bhij_template = [0       2.985;
                         0       2.985;
                         2.985   0];
        bhij = repmat({bhij_template}, 4, 1);

    % Handling time with seasonal scaling (same for all urchins)
        scale = reshape([1 0.9 1.15 1.2], 1,1,4);
        hij = cell(4,1);
        for i = 1:4
            hij{i} = repmat(bhij{i} .* scale, 1, 1, 1, tmax/4);
        end

% Build parameter table:

    Paratable = table(RK, RKstdv, RKbeta, D, muvar, RTk, RKlag, reproWeight, ...
                      g, lambda, rS, c, rD, d, aij, bhij, hij, ...
                      'RowNames', Species);

% Return only selected species
    kelp = Paratable(species, :);

    
end