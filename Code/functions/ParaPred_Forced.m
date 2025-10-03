function pred_forced = ParaPred_Forced(ORSO_data,RR)

% Description:
    % sets predator forcing data for the model using scenario replicates
    % ORSO data is expressed in densities (number of individuals per km^2)
        % need to change to number of individuals or biomass per hectare
    % applies sea otter biomass (mean = 23.02 kg, SD = 4.38 kg)
    % biomass variability reflects a population with a 35% male / 65% female sex ratio

% read the CSV
    otter_data = readtable(ORSO_data, 'VariableNamingRule', 'preserve');

% extract replicate columns (ignore first column: Years)
    replicate_data = otter_data{:, 2:end};

% validate RR
    n_total_reps = size(replicate_data, 2);
    if RR > n_total_reps
        error('RR (%d) exceeds number of available replicates (%d).', RR, n_total_reps);
    end

% select first RR replicates
    replicate_data = replicate_data(:, 1:RR);

% % apply biomass scaling (one value per year per replicate)
%     n_years = size(replicate_data, 1);
% 
% % biomass per individual (kg) with normal variation
%     biomass = 23.02 + 4.38 .* randn(n_years, RR);  % normal distribution (mean ± SD)
% 
% % ensure no negative biomass values (truncate if needed)
%     biomass(biomass < 0) = 0;
% 
% % multiply raw replicate values by biomass per individual
%     scaled_data = (replicate_data ./ 100) .* biomass;
    

% change to otters per hectare (original value is in otters per km^2)
    scaled_data = replicate_data ./ 100;

% repeat each year's row 4 times (for 4 seasons)
    pred_forced = repelem(scaled_data, 4, 1); % final size = (n_years*4 x RR)


end
