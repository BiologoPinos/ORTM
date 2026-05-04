function pred_forced = ParaPred_Forced(ORSO_data,RR)

% Description:
    % sets predator forcing data for the model using scenario replicates
    % ORSO data is expressed in densities (number of individuals per km^2)
        % need to change to number of individuals or biomass per hectare
    % applies sea otter biomass (mean = 32 kg, SD = 4.38 kg)

% Read the CSV
    ORSO_raw_data = readtable(ORSO_data, 'VariableNamingRule', 'preserve');

% Extract data (ignore first column: Years)
    otter_densities = ORSO_raw_data{:, 2:end};

% Validate ORSO replicates
    ORSO_replicates = size(otter_densities, 2);
    if RR > ORSO_replicates
        error('RR (%d) exceeds number of available ORSO replicates (%d).', RR, ORSO_replicates);
    end

% Pair model RR with ORSO replicates
    otter_densities = otter_densities(:, 1:RR);

% Create biomass values (kg per otter) with variation
    biomass = 32 + 4.38 .* randn(size(otter_densities, 1), RR);  % normal distribution (mean ± SD)

% Ensure no negative biomass values (truncate if needed)
    biomass(biomass < 0) = 0;

% Rescale mean otter densities to densities in suitable habitat (Only valid for S6)
    new_otter_densities = otter_densities * (123.96/36.66);
    % ORSO estimates densities to S6 area of 123.96km^2, out of which (roughly) only 
    % 36.66km^2 are actually suitable for otters
    
% Transform otter densities (indv/km^2) to biomass densities (kg/ha)
    scaled_data = (new_otter_densities .* biomass) / 100;

% Replicate Otter annual densities to seasons
    pred_forced = repelem(scaled_data, 4, 1); 

end
