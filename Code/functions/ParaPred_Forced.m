function pred_forced = ParaPred_Forced(ORSO_data,RR)

% Description:
% sets predator forcing data for the model using scenario replicates.
% applies sea otter biomass (29kg ± 5–15kg noise) **once per year**
% and repeats each year's value 4 times (for 4 seasons).

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

% apply biomass scaling (one value per year per replicate)
n_years = size(replicate_data, 1);
noise = 5 + (15 - 5) .* rand(n_years, RR);  % uniform biomass noise
biomass = 29 + noise;                       % size = (n_years x RR)

% Multiply raw replicate values by biomass
scaled_data = replicate_data .* biomass;

% Repeat each year's row 4 times (for 4 seasons)
pred_forced = repelem(scaled_data, 4, 1); % final size = (n_years*4 x RR)


end

% % Description:
% % sets para values (data set) for sea otters portion of the model
% % relevant to ORTM_model_otter.m
% 
% 
% % import the CSV file as a table
% otter_data = readtable(ORSO_data, 'VariableNamingRule', 'preserve');
% 
% % extract column of interest: 'Average Number' (column C)
% average_Number = otter_data.("Average Number"); 
% 
% % replicate each year's data 4 times (winter, spring, summer, autmum)
% replicated_average_number = repelem(average_Number, 4);
% 
% % combine the replicated data into a new table
% otter_pop = table(replicated_average_number, ...
%     'VariableNames', {'Average_Number'});
% 
% % transform table to vector
% otter_vector = table2array(otter_pop);
% 
% % replicates vector times number of replicates 
% pred_forced = repmat(otter_vector,1,RR);
% 
% 
% 
% end