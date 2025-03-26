function pred_forced = ParaPred_Forced_Dataset(ORSO_data,RR)

% Description:
% sets para values (data set) for sea otters portion of the model
% relevant to PredUrchinKelp_ImplicitCC.m


% import the CSV file as a table
otter_data = readtable(ORSO_data, 'VariableNamingRule', 'preserve');

% extract column of interest: 'Average Number' (column C)
average_Number = otter_data.("Average Number"); 

% replicate each year's data 4 times (winter, spring, summer, autmum)
replicated_average_number = repelem(average_Number, 4);

% combine the replicated data into a new table
otter_pop = table(replicated_average_number, ...
    'VariableNames', {'Average_Number'});

% transform table to vector
otter_vector = table2array(otter_pop);

% replicates vector times number of replicates 
pred_forced = repmat(otter_vector,1,RR);



end