function pred = ParaPred_Implicit_v6(tmax)

% Sets para values for kelp portion of the model.
% Relevant to PredUrchinKelp_ImplicitCC_v6.m


% Demographic parameters statistics/data ----------------------
% units in cm or seasonally (3 months)
% where relevant, values have been converted from yearly to seasonally
% references in order of species

% para values follow this order
Species = {'Blue rockfish';'Sheephead'};

% mean recruit size (cm)
% (Nickols et al 2019 | White et al 2021 See IPM_parameters_ChannelIslands.m in github)
RLmu = [7.75; 14.82];

% sd of recruit size (cm)
% (Nickols et al 2019 | White et al 2021 See IPM_parameters_ChannelIslands.m in github)
RLsd = [1.15; 1.83];

% max recruit size (cm)
% (Nickols et al 2019 | White et al 2021 See IPM_parameters_ChannelIslands.m in github)
RLmax = [10; 9];


% instantaneous natural mortality (per 3months) 
% (Dick et al. 2017 Table 51 pg 157 | 
%  Alonzo et al 2004 Stock assessment pg 70)
M = [0.12/4; 0.2/4];

% instantaneous fishing rate (per 3months)
% (Pt Lobos Nickols 2019  | 
%  Alonzo et al 2004 Stock assessment pg 70)
F = [0.19/4; 0.2/4];

% mean size at entry to harvest (cm)
% (Keet et al 2008 & Starr et al 2015 in Nickols et al 2019  | 
%  Alonzo et al 2004 Stock assessment pg 6)
Lfish = [21.03; 30];


% von Bertalanffy growth
    % asymptotic maximum size (cm)
    % (Laidig et al. 2003 in Nickols et al 2019 | 
    %  Hamilton et al 2007 mean of older/unfished data Table 2)
    Linf = [38.15; 57.33];
     
    % instantaneous growth rate (von Bertalanffy growth)(per 3months)
    % (Laidig et al. 2003 in Nickols et al 2019 | 
    %  Alonzo et al 2004 Stock assessment pg 68)
    k = [0.17/4; 0.068/4];

% weight (kg) to length (cm) relationship
    % lenght-weight constant
    % ?? | Alonzo et al 2004
    y = [NaN; 2.7*10^-5];
    
    % length-weigth exponential
    % ?? | Alonzo et al 2004
    z = [NaN; 2.86];

% sd adults size (also sd in size at entry to fishery) (cm)
% (Nickols et al 2019 | Alonzo et al 2004 pg 44)
Lvar = [0.1; 0.11];


% min size for grazing urchins (cm)
% (NaN | Hamilton & Caselle 2015 Selden et al 2017)
Lgraze = [NaN; 20];


% Join in table
Paratable = table(RLmu, RLsd, RLmax, M, F, Linf, k, y, z, Lfish, Lvar, Lgraze, 'RowNames', Species);

% Select species
pred = Paratable('Sheephead',:);


        
% IPM parameters ------------------------
% See White et al (2016) Ecol Apps, specifically appendix S3.

    % number of grids
    pred.meshno = 50; % 100; % 250; % 
    % min mesh size (min length for fish)
    pred.meshmin = 0;
    % max mesh size (larger than any fish is likely to grow)
    pred.meshmax = pred.Linf*2; 
    
    % build mesh (size groupings = lengths)
    pred.x = linspace(pred.meshmin,pred.meshmax,pred.meshno);
    % change in x (mesh/grid size)
    pred.dx = pred.x(2)-pred.x(1);

    % predatory size classes
    pred.Lgraze_ind = find(pred.x>=pred.Lgraze,1);
       

% Settlement/Recrtuiment ------------------------

    % Make pdf vector of settlers/recruits (rho)
    pred.Ro = normpdf(pred.x,pred.RLmu,pred.RLsd);  

    % timing vector of settlement/recrtuiment
    % runs winter, spring, summer, autmum, 
    % Sheephead autum: Warner 1975, Cowen 1990, Cowen et al 1985, Alonzo  et al 2004
    % 0.25 = equal recruitment over the year
    pred.RoT = repmat([0,0,0.05,0.95],1,tmax/4); % repmat([0.25,0.25,0.25,0.25],1,tmax/4); % 


% Open population:
    % mean number of recruits for the year (arbitary/scaling factor)
    % dictates ~mean abundance of unfished population
    % Sheephead: set to PISCO data avg in MPAs in Channel islands in last 10 yrs
    %            = 2.20kg/0.006ha = ~370kg/ha (adults)
    %           (see PISCO_kelp-urchin-sheephead_data.Rmd)
    % 365 w/out noise, 355 w/noise
    pred.Rmean = 355; % 365;%  

    % temporal variability in recruitment (noise)
    % uses normalised (detrended) varaince from PISCO data 
    % (see PISCO_kelp-urchin-sheephead_data.Rmd)
    % set to 0 for no variability
    pred.Rstdv = 0.7188; % 0;%  


end
