% ORTM_model_otter.m
    % June 2025
 
% Authorship: 
    % Andrés Pinos-Sánchez | andres.pinos.sanchez@gmail.com
    % Co-authors: Jess Hopf, Leif Rasmuson, Mark Novak, Will White

% Toolboxes required:
    % Statistics and Machine Learning Toolbox
    % Optimization Toolbox

% Model description:
    % Discrete time
    % Spatially implicit
% This version aims to look at effectiveness of different management


%% TOP-LEVEL STUFF --------------------------

% clear all
    clear   % clear Workspace
    clc     % clear command window
    
% add paths
    addpath('functions\', 'functions\cbrewer')
    addpath('Model outputs\')
    addpath('Model inputs')
    ORSO = './Model inputs';

% set random number
    rng(1) % randomness (stochasticity)

% start stopwatch timer (elapsed time)
    tic 
      

%% PREP MODEL (USER CHOOSES) --------------------------

% Choose your kelp
    % options: 'Giant_kelp' | 'Bull_kelp' | 'Bull_kelp_north' | 'Bull_kelp_south'
    kelp_species = 'Bull_kelp';

% Choose your urchins
    % options: 'Urchins_CA' | 'Urchins_OR'
    urchin_species = 'Urchins_OR';

% Choose your crabs (second prey source)
    % options: 'Dungeness_OR_Norm' | 'Dungeness_OR_LogNorm'
    crab_species = 'Dungeness_OR_LogNorm';

% Choose your sea otters
    % options: 'Table0.csv' | 'Scenario-Norm_SuccesSegment-N3.csv' - Pacific City | 'Scenario-Norm_SuccesSegment-C7.csv' - Newport | 'Scenario-Norm_SuccesSegment-S6.csv' - Port Orford
    otter_population = 'Scenario-Norm_SuccesSegment-N3.csv';    

% Choose your model (this determines which model to run, CA vs OR)
    % options: 'California_Jess' | 'Oregon_Andres'
    model = 'Oregon_Andres';


%% MODEL PARAMETERS --------------------------

% Run-times/time-steps [winter, spring, summer, autumn]
    T1 = 30*4; % for simulated predator runs
    T2 = 40*4; % for kelp-urchin runs

% Number of replicates (RR) - Max 10000 due to ORSO 
    RR =  1000;

% Length of run + buffer
    tmax = T2+100;

% Kelp parameters
    kelp = ParaKelp_Implicit(tmax, kelp_species);
        kelp.mu = 5*10^-5; % 2.5*10^-9; % 3*10^4; % 2.5*10^4; % 1*10^10; % 9*10^5; % kelp.mu;
        kelp.ddD = 0.9;

% Urchin parameters
    urchin = ParaUrchin_Implicit(tmax, urchin_species);
        urchin.RU = 5*10^2; % 8*10^2; % 2*10^2; % 1.5*10^5; % urchin.RU;

% Crab (Dungeness) parameters
    crab = ParaCrab_Implicit(tmax, crab_species);
        % crab.RC = 0; % Potential tuning parameter

% Predator simulated (sheep-head, CA_model)
    pred = ParaPred_Implicit(tmax); % NEED REVISSION

% Predator forced inputs (sea otters, OR_model)
    ORSO_data = fullfile(ORSO, otter_population);
            
    
%% DISTURBANCE (dist) -------------------------- 
    
% Disturbance length
    dist.lngth = 0; % 1; % 2; 

% Disturbance timing 
    if dist.lngth == 0
        dist.yrs = NaN; % no disturbance (default)
    else    
        dist.yrs = (20*4) + repmat(1:4,1,dist.lngth) + repelem(((1:dist.lngth)-1)*4,4);
    end 
       
% How do vital rates change during the disturbance (heatwave) 
    % kelp recruitment
    dist.RK = kelp.RK/7; 
    % Kelp biomass growth reduction during disturbance
    dist.lambda = kelp.lambda.* repmat([1 1 0.5 0.5],1,tmax/4);   
    % urchin grazing increase during disturbance
    dist.hij = repmat(cell2mat(kelp.bhij) .* reshape([1.15 1.05 1.2 1.3],1,1,4), 1, 1, 1, tmax/4);


%% MANAGEMENT SCENARIOS --------------------------

% MPA scenario (remove fishing pressure)🪝:
    pred.F = 0; % 0.1/4;

% Management (mngt) scenarios:  
    mngt_scen =  'none'; % 'none'; % 'restoration'; % 'fishing'; % 'culling'; 
    % get vector values
        mngt = ParaMngt_Implicit(mngt_scen);
        pred.fish = mngt.fish;
        urchin.culling = mngt.culling;
        kelp.restore = mngt.restore;


%% INITIAL CONDITIONS --------------------------

% KELP 🌿 [juvenile, adult, drift] 
    if strcmp(kelp_species, 'Giant_kelp')
        kt0 = [1.17*10^5, 1.17*10^5, 1.17*10^6];  % Giant-kelp
    else
        kt0 = [8*10^4, 8*10^4, 8*10^5];           % Bull-kelp variants
    end

% URCHINS 🟣 [juvenile, hiding, expose]
    ut0 = [0,0,0]; 

% CRABS 🦀 [age_0, age_1-10]
    cft0 = [0, 0, 0, 0, 0, 0];
    cmt0 = [0, 0, 0, 0, 0, 0];
    
    % cft0 = [400, 0, 0, 0, 0, 0];
    % cmt0 = [400, 0, 0, 0, 0, 0];

    % cft0 = [800, 200, 100, 50, 30, 20];
    % cmt0 = [800, 200, 100, 50, 30, 20];
   

%%  PRE-ASSIGN VARIABLES (empty vectors) --------------------------

% Degree (deg) lenght - degree is the last mngt in the matrix
    deglngth = structfun(@numel,mngt);
    deglngth = deglngth(end);

% Kelp biomass 
    kts = NaN(3, T2+1, length(mngt.time), length(mngt.length), deglngth, RR);

% Urchin biomass
    uts = NaN(3, T2+1, length(mngt.time), length(mngt.length), deglngth, RR);

% Crab biomass
    cfts = NaN(6, T2+1, length(mngt.time), length(mngt.length), deglngth, RR);
    cmts = NaN(6, T2+1, length(mngt.time), length(mngt.length), deglngth, RR);

% Average kelp persistence
    kelp_avg = NaN(length(mngt.time), length(mngt.length), deglngth, RR);

% Predator input data
    if strcmp(model,'California_Jess')
        PBE = NaN(T2, length(mngt.time), length(mngt.length), deglngth, RR);
    else
        pred_forced = ParaPred_Forced(ORSO_data, RR);
        buffer = 5*4; % buffer allows kelp–urchin to stabilize
        pred_forced = [zeros(buffer, RR); pred_forced];
    end   
    
    
%% RUN MODELS --------------------------

% Run managements:
    % Run over timing of mngt action
        for h = 1:length(mngt.time)        
            if contains(mngt_scen,'fish'); pred.fishtime = mngt.time(h)+T1; end
            if contains(mngt_scen,'cull'); urchin.culltime = mngt.time(h); end  
            if contains(mngt_scen,'rest'); kelp.resttime = mngt.time(h);  end
    % Run over length of mngt action
        for i = 1:length(mngt.length)
            if contains(mngt_scen,'fish'); pred.fishlgth = mngt.length(i); end
            if contains(mngt_scen,'cull'); urchin.culllgth = mngt.length(i); end  
            if contains(mngt_scen,'rest'); kelp.restlgth = mngt.length(i); end
    % Run over degree of mngt action
        for j = 1:deglngth        
            if contains(mngt_scen,'fish')
                if isscalar(mngt.degreeF); pred.fishF = mngt.degreeF(1);
                else; pred.fishF = mngt.degreeF(j); end 
            end
            if contains(mngt_scen,'cull') 
                if isscalar(mngt.degreeC); urchin.culln = mngt.degreeC(1);
                else; urchin.culln = mngt.degreeC(j); end  
            end
            if contains(mngt_scen,'rest'); kelp.restn = mngt.degreeR(j); end            
         
% Run simulations:
    if strcmp(model,'California_Jess')

        % Simulated predator model (CA version)
        [nt,nb] = run_Predator_Implicit(pred, T1+T2, RR, ones(pred.meshno,1), dist);

        % Predator biomass over kelp–urchin phase
        PBE(:,h,i,j,:) = sum(nb(pred.Lgraze_ind:end,(end-T2+1):end,:),1);

        % Run kelp–urchin model with simulated predator
        [kt2,ut2,~,RK_noise] = run_UrchinKelp_Implicit(kelp, urchin, T2, RR, ...
                               kt0, ut0, squeeze(PBE(:,h,i,j,:)), dist);

        % Set crabs to NaN (not included in CA model)
        cft2 = NaN; 
        cmt2 = NaN;

    elseif strcmp(model,'Oregon_Andres')

        % Forced predator (otter) model with crabs (OR version)
        [kt2,ut2,~,cft2,cmt2,RK_noise] = run_UrchinKelpCrab_Implicit(kelp, urchin, crab, T2, RR, ...
                         kt0, ut0, cft0, cmt0, pred_forced, dist);
    end


%% MODEL OUTPUTS --------------------------

% Outputs
    kts(:,:,h,i,j,:) = kt2;
    uts(:,:,h,i,j,:) = ut2;
    cfts(:,:,h,i,j,:) = cft2;
    cmts(:,:,h,i,j,:) = cmt2;

% Calculate kelp average (mean over last 1 year of run)
    kelp_avg(h,i,j,:) = mean(kt2(2,(end-4):end,:));
          
    
        end    
        end
        end

% Calculate if kelp exists in final time step
    kelp_pres(h,i,j,:) = kt2(2,end-2,:)>0;

% Proportion of sims with kelp persistence over time
    kelp_pers_t = sum(kts(2,:,:,:,:,:)>0,6);

% Read elapsed time from stopwatch (elapsed time)
    toc 
     