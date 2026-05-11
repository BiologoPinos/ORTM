% ORTM_model_otter.m
    % June 2026
 
% Authorship: 
    % Andrés Pinos-Sánchez | andres.pinos.sanchez@gmail.com
    % Co-authors: Jess Hopf, Leif Rasmuson, Mark Novak, Will White

% Toolboxes required:
    % Statistics and Machine Learning Toolbox
    % Optimization Toolbox

% Model description:
    % Discrete time
        % T is at the season level
        % tt is at the daily level (within T)
    % Spatially implicit
        % Roughly 1 hectare

% Instructions for user:
    % User MUST choose the desire conditions in steps 1-4 to set the model
    % User MUSENT change sections 6 and 7, otherwise code might break


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
      

%% 1) PREP MODEL (USER CHOOSES) --------------------------

% Choose your model (this determines which model to run, CA vs OR)
    model = 'Oregon_Andres'; % opts: 'California_Jess' | 'Oregon_Andres'

% Choose your kelp    
    kelp_species = 'Bull_kelp'; % opt: 'Giant_kelp' | 'Bull_kelp' | 'Bull_kelp_north' | 'Bull_kelp_south'

% Choose your urchins
    urchin_species = 'Urchins_OR'; % opt: 'Urchins_CA' | 'Urchins_OR'

% Choose your crabs
    crab_species = 'Dungeness_OR_LogNorm'; % opt: 'Dungeness_OR_Norm' | 'Dungeness_OR_LogNorm'

% Choose your sea otters
    otter_population = 'Scenario-High_SuccessSegment-S6.csv'; % opt: 'Table0.csv' | 'Scenario-High_SuccessSegment-S6.csv'
    

%% 2) MODEL PARAMETERS (Set for "realism") --------------------------

% Run-times/time-steps [winter, spring, summer, autumn]
    T1 = 30*4; % for simulated predator runs
    T2 = 40*4; % for kelp-urchin runs

% Number of replicates (RR) - Max 10000 because of ORSO 
    RR =  10000;

% Length of run + buffer
    tmax = T2+100;

% Kelp parameters
    kelp = ParaKelp_Implicit(tmax, kelp_species);

% Urchin parameters
    urchin = ParaUrchin_Implicit(tmax, urchin_species);
        urchin.phi = 0;
        
% Crab (Dungeness) parameters
    crab = ParaCrab_Implicit(tmax, crab_species);

% Predator simulated (sheep-head, CA_model)
    pred = ParaPred_Implicit(tmax); % NEED REVISSION - CA model not working

% Predator forced inputs (sea otters, OR_model)
    ORSO_data = fullfile(ORSO, otter_population);
    otter_reintro_buffer = 5*4; % 5*4; % 0; % Buffer allows kelp–urchin to stabilize
            

%% 3) INITIAL CONDITIONS (Set for "realism") --------------------------

% KELP 🌿 [juvenile, adult, drift] 
    if strcmp(kelp_species, 'Giant_kelp')
        kt0 = [1.17*10^5, 1.17*10^5, 1.17*10^6];  % Giant-kelp
    else
        kt0 = [8*10^4, 8*10^4, 8*10^5];           % Bull-kelp variants
    end

% URCHINS 🟣 [juvenile, hiding, expose]
    ut0 = [0, 1.92*10^3, 2.15*10^3];
    % ut0 = [0, 0, 0]; 

% CRABS 🦀 [age_0, age_1-10] (COMMON DENSITY (base-case; literature-consistent), for low density x*0.5, for high density x*3)
    cft0 = [0, 753, 753, 0, 0, 587, 587, 0, 0, 555, 555, 0, 0, 501, 501, 0, 0, 391, 391, 0, ...
            0, 331, 331, 0, 0, 243, 243, 0, 0, 203, 203, 0, 0, 165, 165, 0, 0, 142, 142, 0, 0, 118, 118, 0];
    % cft0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

    cmt0 = [0, 753, 753, 0, 0, 587, 587, 0, 0, 555, 555, 0, 0, 501, 501, 0, 0, 391, 391, 0, ...
            0, 241, 241, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    % cmt0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

% Other prey (X)
    xt0 = 5.1253e+03; % 100% of base adult crab biomass T = 81
    % xt0 = 2.5627e+03; % 50% of base adult crab biomass T = 81
    % xt0 = 1.0251e+03; % 20% of base adult crab biomass T = 81
    % xt0 = 0;

    
%% 4) MANAGEMNET & DISTURBANCE -------------------------- 

    % User Note: Modify ParaMngt_Implicit.m script accordingly  
    
% Disturbance length (How long the disturbance will last)
    dist.lngth = 1;  

% Disturbance timing (When the disturbance will happen)
    if dist.lngth == 0
        dist.yrs = NaN; % no disturbance (default)
    else    
        dist.yrs = (7*4) + repmat(1:4,1,dist.lngth) + repelem(((1:dist.lngth)-1)*4,4);
    end 
       
% How do vital rates change during the disturbance (heatwave) 

    % Kelp recruitment
        dist.RK = kelp.RK; % dist.RK = kelp.RK/7; % for heat wave version
    
    % Kelp biomass growth reduction during disturbance
        dist.lambda = kelp.lambda; % dist.lambda = kelp.lambda.* repmat([1 1 0.5 0.5],1,tmax/4); % for heat wave version  
    
    % Urchin grazing increase during disturbance
        dist.hij = kelp.hij; % dist.hij = repmat(cell2mat(kelp.bhij) .* reshape([1.15 1.05 1.2 1.3],1,1,4), 1, 1, 1, tmax/4); % for heat wave version

% Which management scenario to run over?
    mngt_scen = 'cull&rest'; % 'none'; % 'restoration'; % 'culling'; % 'cull&rest';  

% Get vector values
    mngt = ParaMngt_Implicit(mngt_scen);
    pred.fish       = mngt.fish;
    urchin.culling  = mngt.culling;
    kelp.restore    = mngt.restore;
    urchin.season   = mngt.season;
    kelp.season     = mngt.season;
    urchin.strategy = mngt.strategy;
    urchin.time_vec = mngt.time_vec;
    kelp.strategy   = mngt.strategy;
    kelp.time_vec   = mngt.time_vec;


%% 5) PRE-ASSIGN VARIABLES (empty vectors) --------------------------

% Degree (deg) lenght - degree is the last mngt in the matrix
    deglngth = structfun(@numel,mngt);
    deglngth = deglngth(end);

% Kelp biomass 
    kts = NaN(3, T2+1, length(mngt.time), length(mngt.length), deglngth, RR);

% Urchin biomass
    uts = NaN(3, T2+1, length(mngt.time), length(mngt.length), deglngth, RR);

% Crab biomass
    cfts = NaN(44, T2+1, length(mngt.time), length(mngt.length), deglngth, RR);
    cmts = NaN(44, T2+1, length(mngt.time), length(mngt.length), deglngth, RR);

% Average kelp persistence
    kelp_avg = NaN(length(mngt.time), length(mngt.length), deglngth, RR);

% Predator input data
    if strcmp(model,'California_Jess')
        PBE = NaN(T2, length(mngt.time), length(mngt.length), deglngth, RR);
    else
        pred_forced = ParaPred_Forced(ORSO_data, RR);
        buffer = otter_reintro_buffer; % Buffer is set in section 1 
        pred_forced = [zeros(buffer, RR); pred_forced];
    end   

% Prey consumed
    prey_consumed = NaN(4, T2, length(mngt.time), length(mngt.length), deglngth, RR);
    
    
%% 6) RUN MODELS --------------------------

% Run managements:

    % run over timing of mngt action
        for h = 1:length(mngt.time)

            if contains(mngt_scen,'fish')
                pred.fishtime = mngt.time(h)+T1;    end
            if contains(mngt_scen,'cull') 
                urchin.culltime = mngt.time(h);     end  
            if contains(mngt_scen,'rest') 
                kelp.resttime = mngt.time(h);       end
    % run over length of mngt action
        for i = 1:length(mngt.length)

            if contains(mngt_scen,'fish')
                pred.fishlgth = mngt.length(i);     end
            if contains(mngt_scen,'cull') 
                urchin.culllgth = mngt.length(i);   end  
            if contains(mngt_scen,'rest') 
                kelp.restlgth = mngt.length(i);     end

    % Run over degree of mngt action
        for j = 1:deglngth

            if contains(mngt_scen,'fish')
                if isscalar(mngt.degreeF)
                   pred.fishF = mngt.degreeF(1);
                else 
                   pred.fishF = mngt.degreeF(j);     
                end 
            end
            if contains(mngt_scen,'cull') 
                if isscalar(mngt.degreeC)
                    urchin.culln = mngt.degreeC(1);
                else
                    urchin.culln = mngt.degreeC(j);
                end  
            end
            if contains(mngt_scen,'rest') 
                kelp.restn = mngt.degreeR(j);
            end              
         
% Run simulations:
    if strcmp(model,'California_Jess')

        % Simulated predator model (CA version)
        [nt,nb] = run_Predator_Implicit(pred, T1+T2, RR, ones(pred.meshno,1), dist);

        % Predator biomass over kelp–urchin phase
        PBE(:,h,i,j,:) = sum(nb(pred.Lgraze_ind:end,(end-T2+1):end,:),1);

        % Run kelp–urchin model with simulated predator
        [kt2,ut2,~,RK_noise] = run_UrchinKelp_Implicit(kelp, urchin, T2, RR, ...
                               kt0, ut0, squeeze(PBE(:,h,i,j,:)), dist);

        % Set crabs to NaN for CA model
        cft2 = NaN; 
        cmt2 = NaN;

    elseif strcmp(model,'Oregon_Andres')

        % Forced predator (otter) model with crabs (OR version)
        [kt2,ut2,~,cft2,cmt2,RK_noise,prey_consumed2] = run_UrchinKelpCrab_Implicit(kelp, urchin, crab, T2, RR, ...
                         kt0, ut0, cft0, cmt0, xt0, pred_forced, dist);
    end


%% 7) MODEL OUTPUTS --------------------------

% Outputs
    kts(:,:,h,i,j,:)            = kt2;  % kelp biomass
    uts(:,:,h,i,j,:)            = ut2;  % Urchin biomass
    cfts(:,:,h,i,j,:)           = cft2; % Female crab biomass
    cmts(:,:,h,i,j,:)           = cmt2; % Male crab biomass
    prey_consumed(:,:,h,i,j,:)  = prey_consumed2; % Prey consumed

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
     