% Dun_model_otter.m

%% Top-level stuff --------------------------

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
      

%% Model parameters --------------------------

% run-times [winter, spring, summer, autumn]
    % T1 = 30*4; % for simulated predator
    T2 = 60*4; % for kelp-urchin runs 

% number of replicates (RR) - Max 10000 due to ORSO 
    RR =  2; % 10000;

% length of run + buffer
    tmax = T2+100; % 1000;

% crab (Dungeness) parameters
    crab = ParaCrab_Implicit(tmax);

% predator scenario (Choose your sea otters)
    ORSO_data = fullfile(ORSO, 'Table0.csv'); % no otters
    % ORSO_data = fullfile(ORSO, 'Scenario-Norm_SuccesSegment-N3.csv'); % Pacific city


%% Initial conditions --------------------------

% CRABS 
    cft0 = [8*10^5,8*10^4,8*10^4,8*10^4,8*10^4,8*10^4];
    cmt0 = [8*10^5,8*10^4,8*10^4,8*10^4,8*10^4,8*10^4];
   
% pre-assign variables:

    % predator biomass (Choose: forced "otters" or PBE "sheep-head") 
    pred_forced = ParaPred_Forced(ORSO_data, RR);
        % otter reintroduction buffer (allows kelp-urchin to stabilize)
        buffer = 5*4; % 0;
        pred_forced = [zeros(buffer, RR); pred_forced];
    % PBE = NaN(T2, length(mngt.time), length(mngt.length), deglngth, RR);   
    
    
%% Run models --------------------------
    % 2) forced predator [sea otter]
        % [kt2,ut2,~,RK_noise] = run_UrchinKelp_Implicit(kelp, urchin, T2, RR, ...
        %                  kt0, ut0, pred_forced, dist);

        [cft,cmt] = run_Crab_Implicit(crab, T2, RR, ...
            cft0, cmt0, pred_forced);


%% Model outputs --------------------------

% outputs
    % kts(:,:,h,i,j,:) = kt2;
    % uts(:,:,h,i,j,:) = ut2;


% read elapsed time from stopwatch (elapsed time)
    toc 
     