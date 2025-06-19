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
    RR =  10; % 10000;

% length of run + buffer
    tmax = T2+100; % 1000;

% kelp parameters (Giant or Bull inside "ParaKelp_Implicit" func)
    kelp = ParaKelp_Implicit(tmax);
    kelp.mu = 1*10^10; % 9*10^-5; % kelp.mu;

% urchin parameters
    urchin = ParaUrchin_Implicit(tmax);
    urchin.RU = 1.5*10^5; % 3*10^5; % urchin.RU;

% predator scenario (sea otters)
    % ORSO_data = fullfile(ORSO, 'Table0.csv'); % no otters
    % ORSO_data = fullfile(ORSO, 'Scenario-Norm_SuccesSegment-N3.csv'); % Pacific city
    % ORSO_data = fullfile(ORSO, 'Scenario-Norm_SuccesSegment-C7.csv'); % Newport
    ORSO_data = fullfile(ORSO, 'Scenario-Norm_SuccesSegment-S6.csv'); % Port Orford
            
    
%% Disturbance (dist) -------------------------- 
    
% disturbance length
    dist.lngth = 0; % 1; % 2; 

% disturbance timing (20 years pre-run is added to this)
    % non-disturbance scenario:
    dist.yrs = NaN; 
    % disturbance scenario:
    % dist.yrs = (20*4) + repmat(1:4,1,dist.lngth) + repelem(((1:dist.lngth)-1)*4,4);
       
% how do vital rates change during the disturbance (heatwave) 
    % kelp recruitment
    dist.RK = kelp.RK/7; % 0; 
    % change in kelp biomass
    dist.lambda = kelp.lambda .* repmat([1 1 0.5 0.5],1,tmax/4); %* 0.5;   
    % urchin grazing rates
    dist.hij = repmat(cell2mat(kelp.bhij) .* reshape([1.15 1.05 1.2 1.3],1,1,4), 1, 1, 1, tmax/4);


%% Management scenarios --------------------------
% NOTE: run over range of mngt length, timing, degrees

% MPA scenario (remove fishing pressure)🪝:
    pred.F = 0; % 0.1/4;

% management (mngt) scenarios:  
    mngt_scen =  'none'; % 'none'; % 'restoration'; % 'fishing'; % 'culling'; % 'cull&rest'; % 'fish&cull'; % 'fish&rest'; % 'fish&cull'; % 'fish&rest'; % 'fishrestcull'; % 'cull&rest'; %       
        % get vector values
        mngt = ParaMngt_Implicit(mngt_scen);
        pred.fish = mngt.fish;
        urchin.culling = mngt.culling;
        kelp.restore = mngt.restore;


%% Initial conditions --------------------------

% KELP 🌿 (choose between giant & bull kelp & start with high drift)
    % [juvenile, adult, drift] 
    kt0 = [1.17*10^2,1.17*10^2,1.17*10^3]; % Bull-kelp
    % kt0 = [1.17*10^5,1.17*10^5,1.17*10^6]; % Giant-kelp

% URCHINS 🟣 (start low to give system a chance to be in a kelp state)
    % [juvenile, hiding, expose]
    ut0 = [0,0,0]; 
   
% pre-assign variables:

    % degree (deg) lenght - degree is the last mngt in the matrix
    deglngth = structfun(@numel,mngt);
    deglngth = deglngth(end);

    % kelp biomass 
    kts = NaN(3, T2+1, length(mngt.time), length(mngt.length), deglngth,RR);

    % urchin biomass
    uts = kts;

    % average simulation in which kelp survives 
    kelp_avg = NaN(length(mngt.time), length(mngt.length), deglngth, RR);

    % predator biomass (Choose: forced "otters" or PBE "sheep-head") 
    pred_forced = ParaPred_Forced(ORSO_data, RR);
        % otter reintroduction buffer (allows kelp-urchin to stabilize)
        buffer = 5*4; % 0;
        pred_forced = [zeros(buffer, RR); pred_forced];
    % PBE = NaN(T2, length(mngt.time), length(mngt.length), deglngth, RR);   
    
    
%% Run models --------------------------

% run managements:
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
    % run over degree of mngt action
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
                kelp.restn = mngt.degreeR(j);     end            
         
% % run simulated pred model over time (turn on for PBE):
%     % urchin-kelp + simulated predator run (pred)
%         [nt,nb] = run_Predator_Implicit(pred, T1+T2, RR, ones(pred.meshno,1), dist);
%     % biomass of simulated predator predating on urchins
%         PBE(:,h,i,j,:) = sum(nb(pred.Lgraze_ind:end,(end-T2+1):end,:),1);

% run urchin-kelp model section (NOTE: predator biomass needs to be T2 x RR)
    % 1) simulated predator ['Blue rockfish'&'Sheephead']
        % [kt2,ut2,~,RK_noise] = run_UrchinKelp_Implicit(kelp, urchin, T2, RR, ...
        %                   kt0, ut0, squeeze(PBE(:,h,i,j,:)), dist);
    % 2) forced predator [sea otter]
        [kt2,ut2,~,RK_noise] = run_UrchinKelp_Implicit(kelp, urchin, T2, RR, ...
                         kt0, ut0, pred_forced, dist);


%% Model outputs --------------------------

% outputs
    kts(:,:,h,i,j,:) = kt2;
    uts(:,:,h,i,j,:) = ut2;

% calculate kelp average (mean over last 1 year of run)
    kelp_avg(h,i,j,:) = mean(kt2(2,(end-4):end,:));
          
    
        end    
        end
        end

% calculate if kelp exists in final time step
    kelp_pres(h,i,j,:) = kt2(2,end-2,:)>0;

% proportion of sims with kelp persistence over time
    kelp_pers_t = sum(kts(2,:,:,:,:,:)>0,6);

% read elapsed time from stopwatch (elapsed time)
    toc 
     