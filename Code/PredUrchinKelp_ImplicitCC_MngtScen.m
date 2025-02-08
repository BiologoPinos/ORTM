% PredUrchinKelp_ImplicitCC_MngtScen.m
 
% Authorship: 
% Andrés Pinos-Sánchez | andres.pinos.sanchez@gmail.com
% Co-authors: 
    % Jess Hopf (o.g. model developer), 
    % Leif Rasmuson, 
    % Mark Novak,
    % Will White

% Toolboxes required:
    % Statistics and Machine Learning Toolbox
    % Optimization Toolbox

% Model description:
    % Discrete time
    % Spatially implicit
% This version aims to look at effectiveness of different management


%% Top-level stuff --------------------------

% clear all
    clear   %Clear Workspace
    clc     %clear command window
    
% add paths
    addpath('functions\', 'functions\cbrewer')
    addpath('Model outputs\')

% define folder path for otter's abundance data
    ORSO_output_folder = 'C:\Users\pinosa\Documents\Git\OR_Trophic_Model\Code\Model inputs';
% specify otter population
    % ORSO_data = fullfile(ORSO_output_folder, 'Table1.csv'); %this brings ORSO output for "Port Orford"
     ORSO_data = fullfile(ORSO_output_folder, 'Table2.csv'); %this brings ORSO output for "Newport"
    % ORSO_data = fullfile(ORSO_output_folder, 'Table3.csv'); %this brings ORSO output for "Tillamook Bay"

% set random number
    rng(1) %randomness (stochasticity)

%start stopwatch timer (elapsed time)
    tic 
      

%% Model parameters --------------------------

% runtimes [winter, spring, summer, autmum]
    % T1 = 30*4;                % pre for fished predator
    T2 = 40*4; %35*4; %25*4;    % for kelp-urchin runs 

% number of replicates
    RR =  2; %1000; %10000;  

% load in parameter values
    tmax = 1000;
    kelp = ParaKelp_Implicit(tmax);
    urchin = ParaUrchin_Implicit(tmax);
    otter = ParaOtter(ORSO_data,RR); % model predator
    % pred = ParaPred_Implicit(tmax);
        % ignore everything that is "pred" 
        % this is for ['Blue rockfish'&'Sheephead'] pred scenarios


%% Disturbance -------------------------- 
    % no disturbance (heatwave) scenario for this project
    % (ignore this block)

% single acute disturbance 
    % length of disturbance (years)
    dist.lngth = 2; % 1; % 

    % time in model when disturbance happen  (20 years pre run is added to this)
    % 1 & 2 = winter and spring 3 & 4 = summer & autum
        % dist.yrs = (20*4) + repmat(1:4,1,dist.lngth) + repelem(((1:dist.lngth)-1)*4,4); % disturbance is on
        dist.yrs = NaN; % no disturbance scenario

   % How do vital rates change during the disturbance (heatwave) 
        % kelp recruitment
        dist.RK = kelp.RK/7; % 0; 
        % change in kelp biomass
        dist.lambda = kelp.lambda .* repmat([1 1 0.5 0.5],1,tmax/4); %* 0.5;   
        % urchin grazing rates
        dist.hij = repmat(kelp.bhij.*reshape([1.15 1.05 1.2 1.3],1,1,4),1,1,1,tmax/4); % kelp.hij * 1.2; % * 1.4;


%% Run models --------------------------
    % run over range of mngt length, timing, degrees

% management scenarios 
    % MPA (remove fishing) 🪝
        % pred.F = 0; % 0.1/4;
    % which management scenario to run over?
        mngt_scen =  'none'; % 'none'; % 'restoration'; % 'fishing'; % 'culling'; % 'cull&rest'; % 'fish&cull'; % 'fish&rest'; % 'fish&cull'; % 'fish&rest'; % 'fishrestcull'; % 'cull&rest'; %       
    % get vector values
        mngt = ParaMngt_Implicit(mngt_scen);
        pred.fish = mngt.fish;
        urchin.culling = mngt.culling;
        kelp.restore = mngt.restore;

% initial conditions
    % KELP 🌿
        % start with high drift so that it give the system a chance to be in a kelp state
        % [juvenile, adult, drift]
        kt0 = [3.17*10^5,3.17*10^5,3.17*10^6]; % [1.17*10^5,1.17*10^5,1.17*10^6]; % [0,1,0]; % [0,0,0];
    % URCHINS 🟣
        % start with low numbers to give system a chance to be in a kelp state
        % [juvenile, hiding, expose]
        ut0 = [0,0,0]; % [urchin.RU,urchin.RU,0];
    
% pre-assign variables
    deglngth = structfun(@numel,mngt);  %(deg is degree)
    deglngth = deglngth(end); % degree is the last mngt in the matrix
    kts = NaN(3, T2+1, length(mngt.time), length(mngt.length), deglngth,RR);
    uts = kts;
    kelp_avg = NaN(length(mngt.time), length(mngt.length), deglngth, RR);
    % PBE = NaN(T2, length(mngt.time), length(mngt.length), deglngth, RR); % relevant for modeling predators   

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
         
% run pred model over time
    % urchin-kelp + pred run
        % [nt,nb] = run_Predator_Implicit(pred, T1+T2, RR, ones(pred.meshno,1), dist); % relevant for modeling predators
    % biomass predating on urchins
        % PBE(:,h,i,j,:) =
        % sum(nb(pred.Lgraze_ind:end,(end-T2+1):end,:),1); % relevant for modeling predators

% urchin-kelp model section
    % note that input predator biomass needs to be size T2 x RR)
        [kt2,ut2,~,RK_noise] = run_UrchinKelp_Implicit(kelp, urchin, T2, RR, ...
                         kt0, ut0, otter, dist);
        % think about nt_static_forced instead of otter & and create an "if" statement


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
     


% save(['Implicitv6a_MngtScen_v1_',mngt_scen,'_',datestr(now, 'yyyymmdd'),'.mat'],... % '_degreesLRG.mat'],... % 
%       "kelp_pers_t", "T1", "T2", "kelp", "urchin", "pred", 'dist', 'mngt', 'mngt_scen', '-v7.3')

% for baseline no disturbance
    % baseline_persit_overtime = kelp_pers_t(:,:,:,:,1);
    % save(['Implicitv6a_MngtScen_v1_nodist_',datestr(now, 'yyyymmdd'),'.mat'],... % '.mat'],... %
    %  "baseline_persit_overtime", '-v7.3')

% for disturbance no mngt action
    % % baseline
    % baseline_persit_overtime = kelp_pers_t(:,:,:,:,1);
    % 
    % % mean values before disturbance
    % persisting_runs = find(squeeze(prod(kts(2,1:(dist.yrs(1)-1),:,:,:,:)>0,2))>0)';
    % kelpJ_avg_pre = mean(kt2(1,(dist.yrs(1)-40):(dist.yrs(1)-1),persisting_runs),[2,3]);
    % kelpA_avg_pre = mean(kt2(2,(dist.yrs(1)-40):(dist.yrs(1)-1),persisting_runs),[2,3]);
    % urchinA_avg_pre = mean(sum(ut2(2:3,(dist.yrs(1)-40):(dist.yrs(1)-1),persisting_runs),1),[2,3]);
    % save(['Implicitv6a_MngtScen_v1_',mngt_scen,'_',datestr(now, 'yyyymmdd'),'.mat'],... % '.mat'],... %
    %       "kelpJ_avg_pre","kelpA_avg_pre","urchinA_avg_pre", "baseline_persit_overtime", '-v7.3')
    % 


