function [kt,ut,GC,cft,cmt,RK_noise,otterFR,prey_consumed] = run_UrchinKelpCrab_Implicit(kelp, urchin, crab, tmax, RR, kt1, ut1, cft1, cmt1, xt1, nt, dist)

% Description:
    % function running single populations kelp <> urchins model
    % see parent file (Para) for parameter values

% Runs with:
    % ORTM_model_otter.m


%% UNPACK STRUCTURES --------------------------
    % (note: some structures are unpacked later)

% Kelp 🌿
    RKstdv      = kelp.RKstdv;
    RKbeta      = kelp.RKbeta;
    D           = kelp.D;
    muvar       = kelp.muvar;
    RTk         = kelp.RTk;
    reproWeight = kelp.reproWeight;  
    RKlag       = kelp.RKlag;          
    g           = kelp.g;
    rS          = kelp.rS;
    c           = kelp.c;
    d           = kelp.d;
    rD          = kelp.rD;
    aij         = cell2mat(kelp.aij);
    % bhij      = cell2mat(kelp.bhij); % unpacked later
    % hij       = cell2mat(kelp.hij); % unpacked later
    % lambda    = kelp.lambda; % unpacked later

% Urchin 🟣
    RU          = urchin.RU;
    RUstdv      = urchin.RUstdv;
    RUdist      = urchin.RUdist;
    RTU         = urchin.RTU;
    gJ          = urchin.gJ;
    MJU         = urchin.MJU;
    MHU         = urchin.MHU;
    MEU         = urchin.MEU;
    aH          = urchin.aH;
    aE          = urchin.aE;
    bH          = urchin.bH;
    bE          = urchin.bE;
    wu          = urchin.wu;
    phi         = urchin.phi;
    FU          = urchin.FU;
    PLD         = urchin.PLD;
    tau         = urchin.tau;
    w1          = urchin.w1;
    w2          = urchin.w2;
    kmin        = urchin.kmin;

% Crab 🦀
    RC          = crab.RC;
    RCstdv      = crab.RCstdv;
    RCdist      = crab.RCdist;
    RTC         = crab.RTC;
    beta        = crab.beta;
    MRC         = crab.MRC;
    MJC         = crab.MJC;
    MAC         = crab.MAC;
    FC          = crab.FC;
    aC          = crab.aC;
    bC          = crab.bC;
    wc          = crab.wc;
    NCIF        = crab.NCIF{1};
    weight      = crab.weight{1};
    weight      = reshape(weight,44,1,1);

% Other parameters
    dist_yrs = dist.yrs;



%% VECTORS FOR STATE VARIABLES --------------------------

    % Predator
        nt = reshape(nt,[],1,RR); % nt = nt([],1,RR);
    
    % Feeding rate
        otterFR = NaN(4, tmax, RR);
    
    % Prey consumed
        prey_consumed = NaN(4, tmax, RR);
    
    % Kelp
        kt = NaN(3,tmax,RR);
        kt(:,1,:) = repmat(kt1(:),1,1,RR);
        
    % Urchin
        ut = NaN(3,tmax,RR);
        ut(:,1,:) =  repmat(ut1(:),1,1,RR);
    
    % Grazing capacity
        GC = NaN(tmax,1,RR);
    
    % Dungeness crab females
        cft = NaN(44,tmax,RR);
        cft(:,1,:) = repmat(cft1(:),1,1,RR);
    
    % Dungeness crab males
        cmt = NaN(44,tmax,RR);
        cmt(:,1,:) = repmat(cmt1(:),1,1,RR);
    
    % Switching function (urchins), set Psi to 1
        Psi = 1;


%% VECTORS FOR RECRUITMENT (set noise vectors for recruits) --------------------------

    % Kelp recruitment noise
        RK_noise = repelem(max(zeros(tmax/4,1,RR),normrnd(1,RKstdv,tmax/4,1,RR)),4,1,1);
        
    % Urchin recruitment noise
        if strcmp(RUdist{1}, 'normal')
            RU_noise = repelem(max(zeros(tmax/4,1,RR),normrnd(1,RUstdv,tmax/4,1,RR)),4,1,1);
        else 
            RU_noise = repelem(max(zeros(tmax/4,1,RR),lognrnd(1,RUstdv,tmax/4,1,RR)),4,1,1);
        end   
    
    % Crab recruitment noise (CHECK THIS)
        if strcmp(RCdist{1}, 'normal')
            RC_noise = repelem(max(zeros(tmax/4,1,RR), normrnd(1, RCstdv, tmax/4,1,RR)), 4,1,1);
        else
            RC_noise = repelem(max(zeros(tmax/4,1,RR), lognrnd(1, RCstdv, tmax/4,1,RR)), 4,1,1);
        end


%% RUN OVER TIME (seasons first) --------------------------

for t = 1:tmax

    if ismember(t, dist_yrs)        % add disturbance effects if happens
        RK = dist.RK;               % change in kelp recruitment
        lambda = dist.lambda(t);    % change in standing kelp growth
        hij = cell2mat(dist.hij);
        hij = hij(:,:,t);
        % hij = dist.hij(:,:,t);      % change in urchin grazing rates  
    else 
        RK = kelp.RK;
        lambda = kelp.lambda(t);
        hij = cell2mat(kelp.hij);
        hij = hij(:,:,t);
    end


    %% RUN INNER LOOP (daily predation)

        % Number of days within a season (time step)
            ndays = 91;
            
        % Vectors for state variables (daily prey biomass)
        
            % Urchins
                ut_d = NaN(3, ndays+1, RR);
                ut_d(:,1,:) = ut(:,t,:);
        
            % Dungeness crab females
                cft_d = NaN(44, ndays+1, RR);
                cft_d(:,1,:) = cft(:,t,:);
        
            % Dungeness crab males
                cmt_d = NaN(44, ndays+1, RR);
                cmt_d(:,1,:) = cmt(:,t,:);
           
            % Predator (same number)
                nt_season = reshape(nt(t,:,:), 1, 1, RR);

            % Prey vectors (S, A, X)
                S = NaN(1, ndays+1, RR);
                A = NaN(1, ndays+1, RR);
                X = NaN(1, ndays+1, RR);

            % Preference vectors
                p1 = NaN(1, ndays+1, RR);       
                p2 = NaN(1, ndays+1, RR);

            % Feeding rate vectors
                feed_rate_uh = NaN(1, ndays+1, RR);
                feed_rate_ue = NaN(1, ndays+1, RR);
                feed_rate_cf = NaN(1, ndays+1, RR);
                feed_rate_cm = NaN(1, ndays+1, RR);
 
        % Run daily predation loop
        
            for tt = 1:ndays

                % Yodzis preference/switching        

                    % Compute available biomasses per replicate (shape 1 x 1 x RR)
                        S(:,tt,:) = sum(ut_d(2:3, tt, :), 1);   % 1 x 1 x RR -> urchin available biomass
                        A(:,tt,:) = sum(cft_d(9:end, tt, :) + cmt_d(9:end, tt, :) .* weight(9:end), 1); % adult crab available biomass
                        X(:,tt,:) = xt1;

                    % Predator presence mask
                        hasPred = (nt_season > 0); % 1 x 1 x RR
    
                    % Compute switching fractions
                        % p1(:,tt,:) = 1;   % shape 1x1xRR + smallEPS
                        % p2(:,tt,:) = 1;   % shape 1x1xRR + smallEPS
                        p1(:,tt,:) = wu .* (S(:,tt,:) .^ phi) ./ (wu .* (S(:,tt,:) .^ phi) + wc .* (A(:,tt,:) .^ phi) + (1-wu-wc) .* (X(:,tt,:) .^ phi));   % shape 1x1xRR + smallEPS
                        p2(:,tt,:) = wc .* (A(:,tt,:) .^ phi) ./ (wu .* (S(:,tt,:) .^ phi) + wc .* (A(:,tt,:) .^ phi) + (1-wu-wc) .* (X(:,tt,:) .^ phi));   % shape 1x1xRR + smallEPS
                   
                    % No predators → no switching 
                        p1(:,tt,~hasPred) = 0; 
                        p2(:,tt,~hasPred) = 0;

                % Feeding rates
                    
                    % Urchins
                        feed_rate_uh(:,tt,:) = exp(p1(:,tt,:) .* (-Func_TypeII(aH, bH, S(:,tt,:)) .* nt_season(:,:,:))); % hiding
                        feed_rate_ue(:,tt,:) = exp(p1(:,tt,:) .* (-Psi .* Func_TypeII(aE, bE, S(:,tt,:)) .* nt_season(:,:,:))); % exposed

                    % Crabs
                        feed_rate_cf(:,tt,:) = exp(p2(:,tt,:) .*(-Func_TypeII(aC, bC, A(:,tt,:)) .* nt_season(:,:,:))); % females
                        feed_rate_cm(:,tt,:) = exp(p2(:,tt,:) .*(-Func_TypeII(aC, bC, A(:,tt,:)) .* nt_season(:,:,:))); % males

                % Daily predation 
                
                    % Urchins (juveniles unaffected by predation)
                        ut_d(1,tt+1,:) = ut_d(1,tt,:);
                        ut_d(2,tt+1,:) = ut_d(2,tt,:) .* feed_rate_uh(:,tt,:); % hiding
                        ut_d(3,tt+1,:) = ut_d(3,tt,:) .* feed_rate_ue(:,tt,:); % exposed
        
                    % Female crabs (age 0-1 unaffected by predation)
                        cft_d(1:8,tt+1,:) = cft_d(1:8,tt,:);
                        cft_d(9:end,tt+1,:) = cft_d(9:end,tt,:) .* feed_rate_cf(:,tt,:); % female
            
                    % Daily predation on male crabs (age 0-1 unaffected by predation)
                        cmt_d(1:8,tt+1,:) = cmt_d(1:8,tt,:);
                        cmt_d(9:end,tt+1,:) = cmt_d(9:end,tt,:) .* feed_rate_cm(:,tt,:); % male

            end % end daily loop
    
        % Extract post-predation (end-of-season) biomasses
            ut_post     = ut_d(:, end, :);   % 3 x 1 x RR
            cft_post    = cft_d(:, end, :); % 44 x 1 x RR
            cmt_post    = cmt_d(:, end, :); % 44 x 1 x RR

        % Extract the realized sea otter feeding rates
            otterFR(1,t,:) = mean(feed_rate_uh(:,1:ndays,:), 2);   % urchin hiding
            otterFR(2,t,:) = mean(feed_rate_ue(:,1:ndays,:), 2);   % urchin exposed
            otterFR(3,t,:) = mean(feed_rate_cf(:,1:ndays,:), 2);   % female crab
            otterFR(4,t,:) = mean(feed_rate_cm(:,1:ndays,:), 2);   % male crab

        % Realized prey consumed during the season
            prey_consumed(1,t,:) = max(squeeze(ut_d(2,1,:) - ut_d(2,end,:)), 0);                 % hiding urchins
            prey_consumed(2,t,:) = max(squeeze(ut_d(3,1,:) - ut_d(3,end,:)), 0);                 % exposed urchins
            prey_consumed(3,t,:) = max(squeeze(sum(cft_d(9:end,1,:),1) - sum(cft_d(9:end,end,:),1)), 0);  % female crabs
            prey_consumed(4,t,:) = max(squeeze(sum(cmt_d(9:end,1,:),1) - sum(cmt_d(9:end,end,:),1)), 0);  % male crabs


%% RUN CRAB (Dungeness) 🦀 --------------------------

% Survival

    % Females 
       Sf = cat(1, repmat(exp(-MRC),3,1,RR), repmat(exp(-MJC),4,1,RR), repmat(exp(-MAC), 36,1,RR));

    % Males
       Sm = cat(1, repmat(exp(-MRC),3,1,RR), repmat(exp(-MJC),4,1,RR), repmat(exp(-MAC), 8,1,RR), repmat(exp(-MAC - FC(t)), 28,1,RR));

% Projection/Transition matrix
    
    % Females
       Mcf = zeros(44, 44, RR);
       Mcf(2:44, 1:44-1, :) = reshape(eye(44-1), 44-1, 44-1, 1) .* reshape(squeeze(Sf(:,1,:)), 44-1, 1, RR); % 43x43xRR stored in the sub-diagonal of 44x44xRR
    
    % Males
       Mcm = zeros(44, 44, RR);
       Mcm(2:44, 1:44-1, :) = reshape(eye(44-1), 44-1, 44-1, 1) .* reshape(squeeze(Sm(:,1,:)), 44-1, 1, RR);  % 43x43xRR stored in the sub-diagonal of an 44x44xRR

% Recruitment

    % Incoming recruits with timing + noise
        RCnew = zeros(44,1,RR);
        RCnew(1,1,:) = RC .* RTC(t) .* RC_noise(t,:,:);
    
    % Cannibalism (Ricker Style)
        RCnew(1,1,:) = RCnew(1,1,:) .* exp(-beta .* (sum(NCIF .* (cft_post(5:44,:,:) + cmt_post(5:44,:,:)))));

% Advance the adult populations
    cft(:,t+1,:) = pagemtimes(Mcf, cft_post) + (RCnew .* 0.5); % females
    cmt(:,t+1,:) = pagemtimes(Mcm, cmt_post) + (RCnew .* 0.5); % males 

    % Fishery: remove all males age ≥6 (index 25)
    cmt(25:end, t+1, :) = 0;

    % Removal of all elderly (Fall age 10, index 44) crabs (absolute natural mortality)
    cft(44, t+1, :) = 0;
    cmt(44, t+1, :) = 0;



%% RUN URCHINS 🟣 --------------------------
    
% Predation on urchin threshold (set earlier values for delay)
    if t>2                                    % if for last 6 months (2 seasons) kelp biomass is bellow threshold (kmin),                       
        ktdelay = sum(kt(1:2,(t-2:t-1),:));   % then predators stop feeding (and urchin fishing turns off)      
        Psi = sum(ktdelay < kmin) ~= 2;       % Psi is 1 until k < kmin for consecutive time periods
    end
  
% Survival (urchins)

    % Juveniles (juv):
        % density-independent survival
            sJ = exp(-MJU);
        % density-dependent survival (more adults = more survival)
            % sJ = sJ.*(1-exp(-alpha*(ut(2,t)+ut(3,t))));
        % with scaling transition
            % sJ = sJ.*(1-exp(-alpha*(ut(2,t)+ut(3,t))) - 0.5*alpha^2 * exp(-alpha*(ut(2,t)+ut(3,t))) * alphavar);
        % give error if neg survival and set sJ to zero. 
            if sJ<0; sJ = 0; end
        
    % Adults (Hiding and Exposed):
        % type I linear predation (e.g. sheep head)
            % hiding adults 
            % sH = exp(-MH -aH.*nt(t,:,:) - FU);
            % exposed adults
            % sE = exp(-ME -Psi.*(aE.*nt(t,:,:) + FU));
    
        % type II predation (e.g. sea otters)
            % hiding adults 
            sH = exp(-MHU); % exp(-MHU -FU -Func_TypeII(aH,bH,sum(ut(2:3,t,:))).*nt(t,:,:));

            % exposed adults
            sE = exp(-MEU); % exp(-MEU -Psi .* (FU + Func_TypeII(aE,bE,sum(ut(2:3,t,:))).*nt(t,:,:)));

% Proportion being exposed

    % Grazing capacity (ratio of drift kelp to total urchin feeding rate) - using post predation adults 
    GC(t,1,:) = max(kt(3,t,:)./(sum(ut_post(2:3,:,:)).*hij(3,1)),0); % GC(t,1,:) = max(kt(3,t,:)./(sum(ut(2:3,t,:)).*hij(3,1)),0); % pre-predation urchins

    % Grazing behavior switching
    Phi = reshape(Func_Switch(w1, w2, squeeze(GC(t,1,:))),1,1,RR); % Phi = reshape(Func_Switch(w1, w2, squeeze(GC(t,1,:))),1,1,RR); % pre-predation urchins

% Projection/Transition matrix
    Mu = [repmat((1-gJ) * sJ,1,1,RR),  zeros(1,1,RR),   zeros(1,1,RR);
          gJ .* sH .* (1-Phi),         sH .* (1-Phi),   sH .* (1-Phi);
          gJ .* sE .* Phi,             sE .* Phi,       sE .* Phi];

% Incoming recruits with timing + noise
    RUadd = zeros(3,1,RR);
    RUadd(1,1,:) = RU .* RTU(t) .* RU_noise(t,:,:);

% Next yrs numbers
    ut(:,t+1,:) = pagemtimes(Mu, ut_post) + (RUadd .* sJ^tau);

% Urchin culling (mass mortality events)      
    if urchin.culling == "Y"
        if ismember(t, dist_yrs(1) + urchin.culltime + (0:urchin.culllgth-1)) && ismember(rem(t,4), urchin.season)
            % cull all exposed first, then the hiding urchins
                % which reps have more exposed than the cull number
                exp_reps = find(ut(3,t+1,:) > urchin.culln)';
                % ...and which dont
                both_reps = find(ut(3,t+1,:) <= urchin.culln)';
                % for exp_reps, just cull from exposed
                ut(3,t+1,exp_reps) = ut(3,t+1,exp_reps) - urchin.culln;
                % for both_reps, cull all from exposed and rest from hiding
                ut(3,t+1,both_reps) = 0;
                ut(2,t+1,both_reps) = max(ut(2,t+1,both_reps)-(urchin.culln - ut(3,t+1,both_reps)),0);
        end
    end

    
%% RUN KELP 🌿 --------------------------
   
% Kelp restoration (adding juvenile kelp biomass)
    RKr = 0;
    if kelp.restore == "Y"
        if ismember(t, dist_yrs(1) + kelp.resttime + (0:kelp.restlgth-1)) && ismember(rem(t,4), kelp.season)
            RKr = kelp.restn;
        end
    end
    
% Recruitment

    % Per-capita settling spores + timing + noise
        RKset = RK .* RTk(t) .* RK_noise(t,:,:);
        
    % Total incoming settled spores biomass (before DD)  spores from adults + restored juvs
        ksetT = RKset .* sum(reproWeight .* kt(2,max(1,t-RKlag),:)) + RKr;

    % DD survival of young
  
        % Line to avoid NA (if ksetn = 0 & k2n = 0 then the DD function is NaN)
        ksetn = max(ksetT, 1);

        % Adult densities to be used for inter-cohort part of DD
        k2n = max(sum(reproWeight .* kt(2,max(1,t-RKlag),:)), 1);

        % DD functions

            % for Beverton-Holt DD function (note this basic form assumes that the slope at the orgin = 1)
            if D == 1
                sY = 1/(1+ksetn./RKbeta);
    
            else % for mixed and Ricker DD functions
                DD = ( (D-1).*k2n.*exp(RKbeta.*D.*k2n) ) ./ (D.*ksetn.*exp(RKbeta.*D.*k2n) + exp(RKbeta.*k2n) .* ((D-1).*k2n-D.*ksetn)  );
    
                % second derivative of mixed DD function (calculted with matlab solver, see ExploringMixedDD_v0.m) 
                DD2 = (2.*k2n.*exp(D.*RKbeta.*k2n).*(D - 1).*(exp(RKbeta.*k2n).*(D - 1) + RKbeta.*exp(RKbeta.*k2n).*(k2n.*(D - 1) - D.*ksetn) + D^2.*ksetn.*RKbeta.*exp(D.*RKbeta.*k2n)).^2)./(exp(RKbeta.*k2n).*(k2n.*(D - 1) - D.*ksetn) + D.*ksetn.*exp(D.*RKbeta.*k2n)).^3 - (k2n.*exp(D.*RKbeta.*k2n).*(D - 1).*(2.*RKbeta.*exp(RKbeta.*k2n).*(D - 1) + RKbeta.^2.*exp(RKbeta.*k2n).*(k2n.*(D - 1) - D.*ksetn) + D.^3.*ksetn.*RKbeta.^2.*exp(D.*RKbeta.*k2n)))./(exp(RKbeta.*k2n).*(k2n.*(D - 1) - D.*ksetn) + D.*ksetn.*exp(D.*RKbeta.*k2n)).^2 - (2.*exp(D.*RKbeta.*k2n).*(D - 1).*(exp(RKbeta.*k2n).*(D - 1) + RKbeta.*exp(RKbeta.*k2n).*(k2n.*(D - 1) - D.*ksetn) + D.^2.*ksetn.*RKbeta.*exp(D.*RKbeta.*k2n)))./(exp(RKbeta.*k2n).*(k2n.*(D - 1) - D.*ksetn) + D.*ksetn.*exp(D.*RKbeta.*k2n)).^2 + (2.*D.*RKbeta.*exp(D.*RKbeta.*k2n).*(D - 1))./(exp(RKbeta.*k2n).*(k2n.*(D - 1) - D.*ksetn) + D.*ksetn.*exp(D.*RKbeta.*k2n)) + (D.^2.*RKbeta.^2.*k2n.*exp(D.*RKbeta.*k2n).*(D - 1))./(exp(RKbeta.*k2n).*(k2n.*(D - 1) - D.*ksetn) + D.*ksetn.*exp(D.*RKbeta.*k2n)) - (2.*D.*RKbeta.*k2n.*exp(D.*RKbeta.*k2n).*(D - 1).*(exp(RKbeta.*k2n).*(D - 1) + RKbeta.*exp(RKbeta.*k2n).*(k2n.*(D - 1) - D.*ksetn) + D^2.*ksetn.*RKbeta.*exp(D.*RKbeta.*k2n)))./(exp(RKbeta.*k2n).*(k2n.*(D - 1) - D.*ksetn) + D.*ksetn.*exp(D.*RKbeta.*k2n)).^2;
                
                % per capita survival of settlers including scale transition
                sY = DD + 0.5.*DD2.*muvar;
            end            
  
    % Total new biomass of juveniles - post predation urchins
        RKnew = ksetT .* sY .* rS .* exp(-Func_TypeII(aij(1,2),hij(1,2),sum(kt(1:2,t,:))).*ut_post(3,:,:)) .* lambda;
   
% Projection/Transition matrix (bull Kelp) - post predation urchins
    Mk = [zeros(1,1,RR),...
          zeros(1,1,RR),...
          zeros(1,1,RR);

          rS .* g .* (1-c) .* exp(-Func_TypeII(aij(2,2),hij(2,2),sum(kt(1:2,t,:))).*ut_post(3,:,:)) .* lambda,...
          rS .* g .* (1-c) .* exp(-Func_TypeII(aij(2,2),hij(2,2),sum(kt(1:2,t,:))).*ut_post(3,:,:)) .* lambda,...
          zeros(1,1,RR);

          c .* rD .* exp(-Func_TypeII(aij(3,1),hij(3,1),kt(3,t,:)).*ut_post(2,:,:)),...
          c .* rD .* exp(-Func_TypeII(aij(3,1),hij(3,1),kt(3,t,:)).*ut_post(2,:,:)),...
          (1-d) .* rD .* exp(-Func_TypeII(aij(3,1),hij(3,1),kt(3,t,:)).*ut_post(2,:,:))]; 
        
% Next yrs numbers
    kt(:,t+1,:) = pagemtimes(Mk,kt(:,t,:));

% Add juveniles
    kt(1,t+1,:) = RKnew;
    

  
end

end
