function [kt,ut,GC,cft,cmt,RK_noise] = run_UrchinKelpCrab_Implicit(kelp, urchin, crab, tmax, RR, kt1, ut1, cft1, cmt1, nt, dist)

% Description:
    % function running single populations kelp <> urchins model
    % see parent file (Para) for parameter values

% Runs with:
    % ORTM_model_otter.m


%% UNPACK STRUCTURES --------------------------
    % (note: some structures are unpacked later)

% Kelp 🌿
RKstdv      = kelp.RKstdv;
mu          = kelp.mu;
ddD         = kelp.ddD;
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
alpha       = urchin.alpha;
alphavar    = urchin.alphavar;
MHU         = urchin.MHU;
MEU         = urchin.MEU;
PH          = urchin.PH;
PE          = urchin.PE;
HH          = urchin.HH;
HE          = urchin.HE;
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
PC          = crab.PC;
HC          = crab.HC;
NCIF  = crab.NCIF{1};

% Other parameters
dist_yrs = dist.yrs;


%% VECTORS FOR STATE VARIABLES --------------------------

% Predator
    nt = reshape(nt,[],1,RR);
    % nt = nt([],1,RR);

% Kelp
    kt = NaN(3,tmax,RR);
    kt(:,1,:) = repmat(kt1(:),1,1,RR);
    
% Urchin
    ut = NaN(3,tmax,RR);
    ut(:,1,:) =  repmat(ut1(:),1,1,RR);

% Grazing capacity
    GC = NaN(tmax,1,RR);

% Dungeness crab females
    cft = NaN(11,tmax,RR);
    cft(:,1,:) = repmat(cft1(:),1,1,RR);

% Dungeness crab males
    cmt = NaN(11,tmax,RR);
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

% Crab recruitment noise
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
            cft_d = NaN(11, ndays+1, RR);
            cft_d(:,1,:) = cft(:,t,:);
        
            % Dungeness crab males
            cmt_d = NaN(11, ndays+1, RR);
            cmt_d(:,1,:) = cmt(:,t,:);   
    
            % Predator (same number)
            nt_season = reshape(nt(t,:,:), 1, 1, RR); %% CHEKC THIS

        % Run daily predation loop
        
            for tt = 1:ndays
                
                % Daily predation on urchins (juveniles unaffected by predation)
                ut_d(1,tt+1,:) = ut_d(1,tt,:);
                ut_d(2,tt+1,:) = ut_d(2,tt,:) .* exp(-Func_TypeII(PH, HH, sum(ut_d(2:3,tt,:))) .* nt_season(:,:,:)); % hiding
                ut_d(3,tt+1,:) = ut_d(3,tt,:) .* exp(Psi.* -Func_TypeII(PE, HE, sum(ut_d(2:3,tt,:))) .* nt_season(:,:,:)); % exposed
        
                % Daily predation on female crabs (age 0-1 unaffected by predation)
                cft_d(1:2,tt+1,:) = cft_d(1:2,tt,:);
                cft_d(3:end,tt+1,:) = cft_d(3:end,tt,:) .* exp(-Func_TypeII(PC, HC, sum(cft_d(3:end,tt,:) + cmt_d(3:end,tt,:))) .* nt_season(:,:,:));
        
                % Daily predation on male crabs (age 0-1 unaffected by predation)
                cmt_d(1:2,tt+1,:) = cmt_d(1:2,tt,:);
                cmt_d(3:end,tt+1,:) = cmt_d(3:end,tt,:) .* exp(-Func_TypeII(PC, HC, sum(cmt_d(3:end,tt,:) + cft_d(3:end,tt,:))) .* nt_season(:,:,:));
                
        
            end % end daily loop
    
        % Extract post-predation (end-of-season) biomasses
        ut_post = ut_d(:, end, :);   % 3 x 1 x RR
        % ut_post(:,t+1,:) = ut_d(:, end, :);   % 3 x 1 x RR
        cft_post = cft_d(:, end, :); % 11 x 1 x RR
        cmt_post = cmt_d(:, end, :); % 11 x 1 x RR


%% RUN CRAB (Dungeness) 🦀 --------------------------

% Survival

    % Females 
       Sf = cat(1, repmat(exp(-MRC),1,1,RR), repmat(exp(-MJC),1,1,RR), repmat(exp(-MAC), 11-3,1,RR));
    
    % Males
       Sm = cat(1, repmat(exp(-MRC),1,1,RR), repmat(exp(-MJC),1,1,RR), repmat(exp(-MAC), 2,1,RR), repmat(exp(-MAC - FC(t)), 11-5,1,RR));

% Projection/Transition matrix
    
    % Females
        Mcf = zeros(11, 11, RR);
        Mcf(2:11, 1:11-1, :) = reshape(eye(11-1), 11-1, 11-1, 1) .* reshape(squeeze(Sf(:,1,:)), 11-1, 1, RR); 
            % 10x10xRR stored in the sub-diagonal of an 11x11xRR

    % Males
        Mcm = zeros(11, 11, RR);
        Mcm(2:11, 1:11-1, :) = reshape(eye(11-1), 11-1, 11-1, 1) .* reshape(squeeze(Sm(:,1,:)), 11-1, 1, RR); 
            % 10x10xRR stored in the sub-diagonal of an 11x11xRR

% Recruitment

    % Incoming recruits with timing + noise
        RCnew = zeros(11,1,RR);
        RCnew(1,1,:) = RC .* RTC(t) .* RC_noise(t,:,:);
    
    % Cannibalism (Ricker Style)
        RCnew(1,1,:) = RCnew(1,1,:) .* exp(-beta .* (sum(NCIF .* (cft_post(2:11,:,:) + cmt_post(2:11,:,:)))));

% Advance the adult populations
    cft(:,t+1,:) = pagemtimes(Mcf, cft_post) + (RCnew .* 0.5); % females
    cmt(:,t+1,:) = pagemtimes(Mcm, cmt_post) + (RCnew .* 0.5); % males 
        cmt(6:11, t+1, :) = 0; % Fishery: remove all males age ≥5 (set to zero)



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
            % sH = exp(-MH -PH.*nt(t,:,:) - FU);
            % exposed adults
            % sE = exp(-ME -Psi.*(PE.*nt(t,:,:) + FU));
    
        % type II predation (e.g. sea otters)
            % hiding adults 
            sH = exp(-MHU); % exp(-MHU -FU -Func_TypeII(PH,HH,sum(ut(2:3,t,:))).*nt(t,:,:));
            % sH = exp(-MHU - FU); % exp(-MHU -FU -Func_TypeII(PH,HH,sum(ut(2:3,t,:))).*nt(t,:,:));

            % exposed adults
            sE = exp(-MEU); % exp(-MEU -Psi .* (FU + Func_TypeII(PE,HE,sum(ut(2:3,t,:))).*nt(t,:,:)));
            % sE = exp(-MEU - Psi .* FU); % exp(-MEU -Psi .* (FU + Func_TypeII(PE,HE,sum(ut(2:3,t,:))).*nt(t,:,:)));

% Proportion being exposed

    % Grazing capacity (ratio of drift kelp to total urchin feeding rate) - using post predation adults 
    GC(t,1,:) = max(kt(3,t,:)./(sum(ut_post(2:3,:,:)).*hij(3,1)),0);     % use max so that NaN (no drift state) becomes 0

    % Phi = Func_Switch(w1, w2, kt(3,t));
    Phi = reshape(Func_Switch(w1, w2, squeeze(GC(t,1,:))),1,1,RR);

    % % Grazing capacity (ratio of drift kelp to total urchin feeding rate) - using pre predation adults 
    % GC(t,1,:) = max(kt(3,t,:)./(sum(ut(2:3,t,:)).*hij(3,1)),0);     % use max so that NaN (no drift state) becomes 0
    % 
    % % Phi = Func_Switch(w1, w2, kt(3,t));
    % Phi = reshape(Func_Switch(w1, w2, squeeze(GC(t,1,:))),1,1,RR);

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
        % to be used for intra-cohort part of DD
        % if ksetn = 0 & k2n = 0 then the DD function is NaN, and if <1 we get wierd survival rates, 
        % so we will set 1 instead. This wont affect dynamics, since it is then x 0. 
        % NOTE TO SELF: need to make proportional density in 4+ popmodel
        ksetn = max(ksetT, 1);

        % Adult densities to be used for inter-cohort part of DD
        % % k2n = max(kt(2,t,:), 1);
        k2n = max(sum(reproWeight .* kt(2,max(1,t-RKlag),:)), 1);

        % DD functions

            % for Beverton-Holt DD function 
            % (which doesnt use scale transitions, since only settlers)
            % note this basic form assumes that the slope at the orgin = 1
            if ddD == 1
                sY = 1/(1+ksetn./mu);
    
            else % for mixed and ricker DD functions
                DD = ( (ddD-1).*k2n.*exp(mu.*ddD.*k2n) ) ./ (ddD.*ksetn.*exp(mu.*ddD.*k2n) + exp(mu.*k2n) .* ((ddD-1).*k2n-ddD.*ksetn)  );
    
                % second derivative of mixed DD function
                % this has been calculted using the matlab solver, see ExploringMixedDD_v0.m 
                DD2 = (2.*k2n.*exp(ddD.*mu.*k2n).*(ddD - 1).*(exp(mu.*k2n).*(ddD - 1) + mu.*exp(mu.*k2n).*(k2n.*(ddD - 1) - ddD.*ksetn) + ddD^2.*ksetn.*mu.*exp(ddD.*mu.*k2n)).^2)./(exp(mu.*k2n).*(k2n.*(ddD - 1) - ddD.*ksetn) + ddD.*ksetn.*exp(ddD.*mu.*k2n)).^3 - (k2n.*exp(ddD.*mu.*k2n).*(ddD - 1).*(2.*mu.*exp(mu.*k2n).*(ddD - 1) + mu.^2.*exp(mu.*k2n).*(k2n.*(ddD - 1) - ddD.*ksetn) + ddD.^3.*ksetn.*mu.^2.*exp(ddD.*mu.*k2n)))./(exp(mu.*k2n).*(k2n.*(ddD - 1) - ddD.*ksetn) + ddD.*ksetn.*exp(ddD.*mu.*k2n)).^2 - (2.*exp(ddD.*mu.*k2n).*(ddD - 1).*(exp(mu.*k2n).*(ddD - 1) + mu.*exp(mu.*k2n).*(k2n.*(ddD - 1) - ddD.*ksetn) + ddD.^2.*ksetn.*mu.*exp(ddD.*mu.*k2n)))./(exp(mu.*k2n).*(k2n.*(ddD - 1) - ddD.*ksetn) + ddD.*ksetn.*exp(ddD.*mu.*k2n)).^2 + (2.*ddD.*mu.*exp(ddD.*mu.*k2n).*(ddD - 1))./(exp(mu.*k2n).*(k2n.*(ddD - 1) - ddD.*ksetn) + ddD.*ksetn.*exp(ddD.*mu.*k2n)) + (ddD.^2.*mu.^2.*k2n.*exp(ddD.*mu.*k2n).*(ddD - 1))./(exp(mu.*k2n).*(k2n.*(ddD - 1) - ddD.*ksetn) + ddD.*ksetn.*exp(ddD.*mu.*k2n)) - (2.*ddD.*mu.*k2n.*exp(ddD.*mu.*k2n).*(ddD - 1).*(exp(mu.*k2n).*(ddD - 1) + mu.*exp(mu.*k2n).*(k2n.*(ddD - 1) - ddD.*ksetn) + ddD^2.*ksetn.*mu.*exp(ddD.*mu.*k2n)))./(exp(mu.*k2n).*(k2n.*(ddD - 1) - ddD.*ksetn) + ddD.*ksetn.*exp(ddD.*mu.*k2n)).^2;
                
                % per capita survival of settlers including scale transition
                sY = DD + 0.5.*DD2.*muvar;
            end            
  
    % Total new biomass of juveniles - post predation urchins
        RKnew = ksetT .* sY .* rS .* exp(-Func_TypeII(aij(1,2),hij(1,2),sum(kt(1:2,t,:))).*ut_post(3,:,:)) .* lambda;
    % % Total new biomass of juveniles - pre predation urchins
    %     RKnew = ksetT .* sY .* rS .* exp(-Func_TypeII(aij(1,2),hij(1,2),sum(kt(1:2,t,:))).*ut(3,t+1,:)) .* lambda;

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
% % Projection/Transition matrix (bull Kelp)- pre predation urchins
%     Mk = [zeros(1,1,RR),...
%           zeros(1,1,RR),...
%           zeros(1,1,RR);
% 
%           rS .* g .* (1-c) .* exp(-Func_TypeII(aij(2,2),hij(2,2),sum(kt(1:2,t,:))).*ut(3,t+1,:)) .* lambda,...
%           rS .* g .* (1-c) .* exp(-Func_TypeII(aij(2,2),hij(2,2),sum(kt(1:2,t,:))).*ut(3,t+1,:)) .* lambda,...
%           zeros(1,1,RR);
% 
%           c .* rD .* exp(-Func_TypeII(aij(3,1),hij(3,1),kt(3,t,:)).*ut(2,t+1,:)),...
%           c .* rD .* exp(-Func_TypeII(aij(3,1),hij(3,1),kt(3,t,:)).*ut(2,t+1,:)),...
%           (1-d) .* rD .* exp(-Func_TypeII(aij(3,1),hij(3,1),kt(3,t,:)).*ut(2,t+1,:))];            
            
    % if Mk(3,1,:)>1; warning('Drift survival >1'); end
        
% Next yrs numbers
    kt(:,t+1,:) = pagemtimes(Mk,kt(:,t,:));

% Add juveniles
    kt(1,t+1,:) = RKnew;
    

  
end

end
