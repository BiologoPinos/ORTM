function [kt,ut,GC,RK_noise] = run_UrchinKelp_Implicit(kelp, urchin, tmax, RR, kt1, ut1, nt, dist)

% Description:
% function running single populations kelp <> urchins model
% see parent file (Para) for parameter values

% Runs with:
% ORTM_model_otter.m


%% Unpack structures--------------------------
    % note, some parts of structure are not unpacked until later

% KELP
RKstdv = kelp.RKstdv;
mu = kelp.mu;
ddD = kelp.ddD;
muvar = kelp.muvar;
RTk = kelp.RTk;
reproWeight = kelp.reproWeight;  
lag = kelp.lag;          
g = kelp.g;
rS = kelp.rS;
c = kelp.c;
d = kelp.d;
rD = kelp.rD;
aij = cell2mat(kelp.aij);
%bhij = cell2mat(kelp.bhij); % being unpacked later
%hij = cell2mat(kelp.hij); % being unpacked later
%lambda = kelp.lambda; % being unpacked later

% URCHIN
RU = urchin.RU;
RUstdv = urchin.RUstdv;
RUdist = urchin.RUdist;
RTu = urchin.RTu;
gJ = urchin.gJ;
MJ = urchin.MJ;
alpha = urchin.alpha;
alphavar = urchin.alphavar;
MH = urchin.MH;
ME = urchin.ME;
PH = urchin.PH;
PE = urchin.PE;
F = urchin.F;
PLD = urchin.PLD;
tau = urchin.tau;
w1 = urchin.w1;
w2 = urchin.w2;
kmin = urchin.kmin;

dist_yrs = dist.yrs;

    
%% Vectors for state variables--------------------------

    % predators
    nt = reshape(nt,[],1,RR);
    
    % kelp
    kt = NaN(3,tmax,RR);
    kt(:,1,:) = repmat(kt1',1,1,RR);
    
    % urchins
    ut = kt;
    ut(:,1,:) =  repmat(ut1',1,1,RR);  

    % Grazing capacity
    GC = NaN(tmax,1,RR);
    
% set Psi (switching function) to 1
    Psi = 1;

% noice vector for recruitment    
    % set noise vectors for recruits (timeseries of noise)
    % this needs to be done outside the loop becasue we set the same variance for
    % the whole year (= 4 time steps), which is what the repelem() does

    % kelp
        RK_noise = repelem(max(zeros(tmax/4,1,RR),normrnd(1,RKstdv,tmax/4,1,RR)),4,1,1);
    % urchins
        % RU_noise = repelem(max(zeros(tmax/4,1,RR),normrnd(1,RUstdv,tmax/4,1,RR)),4,1,1);
        if strcmp(RUdist{1}, 'normal')
            RU_noise = repelem(max(zeros(tmax/4,1,RR),normrnd(1,RUstdv,tmax/4,1,RR)),4,1,1);
        else 
            RU_noise = repelem(max(zeros(tmax/4,1,RR),lognrnd(1,RUstdv,tmax/4,1,RR)),4,1,1);
        end        


%% Run over time--------------------------

for t = 1:tmax
    % add disturbance effects if happens
    if ismember(t, dist_yrs)
        % kelp recruitment
        RK = dist.RK;
        % change in standing kelp growth
        lambda = dist.lambda(t);
        % urchin grazing rates
        hij = dist.hij(:,:,t); 
    else 
        RK = kelp.RK;
        lambda = kelp.lambda(t);
        hij = cell2mat(kelp.hij);
        hij = hij(:,:,t);  
    end


%% Run urchins🟣--------------------------
    % set earlier values for delay
    % if for the last 6 months (>3months) of kelp biomass are less than the
    % threshold (each year), then predators stop feeding (and urchin fishing turns off)
    % Psi is 1 until  k < kmin for consecutive time periods

if t>2 
    ktdelay = sum(kt(1:2,(t-2:t-1),:));        
    Psi = sum(ktdelay < kmin) ~= 2;
end

% survival
    % juvs
        % density-independent survival
        sJ = exp(-MJ);
        % density-dependent survival (more adults = more survival)
        % sJ = sJ.*(1-exp(-alpha*(ut(2,t)+ut(3,t))));

        % with scaling transition
        % sJ = sJ.*(1-exp(-alpha*(ut(2,t)+ut(3,t))) - 0.5*alpha^2 * exp(-alpha*(ut(2,t)+ut(3,t))) * alphavar);
        % give error if neg survival and set sJ to zero if that is
        % the case (this means that for high variances, such as
        % what we see, the landscape-scale survival will rapidly
        % drop off at moderate-low urchins densities
        % if sJ<0; warning('Urchin survival < 0'); end
        if sJ<0; sJ = 0; end

        % e.g. plot for scaling transition
        % figure
        % hold on
        % plot(1:10000, sJ.*(1-exp(-alpha*(1:10000)))) % basic IDD
        % plot(1:10000, sJ.*(1-exp(-alpha*(1:10000)))- 0.5*alpha^2 * exp(-alpha*(1:10000)) * alphavar % scaling transition
        
    % % type I linear predation (e.g. sheep head)
    %     % hiding adults 
    %     sH = exp(-MH -PH.*nt(t,:,:) - F);
    %     % exposed adults
    %     sE = exp(-ME -Psi.*(PE.*nt(t,:,:) + F));

    % type II predation (e.g. sea otters)
        % hiding adults 
        sH = exp(-MH -F -Func_TypeII(PH,10^3,sum(ut(2:3,t,:))).*nt(t,:,:));
        % exposed adults
        sE = exp(-ME -Psi .* (F + Func_TypeII(PE,10^3,sum(ut(2:3,t,:))).*nt(t,:,:)));

% prop exposed
    % grazing capacity 
    % (ratio of drift kelp to total urchin feeding rate)
    % use max so that NaN (no drift state) becomes 0
    GC(t,:) = max(kt(3,t,:)./(sum(ut(2:3,t,:)).*hij(3,1)),0);   

    % Phi = Func_Switch(w1, w2, kt(3,t));
    Phi = reshape(Func_Switch(w1, w2, GC(t,:)),1,1,RR);

% projection matrix
Mu = [repmat((1-gJ) * sJ,1,1,RR),  zeros(1,1,RR),   zeros(1,1,RR);
      gJ .* sH .* (1-Phi),         sH .* (1-Phi),   sH .* (1-Phi);
      gJ .* sE .* Phi,             sE .* Phi,       sE .* Phi];

% incoming recruits w/ timing + noise
RUadd = zeros(3,1,RR);
RUadd(1,1,:) = RU .* RTu(t) .* RU_noise(t,:,:);

% next yrs numbers
ut(:,t+1,:) = pagemtimes(Mu,ut(:,t,:)) + (RUadd .* sJ^tau);

% urchin culling (mass mortality events)
if urchin.culling == "Y"
    if ismember(t, dist_yrs(1) + urchin.culltime + (0:urchin.culllgth-1))
        % ut(2:3,t+1,:) = max(ut(2:3,t+1,:) - [urchin.culln*0.5; urchin.culln*0.5],[0;0]);
        % cull all exposed first, then the hiding urchins
            % which reps have more exposed than the cull number
            exp_reps = find(ut(3,t+1,:) > urchin.culln)';
            % ...and which dont
            both_reps = find(ut(3,t+1,:) <= urchin.culln)';
            % for exp_reps, just cull from exposed
            ut(3,t+1,exp_reps) = ut(3,t+1,exp_reps) - urchin.culln;
            % for both_reps, cull all from exposed and rest from hiding
            ut(3,t+1,both_reps) = 0;
            ut(2,t+1,both_reps) = max(ut(2,t+1,both_reps)-ut(3,t+1,both_reps),0);
    end
end

    
%% Run kelp 🌿--------------------------
   
% kelp restoration
    % adding juvenile kelp biomass 
    RKr = 0;
    if kelp.restore == "Y"
        if ismember(t, dist_yrs(1) + kelp.resttime + (0:kelp.restlgth-1))
            RK = RK + kelp.restn;
        end
    end
    
% per capita settling spores + timing + noise from adults
RKset = RK .* RTk(t) .* RK_noise(t,:,:);
% % RKset = RK .* RTk(t) .* (sum(reproWeight .* kt(2,max(1,t-lag),:))) .* RK_noise(t,:,:);
        
%% Problem is here (probably)
% total incoming setted sporophyte biomass (before density dependence)
% spores from adults + restored juvs
% % ksetT = RKset .* kt(2,t,:) + RKr;
ksetT = RKset .* sum(reproWeight .* kt(2,max(1,t-lag),:)) + RKr;

     % DD survival of young
            % to be used for intra-cohort part of DD
            % if ksetn = 0 & k2n = 0 then the DD function is NaN, and if <1 we get wierd survival rates, 
            % so we will set 1 instead. This wont affect dynamics, since it is then x 0. 
            % NOTE TO SELF: need to make proportional density in 4+ popmodel
            ksetn = max(ksetT, 1);

            % Adult densities to be used for inter-cohort part of DD
            % NOTE TO SELF: need to make proportional density in 4+ popmodel
            % % k2n = max(kt(2,t,:), 1);
            k2n = max(sum(reproWeight .* kt(2,max(1,t-lag),:)), 1);

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
  
% total new biomass of juveniles
    RKnew = ksetT .* sY .* rS .* exp(-Func_TypeII(aij(1,2),hij(1,2),sum(kt(1:2,t,:))).*ut(3,t+1,:)) .* lambda;

% projection matrix (bull Kelp)
    Mk = [zeros(1,1,RR),...
          zeros(1,1,RR),...
          zeros(1,1,RR);

          rS .* g .* (1-c) .* exp(-Func_TypeII(aij(2,2),hij(2,2),sum(kt(1:2,t,:))).*ut(3,t+1,:)) .* lambda,...
          rS .* g .* (1-c) .* exp(-Func_TypeII(aij(2,2),hij(2,2),sum(kt(1:2,t,:))).*ut(3,t+1,:)) .* lambda,...
          zeros(1,1,RR);

          c .* rD .* exp(-Func_TypeII(aij(3,1),hij(3,1),kt(3,t,:)).*ut(2,t+1,:)),...
          c .* rD .* exp(-Func_TypeII(aij(3,1),hij(3,1),kt(3,t,:)).*ut(2,t+1,:)),...
          (1-d) .* rD .* exp(-Func_TypeII(aij(3,1),hij(3,1),kt(3,t,:)).*ut(2,t+1,:))];            
            
    % if Mk(3,1,:)>1; warning('Drift survival >1'); end
        
% next yrs numbers
    kt(:,t+1,:) = pagemtimes(Mk,kt(:,t,:));
    % This line of code will not work because of
    % kt(:,t+1,:) = pagemtimes(Mk,kt(:,t,:)) + [RKnew;0;0];


% add juveniles
    kt(1,t+1,:) = RKnew;
         
  
end

end
