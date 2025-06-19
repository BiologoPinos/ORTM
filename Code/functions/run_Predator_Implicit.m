function [N_t,B_t] = run_Predator_Implicit(pred, T, RR, nt1, dist)

% runs IMP model for predators

% this version builds on v5 by adding in disturbances, management
% appraoches, and variability in early-life-history stages of kep & urchins

% Relevant to PredUrchinKelp_ImplicitCC_v6.m
% PredUrchinKelp_ImplicitCC_v6_MngtScen_v1.m

% Model run set-up --------------------
    % repilcated rho vector (pdf of recruits)
    Ro = repmat(pred.Ro',1,1,RR);

    % set noise vectors for recruits (timeseries of noise)
    % this needs to be done outside the loop becasue we set the same variance for
    % the whole year (= 4 time steps)
    R_noise = repelem(max(zeros(T/4,1,RR),normrnd(1,pred.Rstdv,T/4,1,RR)),4,1,1);

    % set up matrix to hold model runs 
    % size = y(size classes),t(time),RR(reps)
    N_t = nan(length(pred.x),T,RR);

    % set temporary predator parameters
    predTemp = pred;

    % set intial conditions (size distributions)
    if size(nt1,2) == 1
        N_t(:,1,:) = repmat(nt1,1,RR);
    else
        N_t(:,1,:) = nt1;
    end

    N_R = N_t;
    N_F = N_t;

% Model run ---------------
    for t = 1:T %loop over time steps  

    % Temporarily change fishing pressure
        if pred.fish == "Y"
            if ismember(t, dist.yrs(1) + pred.fishtime + (0:pred.fishlgth-1))
                predTemp.F = pred.fishF;
            else
                predTemp.F = pred.F;
            end
        end

    % GROWTH & SURVIVAL
        % create the kernel that caculates the probability of 
        % growing and moving from size x to y
        % (this is where fishing is included)
        % size = mesh size x mesh size x reps
        % this is determinisic 
        % (if need to add variability to adult survival it happens here, 
        % using vectorised version, it just takes x4 to run)
        kmat = Func_mkkern(predTemp);
              
        % weight by midpoint rule (even weighting)
        kmat = pred.dx .* kmat;

        % replicate
        kmat = repmat(kmat,1,1,RR);


    % RECRUITMENT
        % Open recrutiment
            % Add recrtuiment function (for time of year they recruit) 
            % and variability to recruitment (same for both pops)
              Radd = pred.Ro .* pred.Rmean .* pred.RoT(t).*R_noise(t,:,:);
            % testing - no recruitment
%               Radd = zeros(1,1,RR);

              
    % CALC NEXT TIME STEP
        % keep >= 0 (issue from + var recruits)
        N_t(:,t+1,:) = pagemtimes(kmat,N_t(:,t,:)) + pagetranspose(Radd);
        N_t(:,t+1,:) = max(0,N_t(:,t+1,:)); 

    end

        % biomass (kg)
        B_t = N_t.*repmat(Func_Length2Weight(pred.x',pred.y,pred.z),1,T+1,RR);

end
  