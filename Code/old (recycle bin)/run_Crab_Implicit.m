function [cft,cmt] = run_Crab_Implicit(crab, tmax, RR, cft1, cmt1, nt)
% function [cft,cmt, RK_noise] = run_Crab_Implicit(crab, tmax, RR, cft1, cmt1, nt)


% URCHIN (parameter values)
MJ = crab.MJ;
MA = crab.MA;
F  = crab.F;
P  = crab.P;
H  = crab.H;

%% Vectors for state variables

    % predators
    nt = reshape(nt,[],1,RR);

    % Dungeness crabs females
    cft = NaN(6,tmax+1,RR);
    cft(:,1,:) = repmat(cft1(:),1,1,RR);

    % Dungeness crabs males
    cmt = NaN(6,tmax+1,RR);
    cmt(:,1,:) = repmat(cmt1(:),1,1,RR);

    % % noise vector for recruitment
    %     "Think about this"

% % preset demographic variables:
% Fec_vec = [0, 0, ]; % row vector giving the per capita fecundity of each female age class
% M_infl = []; % row vector giving the per capita cannibalism effect of each male age class
% F_infl = []; % same for females. Probably M_infl and F_infl are the same.


% could also make the Leslie matrices dynamic (changing with time inside the loop) if you
% want to run the model to a steady state and then add otters, etc.

% Initialize values of R, L, cft, cmt

% Need to think about the timing of the annual time step relative to the
% fishing season + reproductive season. Perhaps best to imagine that the
% end of each time step is the close of the fishing season. Settlement &
% recruitment happen over the summer. So Spawning would depend on the
% number of females in the previous time step (t-1). Then advance the male +
% female population, and cannibalism on the new recruits would depend on
% abundance in the *current* time step (t), reflecting the mortality over
% the winter (especially the fishing mortality on males, which is heaviest
% in the winter).

%% Run over time

for t = 1:tmax


%% Run Dungeness

% % reproduction (assumes males are not limiting)
% L(t) = Fec_vec*cft(:,t-1); % matrix operation should produce a scalar
% 
% % could add lognormal variability at this point to reflect ocean conditions
% 
% % post-settlement mortality
% DDmort = M_infl*cmt(:,t-1) + F_infl*cft(:,t-1); % should produce a scalar
% R(t) = L(t)*exp(-DDmort);


% Survival

    % age 0 (natural moratlity + cannibalism)
        sC0 = repmat(exp(-MJ), 1, 1, RR); % density independent

    % age 1 to age +5 female(natural mortality + type II predation)
        sCF6 = exp(-MA -Func_TypeII(P,H,sum(cft(2:6, t, :), 1)).*nt(t,:,:));

    % age 1 to age 3 male (natural mortality + type II predation)
        sCM3 = exp(-MA -Func_TypeII(P,H,sum(cmt(2:4, t, :), 1)).*nt(t,:,:));

    % age 4 to age +5 male (natural moratlity + type II predation + fishery)
        sCM6 = exp(-(MA + F(t) + Func_TypeII(P, H, sum(cmt(5:6, t, :), 1)) .* nt(t, :, :)));

% Projection matrix (Dungeness crab females) (z = zeros(1,1,RR);)
Mcf = [sC0, zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR);

       zeros(1,1,RR), sCF6, zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR);

       zeros(1,1,RR), zeros(1,1,RR), sCF6, zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR);

       zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR), sCF6, zeros(1,1,RR), zeros(1,1,RR);

       zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR), sCF6, zeros(1,1,RR);

       zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR), sCF6;]; % matrix giving the annual survival for each age class (F) 

% projection matrix (Dungeness crab males) (z = zeros(1,1,RR);)
Mcm = [sC0, zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR);

       zeros(1,1,RR), sCM3, zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR);

       zeros(1,1,RR), zeros(1,1,RR), sCM3, zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR);

       zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR), sCM3, zeros(1,1,RR), zeros(1,1,RR);

       zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR), sCM6, zeros(1,1,RR);

       zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR), zeros(1,1,RR), sCM6;]; % matrix giving the annual survival for each age class (M, includes fishing)

% advance the adult populations
cft(:,t+1,:) = pagemtimes(Mcf,cft(:,t,:));
cmt(:,t+1,:) = pagemtimes(Mcm,cmt(:,t,:));

% cft(:,t) = Mcf*cft(:,t-1);
% cmt(:,t) = Mcm*cft(:,t-1);

% % add in new recruits
% cft(1,t) = R(t)/2;
% cmt(1,t) = R(t)/2;



end

end