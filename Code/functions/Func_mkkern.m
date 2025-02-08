function kxy = Func_mkkern(fixparm)
% This function creates the IPM kernel.

% Adapted from Easterling et al. (2001) & Nickols et al. (2019)
% updated from mkkern to allow for vectorisation 
% (i.e., mutiple simultaneous replicates)
% updated for efficency

% RR = replicates

% Set up the integration mesh kernel ----------------------
    % For MLPA monitoring model
    y = fixparm.x;
    %this creates a vector (y) that is equal to x

    [x,y] = meshgrid(fixparm.x,y); 
    % Matlab built in function
    % x is an array (original x by original x) with each
    % row equal to originalvector x
    % y is an array (original y by original y) with each 
    % column equal to original vector y
    % X corresponds to size at time t
    % Y corresponds to size at time t+1


% Populate the kernel ---------------------
    % Define which sizes can be fished:
    isadult = normcdf(x,fixparm.Lfish,fixparm.Lvar); 

    % SURVIVAL PART OF KERNEL
        % natural mortality rate
        % process error
        % if putting this back in, then note that it chews up A LOT of
        % runtime
        % M = max([zeros(1,1,RR),normrnd(fixparm.M,fixparm.PES,1,1,RR)]);
        M = fixparm.M;
    
        % Mortality matrix (natural + fishing if fished age), size x by x by reps
        m = M + (isadult).*fixparm.F; 
        
        % convert mortality rate to survivorship
        % p1 size = length x by length y by reps
        p1 = exp(-m); 

    % GROWTH PART OF KERNEL
        % von Bertalanffy growth
        pmean1 = fixparm.Linf - (fixparm.Linf - x).*exp(-fixparm.k); 
        % add variability
        psig1 = pmean1*fixparm.Lvar; 
    
        % evaluate growth part of kernel
        % p2 size = length x by length y
        p2 = normpdf(y, pmean1, psig1);
        % p2 = repmat(p2,1,1,RR);

    % COMBINE GROWTH AND SURVIVAL
        % make sure no negatives
        p1 = max(0,p1);    
        p2 = max(0,p2);    
         
        % final output matrix
        kxy = p1.*p2;



