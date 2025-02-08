function Mort = Func_TypeII(a,h,prey)
% Type II predation/grazing mortlaity

    % Mort = p.c. mortlaity of resource due to grazing
    % a = attack rate
    % h = 1/handling time (max prey consumed)
    % prey = density of prey

Mort = a./(1+(a.*prey)./h);


end



% Plotting
    % % per predator
    % figure(1)
    % hold on
    % plot(1:100, Func_TypeII(0.9,100,1:100), LineWidth=1)
    % xlabel('Density of resource (kelp or drift)')
    % ylabel('p.c. mortality due to grazing')