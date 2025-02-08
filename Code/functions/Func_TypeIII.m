function Mort = Func_TypeIII(b,c,h,prey)

% Type III 
    % Mort = p.c. mortlaity of resource due to grazing
    % b = 
    % c = 
    % h = 1/handling time
    % prey = density of prey

        Mort = (b.*prey)./(1+(c.*prey)+(b.*prey^2)./h);

end    



% Plotting
%     % per predator
%     figure(1)
%     hold on
%     plot(1:100, TypeIII_func(0.9,2,100,1:100), LineWidth=1)
%     xlabel('Density of resource (kelp or drift)')
%     ylabel('p.c. mortality due to grazing')


