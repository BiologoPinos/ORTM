function Phi = Func_Switch(w1, w2, k)

% Behaviour switching
    % Phi = proportions of urchins exposed garzing kelp
    % w1 = inflection point
    % w2 = slope around inflection point
    % k = density of drift or kelp 

    % v1 & v2 = constants that are solved for 
        % older slower version, requires defining symbols
        % syms v2s 
        % v2 = vpasolve((v2s - 1) * exp((1-v2s)/v2s) == w2 * w1);

    % speedier version 
    v2 = fsolve(@(v2s)(v2s - 1) * exp((1-v2s)/v2s) - w2 * w1,1, optimset('Display','off'));
    v1 = (v2 - 1)/(v2 * w1^v2);

    Phi = exp(-v1.*k.^v2);

end



% Plotting
    % figure
    % hold on 
    % plot(0:0.01:5, Func_Switch(1, 1, 0:0.01:5), LineWidth=1)
    % ylabel('Proportion of urchins grazing standing kelp')
    % xlabel('Ratio of drift to urchin consumptive capacity')