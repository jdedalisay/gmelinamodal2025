function objective = obj(x, exp_freqs)

    Ey = 10^(x(1) + 11); 
    Ex = 10^(x(2) + 11);  
    nuxy = x(3);        
    Gxy = 10^(x(4) + 11); 
    Gyz = 10^(x(5) + 11);  
    Gxz = 10^(x(6) + 11);  

    simulation_freqs = runSimulation([Ey, Ex, nuxy, Gxy, Gyz, Gxz]);

    % mean_exp = mean(exp_freqs);
    % std_exp = std(exp_freqs);
    % max_exp = exp_freqs*1.1;
    % min_exp = exp_freqs*0.9;
    % scaled_exp = (exp_freqs - mean_exp)/std_exp;
    % scaled_max = (max_exp - mean_exp)/std_exp;
    % scaled_min = (min_exp - mean_exp)/std_exp;
    % scaled_sim = (simulation_freqs - mean_exp)/std_exp;


%     % Weight factors (adjust as needed)
%     weight_factors = [1, 1, 1, 1, 1]; % For 5 frequencies
% 
%     % Check for negative frequencies
%     if any(simulation_freqs < 0)
%         objective = 100000; % Penalty for invalid frequencies
%     else
%         % Calculate percentage difference for the first 10 frequencies
%         percent_diff = abs((simulation_freqs - exp_freqs) ./ exp_freqs) * 100;
%         weighted_error = weight_factors .* percent_diff;
% 
%         % Calculate the objective as the root mean square error
%         objective = sqrt(mean(weighted_error.^2));
%     end
% end

% Check for negative frequencies
    if any(simulation_freqs < 0)
        objective = 100000; 
        return;
    end

 
    percent_diff = abs((simulation_freqs - exp_freqs) ./ exp_freqs)*100;
    % percent_diff = abs((scaled_sim - scaled_exp)./ scaled_exp)*100;
    

    maxdiff = max(percent_diff);
    if maxdiff < 10 
        objective = 0; 
    else
    
        objective = sqrt(mean(percent_diff.^2));
    end
    % 
    % ob = zeros(5,1);
    % for n = 1:5
    %     if scaled_sim(n) >= scaled_min(n) && scaled_sim(n) <= scaled_max(n)
    %         ob(n) = 0;
    %     else 
    %         ob(n) = 100;
    %     end
    % end
    % ob = max(ob);
    % if ob == 0
    %     objective = 0;
    % else if ob ~= 0
    %         objective = sqrt(mean(percent_diff.^2));
    % end
            
    
end


% verify the stopping criteria. So what kung objective is around 20%? it
% should continue until <10% nalang
