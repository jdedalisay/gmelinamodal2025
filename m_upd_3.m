clear; clc; 

% exp_freqs = [40.2,61.5, 153, 366.5, 537.7]; %SpecA
% exp_freqs = [38.8, 62.9, 159.5, 384.1, 579]; %SpecB
exp_freqs = [35.8, 62.1, 155.9, 375, 564.4]; %SpecC
% exp_freqs = [30.1, 61.2, 158, 368, 563]; %SpecD
% exp_freqs = [38.1, 62.2, 160.5, 367.5, 548.3]; %SpecF
% exp_freqs = [38.2, 64.1, 162.1, 385, 564]; %SpecH
% exp_freqs = [38.5, 61.9, 156, 360, 563]; %SpecI
% exp_freqs = [39.9, 63.5, 158.4, 387.5, 557.4]; %SpecJ

objective_function = @(x) obj(x, exp_freqs);

samp = 'LHS.xlsx';
lhs_samples = readtable(samp);
slhs = 206;
elhs = 500;
tally = slhs-1;
lhs_samp = lhs_samples(slhs:elhs,2:7);
initial_guesses = table2array(lhs_samp);

% lb = [-3, -3, 0.292, -3, -3, -3];
% ub = [-1, -1, 0.485, -1, -1, -1];
lb = [-3, -4, 0.292, -4, -3, -5]; 
ub = [-0, -2, 0.485, -1, -1, -2];

options = optimoptions(@fmincon, 'UseParallel', false, 'Display', 'iter', 'TolX', 1e-3,'TolFun', 1e-2);
%might be better w/o TolFun
ms = MultiStart('UseParallel', true, 'Display', 'iter');

results = table();
for i = 1:size(initial_guesses, 1)
    initial_guess = initial_guesses(i, :);
    e1 = 10^(initial_guesses(i,1)+11);
    e2 = 10^(initial_guesses(i,2)+11);
    n3 = initial_guesses(i,3);

    check = 1 -n3^2*e1/e2 -0.246^2*343713300/e2 -0.34^2*343713300/e2 -2*n3*0.246*0.34*343713300/e2;
    % disp(['Check: ', num2str(check)]);

    if check <= 0
        disp(['Skipped ', num2str(i+tally), ' because check: ', num2str(check)])
        continue;
    end
        
    problem = createOptimProblem('fmincon', 'objective', objective_function, ...
        'x0', initial_guess, 'lb', lb, 'ub', ub, 'options', options, 'nonlcon', @nu_ex);
    
    [optimized_x, fval, exitflag] = run(ms, problem, 1);
    
    % if isempty(optimized_x)
    %     disp(['Skipped ', num2str(i+tally), ' because optimized_x is empty']);
    %     continue;
    % end
    if exitflag > 0
        disp('Optimization was successful.');
    else
        disp(['Optimization at sample: ', num2str(i+tally), ' failed with exit flag: ', num2str(exitflag)]);
        continue;
    end

    optimized_params = 10.^(optimized_x); 
    optimized_params(1) = 10^(optimized_x(1) + 11);
    optimized_params(2) = 10^(optimized_x(2) + 11);
    optimized_params(3) = optimized_x(3);
    optimized_params(4) = 10^(optimized_x(4) + 11);
    optimized_params(5) = 10^(optimized_x(5) + 11);
    optimized_params(6) = 10^(optimized_x(6) + 11);

    results = [results; table(i, optimized_params(1), optimized_params(2), ...
                               optimized_params(3), optimized_params(4), ...
                               optimized_params(5), optimized_params(6), fval)];
    
    disp(['Initial Guess Row: ', num2str(i+tally)]);
    disp(['Optimized E_y: ', num2str(optimized_params(1)), ' Pa']);
    disp(['Optimized E_x: ', num2str(optimized_params(2)), ' Pa']);
    disp(['Optimized \nu_{xy}: ', num2str(optimized_params(3))]);
    disp(['Optimized G_{xy}: ', num2str(optimized_params(4)), ' Pa']);
    disp(['Optimized G_{yz}: ', num2str(optimized_params(5)), ' Pa']);
    disp(['Optimized G_{xz}: ', num2str(optimized_params(6)), ' Pa']);
    disp(['Objective function value: ', num2str(fval)]);
    disp(['Optimized Ortho Properties: ', num2str(optimized_x)]);
    disp('-----------------------------------');
end
%%
filename = 'LHS_ann.xlsx';
sheet = 'Sheet4'; % or name of the sheet
range = 'D11'; % starting cell
writetable(results, filename, 'Sheet', sheet, 'Range', range);