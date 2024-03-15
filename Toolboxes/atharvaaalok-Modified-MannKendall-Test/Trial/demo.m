%% DEMONSTRATON OF USING THE FUNCTIONS
% Modified_MannKendall_test_optimized - This function uses gpuArrays to optimize sen-slope calculations
clear; clc; close all;


%% GENERATE TIME SERIES TO FIND THE KENDALL TAU VALUE OF

n = 1000;
X = rand(1, n);
t = 1: n;

% Set significance values for hypothesis testing
significance_value_tau = 0.05;
significance_value_ac = 0.05;

N_loop = 10000;

time_taken_MMK = [];
time_taken_MMK_Opt = [];

for i = 1: N_loop
    tic
    % Find the Kendall-tau values and associated quantifiers
    [tau, z_score, p_value, H, time_taken_MMK] = Modified_MannKendall_test(t, X, significance_value_tau, significance_value_ac, time_taken_MMK);
    t1 = toc;
    
    tic
    % The Optimized function uses gpuArrays to improve computation speeds, but requires the Parallel Computing or similar toolboxes
    % Read Details in the function and to find the GPU shift vector length for your system run the codes in Testing_Codes folder
    gpu_shift_critical_size = 550;
    [tau_opt, z_score_opt, p_value_opt, H_opt, time_taken_MMK_Opt] = Modified_MannKendall_test_Optimized(t, X, significance_value_tau, significance_value_ac, gpu_shift_critical_size, time_taken_MMK_Opt);
    t2 = toc;

end


%% CALCULATE AVERAGE

a = mean(time_taken_MMK(:, 1));
b = mean(time_taken_MMK(:, 2));
c = mean(time_taken_MMK(:, 3));
d = mean(time_taken_MMK(:, 4));
e = mean(time_taken_MMK(:, 5));
f = mean(time_taken_MMK(:, 6));

a_opt = mean(time_taken_MMK_Opt(:, 1));
b_opt = mean(time_taken_MMK_Opt(:, 2));
c_opt = mean(time_taken_MMK_Opt(:, 3));
d_opt = mean(time_taken_MMK_Opt(:, 4));
e_opt = mean(time_taken_MMK_Opt(:, 5));
f_opt = mean(time_taken_MMK_Opt(:, 6));

% Print results for the unoptimized Mann-Kendall test
fprintf('Unoptimized and Optimized Comparison\n');
fprintf('a = %f, a_opt = %f\n', a, a_opt);
fprintf('b = %f, b_opt = %f\n', b, b_opt);
fprintf('c = %f, c_opt = %f\n', c, c_opt);
fprintf('d = %f, d_opt = %f\n', d, d_opt);
fprintf('e = %f, e_opt = %f\n', e, e_opt);
fprintf('f = %f, f_opt = %f\n', f, f_opt);
fprintf('\n');


%% PRINT RESULTS

% Print results for the unoptimized Mann-Kendall test
fprintf('Unoptimized Mann-Kendall Function\n');
fprintf('tau \t\t= %.4f\n', tau);
fprintf('z_score \t= %.4f\n', z_score);
fprintf('p_value \t= %.8f\n', p_value);
fprintf('H \t\t\t= %d\n', H);
fprintf('t1 \t\t\t= %.5f\n', t1);
fprintf('\n');

% Print results for the optimized Mann-Kendall test
fprintf('Optimized Mann-Kendall Function\n');
fprintf('tau \t\t= %.4f\n', tau);
fprintf('z_score \t= %.4f\n', z_score);
fprintf('p_value \t= %.8f\n', p_value);
fprintf('H \t\t\t= %d\n', H);
fprintf('t2 \t\t\t= %.5f\n', t2);
fprintf('\n');

