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

tic
% Find the Kendall-tau values and associated quantifiers
[tau, z_score, p_value, H] = Modified_MannKendall_test(t, X, significance_value_tau, significance_value_ac);
t1 = toc;

tic
% The Optimized function uses gpuArrays to improve computation speeds, but requires the Parallel Computing or similar toolboxes
% Read Details in the function and to find the GPU shift vector length for your system run the codes in Testing_Codes folder
gpu_shift_critical_size = 550;
[tau_opt, z_score_opt, p_value_opt, H_opt] = Modified_MannKendall_test_Optimized(t, X, significance_value_tau, significance_value_ac, gpu_shift_critical_size);
t2 = toc;
%% Compare Modified_MannKendall_test with rank correlation 
n = 10001;
for i=1:10
X = rand(1, n);
Xsmoothed=movmean(X,100);
XTrend=X/10+(0:0.0001:1);
t = 1: n;

% Set significance values for hypothesis testing
significance_value_tau = 0.05;
significance_value_ac = 0.1;

% Find the Kendall-tau values and associated quantifiers
[tau, z_score, p_MMK, H] = Modified_MannKendall_test(t, X, significance_value_tau, significance_value_ac);
%[~,p_MMK]=Mann_Kendall_Modified(X,significance_value_tau);
[a,p_K]=corr(t',X','type','Kendall');

[tau, z_score, p_MMK_Smoothed, H] = Modified_MannKendall_test(t, Xsmoothed, significance_value_tau, significance_value_ac);
%[~,p_MMK_Smoothed]=Mann_Kendall_Modified(Xsmoothed,significance_value_tau);
[a,p_K_Smoothed]=corr(t',Xsmoothed','type','Kendall');

[tau, z_score, p_MMK_Trend, H] = Modified_MannKendall_test(t, XTrend, significance_value_tau, significance_value_ac);
%[~,p_MMK_Trend]=Mann_Kendall_Modified(XTrend,significance_value_tau);
[a,p_K_Trend]=corr(t',XTrend','type','Kendall');

figure
subplot(121)
plot(t,X);hold on
plot(t,Xsmoothed)
plot(t,XTrend)
legend('X','Smoothed X','Trend X');
xlabel('Sample');

subplot(122)
bar(1:6,[p_MMK p_K p_MMK_Smoothed p_K_Smoothed p_MMK_Trend p_K_Trend]);
xticks(1:6)
xticklabels({'Modified MannKendall Random','Rank Corr Random','Modified MannKendall Smoothed Random',...
    'Rank Corr Smoothed Random','Modified MannKendall Trend','Rank Corr Trend'});
xtickangle(45);
ylabel('P value')
pause
end


 
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

