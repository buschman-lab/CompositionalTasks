N_sims = 250;
samp_sizes = [25:25:250];
N_resamples = 1000;

wh = waitbar(0, 'running...');
samp_p = NaN*ones(N_sims, length(samp_sizes));
for cur_sim = 1:N_sims,
    samp_dist = zeros(N_resamples, length(samp_sizes));
    raw_data = randn(samp_sizes(end), 1);
    [~, t_p(cur_sim)] = ttest(raw_data, 0, 'Tail', 'right');
   
    raw_data=arrayfun(@(x) mean(raw_data(randsample(1:samp_sizes(end),100)),1),1:samp_sizes(end))';
    for i = 1:length(samp_sizes),
        cur_size = samp_sizes(i);
        for j = 1:N_resamples,
            samp_dist(j, i) = mean(raw_data(randsample(size(raw_data, 1), cur_size, true)));
        end
    end

    samp_p(cur_sim, :) = 1-mean(samp_dist > 0, 1);
    waitbar(cur_sim/N_sims, wh);
end
close(wh);

figure;
plot(t_p, samp_p(:, find(samp_sizes == 250)), 'rx');
hold all;
plot(t_p, samp_p(:, find(samp_sizes == 100)), 'bx');
plot(t_p, samp_p(:, find(samp_sizes == 50)), 'gx');
hold all;
plot([0 1], [0 1], 'k-');
legend('Resampling with 250', 'Resampling with 100', 'Resampling with 50');
ylabel('Bootstrap p-value');
xlabel('t-test p-value');