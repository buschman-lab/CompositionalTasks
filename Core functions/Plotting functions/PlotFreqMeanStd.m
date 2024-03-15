function PlotFreqMeanStd(amp,amp_f,Xlabel,Ylabel,ColInd)
%PLOTMEANSTD Summary of this function goes here
%   Detailed explanation goes here
if nargin<5
    ColInd=3;
end
co = get(gca, 'ColorOrder');
plot(amp_f, nanmean(10*log10(amp), 2), '-', 'Color', co(mod(ColInd-1, size(co, 1))+1, :)); hold on;
plot(amp_f, nanmean(10*log10(amp), 2) + nanstd(10*log10(amp), [], 2), ':', 'Color', co(mod(ColInd-1, size(co, 1))+1, :));
plot(amp_f, nanmean(10*log10(amp), 2) - nanstd(10*log10(amp), [], 2), ':', 'Color', co(mod(ColInd-1, size(co, 1))+1, :));  
xlabel(Xlabel);ylabel(Ylabel)
grid on;
axis tight


end

