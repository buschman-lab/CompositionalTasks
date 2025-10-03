function helperFrequencyAnalysisPlot2(x,y,xlbl,ylbl,ttl,lgnd,xlim)
% Plot helper function for the FrequencyAnalysisExample

% Copyright 2012 The MathWorks, Inc.

plot(x,y)
xlabel(xlbl)
ylabel(ylbl)
if nargin > 4 && ~isempty(ttl)
  title(ttl)
end
if nargin > 5 && ~isempty(lgnd)
  legend(lgnd{:},'Location','southwest')
end
grid on
if nargin > 6 && ~isempty(xlim) 
  ax = axis;
  axis([xlim(1) xlim(2) ax(3:4)])
else
  axis tight
end
