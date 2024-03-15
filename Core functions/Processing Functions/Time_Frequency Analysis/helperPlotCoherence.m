function helperPlotCoherence(wcoh,t,F,coi,xlab,ylab)
%   This helper function is provided solely in support of the featured
%   example wcoherdemo.m
%   This function may be changed or removed in a future release.

Yticks = 2.^(round(log2(min(F))):round(log2(max(F))));
imagesc(t,log2(F),wcoh);
set(gca,'YLim',log2([min(F),max(F)]), ...
    'layer','top', ...
    'YTick',log2(Yticks(:)), ...
    'YTickLabel',num2str(sprintf('%.2f\n',Yticks)), ...
    'layer','top','YDir','normal');
hold on;
plot(t,log2(coi),'w--');
xlabel(xlab);
ylabel(ylab);
hcol = colorbar;
hcol.Label.String = 'Magnitude-Squared Coherence';
title('Wavelet Coherence');
