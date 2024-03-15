function [Sp,h] = GetFig 
%GETMEAFIG gets me a simple figure fast
%   Detailed explanation goes here
FigParams=fig_params;
h=FigParams.RenderFigure(1,[]);
[h{1},Sp]=FigParams.RenderSubplots([1],[1],h{1},1);
end

