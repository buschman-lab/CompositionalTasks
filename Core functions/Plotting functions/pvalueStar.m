function [symb,col]=pvalueStar(pval,varargin)
if ~isempty(varargin);Alpha=varargin{1};else;Alpha=0.05;end

if pval>Alpha
    symb='';
    col=1;
elseif pval<Alpha & pval>=Alpha*0.2
    symb='*';
    col=2;
elseif pval<Alpha*0.2 & pval>=Alpha*0.02
    symb='**';
    col=2;
elseif pval<Alpha*0.02 
    symb='***';
    col=2;
elseif isnan(pval)
    symb='';
    col=1;
end