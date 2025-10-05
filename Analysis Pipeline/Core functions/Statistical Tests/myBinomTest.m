function [pout]=myBinomTest(s,n,p,Sided)

pout=arrayfun(@(x) 1-binocdf(s(x),n(x),p),1:length(n));

