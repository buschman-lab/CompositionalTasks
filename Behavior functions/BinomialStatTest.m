function pval=BinomialStatTest(Ntot,Ncorrect,ChancelLevel)
             pval=arrayfun(@(x) 1-binocdf(Ncorrect(x),Ntot(x),ChancelLevel),1:length(Ntot));
end