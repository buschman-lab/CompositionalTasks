function Selectivity=CalSelIndex(FiringMat)
Nrules=size(FiringMat,2);
SumFiring=sum(FiringMat,2);
MaxFiring=max(FiringMat,[],2);
Selectivity=(Nrules-(SumFiring./MaxFiring))./(Nrules-1);
%Selectivity=Selectivity(~isnan(Selectivity));
end