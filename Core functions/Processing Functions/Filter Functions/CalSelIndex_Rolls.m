function Selectivity=CalSelIndex_Rolls(FiringMat)
N=size(FiringMat,2);
for i=1:size(FiringMat,1)
    FR=FiringMat(i,:);
   Selectivity(i)= (1-((sum(FR,2)/N)^2)/sum((FR.^2)/N))/(1-1/N);
end

 
end