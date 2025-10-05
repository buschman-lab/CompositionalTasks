function p=PopOverlapPermutStatTest(nPop1,nPop2,nOverLap,nTot)

N = nTot;
S = nPop1;
A = nPop2;
SA = nOverLap;
npermute=10000;

sim_SA = [];
for sim = 1:npermute
    sim_S = (randperm(N) <= S);
    sim_A = (randperm(N) <= A);
    sim_SA(sim) = sum(sim_S .* sim_A);
end

p = sum([sim_SA SA] >= SA)/(npermute+1);

end