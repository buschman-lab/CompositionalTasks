
PCATime=mean(SubspaceAnaOpts.SpkCountPeriod,2);
PairAngled_Rule2=cell2mat(arrayfun(@(x) SubspaceAnaResults(nCond).TrialRange(nTrlRng).Rep(x).PCAresults1_Analysis1.PairAngled',1:100,'UniformOutput',0))';
PairAngled_Rule2_1=cell2mat(arrayfun(@(x) SubspaceAnaResults(nCond).TrialRange(nTrlRng).Rep(x).PCAresults1_Analysis2.PairAngled',1:100,'UniformOutput',0))';
PairAngled_Rule2_1_all=cell2mat(arrayfun(@(x) SubspaceAnaResults(nCond).TrialRange(nTrlRng).Rep(x).PCAresults1_Analysis3.PairAngled',1:100,'UniformOutput',0))';

PairAngled_Rule2_3=cell2mat(arrayfun(@(x) SubspaceAnaResults(nCond).TrialRange(nTrlRng).Rep(x).PCAresults2_Analysis2.PairAngled',1:100,'UniformOutput',0))';
PairAngled_Rule2_3_all=cell2mat(arrayfun(@(x) SubspaceAnaResults(nCond).TrialRange(nTrlRng).Rep(x).PCAresults2_Analysis3.PairAngled',1:100,'UniformOutput',0))';

SubspaceComp2=cell2mat(arrayfun(@(x) SubspaceAnaResults(nCond).TrialRange(nTrlRng).Rep(x).PCAresults1_Analysis1.SubspaceCompression2',1:100,'UniformOutput',0))';



nAvg=5;
PCATimeAvg=smoothdata(PCATime',2,'movmean',nAvg);


figure
subplot(121);axis square
hold on
plot(PCATimeAvg,smoothdata(mean(PairAngled_Rule2(1:100,:),1),2,'movmean',nAvg));
plot(PCATimeAvg,smoothdata(mean(PairAngled_Rule2_1,1),2,'movmean',nAvg));
plot(PCATimeAvg,smoothdata(mean(PairAngled_Rule2_3,1),2,'movmean',nAvg));
legend({'PairAngled_Rule2','PairAngled_Rule2_1','PairAngled_Rule2_3'})

subplot(122);axis square
hold on
plot(PCATimeAvg,smoothdata(mean(PairAngled_Rule2,1),2,'movmean',nAvg));
plot(PCATimeAvg,smoothdata(mean(PairAngled_Rule2_1_all,1),2,'movmean',nAvg));
plot(PCATimeAvg,smoothdata(mean(PairAngled_Rule2_3_all,1),2,'movmean',nAvg));
legend({'PairAngled_Rule2','PairAngled_Rule2_1_all','PairAngled_Rule2_3_all'})



yyaxis right
plot(PCATimeAvg,smoothdata(mean(SubspaceComp2,1),2,'movmean',nAvg));