
PCATime=SubspaceAnaOpts.SpkCountPeriod(:,1);
clear Rot_Rule2 Comp_Rule2 Comp_Rule2_tot
for rep=1:SubspaceAnaOpts.Nrep
    Rot_Rule2=cell2mat(cellfun(@(x)  x.Rot,SubspaceAnaResults(nCond).TrialRange(nTrlRng).Rep(rep).AxBresults1_Analysis3,'UniformOutput',0))';
    Rot_Rule2(Rot_Rule2==180)=0;
    SumRot(:,rep)=sum(abs(Rot_Rule2),2)/2;
    Rot_Rule2_Tot(:,:,rep)=Rot_Rule2;
    
    Comp_Rule2(:,:,rep)=cell2mat(cellfun(@(x)  x.Comp,SubspaceAnaResults(nCond).TrialRange(nTrlRng).Rep(rep).AxBresults1_Analysis3,'UniformOutput',0))';
    Comp_Rule2_tot(:,:,rep)=cell2mat(cellfun(@(x)  x.CompTot,SubspaceAnaResults(nCond).TrialRange(nTrlRng).Rep(rep).AxBresults1_Analysis3,'UniformOutput',0))';    
    PCAAng(:,rep)=cell2mat(cellfun(@(x)  x.CosTheta,SubspaceAnaResults(nCond).TrialRange(nTrlRng).Rep(rep).AxBresults1_Analysis1,'UniformOutput',0))';
end

nAvg=10;
PCATimeAvg=smoothdata(PCATime',2,'movmean',nAvg);

Rot_Rule2(Rot_Rule2==180)=0;

figure
clf
subplot(131);axis square
hold on
%  plot(PCATimeAvg,smoothdata(mean((Rot_Rule2_Tot(:,1,:)),3),1,'movmean',nAvg));
%  plot(PCATimeAvg,smoothdata(mean((Rot_Rule2_Tot(:,2,:)),3),1,'movmean',nAvg));
%  plot(PCATimeAvg,smoothdata(mean((Rot_Rule2_Tot(:,3,:)),3),1,'movmean',nAvg));
%  plot(PCATimeAvg,smoothdata(mean((Rot_Rule2_Tot(:,4,:)),3),1,'movmean',nAvg));
 
%plot(PCATimeAvg,smoothdata(mean(SumRot,2),1,'movmean',nAvg));
 plot(PCATimeAvg,smoothdata(rad2deg(circ_mean(deg2rad(SumRot),[],2)),1,'movmean',nAvg));
%plot(PCATimeAvg,smoothdata(mean(PCAAng,2),1,'movmean',nAvg));
 

SumRot
%legend({'Rot1_Rule2','Rot2_Rule2','PairAngled_Rule2_3'})

subplot(132);axis square
hold on
plot(PCATimeAvg,smoothdata(mean(Comp_Rule2(:,1,:),3),1,'movmean',nAvg));
%legend({'PairAngled_Rule2','PairAngled_Rule2_1_all','PairAngled_Rule2_3_all'})

subplot(133);axis square
cla
hold on
plot(PCATimeAvg,smoothdata(mean(Comp_Rule2_tot(:,1,:),[],3),1,'movmean',nAvg));
plot(PCATimeAvg,smoothdata(mean(Comp_Rule2_tot(:,2,:),[],3),1,'movmean',nAvg));
plot(PCATimeAvg,smoothdata(mean(Comp_Rule2_tot(:,3,:),[],3),1,'movmean',nAvg));
plot(PCATimeAvg,smoothdata(mean(Comp_Rule2_tot(:,4,:),[],3),1,'movmean',nAvg));

% plot(PCATimeAvg,smoothdata(mean(Comp_Rule2_tot(:,2,:),3),1,'movmean',nAvg));
% plot(PCATimeAvg,smoothdata(mean(Comp_Rule2_tot(:,3,:),3),1,'movmean',nAvg));
% plot(PCATimeAvg,smoothdata(mean(Comp_Rule2_tot(:,4,:),3),1,'movmean',nAvg));

%legend({'PairAngled_Rule2','PairAngled_Rule2_1_all','PairAngled_Rule2_3_all'})

