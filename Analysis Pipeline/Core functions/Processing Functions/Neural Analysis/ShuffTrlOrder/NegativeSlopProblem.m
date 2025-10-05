ManData=ManipulateData;
AnalysisOpts.Time=-0.6:0.01:0.55;
%PW_method='VCTFPW';
PW_method='3PW';
%ObservedOrg=Observed;
File='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Analysis Pipeline\Input Output Data\Data for Tim\Learning3D_Sh_Co_Ru_AltRule_RB_F_SkTrSq_MS_FEA_FL-75_40_4_50_115AT_ColorCat_D2ShData.mat';
load(File,'Accuracy_Obsv');
ObservedOrg=cell2mat(Accuracy_Obsv);

Observed=ManData.SmoothData(ObservedOrg,15,'SmoothingMethod','movmean','DimSmoothing',2);   
Observed=ManData.SmoothData(Observed,9,'SmoothingMethod','movmean','DimSmoothing',1);  
NTim=size(Observed,2);

NStg=16;%TrlRng(end);
test_data=timetable(datetime(datevec([datenum([2001 01 01]):365:datenum([2000+NStg 01 01])])));
% calculate observed
clear Observed_p Observed_a SenSlope ArimaObserved_a ArimaObserved_p ArimaObserved_slop ArimaResiduals

for t=1:NTim
    test_data.param=(Observed(1:NStg,t));
    try
        results=MK_tempAggr(test_data, 0.001, 'PW_method',PW_method,'alpha_MK',99,'alpha_ak',95, 'alpha_CL',90, 'alpha_Xhomo',99);
        Observed_p(1,t)=results.P;
        Observed_a(1,t)=results.slope;
    catch
        Observed_p(1,t)=nan;
        Observed_a(1,t)=nan;
    end
    SenSlope(1,t)=sens_slope(1:NStg,test_data.param);
end

RepSpace=@(x) strrep(x,'_',' ');
col=copper(NStg);

NPos=find(AnalysisOpts.Time>=0.23,1,'first');
NNeg=find(AnalysisOpts.Time>=0.39,1,'first');
TPos=AnalysisOpts.Time(NPos);
TNeg=AnalysisOpts.Time(NNeg);
ColPos='g';ColNeg='r';

Apos=Observed(1:NStg,NPos);
ANeg=Observed(1:NStg,NNeg);
    
test_data_pos=timetable(datetime(datevec([datenum([2001 01 01]):365:datenum([2000+NStg 01 01])])));
test_data_neg=timetable(datetime(datevec([datenum([2001 01 01]):365:datenum([2000+NStg 01 01])])));

test_data_pos.param=Apos;
test_data_neg.param=ANeg;

results_pos=MK_tempAggr(test_data_pos, 0.01, 'PW_method',PW_method,'alpha_MK',99,'alpha_ak',95, 'alpha_CL',90, 'alpha_Xhomo',99);
results_neg=MK_tempAggr(test_data_neg, 0.01, 'PW_method',PW_method,'alpha_MK',99,'alpha_ak',95, 'alpha_CL',90, 'alpha_Xhomo',99);

 %
 %model = arima(2,0,5); % AR(Lag) model
 ARLag=1; MA=1;
 %model = arima(ARLag,0,MA);%arima('ARLags', ARLag,'MA',MA);
 model = arima('ARLags', ARLag, 'D', 0, 'Constant', 0);%,'MA',2); % AR(Lag) model
 ARIMAtitle=sprintf('ARLag=%i MA=%i',ARLag,MA);
 for t=1:NTim
     try
         [ArimaObserved_a(t),ArimaObserved_p(t),ArimaObserved_slop(t),ArimaResiduals{t}]=...
             ARIMAtrendcorr((Observed(1:NStg,t)),[1:NStg]',model,'Spearman');
     catch me
         fprintf('\n t=%i %s',t, me.message)
         ArimaObserved_a(t)=nan;
         ArimaObserved_p(t)=nan;
         ArimaObserved_slop(t)=nan;
         ArimaResiduals{t}=nan;
     end
 end

%  ArimaObserved_a=cell2mat(ArimaObserved_a);
%  ArimaObserved_p=cell2mat(ArimaObserved_p);
%  ArimaObserved_slop=cell2mat(ArimaObserved_slop);
%  


 figure
 subplot(231)
 plot(AnalysisOpts.Time,Observed_a)
 hold on
 plot(AnalysisOpts.Time,SenSlope)
 legend({PW_method,'Sens Slope'})
 title(PW_method)

 subplot(232)
 plot(AnalysisOpts.Time,Observed_p)
 legend({[PW_method ' pval']})
 title(PW_method)

 subplot(233)
 hold on;arrayfun(@(x) plot(AnalysisOpts.Time,Observed(x,1:NTim),'color',col(x,:)),[1:NStg]);
 xlabel('Time');title('Observed')
 ylabel('Observed')
 yyaxis right
 plot(AnalysisOpts.Time,Observed_a,'b','LineWidth',5)
 ylabel('Slope 3PW')
 plot(TPos,Observed_a(NPos),ColPos,'Marker','*','MarkerSize',20);
 plot(TNeg,Observed_a(NNeg),ColNeg,'Marker','*','MarkerSize',20);


 subplot(234)
 plot(AnalysisOpts.Time,ArimaObserved_a)
 ylabel('Spearman Rho')
 hold on
 
 yyaxis right
 plot(AnalysisOpts.Time,ArimaObserved_slop)
 title(ARIMAtitle)
 xlabel('Time')
 ylabel('Sens slope')
  legend({'ARIMA Spearman rho','Sens Slope'})
  yyaxis left
plot(TPos,ArimaObserved_a(NPos),ColPos,'Marker','*','MarkerSize',20);
 plot(TNeg,ArimaObserved_a(NNeg),ColNeg,'Marker','*','MarkerSize',20);



 subplot(235)
 plot(AnalysisOpts.Time,ArimaObserved_p)
 yyaxis right
 plot(AnalysisOpts.Time,ArimaObserved_slop)
 legend({'ARIMA p','Sens Slope'})
 title(ARIMAtitle)
 xlabel('Time')

 subplot(236)
 plot(1:length(ArimaResiduals{NPos}),ArimaResiduals{NPos},ColPos)
 hold on
 plot(1:length(ArimaResiduals{NPos}),ArimaResiduals{NNeg},ColNeg)
 legend({'residuals pos slop','residuals neg slop'})
 title(ARIMAtitle)