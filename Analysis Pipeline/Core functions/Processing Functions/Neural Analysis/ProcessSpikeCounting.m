function [TrialSpkAnalysis]=ProcessSpikeCounting(SpikeTimesData,CriticalTimes)
%%%this also builds a matrix through time
global AnalysisOpts
TrialSpkAnalysis.CriticalTimes=CriticalTimes;
ShowPlot=0;
SpkParams=AnalysisOpts.SpkParams;
SpkParams.numbins=(SpkParams.PeriodLength-SpkParams.Bin)/SpkParams.BinShift;
%%%
SpkParams.NumTrial=size(SpikeTimesData,2);
 
AllSpikeTimesData=zeros(SpkParams.NumTrial,floor(AnalysisOpts.SpkParams.PeriodLength*1000));

%%% find when the trial is starting and Ending
StartTrialInd=strcmp(AnalysisOpts.TrialTimesFields,'START_TRIAL');
StopTrialInd=strcmp(AnalysisOpts.TrialTimesFields,'END_TRIAL');

%%% find when we want to count the number of spikes related to the task 
StartCountInd=strcmp(AnalysisOpts.TrialTimesFields,AnalysisOpts.SpkParams.StartFieldName);
StopCountInd=strcmp(AnalysisOpts.TrialTimesFields,AnalysisOpts.SpkParams.StopFiledName);

RewardInd=strcmp(AnalysisOpts.TrialTimesFields,'REWARD');

for Tri=1:SpkParams.NumTrial   
      
      SpkParams.SpkCountStr  =CriticalTimes(Tri,StartCountInd)-SpkParams.RuleSelectvityPeriod; %%% count 
      SpkParams.SpkCountStp  =CriticalTimes(Tri,StartCountInd);     %%%CriticalTimes(Tri,StopCountInd);
      
      SpkParams.SpkCountStr2  =CriticalTimes(Tri,StartCountInd); %%% count 
      SpkParams.SpkCountStp2  =CriticalTimes(Tri,StartCountInd)+SpkParams.RuleSelectvityPeriod;     %%%CriticalTimes(Tri,StopCountInd);
      
      SpkParams.RewardTime   =CriticalTimes(Tri,RewardInd);
      SpkParams.TimePCAStart =CriticalTimes(Tri,StartCountInd)-AnalysisOpts.SpkParams.BaselineDelay;  %% start at the time of fixation
      SpkParams.StartTrial   =SpkParams.TimePCAStart;
      SpkParams.EndTrial     =SpkParams.StartTrial+AnalysisOpts.SpkParams.PeriodLength;
      SpkParams.TrialDuration=SpkParams.EndTrial-SpkParams.StartTrial;
      SpkParams.SpkCountDur  =SpkParams.SpkCountStp-SpkParams.SpkCountStr;

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      SpkCount=find(SpikeTimesData{Tri}>=SpkParams.SpkCountStr & SpikeTimesData{Tri}<=SpkParams.SpkCountStp); 
      TrialSpkAnalysis.SpikeCountData.NumSpk_befStim(Tri)=length(SpkCount);
      
      SpkCount=find(SpikeTimesData{Tri}>=SpkParams.SpkCountStr2 & SpikeTimesData{Tri}<=SpkParams.SpkCountStp2); 
      TrialSpkAnalysis.SpikeCountData.NumSpk_aftStim(Tri)=length(SpkCount);
      
      %%% find baselinefiring at the end of 1sec 
      %SpkCountBaseline=find(SpikeTimesData{Tri}{Neu}>=(1-SpkParams.SpkCountDur) & SpikeTimesData{Tri}{Neu}<=1); 
      %SpikeData{Tri}.NumSpkBaslin(Neu)=length(SpkCountBaseline);
      
          for Nbin=1:SpkParams.numbins
              %   SpkParams.TimePCAStart=SpkParams.SpkCountStr;
                  Stp=SpkParams.TimePCAStart+SpkParams.Bin+SpkParams.BinShift*(Nbin-1);
                  Str=Stp-SpkParams.Bin;
                  TrialSpkAnalysis.PSTH.SpkCountBin(Tri,Nbin)=length(find(SpikeTimesData{Tri}>=Str & SpikeTimesData{Tri}<=Stp));
                  TrialSpkAnalysis.PSTH.SpkCountTim(Nbin)=Str;
          end   
    
      %%% Plot raster  
      ValidSpikeTimes=SpikeTimesData{Tri}(SpikeTimesData{Tri}>=SpkParams.StartTrial & SpikeTimesData{Tri}<=SpkParams.EndTrial)-SpkParams.StartTrial;
      SpkTimIndCh    =floor(ValidSpikeTimes*1000);
      SpkTimIndCh(SpkTimIndCh==0)=1;
      TrialSpkAnalysis.AllSpikeTimesData(Tri,SpkTimIndCh)=1;
end
TrialSpkAnalysis.SpkParams=SpkParams;
%%% now plot the results for each neuron
if ShowPlot
     for Neu=1:NumNeu 
             %%%raster plot
           %  subplot(4,ceil(NumNeu/4),Neu);
          %  figure(90000);cla
           % rasterPlot(AllSpikeTimesData(:,:,Neu),0.001,SpkParams,Col);drawnow  
            %%% psth plot
            figure(80000);cla
          %  set(gca,'units','normalized','position',[0 0 1 1])     
         %   subplot(4,ceil(NumNeu/4),Neu);hold on
            thisPSTH=squeeze(mean(PSTH.SpkCountBin(:,Neu,:),1));   
            plot(PSTH.SpkCountTim,thisPSTH/SpkParams.Bin,'Color',Col(Neu,:),'LineWidth',2);hold on
            plot([SpkParams.SpkCountStr SpkParams.SpkCountStr],[0 max(thisPSTH(:))],'r','linewidth',5)
            plot([SpkParams.SpkCountStp SpkParams.SpkCountStp],[0 max(thisPSTH(:))],'b','linewidth',5)
            plot([SpkParams.RewardTime  SpkParams.RewardTime ],[0 max(thisPSTH(:))],'k','linewidth',5)
            title([num2str(Neurons(Neu))])
            set(gcf,'Name','PSTH for Visual or Electrical Response')
            pause
     end
end
    

 