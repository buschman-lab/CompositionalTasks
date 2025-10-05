function Perf=ProcessBhvData(bhv) % processes behavioral data for each day only 
%@bhv recorded behavior file

clearvars -global -except AnalysisOpts 
%% define varibales we need first
global AnalysisOpts LastTrial LastTrialMat Ncol LastTrialSample N PerfFit PerfFitTot TrialInfoBlk N_Ind LastTrialMat_Ind LastTrial_Ind 

%% load relevant classes
FigParams = fig_params; %load figure class
bhvAnaFunc=BhvAnalysisFuncs;

%% define some options of the function
AnalysisOpts.BhvAna.Kidname=bhv.Monkey;
AnalysisOpts.BhvAna.Nback=15;
AnalysisOpts.BhvAna.NbackSamp=5;
AnalysisOpts.BhvAna.FileMonth='';
AnalysisOpts.BhvAna.Perc_nonResTH=0.2;

%% perepare important functions and informaton
if ispc
    save_dir = ['Z:\Projects\Rule_Representation\Behavior AnalysisCode\Analysis_Code\Results\' AnalysisOpts.BhvAna.Kidname ' Performance\'];
 else
    save_dir =['/Volumes/buschman/Projects/Rule_Representation/Behavior AnalysisCode/Analysis_Code/Results/' AnalysisOpts.BhvAna.Kidname ' Performance/'];
 end
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%% Do some initial anaysis on the files
[Trials,SwitchTrials,Ltrials,Ltrialsthisblk,TrialOpts]=ExtractFileInformation(bhv);

%% Initialize variables 
Ncol=zeros(1,3);
AllConds=[Trials(SwitchTrials(1:end-1)).Condition];
NumCond=arrayfun(@(x) sum(AllConds==x),1:3);
PerfFitTot.L=[];PerfFitTot.k=[];PerfFitTot.x0=[];PerfFitTot.Ntrl=[];ConditionTot=[];TrialInfoBlk=[];
LastTrialSample{1}=zeros(1,4);LastTrialSample{2}=zeros(1,4);LastTrialSample{3}=zeros(1,4); 
LastTrial=zeros(1,3);LastTrialMat{1}=[1];LastTrialMat{2}=[1];LastTrialMat{3}=[1];N=ones(1,3);ChunkTrials=[];
N_Ind=ones(1,3); LastTrialMat_Ind{1}=[1];LastTrialMat_Ind{2}=[1];LastTrialMat_Ind{3}=[1]; 
LastTrial_Ind=zeros(1,3);
NSwitch=0;CorrectMatTot{1}=[];CorrectMatTot{2}=[];CorrectMatTot{3}=[];
CC=ones(1,3); % count the PSM blocks
CCA=ones(1,3); % count all of the blocks

%% initialize plots
fh=[]; %% total performance plot
fh_ind=[]; % ind perf plot
fs=cell(1,3); %% psychometics plot
pf=figure(bhv.FigNum);%figure(1);  %% main figure
NumTrials2LookFlora=[100 20 100]; % in case we are using trial2look that Flora has used

%% Start analyzing the data 
for i=1:length(SwitchTrials)-1
    TrialInterval=SwitchTrials(i):SwitchTrials(i+1)-1;
    
    if ~isempty(TrialInterval)
        PlotInterval=SwitchTrials(i):SwitchTrials(i+1)- AnalysisOpts.BhvAna.Nback;
        Condition(i)=unique([Trials(TrialInterval).Condition]);
        ThisCond=Condition(i);       
        BlockCorrIncorrMat=double(([Trials(TrialInterval).StopCondition] == 1 | [Trials(TrialInterval).StopCondition] ==-1) & [Trials(TrialInterval).Condition]==Condition(i));
        % grab number of trials to look
        NumTrials2Look=[Trials(TrialInterval).NumTrail2Look];
        NumTrials2Look=NumTrials2Look(end);
        %NumTrials2Look=NumTrials2LookFlora(ThisCond); % this is to match flora's code
        
        if sum(BlockCorrIncorrMat)>=NumTrials2Look 
            %% charctrize winstay lose switch performance
            I_WinStayLoseSwitch=CalWinStayLoseSwitchVer1(Condition(i),NSwitch,Trials,AnalysisOpts.BhvAna.Nback,TrialInterval,fh);
           
            %% Caculate the overall performance in each block and get all trials for the block
            [fh,Perf.OverAllPerf{ThisCond}{CCA(ThisCond)},Perf.AllBlkTrls{ThisCond}{CCA(ThisCond)}]=...
                CalculateOverallPerformance(Condition(i),NSwitch,Trials,AnalysisOpts.BhvAna.Nback,TrialInterval,I_WinStayLoseSwitch,fh);
           
            %% calculate performance for individual morph levels
            fh_ind=CalculateOverallPerformanceIndMorphVer1(Condition(i),NSwitch,Trials,NumTrials2Look,TrialInterval,fh_ind);            
            ConditionTot=[ConditionTot Condition(i)]; % get the sequence of conditions that are included 
           
            %% get last NumTrials2Look trials so we can look at 
            Head=TrialInterval(end);       %% calclulate where we have total correct incorrect that is equal to NumTrial2Look
            ConvCorrIncorr=conv(BlockCorrIncorrMat,ones(1,length(BlockCorrIncorrMat)));
            IndLook=find(ConvCorrIncorr==NumTrials2Look);
            IndLook=length(ConvCorrIncorr)-IndLook(end);
            TrialRange=Head-IndLook:Head;
            ChunkTrials=[Trials(TrialRange)]; % get this last chunk of trials
            
            %% run individual sample and psychometric curve analysis on the data 
            if ~(i==(length(SwitchTrials)-1) & AnalysisOpts.BhvAna.IgnoreLastBlk) % if we dont' want to include the last block
                [Perf.PsycPerf{ThisCond}(CC(ThisCond),:),Perf.PsycPerfLvl,Perf.PsycCount{ThisCond}(CC(ThisCond),:),fs]=PlotPsychoMetricCurveIndBlockSummery(bhv,ChunkTrials,ThisCond,0,[NumCond(ThisCond) CC(ThisCond)],fs);
                Perf.IndSampInfo{ThisCond}{CC(ThisCond)}=bhvAnaFunc.CalIndividualSampleInfo(bhv,ChunkTrials,ThisCond);
                % get info for all trials of this block
                Perf.IndSampInfoAllTrls{ThisCond}{CC(ThisCond)}=bhvAnaFunc.CalIndividualSampleInfo(bhv,Trials(TrialInterval),ThisCond);
                CC(ThisCond)=CC(ThisCond)+1;
            end
            
            CCA(ThisCond)=CCA(ThisCond)+1;
        end
    end
    NSwitch=NSwitch+1;
end

%% report extra varibles
AttemptedTrls=[Trials.StopCondition]==1 | [Trials.StopCondition]==-1;
Perf.TotAttemptedTrls=sum(AttemptedTrls);
Perf.TotCorrectTrls=sum([Trials.StopCondition]==1);
Perf.RewardPulse=bhv.opts.RewardPulse;
Perf.NumRewards.Switch=bhv.NRewardChange;
Perf.ConditionTot=ConditionTot;
% [BlockNums,BlockIden]=SegregateBhvBlocks(ConditionTot);
% PlotBlockPerf(BlockNums,ConditionTot,BlockIden,TrialOpts)
% AnalyzeTrialDist(bhv)

%% save figures and variables
if ~isempty(fh)
    %% Now save all of the figures and append them into a PDF
    [~,SaveFileName]= fileparts(bhv.SaveFilename);
    FigParams.SaveFigSeries(SaveFileName,save_dir,[{fh_ind} fs {pf}])    
    close(fh)
else
    Perf.OverAllPerf=[];
    Perf.IndSampInfo=[];
    Perf.PsycPerf=[];
    Perf.PsycPerfLvl=[];
    Perf.PsycCount=[];
    Perf.TotAttemptedTrls
end








