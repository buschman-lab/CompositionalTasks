function Perf=PreprocessBhvData(Animal) % preprocesses behavioral data

global AnalysisOpts
bhvAna=BhvAnalysisFuncs;
if strcmp(Animal,'ALL')
    AnimalSet={'Silas','Chico'};
    AnimalNameSave='ALL';
else
    AnimalSet={Animal};
    AnimalNameSave=Animal;
end
PerfTot=[];
for Animal=AnimalSet
    Animal=Animal{1};
    %% define analysis paramters here
    opts.InitialDate= AnalysisOpts.BhvAna.([Animal '_InitialDate']) ;%only dates bigger than this are chosen
    opts.EndDate    = AnalysisOpts.BhvAna.([Animal '_EndDate']) ; %only dates smaller than this are chosen
    if sum(opts.InitialDate==opts.EndDate)==3;opts.PltSingleDay=1;else;opts.PltSingleDay=0;end % are we plotting single days
    opts.Kidname=Animal;
    opts.PATH=['Z:\Projects\Rule_Representation\Data\' opts.Kidname '_Recording\Useable Behavioral Data\'];
    opts.CodePath=pwd;
    
    if ispc
        save_dir=['Z:\Projects\Rule_Representation\Behavior AnalysisCode\Analysis_Code\Results\' opts.Kidname ' Performance\'];
    else
        save_dir=['/Volumes/buschman/Projects/Rule_Representation/Behavior AnalysisCode/Analysis_Code/Results/' opts.Kidname ' Performance/'];
    end
    
    %% load behavior files
    cd(opts.PATH)
    warning off
    if ~opts.PltSingleDay
        BHVs=ChooseBHVfilebyDate(opts.InitialDate,opts.EndDate,opts);
    else
        BHVs=ChooseBHVfileSpecificDate(opts.Kidname,opts,opts.InitialDate(1),opts.InitialDate(2),opts.InitialDate(3),0) ;
    end
    
    %% preprocess each day's data
    Perf=cell(1,length(BHVs));
    for i=1:length(BHVs)
        fprintf('\nprocesessing date %s',BHVs{i}.bhv.Date)
        Perf{i}=ProcessBhvData(BHVs{i}.bhv);
        close all;
    end
    PerfTot=[PerfTot Perf];
end

[AllPSMPerf,AllTrlPerf,IndSamp,AllTrlCount,AllTrlCountDay,NBlocksDay,RewardPulse,NCorrectTrl,NumRewards,AllSeqHist] =...
    bhvAna.CancatinateInfoDays(PerfTot);
Perf=PerfTot;
%% save data from this analysis
save([AnalysisOpts.BhvMdlPath AnimalNameSave '_BhvData_' char(datetime(opts.InitialDate)) '_till_' char(datetime(opts.EndDate)) '.mat'],...
    'Perf','opts','AnalysisOpts','AllPSMPerf','AllTrlPerf','IndSamp','AllTrlCount','AllTrlCountDay','NBlocksDay',...
    'RewardPulse','NCorrectTrl','NumRewards','AllSeqHist','-v7.3')

end