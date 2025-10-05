%%%  grabs basic information about a behavioral file 
function [Trials,SwitchTrials,Ltrials,Ltrialsthisblk,TrialOpts]= ExtractFileInformation(bhv)

% initialze variables
SwitchTrials=[];Trials=[];Ltrials=0;
TrialOpts.BiasSoftmax     =[]; 
TrialOpts.RewardPulse     =[]; 
TrialOpts.ITITime         =[]; 
TrialOpts.IncorrectTimeout=[]; 
TrialOpts.Date=[];

if isfield(bhv,'Trials')
    Trials=[Trials, bhv.Trials];
    SwitchTrials=[SwitchTrials, bhv.SwitchTrialsInd+Ltrials];
    Ltrials=length(Trials);
    Ltrialsthisblk=length(bhv.Trials);
    SwitchTrials=[SwitchTrials Ltrials];
     %%%% what are options of training
     if isfield(bhv,'BiasSoftmax')
         TrialOpts.BiasSoftmax=[TrialOpts.BiasSoftmax bhv.BiasSoftmax*ones(1,Ltrialsthisblk)]; 
     else 
         TrialOpts.BiasSoftmax=[TrialOpts.BiasSoftmax NaN];
     end
    TrialOpts.RewardPulse=[TrialOpts.RewardPulse bhv.opts.RewardPulse*ones(1,Ltrialsthisblk)]; 
    TrialOpts.ITITime    =[TrialOpts.ITITime bhv.opts.ITITime*ones(1,Ltrialsthisblk)]; 
    TrialOpts.IncorrectTimeout=[TrialOpts.IncorrectTimeout bhv.opts.IncorrectTimeout*ones(1,Ltrialsthisblk)]; 
    TrialOpts.Date=[TrialOpts.Date str2double(bhv.Date)*ones(1,Ltrialsthisblk)];

end
