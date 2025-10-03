function RuleBlockTrials=SegregateBlocks(TrialTimes)
% finds trials for each rule and blocks
global AnalysisOpts
RuleInd=6;%strcmp(AnalysisOpts.TrialTimesFields,'CONDITION_NUM_OFFSET');%% we don't use this here but still good to have it
Ntrials=size(TrialTimes,1);
TrialRule=TrialTimes(:,RuleInd)';
temp=find(diff(TrialRule)~=0);
if temp(end)~=Ntrials
    SwitchTrials=[0 temp Ntrials];
else 
    SwitchTrials=[0 temp];
end
N=ones(1,3);
for i=1:length(SwitchTrials)-1            
   TrialInterval=SwitchTrials(i)+1:SwitchTrials(i+1);
   Rule=unique([TrialRule(TrialInterval)]);
   RuleBlockTrials.Rule{Rule}(N(Rule),:)=[SwitchTrials(i)+1 SwitchTrials(i+1) i];
   N(Rule)=N(Rule)+1;
   RuleBlockTrials.Seq(i,:)=[SwitchTrials(i)+1 SwitchTrials(i+1) Rule];
end
   
   