function FindRewardPatterns

DIR=dir;

StartDate = datetime(2019,01,01);N=1;
 for HH=1:length(DIR)
     [~,~,Ext]=fileparts(DIR(HH).name);
            if DIR(HH).date>StartDate && strcmp(Ext,'.mat') && ~isempty(findstr(DIR(HH).name,'Chico'))
                 disp(['Loading ....' char(DIR(HH).name)])
                BHV{N}=load(DIR(HH).name,'bhv');  
                N=N+1;
            end
 end
 close all
Dates= cellfun(@(x) x.bhv.Date,BHV,'UniformOutput',0);
Stats.ITI=cellfun(@(x) x.bhv.opts.ITITime,BHV);
Stats.NumRewards=cellfun(@(x) x.bhv.opts.NumRewards,BHV);
Stats.RewardPulse=cellfun(@(x) x.bhv.opts.RewardPulse,BHV);
Stats.IncorrectTimeout=cellfun(@(x) x.bhv.opts.IncorrectTimeout,BHV);
Stats.TotalCorrTrls=cellfun(@(x) sum([x.bhv.Trials.StopCondition] == 1),BHV);
Stats.TotalTrls=cellfun(@(x) sum([x.bhv.Trials.StopCondition] == 1 | [x.bhv.Trials.StopCondition] == -1 ),BHV);

save stats Stats Dates StartDate
%%% now find patterns of best reward 


%X=[Stats.ITI' Stats.NumRewards' Stats.RewardPulse' [Stats.RewardPulse.*Stats.NumRewards]' Stats.IncorrectTimeout'];
X=[Stats.ITI' Stats.NumRewards' Stats.RewardPulse'  Stats.IncorrectTimeout'];

