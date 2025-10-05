function [Trials,SwitchTrials,TrialOpts]=GrabBehavioralFiles

global Anaopts

Trials=[];Ltrials=0;SwitchTrials=[];

if 1
[file_list, folder_list] = GrabFiles(['bhv' Anaopts.FileMonth], 1);
 

data_all = cell(1,numel(file_list)); 
%%%%%
TrialOpts.BiasSoftmax     =[]; 
TrialOpts.RewardPulse     =[]; 
TrialOpts.ITITime         =[]; 
TrialOpts.IncorrectTimeout=[]; 
TrialOpts.Date=[]; 

%%%%%

for cur_block = 1:numel(file_list)
  
   disp(['Loading ....' file_list{cur_block} ])
    
   data_all{cur_block} = load(file_list{cur_block},'bhv');
   close all
   bhv=data_all{cur_block}.bhv;
   if isfield(bhv,'Trials')
         Trials=[Trials bhv.Trials];
         SwitchTrials=[SwitchTrials bhv.SwitchTrialsInd+Ltrials];
         Ltrials=length(Trials);
         Ltrialsthisblk=length(bhv.Trials);
         %%%%
         TrialOpts.BiasSoftmax=[TrialOpts.BiasSoftmax bhv.BiasSoftmax*ones(1,Ltrialsthisblk)]; 
         TrialOpts.RewardPulse=[TrialOpts.RewardPulse bhv.opts.RewardPulse*ones(1,Ltrialsthisblk)]; 
         TrialOpts.ITITime    =[TrialOpts.ITITime bhv.opts.ITITime*ones(1,Ltrialsthisblk)]; 
         TrialOpts.IncorrectTimeout=[TrialOpts.IncorrectTimeout bhv.opts.IncorrectTimeout*ones(1,Ltrialsthisblk)]; 
         TrialOpts.Date=[TrialOpts.Date str2double(bhv.Date)*ones(1,Ltrialsthisblk)];
   
   end
      
end

SwitchTrials=[SwitchTrials Ltrials];
end
save('playdata1months','Trials','SwitchTrials','TrialOpts','-v7.3')
%%%%%%%%%%%% Another option to open files 
if 0

PATH=['/Volumes/buschman/Projects/Rule Representation/Data/' Anaopts.Kidname '/'];
cd(PATH)


         ip=1;BB=[zeros(1,11)]; %[0 1  0  0  1  2  3   0]%[0  1  2 0  1  0   1  2 6];%
        for HH= [15 16 17 20 21 22 23 24 28 29 30]% 18 19 21 22]%[15 15 11 19 19 19 19 20]%[18 21 21 22 22 23 23 23 23   ]%
            A=num2str(padarray(HH,[0 2-length(num2str(HH))],0,'pre'));
            HH=A(~isspace(A));
            
            Block=BB(ip);           
            B=num2str(padarray(Block,[0 2-length(num2str(Block))],0,'pre'));
            BC=B(~isspace(B));
            
            FileNames=[Anaopts.Kidname '_1905' num2str(HH) '_' num2str(BC) '_bhv.mat'];
            disp(['Loading ....' FileNames ])
            load(FileNames,'bhv');
            disp(['Incorrect timeout=' num2str(bhv.opts.IncorrectTimeout)])
            if isfield(bhv,'Trials')
                Trials=[Trials bhv.Trials];
                SwitchTrials=[SwitchTrials bhv.SwitchTrialsInd+Ltrials];
                Ltrials=length(Trials);
            end
            ip=ip+1;
        end
end
