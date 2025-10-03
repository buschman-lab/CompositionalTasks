%%% Choosing Files Functions
%%% This is example cose of grabbing file s


%% Example 1
[Trials,SwitchTrials,TrialOpts]=GrabBehavioralFiles;
  close all  ;drawnow
 
%% example 2

if 1
    
    TrialOpts.BiasSoftmax     =[]; 
    TrialOpts.RewardPulse     =[]; 
    TrialOpts.ITITime         =[]; 
    TrialOpts.IncorrectTimeout=[]; 
    TrialOpts.Date=[];
    
        ip=1;BB=[0   ] %[0 1  0  0  1  2  3   0]%[0  1  2 0  1  0   1  2 6];%
        for HH=[12  ]%[15 15 11 19 19 19 19 20]%[18 21 21 22 22 23 23 23 23   ]%
            A=num2str(padarray(HH,[0 2-length(num2str(HH))],0,'pre'))
            HH=A(~isspace(A));
            
            Block=BB(ip);           
            B=num2str(padarray(Block,[0 2-length(num2str(Block))],0,'pre'))
            BC=B(~isspace(B));
            
            FileNames=[Anaopts.Kidname '_1909' num2str(HH) '_' num2str(BC) '_bhv.mat'];
            disp(['Loading ....' FileNames ])
            load([PATH FileNames],'bhv');
           
            ip=ip+1;
        end
end

%% example 3  Choose based on date
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

