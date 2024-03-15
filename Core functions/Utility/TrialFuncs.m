classdef TrialFuncs < ManipulateData
    %TRIALFUNCS Summary of this class goes here
    %   Detailed explanation goes here

    properties
        AnalysisType
        Fs   % sampling frequency of data
        WaveletDownSampleFactor  % how much downsampel for wavelet
        Rule_2look=[1,2,3];
        Nblk_2Look
        MaxNTrial   % Maximum number of trials for each block
        NTrial2Swtch % Number of trials to switch that we want to look at
        MinTrlsRuleBlk=[100 40 100]; % minimum number of trials a block needs to have in order to be included in the analysis (taking flora's numbers)
        UseDataPointer=0; % if this is on the data pointer will be used instead of ral variable
        DataPointerVar=''; % what is the varibale that the pointer is using
    end
    properties (Access=private)
        ManData=ManipulateData;
        bhvAna=BhvAnalysisFuncs;
    end

    methods
        function obj = TrialFuncs(varargin)
            if nargin~=0 % initialize vars
                obj=obj.ParseParams(varargin) ; %%Process optional inputs
            end
        end
        function obj=ParseParams(obj,InputArgs)
            %Process optional inputs
            if mod(length(InputArgs), 2) ~= 0, error('Must pass key/value pairs for options.'); end
            for i = 1:2:length(InputArgs)
                try
                    obj.(InputArgs{i}) = InputArgs{i+1};
                catch
                    error('Couldn''t set option ''%s''.', InputArgs{2*i-1});
                end
            end
        end
        function Trial_Data=GrabDataTimeTrial(obj,data,TrialTimes,DesiredParameter,varargin)

            %%% grabs trial data from a straem of data given the start stop spec and
            %%% full data stream and timing of trials
            %%% AnalysisOpts contains desired paramters
            %%% input:
            %%% 1- data: raw data that needs to be chopped based on timing and trial
            %%% 2- TrialTimes: trial times matrix that contains timings of trials
            %%% 3- DesiredParameter: Parameters that we desire to be inspected such as
            %%% frequency band, leave empty if none
            %%% output: Trial/Time based chpped data in the desired paramter range

            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            if isempty(DesiredParameter);DesiredParameter=1;end
            %set up data pointer

            if strcmpi(AnalysisOpts.Project,'Learning attentional templates')
                StartCountInd=strcmp(AnalysisOpts.TrialTimesFields,AnalysisOpts.TrialTiming.StartFieldName); %% time to start chopping
                %%% define time axis
                if ~obj.UseDataPointer
                    Time=(0:1/obj.Fs:(length(data)-1)/obj.Fs);
                else
                    eval(['Npoints=size(data.' obj.DataPointerVar ',2);'])
                    Time=(0:1/obj.Fs:(Npoints-1)/obj.Fs);
                end
                Trial_Data=[];

                % get the trail data
                for Tri=1:size(TrialTimes,1)  %loop through trials
                    fprintf('\n*** Grabing data from trial: %s ***\n ', num2str(Tri));
                    StartTime  =TrialTimes(Tri,StartCountInd)-AnalysisOpts.TrialTiming.BaselineDelay;  %% start of chunk time
                    StopTime   =StartTime+AnalysisOpts.TrialTiming.PeriodLength;   %% stop of chunk time
                    Duration   =StopTime-StartTime;   % duration of chunk

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run your analysis on each
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% trial

                    if    strcmpi(obj.AnalysisType,'LFP_Trial')
                        for DesPar=1:size(DesiredParameter,1) %% now loop on the desired paramters you want such as frequncy band

                            if obj.UseDataPointer
                                eval(['RawData=data.' obj.DataPointerVar '(DesPar,Time>=StartTime & Time<StopTime);']);
                            else
                                RawData=data(DesPar,Time>=StartTime & Time<StopTime);
                            end
                            RawData=RawData';
                            %% apply downsampling analysis we want
                            [RawData] = ApplyNotchFilter(RawData,AnalysisOpts.NeuralAnalysis.SamplingFrequency); % notch filter the data
                            RawData = ApplyLFPLowPass(RawData,AnalysisOpts.NeuralAnalysis.SamplingFrequency);  % apply low pass filter to only look at lfp
                            Trial_Data(:,Tri,DesPar)=downsample(RawData,obj.Fs/AnalysisOpts.LFPParams.FilterOpts.TargetLFPSamplingFrequency);        % downsample the data to a reasonable range
                        end

                    elseif strcmpi(obj.AnalysisType,'Wavelet_Trial')
                        if obj.UseDataPointer
                            eval(['TFData=data.' obj.DataPointerVar '(:,Time>=StartTime & Time<StopTime);']);
                        else
                            TFData=data(:,Time>=StartTime & Time<StopTime);
                        end

                        TFData=data(:,Time>=StartTime & Time<StopTime);
                        Trial_Data(:,:,Tri)=transpose(downsample(TFData',obj.WaveletDownSampleFactor));        % downsample the data to a reasonable range

                    elseif strcmpi(obj.AnalysisType,'Wavelet_Trial_H')  % reorganize H matrix
                        if isnan(StartTime) || isnan(StopTime)
                            TFData=nan*ones(size(data,1),obj.Fs*AnalysisOpts.TrialTiming.PeriodLength);
                        else
                            if obj.UseDataPointer
                                eval(['TFData=data.' obj.DataPointerVar '(:,find(Time>=StartTime & Time<StopTime));']);
                            else
                                TFData=data(:,Time>=StartTime & Time<StopTime);
                            end
                        end
                        Trial_Data=[Trial_Data transpose(downsample(TFData',obj.WaveletDownSampleFactor))];        % downsample the data to a reasonable range
                    end
                end

            elseif strcmpi(AnalysisOpts.Project,'Rule Representation')
                %%
                %%% find when the trial is starting and Ending
                StartTrialInd=strcmp(AnalysisOpts.TrialTimesFields,'START_TRIAL'); %% we don't use this here but still good to have it
                StopTrialInd=strcmp(AnalysisOpts.TrialTimesFields,'END_TRIAL');%% we don't use this here but still good to have it

                %%% find the indices related to the times we want to look int the task
                StartCountInd=strcmp(AnalysisOpts.TrialTimesFields,AnalysisOpts.TrialTiming.StartFieldName); %% time to start chopping
                %StopCountInd=strcmp(AnalysisOpts.TrialTimesFields,AnalysisOpts.LFPParams.StopFiledName); %% time to stop chopping
                % we are not using the stop time because we want to center analysis on a
                % specific event
                RewardInd=strcmp(AnalysisOpts.TrialTimesFields,'CORRECT_TRIAL');  %% reward dosn't exist but we still put it here
                %%% define time axis
                if ~obj.UseDataPointer
                    Time=(0:1/obj.Fs:(length(data)-1)/obj.Fs);
                else
                    eval(['Npoints=size(data.' obj.DataPointerVar ',2);'])
                    Time=(0:1/obj.Fs:(Npoints-1)/obj.Fs);
                end
                Trial_Data=[];

                % get the trail data
                for Tri=1:size(TrialTimes,1)  %loop through trials
                    fprintf('\n*** Grabing data from trial: %s ***\n ', num2str(Tri));
                    RewardTime =TrialTimes(Tri,RewardInd);  %% reward time
                    StartTime  =TrialTimes(Tri,StartCountInd)-AnalysisOpts.TrialTiming.BaselineDelay;  %% start of chunk time
                    StopTime   =StartTime+AnalysisOpts.TrialTiming.PeriodLength;   %% stop of chunk time
                    Duration   =StopTime-StartTime;   % duration of chunk

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run your analysis on each
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% trial

                    if    strcmpi(obj.AnalysisType,'LFP_Trial')
                        for DesPar=1:size(DesiredParameter,1) %% now loop on the desired paramters you want such as frequncy band
                            RawData=data(DesPar,Time>=StartTime & Time<StopTime)';
                            %% apply downsampling analysis we want
                            [RawData] = ApplyNotchFilter(RawData,AnalysisOpts.NeuralAnalysis.SamplingFrequency); % notch filter the data
                            RawData = ApplyLFPLowPass(RawData,AnalysisOpts.NeuralAnalysis.SamplingFrequency);  % apply low pass filter to only look at lfp
                            Trial_Data(:,Tri,DesPar)=downsample(RawData,obj.Fs/AnalysisOpts.LFPParams.FilterOpts.TargetLFPSamplingFrequency);        % downsample the data to a reasonable range
                        end

                    elseif strcmpi(obj.AnalysisType,'Wavelet_Trial')
                        TFData=data(:,Time>=StartTime & Time<StopTime);
                        Trial_Data(:,:,Tri)=transpose(downsample(TFData',obj.WaveletDownSampleFactor));        % downsample the data to a reasonable range

                    elseif strcmpi(obj.AnalysisType,'Wavelet_Trial_H')  % reorganize H matrix
                        TFData=data(:,Time>=StartTime & Time<StopTime);
                        Trial_Data=[Trial_Data transpose(downsample(TFData',obj.WaveletDownSampleFactor))];        % downsample the data to a reasonable range
                    end
                end
            end
        end
        function BhvModelInfo=GetBhvModelInfo(obj,RuleBlockTrials,TrialTimes,varargin)
            % retrieves model information if we are using it
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs

            if AnalysisOpts.Bhv.UseModel % incorporate Bhv model in all of our block spec calculations
                % this means to align everything to the switch
                % trial
                ModelFileNameTemplate='Hard_switch_6_channels_surprise_precision_color_pref.mat';%'Hard_switch_12_channels_surprise_color_pref.mat'
                ModelFile=[AnalysisOpts.RecDatePath 'bhv' AnalysisOpts.FS 'Model' AnalysisOpts.FS ModelFileNameTemplate ];
                [ModelFileFold]=fileparts(ModelFile);
                if ~exist(ModelFile,'file')
                    fprintf('\nModel file does not exist... trying to copy it from source');
                    % check if the model exists at all?
                    OrgModelFile=[AnalysisOpts.BhvModelPath ModelFileNameTemplate];
                    if exist(OrgModelFile,'file')
                        if ~exist(ModelFileFold,'file')
                            mkdir(ModelFileFold);
                        end
                        copyfile(OrgModelFile,ModelFile)
                        disp('\nModel coppied from source...');
                    else
                        disp('\nModel source does not exist...');
                        BhvModelInfo=[];
                        return;
                    end
                end
                ModelInfo=load(ModelFile,'Model_predictions');
                %% do some preprocessing on the data
                %% down sample Value for Choice to 4 bins
                inG=ModelInfo.Model_predictions.model_input.inG;
                N_channels=4;
                ModelInfo.Model_predictions.model_outputs.down_sampled_Value_for_choice4=NaN(N_channels,inG.N_forward_trials,inG.N_learning);
                for i=1:inG.N_forward_trials
                    for j=1:inG.N_learning
                        if ~isnan(inG.X_learning(i,11,j))
                            for n=1:N_channels
                                ModelInfo.Model_predictions.model_outputs.down_sampled_Value_for_choice4(n,i,j)=mean(ModelInfo.Model_predictions.model_outputs.Value_for_choice(1+(n-1)*round(size(ModelInfo.Model_predictions.model_outputs.Value_for_choice,1)/N_channels):n*round(size(ModelInfo.Model_predictions.model_outputs.Value_for_choice,1)/N_channels),i,j),1);
                            end
                        end
                    end
                end
                %% Copy data into a strcture similar to what we have
                ModelOutputs=ModelInfo.Model_predictions.model_outputs;
                BhvModelInfo.model_outputs.NBlocks=size(ModelOutputs.RPE,2);
                BhvModelInfo.model_outputs.BlockLength=arrayfun(@(x) RuleBlockTrials.Seq(x,2)-RuleBlockTrials.Seq(x,1)+1,1:BhvModelInfo.model_outputs.NBlocks);
                % define a function to grab the data
                GrabData2=@(Y,Z) arrayfun(@(x) Y.(Z)(1:BhvModelInfo.model_outputs.BlockLength(x),x),1:BhvModelInfo.model_outputs.NBlocks,'UniformOutput',0);
                GrabData3=@(Y,Z) arrayfun(@(x) Y.(Z)(:,1:BhvModelInfo.model_outputs.BlockLength(x),x),1:BhvModelInfo.model_outputs.NBlocks,'UniformOutput',0);
                %% copy all of the data from model output into a new structure that matches what we have
                Model_output_fields=fieldnames(ModelInfo.Model_predictions.model_outputs);
                for f=1:length(Model_output_fields)
                    Siz1=size(ModelOutputs.(Model_output_fields{f}),1);
                    Siz2=size(ModelOutputs.(Model_output_fields{f}),2);
                    Siz3=size(ModelOutputs.(Model_output_fields{f}),3);
                    if  Siz3>1
                        %     BhvModelInfo.model_outputs.precision_ind=arrayfun(@(x) ModelOutputs.precision_ind(1:BhvModelInfo.BlockLength(x),x),1:BhvModelInfo.NBlocks,'UniformOutput',0);
                        BhvModelInfo.model_outputs.(Model_output_fields{f})= GrabData3(ModelOutputs,Model_output_fields{f});
                    elseif  Siz1>1 && Siz2>1
                        BhvModelInfo.model_outputs.(Model_output_fields{f})= GrabData2(ModelOutputs,Model_output_fields{f});
                    elseif (Siz1==1 && Siz2>1) || (Siz1>1 && Siz2==1)
                        BhvModelInfo.model_outputs.(Model_output_fields{f})= ModelOutputs.(Model_output_fields{f});
                    end
                end

                %% Align all of the blocks to the switch trial
                BhvModelInfo.BlockSpec.TrialOrder=-obj.NTrial2Swtch:obj.MaxNTrial;
                GrabTrialData=@(x,y)  [nan*ones(length(find(x<1)),size(y,2)) ;y(x(x>=1),:)];
                GrabTrialData2=@(x,y)  [nan*ones(size(y,1),length(find(x<1))) y(:,x(x>=1))];
                for b=1:BhvModelInfo.model_outputs.NBlocks
                    % if there is a switch in this block then align to that
                    if BhvModelInfo.model_outputs.Total_Switch(b)==1
                        SwitchTrl=BhvModelInfo.model_outputs.Total_Switch(b)*BhvModelInfo.model_outputs.Total_Switch_when(b)-1;
                    else
                        SwitchTrl=0;
                    end
                    BhvModelInfo.BlockSpec.ThisBlkTrials{b}=RuleBlockTrials.Rule{b}(1)+SwitchTrl+BhvModelInfo.BlockSpec.TrialOrder;
                    % check if trial zero is really switch
                    % trial
                    %                         SwitchIndTrlTimes=find(diff(TrialTimes(BlockSpec.ThisBlkTrials{1},6)==1))+1;
                    %                         SwitchIndTrlOrder=find(BlockSpec.TrialOrder==0);
                    %                         if SwitchIndTrlTimes~=SwitchIndTrlOrder
                    %                             error('Trouble finding switch trial')
                    %                         end
                    BhvModelInfo.BlockSpec.ThisBlkTrialTimes{b}=GrabTrialData(BhvModelInfo.BlockSpec.ThisBlkTrials{b},TrialTimes);
                    BhvModelInfo.BlockSpec.Rule(b)=b;
                    BhvModelInfo.BlockSpec.ThisBlkOrder(b)=RuleBlockTrials.Rule{b}(3);
                    BhvModelInfo.BlockSpec.ThisBlkTrialNum(b)=length(BhvModelInfo.BlockSpec.ThisBlkTrials{b});
                    %    end
                end
                %% align all of the model outputs to the switch trial
                WantedFields={'RPE','Value_for_choice','down_sampled_Value_for_choice4','mean_ind_bin','precision_ind','mean_ind'};
                for f=WantedFields%1:length(Model_output_fields)
                    if size(BhvModelInfo.model_outputs.(f{1}){1},2)>1
                        ThisVar=cell2mat(BhvModelInfo.model_outputs.(f{1}));
                        for b=1:BhvModelInfo.model_outputs.NBlocks
                            BhvModelInfo.BlockSpec.(f{1}){b}=GrabTrialData2(BhvModelInfo.BlockSpec.ThisBlkTrials{b},ThisVar);
                        end
                    else
                        ThisVar=cell2mat(BhvModelInfo.model_outputs.(f{1})');
                        for b=1:BhvModelInfo.model_outputs.NBlocks
                            BhvModelInfo.BlockSpec.(f{1}){b}=GrabTrialData(BhvModelInfo.BlockSpec.ThisBlkTrials{b},ThisVar);
                        end
                    end
                end
            end
        end
        function BlockSpec=GetBlockTrialInfo(obj,RuleBlockTrials,TrialTimes,varargin)
            % gets the information and timing for each block and its
            % trials
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs

            if strcmpi(AnalysisOpts.Project,'Learning attentional templates')
                % get the number of trials for each block and remove the ones
                % don't have maxntrials
                if ~isempty(obj.MaxNTrial)
                    NTrlsBlk=RuleBlockTrials.Seq(:,2)-RuleBlockTrials.Seq(:,1);
                    WantedBlk=find(NTrlsBlk>=obj.MaxNTrial);
                    RuleBlockTrials.Rule=RuleBlockTrials.Rule(WantedBlk);
                    RuleBlockTrials.Seq=RuleBlockTrials.Seq(WantedBlk,:);
                end
                NBlks=length(RuleBlockTrials.Rule);
                if strcmpi(obj.Nblk_2Look,'MAX');obj.Nblk_2Look=NBlks;end %
                ThisBlks=min(NBlks,obj.Nblk_2Look);
                if ~(isempty(obj.MaxNTrial) || isempty(obj.NTrial2Swtch));ThisBlks=ThisBlks-1;end
                for b=1:ThisBlks
                    % if we don't have max num trials then set it to max
                    if isempty(obj.MaxNTrial)     % if we want to look at all of the block
                        BlockSpec.ThisBlkTrials{b}=...
                            RuleBlockTrials.Rule{b}(1):RuleBlockTrials.Rule{b}(2);
                        BlockSpec.TrialOrder=0:length(BlockSpec.ThisBlkTrials{b})-1;
                    elseif ~isempty(obj.MaxNTrial) && isempty(obj.NTrial2Swtch) % if we only want to look at after switch
                        BlockSpec.ThisBlkTrials{b}=...
                            RuleBlockTrials.Rule{b}(1):RuleBlockTrials.Rule{b}(1)+obj.MaxNTrial-1;
                        BlockSpec.TrialOrder=0:obj.MaxNTrial-1;
                    elseif ~(isempty(obj.MaxNTrial) || isempty(obj.NTrial2Swtch)) %&& ~AnalysisOpts.Bhv.UseModel % if we want to look before and after switch
                        BlockSpec.ThisBlkTrials{b}=...
                            [RuleBlockTrials.Rule{b}(2)-obj.NTrial2Swtch+1:RuleBlockTrials.Rule{b}(2) RuleBlockTrials.Rule{b+1}(1):RuleBlockTrials.Rule{b+1}(1)+obj.MaxNTrial-1];
                        BlockSpec.TrialOrder=-obj.NTrial2Swtch:obj.MaxNTrial-1;
                        % check if trial zero is really switch
                        % trial
                        SwitchIndTrlTimes=find(diff(TrialTimes(BlockSpec.ThisBlkTrials{1},6)==1))+1;
                        SwitchIndTrlOrder=find(BlockSpec.TrialOrder==0);
                        if SwitchIndTrlTimes~=SwitchIndTrlOrder
                            error('Trouble finding switch trial')
                        end

                    end
                    BlockSpec.ThisBlkTrialTimes{b}=TrialTimes(BlockSpec.ThisBlkTrials{b},:);
                    BlockSpec.Rule(b)=b;
                    BlockSpec.ThisBlkOrder(b)=RuleBlockTrials.Rule{b}(3);
                    BlockSpec.ThisBlkTrialNum(b)=length(BlockSpec.ThisBlkTrials{b});
                end
            elseif strcmpi(AnalysisOpts.Project,'Rule Representation')
                if strcmp(AnalysisOpts.RecDate,'061518')
                    obj.MinTrlsRuleBlk=[70 40 70];
                    warning('Changing MinTrlsRuleBlk for 061518 recording...')
                elseif strcmp(AnalysisOpts.RecDate,'073021')
                    obj.MinTrlsRuleBlk=[100 35 100];
                    warning('Changing MinTrlsRuleBlk for 073021 recording...')
                end
                % get trails after the switch, before the switch and around
                % the switch
                NonTrl=@(x) find(x<1); % get a function to get out the negative and zero trials numbers
                Trls=@(x) find(x>0);
                OrgTrlTimes=@(x,y) [nan*ones(length(NonTrl(x)),size(y,2)) ; y(x(Trls(x)),:)];
                %% loop on the rule
                for rule=obj.Rule_2look
                    for b=1:size(RuleBlockTrials.Rule{rule},1)
                        BlkCnt=RuleBlockTrials.Rule{rule}(b,3);
                        % number of trials in this block
                        NTrlsBlk=RuleBlockTrials.Rule{rule}(b,2)-RuleBlockTrials.Rule{rule}(b,1)+1;
                        BlockSpec.ThisBlkTrialNum(BlkCnt)=NTrlsBlk;

                        if NTrlsBlk<obj.MaxNTrial | NTrlsBlk<obj.MinTrlsRuleBlk(rule)
                            BlockSpec.ThisBlkTrials.FromSwitch{BlkCnt}=[];BlockSpec.ThisBlkTrials.ToSwitch{BlkCnt}=[];BlockSpec.ThisBlkTrials.AroundSwitch{BlkCnt}=[];
                            BlockSpec.BlockIncluded(BlkCnt)=false;
                        else
                            % if we don't have max num trials then set it to max
                            if isempty(obj.MaxNTrial)
                                BlockSpec.ThisBlkTrials.FromSwitch{BlkCnt}=...
                                    RuleBlockTrials.Rule{rule}(b,1):RuleBlockTrials.Rule{rule}(b,2);
                            else
                                % Trials forward from switch
                                BlockSpec.ThisBlkTrials.FromSwitch{BlkCnt}=...
                                    RuleBlockTrials.Rule{rule}(b,1):RuleBlockTrials.Rule{rule}(b,1)+obj.MaxNTrial-1;
                                % Trials Backward to switch
                                BlockSpec.ThisBlkTrials.ToSwitch{BlkCnt}=...
                                    RuleBlockTrials.Rule{rule}(b,2)-obj.NTrial2Swtch+1:RuleBlockTrials.Rule{rule}(b,2);
                                % Trials around the switch
                                BlockSpec.ThisBlkTrials.AroundSwitch{BlkCnt}=...
                                    RuleBlockTrials.Rule{rule}(b,1)-obj.NTrial2Swtch:RuleBlockTrials.Rule{rule}(b,1)+obj.MaxNTrial-1;
                                % take all of the trials
                                BlockSpec.ThisBlkTrials.AllTrials{BlkCnt}=...
                                    RuleBlockTrials.Rule{rule}(b,1):RuleBlockTrials.Rule{rule}(b,2); % take all trials of this block
                            end

                            BlockSpec.ThisBlkTrialTimes.FromSwitch{BlkCnt}  =OrgTrlTimes(BlockSpec.ThisBlkTrials.FromSwitch{BlkCnt},TrialTimes);
                            BlockSpec.ThisBlkTrialTimes.ToSwitch{BlkCnt}    =OrgTrlTimes(BlockSpec.ThisBlkTrials.ToSwitch{BlkCnt},TrialTimes);
                            BlockSpec.ThisBlkTrialTimes.AroundSwitch{BlkCnt}=OrgTrlTimes(BlockSpec.ThisBlkTrials.AroundSwitch{BlkCnt},TrialTimes);
                            BlockSpec.ThisBlkTrialTimes.AllTrials{BlkCnt}   =OrgTrlTimes(BlockSpec.ThisBlkTrials.AllTrials{BlkCnt},TrialTimes);

                            BlockSpec.Rule(BlkCnt)=rule;
                            BlockSpec.ThisBlkOrder(BlkCnt)=RuleBlockTrials.Rule{rule}(b,3);
                            BlockSpec.BlockIncluded(BlkCnt)=true;
                            BlockSpec.SeqHist(BlkCnt)=RuleBlockTrials.SeqHist(BlkCnt); % history of rule sequence
                            BlockSpec.TrialNum{BlkCnt}=1:NTrlsBlk;
                        end
                    end
                    BlockSpec.TrialOrder(rule).AroundSwitch=-obj.NTrial2Swtch:obj.MaxNTrial-1;
                    BlockSpec.TrialOrder(rule).ToSwitch=0:obj.MaxNTrial-1;
                    BlockSpec.TrialOrder(rule).FromSwitch=-obj.NTrial2Swtch:-1;
                    BlockSpec.TrialOrder(rule).AllTrials=nan;
                end

                %% if you want to keep the sequence of blocks as they appear in the experiment 
%                 BlockSequence=RuleBlockTrials.Seq(:,3)';RuleCnt=[0 0 0];
%                 %for rule=obj.Rule_2look
%                 %  for b=1:size(RuleBlockTrials.Rule{rule},1)
%                 for BlkSeq=1:length(BlockSequence)
%                     rule=BlockSequence(BlkSeq);
%                     RuleCnt(rule)=RuleCnt(rule)+1;
%                     b=RuleCnt(rule);
%                     BlkCnt=RuleBlockTrials.Rule{rule}(b,3);
%                     % number of trials in this block
%                     NTrlsBlk=RuleBlockTrials.Rule{rule}(b,2)-RuleBlockTrials.Rule{rule}(b,1)+1;
%                     BlockSpec.ThisBlkTrialNum(BlkCnt)=NTrlsBlk;
% 
%                     if NTrlsBlk<obj.MaxNTrial | NTrlsBlk<obj.MinTrlsRuleBlk(rule)
%                         BlockSpec.ThisBlkTrials.FromSwitch{BlkCnt}=[];BlockSpec.ThisBlkTrials.ToSwitch{BlkCnt}=[];BlockSpec.ThisBlkTrials.AroundSwitch{BlkCnt}=[];
%                         BlockSpec.BlockIncluded(BlkCnt)=false;
%                     else
%                         % if we don't have max num trials then set it to max
%                         if isempty(obj.MaxNTrial)
%                             BlockSpec.ThisBlkTrials.FromSwitch{BlkCnt}=...
%                                 RuleBlockTrials.Rule{rule}(b,1):RuleBlockTrials.Rule{rule}(b,2);
%                         else
%                             % Trials forward from switch
%                             BlockSpec.ThisBlkTrials.FromSwitch{BlkCnt}=...
%                                 RuleBlockTrials.Rule{rule}(b,1):RuleBlockTrials.Rule{rule}(b,1)+obj.MaxNTrial-1;
%                             % Trials Backward to switch
%                             BlockSpec.ThisBlkTrials.ToSwitch{BlkCnt}=...
%                                 RuleBlockTrials.Rule{rule}(b,2)-obj.NTrial2Swtch+1:RuleBlockTrials.Rule{rule}(b,2);
%                             % Trials around the switch
%                             BlockSpec.ThisBlkTrials.AroundSwitch{BlkCnt}=...
%                                 RuleBlockTrials.Rule{rule}(b,1)-obj.NTrial2Swtch:RuleBlockTrials.Rule{rule}(b,1)+obj.MaxNTrial-1;
%                             % take all of the trials
%                             BlockSpec.ThisBlkTrials.AllTrials{BlkCnt}=...
%                                 RuleBlockTrials.Rule{rule}(b,1):RuleBlockTrials.Rule{rule}(b,2); % take all trials of this block
%                         end
% 
%                         BlockSpec.ThisBlkTrialTimes.FromSwitch{BlkCnt}  =OrgTrlTimes(BlockSpec.ThisBlkTrials.FromSwitch{BlkCnt},TrialTimes);
%                         BlockSpec.ThisBlkTrialTimes.ToSwitch{BlkCnt}    =OrgTrlTimes(BlockSpec.ThisBlkTrials.ToSwitch{BlkCnt},TrialTimes);
%                         BlockSpec.ThisBlkTrialTimes.AroundSwitch{BlkCnt}=OrgTrlTimes(BlockSpec.ThisBlkTrials.AroundSwitch{BlkCnt},TrialTimes);
%                         BlockSpec.ThisBlkTrialTimes.AllTrials{BlkCnt}   =OrgTrlTimes(BlockSpec.ThisBlkTrials.AllTrials{BlkCnt},TrialTimes);
% 
%                         BlockSpec.Rule(BlkCnt)=rule;
%                         BlockSpec.ThisBlkOrder(BlkCnt)=RuleBlockTrials.Rule{rule}(b,3);
%                         BlockSpec.BlockIncluded(BlkCnt)=true;
%                         BlockSpec.SeqHist(BlkCnt)=RuleBlockTrials.SeqHist(BlkCnt); % history of rule sequence
%                         BlockSpec.TrialNum{BlkCnt}=1:NTrlsBlk;
%                     end
% 
%                     BlockSpec.TrialOrder(rule).AroundSwitch=-obj.NTrial2Swtch:obj.MaxNTrial-1;
%                     BlockSpec.TrialOrder(rule).ToSwitch=0:obj.MaxNTrial-1;
%                     BlockSpec.TrialOrder(rule).FromSwitch=-obj.NTrial2Swtch:-1;
%                     BlockSpec.TrialOrder(rule).AllTrials=nan;
%                 end
                if ~sum(BlockSpec.ThisBlkTrialNum)==size(TrialTimes,1)
                    error('TrialTimes and Block Specs are not matching');
                end
            end
        end
        function RuleBlockTrials=SegregateBlocks(obj,TrialTimes)
            % finds trials for each rule and blocks
            global AnalysisOpts
            if strcmpi(AnalysisOpts.Project,'Learning attentional templates')
                BlockInd=strcmp(AnalysisOpts.TrialTimesFields,'BlockNo');%% we don't use this here but still good to have it
                Ntrials=size(TrialTimes,1);
                TrialBlock=TrialTimes(:,BlockInd)';
                temp=find(diff(TrialBlock)~=0);
                if temp(end)~=Ntrials
                    SwitchTrials=[0 temp Ntrials];
                else
                    SwitchTrials=[0 temp];
                end
                for i=1:length(SwitchTrials)-1
                    TrialInterval=SwitchTrials(i)+1:SwitchTrials(i+1);
                    Block=unique(TrialBlock(TrialInterval));
                    RuleBlockTrials.Rule{Block}=[SwitchTrials(i)+1 SwitchTrials(i+1) i];
                    RuleBlockTrials.Seq(i,:)=[SwitchTrials(i)+1 SwitchTrials(i+1) Block];
                end


            elseif strcmpi(AnalysisOpts.Project,'Rule Representation')
                RuleInd=strcmp(AnalysisOpts.TrialTimesFields,'CONDITION_NUM_OFFSET');%% we don't use this here but still good to have it
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
                % determine if block with the same axis before this block(2 blocks before) have the same rule(0) or not(1)
                RuleSeq=RuleBlockTrials.Seq(:,3)';
                RuleBlockTrials.SeqHist=nan(1,size(RuleBlockTrials.Seq,1));
                RuleBlockTrials.SeqHist(1:2:end)=[false logical(diff(RuleSeq(1:2:end)))];
                RuleBlockTrials.SeqHist(2:2:end)=[false logical(diff(RuleSeq(2:2:end)))];
            end
        end
        function [TrialTimes]=FindWantedTrialinArea(obj)
            %%% Finds a specific trial types in a specific area in a recording
            %%% Inputs:
            %%% TrialType: can take 'ALL', 'CORRECT', 'INCORRECT'
            %%% Area: Cell that contains the name of the Areas we are intrested in
            global AnalysisOpts
            if strcmpi(AnalysisOpts.Project,'Learning attentional templates')
                [TrialTimes]=obj.ExtractTrialTimes;

            elseif strcmpi(AnalysisOpts.Project,'Rule Representation')
                if isempty(AnalysisOpts.RecDate) % change the recdate to first recoridng day and give  awarning 
                    AnalysisOpts.RecDate=AnalysisOpts.DateSet{1};
                    warning('Recdate was empty. it was changed to %s',AnalysisOpts.RecDate)                    
                end
                %%%  get the trial times
                load([AnalysisOpts.RecDatePath 'Rec' AnalysisOpts.RecDate '_TimingData.mat' ],'TrialInfo','AdjTiming');
                %% combine ajusted timings with digital timings
                mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);%This works IF AND ONLY IF there are no common fields (duplicate fields) in the 2 structures.
                CombinedTrialInfo=mergestructs(TrialInfo,AdjTiming);
                clear TrialInfo AdjTiming
                % replace adjested timings in the strcuture
                for ntrl=1:length(CombinedTrialInfo) % unfortunately I am using for here but I might find a more elegant solution later
                    CombinedTrialInfo(ntrl).SAMPLE_ON=CombinedTrialInfo(ntrl).SAMPLE_ON_Photdiode; % replace Sample on with photodiode value
                    CombinedTrialInfo(ntrl).SAMPLE_OFF=CombinedTrialInfo(ntrl).SAMPLE_OFF_Photdiode; % replace sample off with photodiode value
                    CombinedTrialInfo(ntrl).FIXATE_ACQUIRED=CombinedTrialInfo(ntrl).FIXATE_ACQUIRED(1); % the first fix acuiqred is for sample on the second is for response
                end
                %%% Extract timing and information for each trial
                [TrialTimes]=obj.ExtractTrialTimes(CombinedTrialInfo);
                %%% from these trials take the ones with this specific rule and cut the
                %%% trials we don't want
                if ~strcmp(AnalysisOpts.Trial.TrialsToAnalyz,'ALL')
                    RuleFieldInd=strcmp(AnalysisOpts.TrialTimesFields,'CONDITION_NUM_OFFSET');
                    RuleTrialsInd=arrayfun(@(x) find(TrialTimes(:,RuleFieldInd)==x),AnalysisOpts.Trial.Rule,'UniformOutput',0);
                    TrialTimes=TrialTimes(RuleTrialsInd{1}(1:AnalysisOpts.Trial.TrialsToAnalyz),:);  % cut this rul trials
                end
            end

        end
        function [WantedTrials,WantedTrialTimes,WantedStimMap]=FindWantedTrialType(~,TrialTimes,Rule,Category,MorphLevel,TrialType,varargin)%gets wanted trial type in TrialTimes matrix
            % Finds a specific trial
            %BlockSpecs: TrialTimes mat
            %Rule:Rule [1 2 3]
            %Category:{'ShapeCategory','ColorCategory'}; Shape: 1 bunny 2 Tee /1 Red 2 Green Examle {[1],[2]}
            %or {'All',[2]}
            %MoprhLevel: {'ShapeML' 'ColorML'}; Example{[0:200],[0:100]} or
            %{'ALL',[0:100]

            % TrialType: can take 'ALL', 'CORRECT', 'INCORRECT'
            global AnalysisOpts

            if isempty(TrialType);TrialType='ALL';end
            if isempty(Rule);Rule=[1 2 3];end
            for i=find(strcmp(Category,'ALL'));Category{i}=[1 2 -1];end
            for i=find(strcmp(MorphLevel,'ALL'));MorphLevel{i}=[0:5:200];end

            TrialTypeInd=strcmp(AnalysisOpts.TrialTimesFields,'CORRECT_TRIAL');
            RuleInd=strcmp(AnalysisOpts.TrialTimesFields,'CONDITION_NUM_OFFSET');
            StimNumInd=strcmp(AnalysisOpts.TrialTimesFields,'SAMPLE_NUM_OFFSET');
            StimNum=TrialTimes(:,StimNumInd);

            TrialTimesStimMap=nan*ones(length(StimNum),4);
            GoodTrl=find(~isnan(StimNum'));
            TrialTimesStimMap(GoodTrl,:)=AnalysisOpts.StimMap.cuemapshape(StimNum(GoodTrl),2:end);
            TrialTimesRule=TrialTimes(:,RuleInd)';
            TrialTimesTrialType=~isnan(TrialTimes(:,TrialTypeInd))';

            % now get what we want
            ShapeMLStimMapInd=strcmp(AnalysisOpts.StimMap.FieldNames,'ShapeML');
            ColorMLStimMapInd=strcmp(AnalysisOpts.StimMap.FieldNames,'ColorML');
            ShapeCatStimMapInd=strcmp(AnalysisOpts.StimMap.FieldNames,'ShapeCat');
            ColorCatStimMapInd=strcmp(AnalysisOpts.StimMap.FieldNames,'ColorCat');

            RuleTrls=arrayfun(@(x) sum(x==Rule)>0,TrialTimesRule);
            MorphLevelTrlsShape=arrayfun(@(x) sum(x==MorphLevel{1})>0,TrialTimesStimMap(:,ShapeMLStimMapInd))';
            MorphLevelTrlsColor=arrayfun(@(x) sum(x==MorphLevel{2})>0,TrialTimesStimMap(:,ColorMLStimMapInd))';
            CategoryTrlShape   =arrayfun(@(x) sum(x==Category{1})>0,TrialTimesStimMap(:,ShapeCatStimMapInd))';
            CategoryTrlColor   =arrayfun(@(x) sum(x==Category{2})>0,TrialTimesStimMap(:,ColorCatStimMapInd))';
            if strcmp(TrialType,'CORRECT');TrialTypeTrls=TrialTimesTrialType==1;
            elseif strcmp(TrialType,'INCORRECT');TrialTypeTrls=TrialTimesTrialType~=1;
            elseif strcmp(TrialType,'ALL');TrialTypeTrls=ones(1,length(TrialTimesTrialType));
            end

            WantedTrials=RuleTrls & MorphLevelTrlsShape & MorphLevelTrlsColor & CategoryTrlShape & CategoryTrlColor & TrialTypeTrls;
            WantedTrialTimes=TrialTimes(WantedTrials,:);
            WantedStimMap=TrialTimesStimMap(WantedTrials,:);
        end
        function [NewBlockSpec]=FindWantedTrialTypeBlockSpec(obj,BlockSpec,Rule,Category,MorphLevel,TrialType,varargin)% finds trial type in BlockSpec
            % gets the wanted trials from each block of the blockSpec
            % matrix
            global AnalysisOpts
            if isempty(BlockSpec)% then load it from the recording;
                [TrialTimes]=obj.FindWantedTrialinArea;   % retreive channel information of this rec
                [RuleBlockTrials]=obj.SegregateBlocks(TrialTimes); % get the sequence of trials for each block
                BlockSpec=obj.GetBlockTrialInfo(RuleBlockTrials,TrialTimes,'MaxNTrial',AnalysisOpts.Bhv.MaxNTrial,'NTrial2Swtch',AnalysisOpts.Bhv.NTrial2Swtch,'Nblk_2Look','MAX'); % get blockspecs
            end
            AnalysisOpts.NTrlsRec=sum(BlockSpec.ThisBlkTrialNum(BlockSpec.BlockIncluded));

            if isempty(Rule);Rule=[1 2 3];end
            Cnt=1;
            for R=Rule
                RuleInd=find(BlockSpec.Rule==R);
                for b=RuleInd
                    for StrName=AnalysisOpts.BlockSpecsFeilds
                        StrName=StrName{1};
                        if strcmp(StrName,'AroundSwitch');ThisBlkRule=[1 2 3];else;ThisBlkRule=R;end
                        [temp,NewBlockSpec.ThisBlkTrialTimes.(StrName){Cnt},NewBlockSpec.WantedStimMap.(StrName){Cnt}]=...
                            obj.FindWantedTrialType(BlockSpec.ThisBlkTrialTimes.(StrName){b},ThisBlkRule,Category,MorphLevel,TrialType);
                        NewBlockSpec.ThisBlkTrials.(StrName){Cnt}=BlockSpec.ThisBlkTrials.(StrName){b}(temp);

                        % get information from differnt types of behavioral models
                        for bhvmdl=AnalysisOpts.Bhv.BhvModelTypes
                            NewBlockSpec.([bhvmdl{1} 'BhvMdl']).(StrName){Cnt}=obj.GetRuleRepresentationBhvMdl(BlockSpec.ThisBlkTrials.(StrName){b}(temp),bhvmdl{1});
                        end

                    end
                    NewBlockSpec.ThisBlkTrialNum(Cnt)=BlockSpec.ThisBlkTrialNum(b);
                    NewBlockSpec.Rule(Cnt)=BlockSpec.Rule(b);
                    NewBlockSpec.ThisBlkOrder(Cnt)=BlockSpec.ThisBlkOrder(b);
                    NewBlockSpec.SeqHist(Cnt)=BlockSpec.SeqHist(b);
                    Cnt=Cnt+1;
                end
            end
            NewBlockSpec.TrialOrder=BlockSpec.TrialOrder;
            %% if you want to have sequence of blocks are the appear in the experiment 
%             if isempty(Rule);Rule=[1 2 3];end
%             Cnt=1;
%             for b=BlockSpec.ThisBlkOrder
%              R=BlockSpec.Rule(b);
%                     for StrName=AnalysisOpts.BlockSpecsFeilds
%                         StrName=StrName{1};
%                         if strcmp(StrName,'AroundSwitch');ThisBlkRule=[1 2 3];else;ThisBlkRule=R;end
%                         [temp,NewBlockSpec.ThisBlkTrialTimes.(StrName){Cnt},NewBlockSpec.WantedStimMap.(StrName){Cnt}]=...
%                             obj.FindWantedTrialType(BlockSpec.ThisBlkTrialTimes.(StrName){b},ThisBlkRule,Category,MorphLevel,TrialType);
%                         NewBlockSpec.ThisBlkTrials.(StrName){Cnt}=BlockSpec.ThisBlkTrials.(StrName){b}(temp);
% 
%                         % get information from differnt types of behavioral models
%                         for bhvmdl=AnalysisOpts.Bhv.BhvModelTypes
%                             NewBlockSpec.([bhvmdl{1} 'BhvMdl']).(StrName){Cnt}=obj.GetRuleRepresentationBhvMdl(BlockSpec.ThisBlkTrials.(StrName){b}(temp),bhvmdl{1});
%                         end
% 
%                     end
%                     NewBlockSpec.ThisBlkTrialNum(Cnt)=BlockSpec.ThisBlkTrialNum(b);
%                     NewBlockSpec.Rule(Cnt)=BlockSpec.Rule(b);
%                     NewBlockSpec.ThisBlkOrder(Cnt)=BlockSpec.ThisBlkOrder(b);
%                     NewBlockSpec.SeqHist(Cnt)=BlockSpec.SeqHist(b);
%                     Cnt=Cnt+1;
%              end
        end
        function out=GetRuleRepresentationBhvMdl(obj,TrialNum,ModelType,varargin)
            global AnalysisOpts
            obj=obj.ParseParams(varargin);

            % load the data from bhv model first
            LoadFileName=[AnalysisOpts.BhvMdlPath ModelType AnalysisOpts.FS AnalysisOpts.Animal AnalysisOpts.FS AnalysisOpts.Animal '_' AnalysisOpts.RecDate(5:6) AnalysisOpts.RecDate(1:4) '_00_bhv.mat'  ];

            BhvMdl=load(LoadFileName);
            % if the number of trials in this model doesn't match the total
            % number of trials then give a warnign and return NaN
            if AnalysisOpts.NTrlsRec~=size(BhvMdl.Rewards,2)
                warning('Bhv model does not match the behavior... aborting')
                AnalysisOpts.NonMatchingRecs=[AnalysisOpts.NonMatchingRecs,{AnalysisOpts.RecDate}];
                out=NaN;
                return
            end

            % get all of the data from fields we care about or all of the
            % feilds
            BhvMdlFields=fieldnames(BhvMdl);
            out=cellfun(@(x) obj.ManData.GetTrlsFromData(BhvMdl.(x),TrialNum),BhvMdlFields,'UniformOutput',0);
            out=cell2struct(out,BhvMdlFields,1); % convert back to a structure

        end
        function [TrialSpikeData,TimingStimData]=GrabWantedSpksWantedTrial(obj,WantedBlockSpec,SpikeData,varargin) % grabs spikes from wanted trials in blockspec
            global AnalysisOpts
            obj=obj.ParseParams(varargin);
            if isempty(SpikeData);TrialSpikeData=[];TimingStimData=[];return;end
            % loop on the cells
            NNeu=length(SpikeData);
            % TrialSpikeTime.(FieldName).Neuron.Block.Trial
            StrName=AnalysisOpts.TrlSpkTimeFieldName;
            RewardInd=strcmp(AnalysisOpts.TrialTimesFields,'REWARD_START');
            for Neu=1:NNeu % loop on channels
                %% get all of the information about this channel in one structure to make life easier!
                TrialSpikeData(Neu).TrialSpikeTime=cellfun(@(x) SpikeData(Neu).TrialSpikeTime(x),WantedBlockSpec.ThisBlkTrials.(StrName),'UniformOutput',0);
                TrialSpikeData(Neu).TrlSpkTimeFieldName=StrName;
                TrialSpikeData(Neu).Ch=SpikeData(Neu).Ch;
                TrialSpikeData(Neu).Cluster=SpikeData(Neu).Cluster;
                TrialSpikeData(Neu).ChannelArea=SpikeData(Neu).ChannelArea;
                TrialSpikeData(Neu).RecDate=SpikeData(Neu).RecDate;
                TrialSpikeData(Neu).Specs=SpikeData(Neu).Specs;
                %% get the information about timings and stimulus
                % information about timing and stimulus are here
                NBlks=length(WantedBlockSpec.ThisBlkTrials.(StrName));
                % TimingStimData(Neu).CriticalTimes=arrayfun(@(x) [SpikeData(Neu).CriticalTimes(WantedBlockSpec.ThisBlkTrials.(StrName){x},:) ],...
                %     1:NBlks,'UniformOutput',0);
                TimingStimData(Neu).CriticalTimes=arrayfun(@(x) obj.ExtractCriticalTimes(WantedBlockSpec.ThisBlkTrialTimes.(StrName){x}),1:NBlks,'UniformOutput',0);
                TimingStimData(Neu).TrialTimes=arrayfun(@(x) WantedBlockSpec.ThisBlkTrialTimes.(StrName){x},1:NBlks,'UniformOutput',0);
                TimingStimData(Neu).StimulusInfo=arrayfun(@(x) WantedBlockSpec.WantedStimMap.(StrName){x},1:NBlks,'UniformOutput',0);
                TimingStimData(Neu).Rule=WantedBlockSpec.Rule;
                TimingStimData(Neu).TrialOrder=WantedBlockSpec.TrialOrder.(StrName);
                TimingStimData(Neu).ThisBlkOrder=WantedBlockSpec.ThisBlkOrder;
                TimingStimData(Neu).SeqHist=WantedBlockSpec.SeqHist;
                TimingStimData(Neu).ThisBlkTrialNum=WantedBlockSpec.ThisBlkTrialNum;
                % calculate block performance and fit a logistic function
                TimingStimData(Neu).BlkPerf=cellfun(@(x) obj.bhvAna.CalculateBhvPerf(double(~isnan(x(:,RewardInd)))',AnalysisOpts.BhvAna.NTrlBlkPerf),TimingStimData(Neu).CriticalTimes,'UniformOutput',0);
%                 % calculte block performance for different lengths as well 
%                 for L=[5 10 40 50 60]
%                     TimingStimData(Neu).(['BlkPerf' num2str(L)])=cellfun(@(x) obj.bhvAna.CalculateBhvPerf(double(~isnan(x(:,RewardInd)))',L),TimingStimData(Neu).CriticalTimes,'UniformOutput',0);
%                 end
                
                % get the data from the behavioral model now
                for bhvmdl=AnalysisOpts.Bhv.BhvModelTypes
                    TimingStimData(Neu).([bhvmdl{1} 'BhvMdl'])=arrayfun(@(x) WantedBlockSpec.([bhvmdl{1} 'BhvMdl']).(StrName){x},1:NBlks,'UniformOutput',0);
                end
            end
        end
        function [cuemapshape,cue,FieldNames]=GenerateStimulusMapping(obj) % Generate stimuls mapping only for rule representation project
            global AnalysisOpts

            %% Define Cue Images
            srgb = iccfind(iccroot,'srgb');
            srgb = srgb{1};
            dell = iccfind(iccroot,'rgb');
            dell = dell{1};
            rgb2correct = makecform('icc',srgb,dell);
            Prop=[0];k=1;
            for Colprop=0:10:200
                IMAGEPATH=[AnalysisOpts.CuesPath 'Color_' num2str(Colprop/100) '/'];
                for p=1:length(Prop)
                    for m=0:10:200
                        cue(k).FileName = [IMAGEPATH 'SphericalMorph_P_' num2str(m/100,2) '_Prop_' num2str(Prop(p)) '_col_' num2str(Colprop/100)  '.png'];
                        cue(k).Array = imread(cue(k).FileName,'PNG');
                        %          %%%%%%%Stimulus Correciton Code
                        cue(k).Array = applycform(cue(k).Array,rgb2correct);
                        cuemapshape(k,:)=[p Colprop (m) 0 0];
                        k=k+1;
                    end
                end
            end
            Obj1Ind=find((cuemapshape(:,3)>=0 & cuemapshape(:,3)<50) | (cuemapshape(:,3)>150 & cuemapshape(:,3)<=200));
            Obj2Ind=find((cuemapshape(:,3)> 50 & cuemapshape(:,3)<150));
            RandObj=find((cuemapshape(:,3)==50) | (cuemapshape(:,3)==150));

            cuemapshape(Obj1Ind,4)=1;
            cuemapshape(Obj2Ind,4)=2;
            cuemapshape(RandObj,4)=-1;

            Col1Ind=find((cuemapshape(:,2)>=0 & cuemapshape(:,2)<50) | (cuemapshape(:,2)>150 & cuemapshape(:,2)<=200));
            Col2Ind=find((cuemapshape(:,2)> 50 & cuemapshape(:,2)<150));
            RandObjCol=find((cuemapshape(:,2)==50) | (cuemapshape(:,2)==150));

            cuemapshape(Col1Ind,5)=1;
            cuemapshape(Col2Ind,5)=2;
            cuemapshape(RandObjCol,5)=-1;
            FieldNames={'ColorML','ShapeML','ShapeCat','ColorCat'};
        end
        function [TrialTimes]=ExtractTrialTimes(obj,varargin)

            %%TrialType=0 All trials,1 only correct, 2 only incorrect
            global AnalysisOpts
            if strcmpi(AnalysisOpts.Project,'Learning attentional templates')
                % creates trial times matrix for Learning Attentional Templates project
                % CutTrialsFile is the file that was reconstrcuted from timings of NEV
                % file during recroding
                % BhvFile behavioral file. We need both these files to generate the
                % TrialTimes matrix
                CutTrialsFile=[AnalysisOpts.RecDatePath 'CutTrials' AnalysisOpts.FS 'CutTrials_' AnalysisOpts.Animal '_' AnalysisOpts.RecDate '_00_bhv.mat'];
                load(CutTrialsFile,'cutTrials');
                BhvFile=[AnalysisOpts.RecDatePath 'bhv' AnalysisOpts.FS AnalysisOpts.Animal '_' AnalysisOpts.RecDate '_00_bhv.mat'];
                load(BhvFile,'trials');

                %% get some params
                %%% fixed field names
                EphysCodes.Name={'FIX_ON','TARGET_ON','HOLD_TARGET_ON','GIVE_REWARD'}; % these are the codes that have been used in ephys cutTrials.Timing.NEV
                EphysCodes.Val=[2,12,21,22];

                %%% define Charactristics of trial that you are intrested in first
                TrlFeature={'FIX_ON','TARGET_ON','HOLD_TARGET_ON','GIVE_REWARD','RESPONSE_ON','BlockNo',...
                    'BestColor','ResponseColor','StopCondition','ResponseError','NumRewards'};

                ResponseOnInd=(strcmpi(TrlFeature,'RESPONSE_ON'));
                BlockNoInd=(strcmpi(TrlFeature,'BlockNo'));
                BestColorInd=(strcmpi(TrlFeature,'BestColor'));
                StopConditionInd=(strcmpi(TrlFeature,'StopCondition'));
                ResponseErrorInd=(strcmpi(TrlFeature,'ResponseError'));
                NumRewardsInd=(strcmpi(TrlFeature,'NumRewards'));
                ResponseColorInd=(strcmpi(TrlFeature,'ResponseColor'));

                Nfeatur=length(TrlFeature);
                Ntrls=length(cutTrials);

                TrialTimes=nan*ones(Ntrls,Nfeatur);
                for trl=1:Ntrls
                    %% look at NEV file and get required information
                    StageSequence=cutTrials(trl).Timing.NEV.StageSequence; % get first 4 timings
                    if ~isnan(StageSequence)
                        TrialTimes(trl,arrayfun(@(x) find(EphysCodes.Val==x),StageSequence))=cutTrials(trl).Timing.NEV.TimeStampSec;

                        TrialTimes(trl,ResponseOnInd) =cutTrials(trl).Timing.NEV.RESPONSE_ON_Sec; % get response on
                        %% look at BHV file and get condition details

                        TrialTimes(trl,BlockNoInd)  =trials(trl).Conditions.BlockNo; % block number
                        TrialTimes(trl,BestColorInd)=trials(trl).Conditions.BestColor; % best color
                        TrialTimes(trl,StopConditionInd)=trials(trl).Behavior.StopCondition; % stop condition
                        TrialTimes(trl,NumRewardsInd)=trials(trl).Reward.NumRewards; % num rewards

                        % find the deg of response error
                        if ~isnan(trials(trl).Behavior.ColorResponse)
                            Response=trials(trl).Behavior.ColorResponse;
                            BestTargInd=trials(trl).Conditions.BestTargetIndex;
                            %  TrialTimes(trl,ResponseErrorInd)=trials(trl).Conditions.TargetColorDistance(Response);
                            TrialTimes(trl,ResponseErrorInd)=trials(trl).Conditions.TargetColorAngles(Response)-trials(trl).Conditions.TargetColorAngles(BestTargInd);
                            TrialTimes(trl,ResponseColorInd)=trials(trl).Conditions.TargetColorAngles(Response);
                        end
                    end
                end
                % remove trials where there is no response
                UseableTrilInd=~isnan(TrialTimes(:,NumRewardsInd));
                TrialTimes=TrialTimes(UseableTrilInd,:);
                %
                AnalysisOpts.TrialTimesFields=TrlFeature;

                %%% Determine if this field has any time information
                AnalysisOpts.FieldIden=zeros(1,Nfeatur);
                AnalysisOpts.FieldIden(1:5)=1; % field 1 to 5 have timing info

            elseif strcmpi(AnalysisOpts.Project,'Rule Representation')
                TrialInfoMat=varargin{1};
                %%% fixed field names
                StartEncodeInd      ='START_TRIAL';
                StopEncodeInd       ='END_TRIAL';
                RuleInd             ='CONDITION_NUM_OFFSET';
                CorrectTrialInd     ='CORRECT_TRIAL';
                InCorrectTrialInd   ='INCORRECT_TRIAL';

                TrialType=AnalysisOpts.Trial.TrialType;

                %%% define Charactristics of trial that you are intrested in first

                %TrlFeature={'FIXATE_ACQUIRED','SAMPLE_ON','REWARD','CONDITION_NUM_OFFSET','SAMPLE_NUM_OFFSET','RULE_SWITCH_ON','RESLOC_NUM_OFFSET'};
                %TrlFeature={'FIXATE_ACQUIRED','SAMPLE_ON','CORRECT_TRIAL','CONDITION_NUM_OFFSET','SAMPLE_NUM_OFFSET','RULE_SWITCH_ON','RESLOC_NUM_OFFSET'};
                TrlFeature={'FIXATE_ACQUIRED','SAMPLE_ON','SACCADE_START','SACCADE_STOP','CORRECT_TRIAL','REWARD_START','REWARD_END','RULE_SWITCH_ON','RULE_SWITCH_OFF','RT_EyeTracker','CONDITION_NUM_OFFSET','SAMPLE_NUM_OFFSET','RESLOC_NUM_OFFSET'};
                FeatureFieldIden=[1,1,1,1,1,1,1,1,1,0,0,0,0]; % is this a field that is timing from start of the exp or not
                NFeatures=size(TrlFeature,2);
                m=1;

                for tri=1:size(TrialInfoMat,2)
                    temp=[];
                    if strcmp(TrialType,'ALL')   %all

                        if ( ~isnan(TrialInfoMat(tri).(StartEncodeInd)) & ~isnan(TrialInfoMat(tri).(StopEncodeInd))) & (~isnan(TrialInfoMat(tri).(InCorrectTrialInd)) | ~isnan(TrialInfoMat(tri).(CorrectTrialInd)))


                            for f=1:NFeatures
                                temp=[temp obj.ManData.ReplaceEmptywithNaN(TrialInfoMat(tri).(TrlFeature{f}))];
                            end

                            TrialTimes(m,:)=[TrialInfoMat(tri).(StartEncodeInd) TrialInfoMat(tri).(StopEncodeInd) temp];

                            m=m+1;
                        end

                    elseif strcmp(TrialType,'CORRECT')  % only correct

                        if  ~isnan(TrialInfoMat(tri).(StartEncodeInd)) & ~isnan(TrialInfoMat(tri).(StopEncodeInd)) & ~isnan(TrialInfoMat(tri).(CorrectTrialInd))

                            for f=1:NFeatures
                                temp=[temp obj.ManData.ReplaceEmptywithNaN(TrialInfoMat(tri).(TrlFeature{f}))];
                            end

                            TrialTimes(m,:)=[TrialInfoMat(tri).(StartEncodeInd) TrialInfoMat(tri).(StopEncodeInd) temp];
                            m=m+1;
                        end

                    elseif strcmp(TrialType,'INCORRECT')   %only incorrect

                        if ~isnan(TrialInfoMat(tri).(StartEncodeInd)) & ~isnan(TrialInfoMat(tri).(StopEncodeInd))  & ~isnan(TrialInfoMat(tri).(InCorrectTrialInd))

                            for f=1:NFeatures
                                temp=[temp obj.ManData.ReplaceEmptywithNaN(TrialInfoMat(tri).(TrlFeature{f}))];
                            end

                            TrialTimes(m,:)=[TrialInfoMat(tri).(StartEncodeInd) TrialInfoMat(tri).(StopEncodeInd) temp];
                            m=m+1;
                        end
                    end

                end

                AnalysisOpts.TrialTimesFields=[{StartEncodeInd,StopEncodeInd} TrlFeature];
                AnalysisOpts.FieldIden=logical([1 1 FeatureFieldIden]);
                %%% Determine if this field has any time
                %%% information(replaced by the code above)
                %                 NonTimeFields=strfind(AnalysisOpts.TrialTimesFields,'OFFSET');
                %                 for i=1:size(AnalysisOpts.TrialTimesFields,2)
                %                     if isempty(NonTimeFields{i})
                %                         AnalysisOpts.FieldIden(i)=1;
                %                     else
                %                         AnalysisOpts.FieldIden(i)=0;
                %                     end
                %                 end
            end
        end
        function CriticalTimes=ExtractCriticalTimes(~,TrialTimes)
            global AnalysisOpts
            CriticalTimes=[TrialTimes(:,AnalysisOpts.FieldIden)-TrialTimes(:,1) TrialTimes(:,~AnalysisOpts.FieldIden)];
        end
        function [ChannelInfo,ChannelArea]=FindWantedChsArea(obj)
            %%%  get the information about the area of each channel

            global AnalysisOpts
            if strcmpi(AnalysisOpts.Project,'Learning attentional templates')
                % load NS directory first
                load([AnalysisOpts.NS6DirPath 'NS6Directory_AC_corrected.mat']); % load NS6 data
                % now load the electode information of rec day too % not
                % looking at delayed succade data for now
                RawDataFile=[AnalysisOpts.RecDatePath 'LearnAttTemplate_' AnalysisOpts.Animal '_' AnalysisOpts.RecDate '002ns3.mat' ];
                if ~exist(RawDataFile,'file')
                    ReadDataFile(AnalysisOpts.RecDatePath,AnalysisOpts.Animal,AnalysisOpts.RecDate,'2','ns3')
                end
                load(RawDataFile,'MetaTags');

                FileNameTarg=['LearnAttTemplate_' AnalysisOpts.Animal '_' AnalysisOpts.RecDate '002.ns6'];
                % now find this file info in NS6
                FileInd=arrayfun(@(x) strcmpi(ns6directory(x).FileName,FileNameTarg),1:length(ns6directory));
                FileInfo=ns6directory(FileInd); % information about this recording day

                r=1; %% loop in the wanted areas
                for Ar=AnalysisOpts.Area_2look
                    ChDig=logical(zeros(1,256));
                    Ar=Ar{1};
                    ChannelInfo(r).AreaName= Ar;
                    if isfield(FileInfo,Ar)
                        ChDig(FileInfo.(Ar))=1;
                        ChannelInfo(r).AreaChs = ChDig;
                        ChannelInfo(r).ToUseChs= ChannelInfo(r).AreaChs;
                        ChannelInfo(r).ChannelID=MetaTags.ChannelID; %[arrayfun(@(x) find(MetaTags.ChannelID==x),find(ChDig))]; % correspoding channel in the main rec file
                    else
                        ChannelInfo(r).AreaChs=ChDig;ChannelInfo(r).ToUseChs=ChDig;ChannelInfo(r).ChannelID=MetaTags.ChannelID;
                    end
                    r=r+1;
                end
                %come up with a matrix
                A=arrayfun(@(x) x*ChannelInfo(x).AreaChs',1:length(AnalysisOpts.Area_2look),'UniformOutput',0);
                ChannelArea=sum(cell2mat(A),2);

            elseif strcmpi(AnalysisOpts.Project,'Rule Representation')
                %%% the channels that we are sure there is a signal in there by looking at
                %%% manual check variable "Check'

                %               SummeryFile=['Rec' AnalysisOpts.RecDate '_TH_' num2str(AnalysisOpts.NeuralAnalysis.StdTH2look) '_SummeryData.mat'];
                if str2num(AnalysisOpts.RecDate(end-1:end))==18 %if we are processing Silas data 
                    SummeryFile=['Rec' AnalysisOpts.RecDate '_SummeryData_ManuSpkSrt.mat'];
                elseif str2num(AnalysisOpts.RecDate(end-1:end))==21 % if we are processing Chico data 
                    SummeryFile=['Rec' AnalysisOpts.RecDate '_SummeryDataV2_ManuSpkSrt.mat'];
                end

               
                % for check (1:Sing Unit/2:Multi/3:Low Quality/4:Not Fulltime/0:Not Useable)
                load([AnalysisOpts.RecDatePath SummeryFile],'AllChannelsSpec','Check')
                r=1; %% loop in the wanted areas
                for Ar=AnalysisOpts.Area_2look
                    Ar=Ar{1};
                    ChannelInfo(r).AreaName=Ar;
                    ChannelInfo(r).AreaChs=arrayfun(@(x) ~isempty(strfind(lower(AllChannelsSpec(x).AreaName),lower(Ar))),1:length(Check)); %% finds this areas Channels
                    %% (1:Sing Unit/2:Multi/3:Low Quality/4:Not Fulltime/0:Not Useable)
                    % find valid single neuron channels
                    ChannelInfo(r).ValidChs=arrayfun(@(x) sum(cell2mat(Check(x).Valid)==1),1:length(Check));
                    ChannelInfo(r).ValidChsClustNum=arrayfun(@(x) find(cell2mat(Check(x).Valid)==1),1:length(Check),'UniformOutput',0);
                    % find valid multi units
                    ChannelInfo(r).ValidMultiChs=arrayfun(@(x) sum(cell2mat(Check(x).Valid)==2),1:length(Check));
                    ChannelInfo(r).ValidMultiChsClustNum=arrayfun(@(x) find(cell2mat(Check(x).Valid)==2),1:length(Check),'UniformOutput',0);
                    % find lowquality multiunits
                    ChannelInfo(r).ValidLowQuaMultiChs=arrayfun(@(x) sum(cell2mat(Check(x).Valid)==3),1:length(Check));
                    ChannelInfo(r).ValidLowQuaMultiChsClustNum=arrayfun(@(x) find(cell2mat(Check(x).Valid)==3),1:length(Check),'UniformOutput',0);
                    % find not full time units
                    ChannelInfo(r).ValidNotFulTimChs=arrayfun(@(x) sum(cell2mat(Check(x).Valid)==4),1:length(Check));
                    ChannelInfo(r).ValidNotFulTimChsClustNum=arrayfun(@(x) find(cell2mat(Check(x).Valid)==4),1:length(Check),'UniformOutput',0);
                    % usable channels
                    ChannelInfo(r).ToUseChs=ChannelInfo(r).ValidChs & ChannelInfo(r).AreaChs;
                    r=r+1;
                end
                %come up with a matrix
                A=arrayfun(@(x) x*ChannelInfo(x).AreaChs',1:length(AnalysisOpts.Area_2look),'UniformOutput',0);
                ChannelArea=sum(cell2mat(A),2);
            end
        end
        function ChsSet=FindWantedChSet(obj,ChannelInfo,varargin) % finds the channels specified in the AnalysisOpts
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs

            if isempty(AnalysisOpts.Ch_2look) & isempty(AnalysisOpts.File_2look)
                ChsSet=arrayfun(@(x)  find(ChannelInfo(x).ToUseChs),1:length(AnalysisOpts.Area_2look),'UniformOutput',0);  %%% loop thourgh wanted Channels
                ChsSet=cell2mat(ChsSet);
            elseif ~isempty(AnalysisOpts.File_2look)
                [~,~,loadFile2lookName]=GenerateFileName(FS,AnalysisOpts.DataSavePath,AnalysisOpts.AnalysisPathName,[],[],[],[],'SelfName',1,'SelfNameTxt',AnalysisOpts.File_2look);
                load(loadFile2lookName,'AnalysisData')  % replace analysis data with this
                ChsSet=AnalysisData.Ch;
            else
                % ChsSet=arrayfun(@(x)  find(ChannelInfo(x).ToUseChs),1:length(AnalysisOpts.Area_2look),'UniformOutput',0);  %%% loop thourgh wanted Channels
                % ChsSet=cell2mat(ChsSet);
                ChsSet=1:112;%AnalysisOpts.Ch_2look;
            end
        end
        function [ChNumAll,ClustNumAll,AreaIndAll,QualIndAll]=FindWantedSpkChSet(obj,SpkSrtQual,Area,varargin)
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            [ChannelInfo]=obj.FindWantedChsArea();
            if iscell(Area) | ischar(Area)
                AreaInd=cellfun(@(y) find(arrayfun(@(x) contains(ChannelInfo(x).AreaName,y),1:length(ChannelInfo))),Area);
            else
                AreaInd=Area;
            end
            ChNumAll=[];ClustNumAll=[];AreaIndAll=[];QualIndAll=[];
            for ArInd=AreaInd
                for QualInd=SpkSrtQual
                    switch QualInd
                        case 1
                            ChNum=cell2mat(arrayfun(@(x) x*ones(1,ChannelInfo(ArInd).ValidChs(x)),find(ChannelInfo(ArInd).ValidChs & ChannelInfo(ArInd).AreaChs),'UniformOutput',false));
                            ClustNum=cell2mat(arrayfun(@(x) cell2mat(ChannelInfo(ArInd).ValidChsClustNum(x)),find(ChannelInfo(ArInd).ValidChs & ChannelInfo(ArInd).AreaChs),'UniformOutput',false));
                        case 2
                            ChNum=cell2mat(arrayfun(@(x) x*ones(1,ChannelInfo(ArInd).ValidMultiChs(x)),find(ChannelInfo(ArInd).ValidMultiChs & ChannelInfo(ArInd).AreaChs),'UniformOutput',false));
                            ClustNum=cell2mat(arrayfun(@(x) cell2mat(ChannelInfo(ArInd).ValidMultiChsClustNum(x)),find(ChannelInfo(ArInd).ValidMultiChs & ChannelInfo(ArInd).AreaChs),'UniformOutput',false));
                        case 3
                            ChNum=cell2mat(arrayfun(@(x) x*ones(1,ChannelInfo(ArInd).ValidLowQuaMultiChs(x)),find(ChannelInfo(ArInd).ValidLowQuaMultiChs & ChannelInfo(ArInd).AreaChs),'UniformOutput',false));
                            ClustNum=cell2mat(arrayfun(@(x) cell2mat(ChannelInfo(ArInd).ValidLowQuaMultiChsClustNum(x)),find(ChannelInfo(ArInd).ValidChs & ChannelInfo(ArInd).AreaChs),'UniformOutput',false));
                        case 4
                            ChNum=cell2mat(arrayfun(@(x) x*ones(1,ChannelInfo(ArInd).ValidNotFulTimChs(x)),find(ChannelInfo(ArInd).ValidNotFulTimChs & ChannelInfo(ArInd).AreaChs),'UniformOutput',false));
                            ClustNum=cell2mat(arrayfun(@(x) cell2mat(ChannelInfo(ArInd).ValidNotFulTimsClustNum(x)),find(ChannelInfo(ArInd).ValidChs & ChannelInfo(ArInd).AreaChs),'UniformOutput',false));
                    end
                    ChNumAll=cat(2,ChNumAll,ChNum);
                    ClustNumAll=cat(2,ClustNumAll,ClustNum);
                    AreaIndAll=cat(2,AreaIndAll,ones(1,length(ChNum))*ArInd);
                    QualIndAll=cat(2,QualIndAll,ones(1,length(ChNum))*QualInd);
                end
            end
        end
        function [StrTim,StpTim]=CalStrStpRec(obj,BlockSpec,RawData,Fs,varargin) % calculates Start and stop of recording
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs

            if ~isempty(AnalysisOpts.PreProcess.MaxNTrial) % trial based analysis
                % grab the whole trial data % havent set up loop analysis yet
                % %%%%%%%%
                [~,MinOrderBlk]=min(BlockSpec.ThisBlkOrder);
                [~,MaxOrderBlk]=max(BlockSpec.ThisBlkOrder);
                StrTim=BlockSpec.ThisBlkTrialTimes{MinOrderBlk}(1,strcmp(AnalysisOpts.TrialTimesFields,'START_TRIAL'));
                StpTim=BlockSpec.ThisBlkTrialTimes{MaxOrderBlk}(end,strcmp(AnalysisOpts.TrialTimesFields,'END_TRIAL'));
            elseif ~isempty(AnalysisOpts.PreProcess.MaxTime) % Time based analysis
                % grab time data
                if strcmpi(AnalysisOpts.PreProcess.MaxTime,'ALL')
                    Timetot=size(RawData.data,2)/Fs; % how many secs of data
                    AnalysisOpts.PreProcess.StrTime=0; % reset start time
                else
                    Timetot= AnalysisOpts.PreProcess.MaxTime;
                end
                StrTim=AnalysisOpts.PreProcess.StrTime;StpTim=AnalysisOpts.PreProcess.StrTime+Timetot; % in secc
            end
        end
        function AreaNum=FindChArea(~,Ch,ChannelInfo) % finds the area of this channel
            AreaNum=find(arrayfun(@(x) sum(ChannelInfo(x).AreaChs(Ch)==1),1:length(ChannelInfo))); % find which area this channel belongs
        end
        function [TrialTimes,RuleBlockTrials,ChannelInfo,ChannelArea,ChsSet,BlockSpec,StimMap]=InitializeTrialFuncs(obj,varargin)  % run the trial functions that are common to all of the fucntions

            global AnalysisOpts AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            CurrentAnimal=AnalysisOpts.Animal;
            CurrentRecDate=AnalysisOpts.RecDate;
            if strcmpi(AnalysisOpts.RecDate,'ALL') | isempty(CurrentAnimal) | isempty(CurrentRecDate)
                % AnalysisOpts.Animal=AnalysisOpts.AnimalSet{1};
                %  AnalysisOpts.RecDate=AnalysisOpts.DateSet{1};
                [TrialTimes]=[];[RuleBlockTrials]=[]; [ChannelInfo]=[];ChannelArea=[];ChsSet=[];BlockSpec=[];
            else
                [TrialTimes]=obj.FindWantedTrialinArea;   % retreive channel information of this rec
                [RuleBlockTrials]=obj.SegregateBlocks(TrialTimes); % get the sequence of trials for each block
                [ChannelInfo,ChannelArea]=obj.FindWantedChsArea;  % find channels we want from this recording
                ChsSet=obj.FindWantedChSet(ChannelInfo); % find desired channels
                BlockSpec=obj.GetBlockTrialInfo(RuleBlockTrials,TrialTimes,'MaxNTrial',AnalysisOpts.Bhv.MaxNTrial,'NTrial2Swtch',AnalysisOpts.Bhv.NTrial2Swtch,'Nblk_2Look','MAX'); % get blockspecs
            end
            AnalysisOpts.Animal=CurrentAnimal;
            AnalysisOpts.RecDate=CurrentRecDate;
            
            %  BlockSpec.BhvModelInfo=obj.GetBhvModelInfo(RuleBlockTrials,TrialTimes,'MaxNTrial',AnalysisOpts.Bhv.MaxNTrial,'NTrial2Swtch',AnalysisOpts.Bhv.NTrial2Swtch);
            %             else
            %                 TrialTimes=[];RuleBlockTrials=[];ChannelInfo=[];ChannelArea=[];ChsSet=[];BlockSpec=[];StimMap=[];
            %             end
            StimMapFile=[AnalysisOpts.DataPath 'StimulusMapping.mat'  ];
            if exist(StimMapFile,'file')
                load(StimMapFile,'StimMap');
            else
                [StimMap.cuemapshape,StimMap.cue,StimMap.FieldNames]=obj.GenerateStimulusMapping;
                save(StimMapFile,'StimMap');
            end
            % copy some of this data into AnalysisData
            obj.ManData.CopyVars2AnalysisData('TrialTimes',TrialTimes,'ChannelInfo',ChannelInfo,'ChannelArea',ChannelArea,'ChsSet',ChsSet,'StimMap',StimMap);
            AnalysisData.Ch=ChannelArea~=0;
            AnalysisData.ChArea=ChannelArea; % get the area of each Channel
            AnalysisOpts.StimMap=StimMap;          
        end
        function ParseChannelInformation(obj,Chns,RecSpecs,varargin) % parses all of the channel information and their recspecs into a coherent matrix for further use
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            k=1;
            Nrecs=size(Chns,1);
            for rec=1:Nrecs
                Nchs=length(Chns{rec,1});
                for nch=1:Nchs
                    %fill in the information about each channel
                    AllChannelsInfo(k).ChNum=Chns{rec,1}(nch);
                    AllChannelsInfo(k).ChClust=Chns{rec,2}(nch);
                    AllChannelsInfo(k).SpkSrtQual=Chns{rec,4}(nch);
                    AllChannelsInfo(k).ChAreaNum=Chns{rec,3}(nch);
                    AllChannelsInfo(k).ChAreaName=AnalysisOpts.AreaNames{Chns{rec,3}(nch)};
                    AllChannelsInfo(k).RecDate=RecSpecs(rec).RecDate;
                    AllChannelsInfo(k).Animal=RecSpecs(rec).Animal;
                    AllChannelsInfo(k).RecChNum=nch; % which channel number in this recording
                    k=k+1;
                end
            end
            % now if you have specifiv channels you want to look at only
            % take those from all of these channels 
            if ~isempty(AnalysisOpts.Ch_2look);AllChannelsInfo=AllChannelsInfo(AnalysisOpts.Ch_2look);end
            % update informations about channels
            AnalysisOpts.AllChannelsInfo=AllChannelsInfo;
            AnalysisOpts.Ch_2look_ChNum=AnalysisOpts.Ch_2look;
            AnalysisOpts.Ch_2look=[AllChannelsInfo.ChNum]; % update channels to look
            AnalysisOpts.nCh_2look=[1:length(AnalysisOpts.AllChannelsInfo)];
            if isempty(AnalysisOpts.Ch_2look_ChNum);AnalysisOpts.Ch_2look_ChNum=AnalysisOpts.nCh_2look;end
            AnalysisOpts.Ch_2look_Clust=[AllChannelsInfo.ChClust];
            AnalysisOpts.Ch_2look_AreaNum=[AllChannelsInfo.ChAreaNum];
            AnalysisOpts.Ch_2look_AreaName=arrayfun(@(x) AllChannelsInfo(x).ChAreaName,AnalysisOpts.nCh_2look,'UniformOutput',0);
            AnalysisOpts.Ch_2look_RecDate=arrayfun(@(x) AllChannelsInfo(x).RecDate,AnalysisOpts.nCh_2look,'UniformOutput',0);
            AnalysisOpts.Ch_2look_Animal=arrayfun(@(x) AllChannelsInfo(x).Animal,AnalysisOpts.nCh_2look,'UniformOutput',0);
            AnalysisOpts.Ch_2look_RecChNum=[AllChannelsInfo.RecChNum];
            % update Current channel to the last channel
            obj.UpdateCurrentCh(AnalysisOpts.nCh_2look(end));

            AnalysisOpts.FieldNames=fieldnames(AnalysisOpts);
            Ch2LookInd=contains(AnalysisOpts.FieldNames,'Ch_2look');
            for Field=AnalysisOpts.FieldNames(Ch2LookInd)'
                AnalysisOpts.Ch2lookRep.(Field{1})=AnalysisOpts.(Field{1});
            end
        end
        function ParseChannelInformationTrialSpikeData(obj,TrialSpikeData,varargin) % parses all of the channel information and their recspecs into a coherent matrix for further use using TrialSpikeData instead of RecSpecs
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            k=1;
            Nchs=length(TrialSpikeData{1});            
            for nch=1:Nchs
                %fill in the information about each channel
                AllChannelsInfo(k).ChNum=TrialSpikeData{1}(nch).Ch;
                AllChannelsInfo(k).ChClust=TrialSpikeData{1}(nch).Cluster;
                AllChannelsInfo(k).SpkSrtQual=TrialSpikeData{1}(nch).Specs.SpkQuality;
                AllChannelsInfo(k).ChAreaNum=TrialSpikeData{1}(nch).ChannelArea;
                AllChannelsInfo(k).ChAreaName=AnalysisOpts.AreaNames{TrialSpikeData{1}(nch).ChannelArea};
                AllChannelsInfo(k).RecDate=TrialSpikeData{1}(nch).RecDate;
                AllChannelsInfo(k).Animal=obj.ManData.DetermineDateAnimal(TrialSpikeData{1}(nch).RecDate);
                k=k+1;
            end
            % update informations about channels
            AnalysisOpts.AllChannelsInfo=AllChannelsInfo;
            AnalysisOpts.Ch_2look_ChNum=AnalysisOpts.Ch_2look;
            AnalysisOpts.Ch_2look=[AllChannelsInfo.ChNum]; % update channels to look
            AnalysisOpts.nCh_2look=[1:length(AnalysisOpts.AllChannelsInfo)];
            AnalysisOpts.Ch_2look_Clust=[AllChannelsInfo.ChClust];
            AnalysisOpts.Ch_2look_AreaNum=[AllChannelsInfo.ChAreaNum];
            AnalysisOpts.Ch_2look_AreaName=arrayfun(@(x) AllChannelsInfo(x).ChAreaName,AnalysisOpts.nCh_2look,'UniformOutput',0);
            AnalysisOpts.Ch_2look_RecDate=arrayfun(@(x) AllChannelsInfo(x).RecDate,AnalysisOpts.nCh_2look,'UniformOutput',0);
            AnalysisOpts.Ch_2look_Animal=arrayfun(@(x) AllChannelsInfo(x).Animal,AnalysisOpts.nCh_2look,'UniformOutput',0);
            % if we have not specified any channel that means we want all
            % of the channels 
            if isempty(AnalysisOpts.Ch_2look_ChNum);AnalysisOpts.Ch_2look_ChNum=AnalysisOpts.nCh_2look;end
            % update Current channel to the last channel
            obj.UpdateCurrentCh(AnalysisOpts.nCh_2look(end));

            AnalysisOpts.FieldNames=fieldnames(AnalysisOpts);
            Ch2LookInd=contains(AnalysisOpts.FieldNames,'Ch_2look');
            for Field=AnalysisOpts.FieldNames(Ch2LookInd)'
                AnalysisOpts.Ch2lookRep.(Field{1})=AnalysisOpts.(Field{1});
            end
        end

        function UpdateCurrentCh(~,ChNum) % updates the current channel
            global AnalysisOpts
            if length(ChNum)>1 | isempty(ChNum)
                AnalysisOpts.CurrentCh='';
                AnalysisOpts.CurrentChClust='';
                AnalysisOpts.CurrentCh_AreaNum='';
                AnalysisOpts.CurrentCh_AreaName= '';
                AnalysisOpts.CurrentCh_RecDate='';
                AnalysisOpts.CurrentCh_Animal='';
                AnalysisOpts.CurrentCh_ChNum='';
            else
                AnalysisOpts.CurrentCh=AnalysisOpts.Ch_2look(ChNum);
                AnalysisOpts.CurrentChClust=AnalysisOpts.Ch_2look_Clust(ChNum);
                AnalysisOpts.CurrentCh_AreaNum= AnalysisOpts.Ch_2look_AreaNum(ChNum);
                AnalysisOpts.CurrentCh_AreaName= AnalysisOpts.Ch_2look_AreaName{ChNum};
                AnalysisOpts.CurrentCh_RecDate= AnalysisOpts.Ch_2look_RecDate{ChNum};
                AnalysisOpts.CurrentCh_Animal= AnalysisOpts.Ch_2look_Animal{ChNum};
                AnalysisOpts.CurrentCh_ChNum=AnalysisOpts.Ch_2look_ChNum(ChNum);              
            end

        end
        function UpdateCh_2look(~,ChNum) %removes channels from list of channels
            global AnalysisOpts
            AnalysisOpts.Ch_2look=AnalysisOpts.Ch_2look(ChNum);
            AnalysisOpts.Ch_2look_Clust=AnalysisOpts.Ch_2look_Clust(ChNum);
            AnalysisOpts.Ch_2look_AreaNum= AnalysisOpts.Ch_2look_AreaNum(ChNum);
            AnalysisOpts.Ch_2look_AreaName= AnalysisOpts.Ch_2look_AreaName(ChNum);
            AnalysisOpts.Ch_2look_RecDate= AnalysisOpts.Ch_2look_RecDate(ChNum);
            AnalysisOpts.Ch_2look_Animal= AnalysisOpts.Ch_2look_Animal(ChNum);
            AnalysisOpts.nCh_2look=1:length(AnalysisOpts.Ch_2look);
            AnalysisOpts.Ch_2look_ChNum=AnalysisOpts.Ch_2look_ChNum(ChNum); 
        end
        function RevertCh_2look(~) % reverts channel to look into it's original state
            global AnalysisOpts

            FieldNames=fieldnames(AnalysisOpts.Ch2lookRep)';
            for Field=FieldNames
                AnalysisOpts.(Field{1})=AnalysisOpts.Ch2lookRep.(Field{1});
            end
        end
        function [ChNum,AreaNames]=GrabNeuronswithArea(~,AreaNum)
            global AnalysisOpts
            if sum(AreaNum)==0 | isnan(AreaNum) | isempty(AreaNum)
                ChNum=[];AreaNames='';return;
            end

            ChNum=cell2mat(arrayfun(@(x) find(AnalysisOpts.Ch_2look_AreaNum==x),AreaNum,'UniformOutput',0));
            AreaNames=AnalysisOpts.AreaNames(AreaNum);
            if sum(AreaNum==1:5)==5;AreaNames='All';end
        end
        function [Chns,PSTHRaster,TrialSpikeData,ChAreas,RecSpecs]=GetSpecificChData(obj,DataFileName)
            global AnalysisOpts
            
            load(DataFileName,'Chns','RecSpecs','TrialSpikeData',['PSTHRaster_' num2str(AnalysisOpts.PopulationAna.PSTHbin) 'ms']);           
            if ~isempty(AnalysisOpts.Ch_2look)
                ChNum=AnalysisOpts.Ch_2look;
                ChNumAll=cell2mat(cellfun(@(x) x,Chns(:,1)','UniformOutput',0));
                ClustNumAll=cell2mat(cellfun(@(x) x,Chns(:,2)','UniformOutput',0));
                AreaIndAll=cell2mat(cellfun(@(x) x,Chns(:,3)','UniformOutput',0));
                QualIndAll=cell2mat(cellfun(@(x) x,Chns(:,4)','UniformOutput',0));
                Chns=[{ChNumAll(ChNum)} {ClustNumAll(ChNum)} {AreaIndAll(ChNum)} {QualIndAll(ChNum)}];
                eval(['PSTHRaster=PSTHRaster_' num2str(AnalysisOpts.PopulationAna.PSTHbin) 'ms(ChNum);']);
                TrialSpikeData=[{TrialSpikeData{1}(ChNum)} {TrialSpikeData{2}(ChNum)}];             
            else 
                eval(['PSTHRaster=PSTHRaster_' num2str(AnalysisOpts.PopulationAna.PSTHbin) 'ms;']);
            end
            ChAreas=unique(cell2mat(Chns(:,3)')); % find the areas of the channels we have
        end
        function [ThisDateNum,RecChns,ThisRecChstxt,NChs]=GetRecordingChannels(obj,DateNum,AreaNum,ChNum,UpdataChInfoFlag)
            % @UpdataChInfoFlag are we also updating our Ch2Look
            % information matrix with these channels
            global AnalysisOpts
            
            AreaName=AnalysisOpts.AreaNames(AreaNum);
            load('RecNums.mat','SilasRecs','ChicoRecs');
            if strcmp(DateNum,'Silas')
                ThisDateNum=SilasRecs;
                AnalysisOpts.Animal='Silas';
            elseif strcmp(DateNum,'Chico')
                ThisDateNum=ChicoRecs;
                AnalysisOpts.Animal='Chico';
            elseif strcmp(DateNum,'ALL')
                ThisDateNum=0;
                AnalysisOpts.Animal='ALL';
            else
                ThisDateNum=DateNum;
            end
            % get the channels for this recording
            DataFileName=[AnalysisOpts.DataSavePath 'Core Data' filesep 'SpikingData' filesep AnalysisOpts.Animal '_' AnalysisOpts.AreaNames{AreaNum} '_' AnalysisOpts.TrlSpkTimeFieldName '_' AnalysisOpts.SpkCntStartFieldName '.mat'];
            if  ~sum(strcmp(DateNum,{'Silas','Chico','ALL'}))
                [Chns,RecSpecs]=obj.ManData.RunFuncOnRec(TrialFunc,'FindWantedSpkChSet',ThisDateNum,{[AnalysisOpts.SpikeQuality2Look],AreaName},3);
                RecChns=cell2mat(Chns(:));
            elseif sum(strcmp(DateNum,{'Silas','Chico','ALL'})) % load channels from these rec
                if UpdataChInfoFlag
                   load(DataFileName,'RecSpecs');
                end
                load(DataFileName,'Chns');
                ChNumAll=cell2mat(cellfun(@(x) x,Chns(:,1)','UniformOutput',0));
                RecChns=1:length(ChNumAll);                                
            end
            if ~isempty(ChNum)
                RecChns=RecChns(ChNum);
            end
            ThisRecChstxt=obj.ManData.ConvMat2Char(RecChns);
            NChs=length(RecChns);
            if UpdataChInfoFlag % we are updating Ch2Look info as well
                AnalysisOpts.Ch_2look=[]; % initialize Ch_2Look
                obj.ParseChannelInformation(Chns,RecSpecs);
            end
        end

    end

end
