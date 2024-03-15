
classdef HelperFuncs
    %TOOLBOXHELPERS Summary of this class goes here
    %   Detailed explanation goes here
    % Sina Tafazoli
    properties
        UseOldData=0; % if using old data with 3hz freq spacing(this is redundanct and she be removed later)
    end
    properties (Access=private)
        ManData=ManipulateData;
    end
    methods

        function  obj=HelperFuncs(varargin)
            if nargin~=0 % initialize vars
                obj=obj.ParseParams(varargin) ; %%Process optional inputs
            end
        end
        function  obj=ParseParams(obj,InputArgs)
            %Process optional inputs
            if mod(length(InputArgs), 2) ~= 0, error('Must pass key/value pairs for options.'); end
            for i = 1:2:length(InputArgs)
                try
                    obj.(InputArgs{i}) = InputArgs{i+1};
                catch
                    error('Couldn''t set option ''%s''.', InputArgs{2*i-1});
                end
            end
        end  % parse parameters
        %% data processing
        function   [DataSpace,ParamSpace,BlockSpec]=getEphysData(obj,Ch,varargin)
            % loads preprocessed Ephysdata
            global AnalysisOpts AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            fprintf('\nloading wavelet data from channel %i',Ch)
            DataFileName=[AnalysisOpts.DataSavePath AnalysisOpts.FS 'Wavelet_Data' AnalysisOpts.FS 'Rec' AnalysisOpts.RecDate 'WaveletData_Channel_' num2str(Ch) 'Width' ...
                num2str(AnalysisOpts.MotifAnalysis.WaveletWidth) '.mat'];
            load(DataFileName,'DataSpace','ParamSpace','BlockSpec');

            if strcmpi(AnalysisOpts.MotifAnalysis.CwtFreqScale,'LinScale')
                AnalysisData.cwt=DataSpace.Wavelet_Linear;
                AnalysisData.cwt_f=DataSpace.f_Linear';
            elseif strcmpi(AnalysisOpts.MotifAnalysis.CwtFreqScale,'LogScale')
                AnalysisData.cwt=DataSpace.Wavelet_Log;
                AnalysisData.cwt_f=DataSpace.f_Log';
            end
            AnalysisData.RawData=DataSpace.RawData; % get the raw data
            % if ~isempty(BlockSpec)
            %    AnalysisData.TrialTimes=BlockSpec.ThisBlkTrialTimes{1}-BlockSpec.ThisBlkTrialTimes{1}(1,strcmp(AnalysisOpts.TrialTimesFields,'START_TRIAL'));
            % else
            AnalysisData.TrialTimes=[];
            % end
            AnalysisOpts.NeuralAnalysis.SamplingFrequency=ParamSpace.Fs;
            %  AnalysisOpts.LFPParams.FilterOpts.TargetLFPSamplingFrequency=ParamSpace.FsLFP;
            AnalysisOpts.NeuralAnalysis.FsWave=ParamSpace.FsWave;
        end
        function   [RawEphysData,Fs,DataPointer]=getRawEphysData(obj,Ch,UseDataPointer,varargin)
            % loads preprocessed Ephysdata( only output the RawData)
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            DataPointer=cell(1,length(Ch));
            ExtraTxt=['Width' num2str(AnalysisOpts.MotifAnalysis.WaveletWidth)];
            for ch=Ch % loop on chs
                fprintf('\nloading raw data from channel %i',ch)
                AnalysisOpts.CurrentCh=ch; % change the current channel tho this one now
                if ~UseDataPointer
                    % load(DataFileName,'RawData');
                    RawEphysData{find(Ch==ch)}=obj.ManData.LoadVar('Wavelet_Data','RawData',ExtraTxt,UseDataPointer); % get the raw data
                    DataPointer=[];
                else
                    DataPointer{find(Ch==ch)}=obj.ManData.LoadVar('Wavelet_Data','RawData',ExtraTxt,UseDataPointer);
                    RawEphysData=[];
                end
                if ch==Ch(1)
                    ParamSpace=obj.ManData.LoadVar('Wavelet_Data','ParamSpace',ExtraTxt,UseDataPointer); % get the raw data
                    Fs=ParamSpace.FsLFP;
                end

            end
        end
        function   [WaveletData,f_Linear,FsWave,DataPointer]=getWaveletData(obj,Ch,UseDataPointer,varargin)
            % loads wavelet data from file
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            DataPointer=cell(1,length(Ch));
            if obj.UseOldData
                ExtraTxt=['Width' num2str(AnalysisOpts.MotifAnalysis.WaveletWidth) 'old'];
            else
                ExtraTxt=['Width' num2str(AnalysisOpts.MotifAnalysis.WaveletWidth)];
            end
            for ch=Ch % loop on chs
                fprintf('\nloading wavelet data from channel %i',ch)
                AnalysisOpts.CurrentCh=ch; % change the current channel tho this one now
                if ~UseDataPointer
                    WaveletData{find(Ch==ch)}=obj.ManData.LoadVar('Wavelet_Data','Wavelet_Linear',ExtraTxt,UseDataPointer);
                    f_Linear=obj.ManData.LoadVar('Wavelet_Data','f_Linear',ExtraTxt,UseDataPointer);
                    ParamSpace=obj.ManData.LoadVar('Wavelet_Data','ParamSpace',ExtraTxt,UseDataPointer);
                    FsWave=ParamSpace.FsWave;
                    DataPointer=[];
                else
                    DataPointer{find(Ch==ch)}=obj.ManData.LoadVar('Wavelet_Data','Wavelet_Linear',ExtraTxt,UseDataPointer);
                    RawEphysData=[];
                    WaveletData=[];
                    f_Linear=[];
                    FsWave=[];
                end
            end
        end
        function [SpikeData]=getSpikeData(obj,Chs,Cluster,ChannelArea,varargin)
            % loads  spike data for Chs we want
            global AnalysisOpts 
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            if isempty(Chs)
                SpikeData=[];
                return;
            end
           % GenerateClusterPath(AnalysisOpts.AnalysisType) % update the path we need for this recording
            for i=1:length(Chs) % loop on channels
                Ch=Chs(i);
                Clust=Cluster(i);
                % find the corresponding channel and fill the in the
                % threshold
                fstruct = dir([AnalysisOpts.RecDatePath 'Neuron_' num2str(Ch) filesep 'TrialData*ManuSpkSrt.mat']);
                if ~isempty(fstruct)
                    FileSpecs=SparseFileName([fstruct.folder filesep fstruct.name]);
                    StdTH2look=FileSpecs.Th;
                    SpikeDataFile=[fstruct(1).folder filesep 'TrialData_Rec_' AnalysisOpts.RecDate '_Channel_' num2str(Ch) '_Cluster_' num2str(Clust) '_TH' num2str(StdTH2look) AnalysisOpts.ManualSpkSortingTxt '.mat'];
                       % get to the proper file and load the variables
                       % inside( regenerating critical times from trial
                       % times again)
                  %  SpikeData(i).CriticalTimes    =obj.ManData.LoadVar('','CriticalTimes',SpikeDataFile,0,'SelfName',1); 
                    SpikeData(i).SpikeShape       =mean(obj.ManData.LoadVar('','SpikeShape',SpikeDataFile,0,'SelfName',1),1);% taking the mean to save the space
                    SpikeData(i).TrialSpikeTime   =obj.ManData.LoadVar('','TrialSpikeTime',SpikeDataFile,0,'SelfName',1);
                    SpikeData(i).Specs            =obj.ManData.LoadVar('','Specs',SpikeDataFile,0,'SelfName',1);
                    SpikeData(i).SpikeTime        =obj.ManData.LoadVar('','SpikeTime',SpikeDataFile,0,'SelfName',1);
                    % fill in all of the info about this channel
                    SpikeData(i).Ch=Ch;
                    SpikeData(i).ChannelArea=ChannelArea(i);
                    SpikeData(i).RecDate    =AnalysisOpts.RecDate;
                    SpikeData(i).Specs.FileName   =SpikeDataFile;
                    SpikeData(i).Cluster=Clust;
                     fprintf('\nReading Spikes from Rec %s Channel %i Cluster %i ...',AnalysisOpts.RecDate,Ch,Clust);
                else
                    fprintf('\nSpike data for Rec %s Channel %i Cluster %i does not exist...',AnalysisOpts.RecDate,Ch,Clust);
                    SpikeData{i}=[];
                end
            end
        end
        %% AUX and cluster funcs
        function [TimeFreq,FigPrms,ManData,TrialFunc,FilterFunc,NeuAnaFunc,CrosFreqCop]=LoadAnalysisClasses(~)% loads all of the classes we have so far
            global AnalysisOpts 
            TimeFreq=TimeFreqAnalysis; % get the class of Time FreqAnalysis
            FigPrms=fig_params;        % fig param calss
            ManData=ManipulateData;    % manipulate data class
            TrialFunc=TrialFuncs;      % getting trials from data (this is private to each project)
            FilterFunc=FilterFuncs;    % Filtering functions 
            NeuAnaFunc=NeuralAnalysisFuncs;% Neural Analysis functions
            CrosFreqCop=CrossFreqCopling; % Cross Frequency Coupling Analysis 
         end
        function [Fs,FsLFP,FsWave,FsWaveTarg]=SetupGeneralFuncs(obj,varargin)
            global AnalysisOpts 
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % add a string at the end of the file in case we are
            % rereferenceing
            AnalysisOpts.ExtraStr=AnalysisOpts.ReRefChannelsStr{AnalysisOpts.ReRefChannels+1};
            %%manage plots
            if ~AnalysisOpts.ShowPlot
                set(0,'DefaultFigureVisible','off');   % don't show any plot if we don't want it
            else
                set(0,'DefaultFigureVisible','on');   % don't show any plot if we don't want it
            end

            % load frequencies into single var so it is easier to use them
            Fs=AnalysisOpts.NeuralAnalysis.SamplingFrequency;   % take Fs so it's a smaller string
            FsLFP=AnalysisOpts.MotifAnalysis.FsLFP ;% take FsLFP so it's a smaller string
            FsWave=AnalysisOpts.MotifAnalysis.FsWave; %  frequency of Wavelet
            FsWaveTarg=AnalysisOpts.MotifAnalysis.FsWaveTarg; % target frequency for Wavelet
        end

    end
end

