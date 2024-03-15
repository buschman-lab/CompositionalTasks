% a prototype fucntion to run functions we want 
function RunFullFunconClust(RunonCluster,DateNum,ChNum,AreaNum,PairNum,varargin)
global AnalysisOpts
% load clusters
ClusterFunc=ClusterFuncs;
SetAnalysisOptions(varargin) % set all of the analysis functions
SpockPass='Poqttirwpir2084';
cd('Z:\Users\Sina\Rule Representation Project\Server Analysis forlders\Cluster Analsis General Functions\');
for ThisDateNum=DateNum
    
    if ThisDateNum==0;ThisDateNum='ALL';end
    if strcmp(ThisDateNum,'ALL')
        ThisRecs=AnalysisOpts.DateSet_2look;
    else
        ThisRecs=ThisDateNum;
    end
    if exist('RecDateInfo.mat','file')
        load('RecDateInfo.mat');
    else
        %get all of the recordings information
        for Rec=1:length(ThisRecs)
            [TrialTimes{Rec},RuleBlockTrials{Rec},ChannelInfo{Rec},ChannelArea{Rec},ChsSet{Rec}]=GetRecordingInfo(RunonCluster,ThisRecs(Rec),[],6,[],varargin);
            ClusterChannels{Rec}=ClusterFunc.ConvMat2Char(ChsSet{Rec});
        end
    end
    ThisRecstxt=ClusterFunc.ConvMat2Char(ThisRecs);
%     for Rec=1:length(ThisRecs)
%         % Prepare ephysdata
%         [response_PrepareEphysData,job_id_PrepareEphysData,script_name_PrepareEphysData]=ClusterFunc.RunFunctiononCluster('Run_PrepareEphysData','PrepareEphysData',...
%             {1,ThisRecs(Rec),'$SLURM_ARRAY_TASK_ID','[]','[]','ReRefChannels',AnalysisOpts.ReRefChannels,'MotifAnalysis.L_ms',AnalysisOpts.MotifAnalysis.L_ms},...
%             {'%i','%i','%s','%s','%s',"'%s'",'%i',"'%s'",'%i'},...
%             'RunWithDepenency',0,'Array',ClusterChannels{Rec},'Time',700,'password',SpockPass,'RunThisFunc',0);
%        
%     end
     % discover motifs
    [response,job_id ,script_name]=ClusterFunc.RunFunctiononCluster('Run_AnyFunc','Motif_Cluster',...
        {1,ThisDateNum,'[]',6,'$SLURM_ARRAY_TASK_ID','ReRefChannels',AnalysisOpts.ReRefChannels,...
        'MotifAnalysis.L_ms',AnalysisOpts.MotifAnalysis.L_ms,'MotifAnalysis.ProcessStep',3,'MotifAnalysis.ClustringMethod',AnalysisOpts.MotifAnalysis.ClustringMethod},...
        {'%i',"'%s'",'%s','%i','%s',"'%s'",'%i',"'%s'",'%i',"'%s'",'%i',"'%s'","'%s'"},...
        'RunWithDepenency',0,'job_id_dep',[],'Array','1','Time',200,'password',...
        SpockPass,'RunThisFunc',0);
    
    [response,job_id ,script_name]=ClusterFunc.RunFunctiononCluster('Run_AnyFunc','Motif_SelectivityAnalysis',...
        {1,ThisDateNum,'[]',6,'$SLURM_ARRAY_TASK_ID','ReRefChannels',AnalysisOpts.ReRefChannels,...
        'MotifAnalysis.L_ms',AnalysisOpts.MotifAnalysis.L_ms,'MotifAnalysis.ProcessStep',2,'MotifAnalysis.ClustringMethod',AnalysisOpts.MotifAnalysis.ClustringMethod},...
        {'%i',"'%s'",'%s','%i','%s',"'%s'",'%i',"'%s'",'%i',"'%s'",'%i',"'%s'","'%s'"},...
        'RunWithDepenency',0,'job_id_dep',[],'Array','1','Time',200,'password',...
        SpockPass,'RunThisFunc',1);
    

end
end
