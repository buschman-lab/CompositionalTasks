function [response] = ExampleClusterFunction
%EXAMPLECLUSTERFUNCTION Summary of this function goes here
%   Detailed explanation goes here
 response=RunFunctiononCluster('Run1','Motif_CoherenceAnalysis',...
     {1,1,'[]',6,'$SLURM_ARRAY_TASK_ID','ReRefChannels',0},{'%i','%i','%s','%i','%s',"'%s'",'%i'},...
    sprintf('%i-%i',1,5),50,200);
end

