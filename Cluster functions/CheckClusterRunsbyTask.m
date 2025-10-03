Area={'PFC','Striatum','IT','FEF','LIP'};
Condition={'PopulationAnalysis_ConcatClassifier','PopulationAnalysis_RunClassifier_Shuff'};
Task={'3D Color Category Response Xgen BalInCongV2_MS','3D Color Category Discrimination Xgen_MS','3D Shape Category Discrimination Xgen_MS'};
RunTime=[1000,140];
removespaces=@(x) regexprep(x, '[^a-zA-Z]', '');
for c=1:length(Condition)
    for t=1:length(Task)
        for i=1:5
            Ind=find(strcmp({ClusterStatus.JobName},Condition{c}) & ... %PopulationAnalysis_RunClassifier_Shuff
                strcmp({ClusterStatus.TaskName},Task{t}) & ...%'3D Color Category Discrimination Xgen_MS') & ...
                strcmp({ClusterStatus.ExtraInfo},[Area{i} '-100-SAMP-3-' num2str(RunTime(c))]));
            if ~isempty(Ind)
                if ~strcmp(ClusterStatus(Ind(end)).Message,'END')
                    ConditionCheck.(removespaces(Area{i})).(removespaces(Task{t})).(removespaces(Condition{c}))=0;
                    fprintf('\n%s for %s for %s has error ',Condition{c},Task{t},Area{i});
                    arrayfun(@(x) OutputFileClust(ClusterStatus(x).SLURM_ARRAY_JOB_ID,ClusterStatus(x).SLURM_ARRAY_TASK_ID),Ind(end))
                else
                     ConditionCheck.(removespaces(Area{i})).(removespaces(Task{t})).(removespaces(Condition{c}))=1;
                end
            else
                ConditionCheck.(removespaces(Area{i})).(removespaces(Task{t})).(removespaces(Condition{c}))=-1;
                fprintf('\n%s for %s for %s is missing ',Condition{c},Task{t},Area{i});
            end
        end
    end
end
%arrayfun(@(x) OutputFileClust(ClusterStatus(x).SLURM_ARRAY_JOB_ID,ClusterStatus(x).SLURM_ARRAY_TASK_ID),Ind)


IDs=find([ClusterStatus.SLURM_ARRAY_JOB_ID]==34257062);
 problematic=[];
for id=IDs
    if ~(strcmp(ClusterStatus(id).Message,'END') | strcmp(ClusterStatus(id).Message,'STR'))
        problematic=[problematic ClusterStatus(id).SLURM_ARRAY_TASK_ID];
    end
    
end