function GetSpecsofThisRun(JobSets,JobID,ArrayID)

ind=(find([JobSets.SLURM_ARRAY_JOB_ID]==JobID & [JobSets.SLURM_ARRAY_TASK_ID]==ArrayID,1,'last'));

if isempty(ind);ind=find([JobSets.SLURM_ARRAY_JOB_ID]==JobID ,1,'last');end

JobSets(ind)


end