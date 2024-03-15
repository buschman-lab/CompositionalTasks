function  OutputFileClust(SLURM_ARRAY_JOB_ID,SLURM_ARRAY_TASK_ID) % opens output file from cluster

CL=ClusterFuncs;
if ~exist('SLURM_ARRAY_TASK_ID','var') & ischar(SLURM_ARRAY_JOB_ID);SLURM_ARRAY_TASK_ID=[];end
CL.OutputFileClust(SLURM_ARRAY_JOB_ID,SLURM_ARRAY_TASK_ID)
end

