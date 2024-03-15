% cancels jobs on cluster
function CancelClusterJobs(JobID)
CL=ClusterFuncs;
if exist('JobID','var') % cancel specific job
    CL.CancelJobs(JobID)
else
    CL.CancelJobs % cancel all jobs
end
