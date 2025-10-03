% gets a rough report on what is going on in the cluster. 
% @RefDate and @RefTime are relative date/time with respect to now
% e.g RefDate=1 RefTime=1 means yesterday 1 hours ago from now
% they need to be defined outside of here
CL=ClusterFuncs;
CL.RunningStatusCluster;
Y=input('\nDo you want to open error Files(y)','s');
if strcmp(Y,'y');CL.ShowErrorRuns([],1);else;CL.ShowErrorRuns([],0);end
[ClusterStatus,SubmittedJobs]=CL.ClusterStatusReport;
if exist('Ref','var')
    RefDate=Ref(1);RefTime=Ref(2);
    ClusterStatus=CL.LimitDateTimeClusterStatus(ClusterStatus,RefDate,RefTime);
end
openvar('ClusterStatus');
openvar('SubmittedJobs');

