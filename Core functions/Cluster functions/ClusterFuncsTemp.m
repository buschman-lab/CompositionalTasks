classdef ClusterFuncsTemp
    %CLUSTERFUNCS Functions to run function on cluster 
    % from https://npcdocs.princeton.edu/index.php/SpockResources
    % Resource Maximums
    % Current Resource Maximums on Spock
    % CPUs	Memory
    % 24 Cores	250 GB
    % 
    % To access 250GB per node, you need to ask for all 24 cores on the node. By default users are allocated 10GB per core, and you can request to 12GB per core.
    % 
    % -t, --time=
    % The last one is a bit looser, it has no impact on your overall scheduling priority but will impact the priority of the specific task you are submitting, with longer tasks having lower priority overall. However if you exceed the amount of time you ask for your job will be killed by the scheduler. While there are a number of ways specify time listed in the manual or local scheduling requires you to specify it in minutes like this -t 120 to ask for two hours for example. You can ask for up to 4 days (5760 minutes) for any task as of this writing.
    % 
    % Current time limits for jobs SPOCK
    % test <= 5 min	
    % short <= 4 hours	240 minutes
    % long <= 48 hours	2880 minutes
    % very-long <= 4 days	5760 minutes
    properties
       %  script_path='Z:\Users\Sina\Rule Representation Project\Server Analysis forlders\Cluster Analsis General Functions\'
       %  script_path='/Volumes/buschman/Users/Sina/Rule Representation Project/Server Analysis forlders/Cluster Analsis General Functions/'
       %  function_path_PC='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Analysis Pipeline\Code\Pipeline\Primary Analysis\'
       %  clust_status_path_PC='Z:\Projects\Rule_Representation\ElecPhys_Analysis\ClusterStatus\';
       %  function_path_PC='/Volumes/buschman/Projects/Rule_Representation/ElecPhys_Analysis/Rule Representation Project/Analysis Pipeline/Code/Pipeline/Primary Analysis/'
       %  clust_status_path_PC='/Volumes/buschman/Projects/Rule_Representation/ElecPhys_Analysis/ClusterStatus/';
      
        local_bucket='/jukebox/buschman/Users/Sina/Rule\ Representation\ Project/Server\ Analysis\ forlders/Cluster\ Analsis\ General\ Functions/'
        function_path='/jukebox/buschman/Projects/Rule_Representation/ElecPhys_Analysis/Rule Representation Project/Analysis Pipeline/Code/Pipeline/Primary Analysis/'
        clust_status_path='/jukebox/buschman/Projects/Rule_Representation/ElecPhys_Analysis/ClusterStatus/';
        Array='1';
        Memory=4;
        Time=700;
        Ncores=10; %--cpus-per-task 
        RunWithDepenency=0; % should we run this just with dependency 
        DepenencyType='afterok';%dependency type can be after, afterany  % check https://hpc.nih.gov/docs/job_dependencies.html
        job_id_dep=''; % what is the job that we are depending on
        password='';
        username='tafazoli';
        RunThisFunc=1; % should we run this function?
        ExcludeRedshirt=0; % should we exclude redshirt?
        CompResource='spock'; % what computational rescource we are using spock or della
       
    end
    
    methods
        function obj = ClusterFuncs(varargin)
           if nargin~=0 % initialize vars
                obj=obj.ParseParams(varargin) ; %%Process optional inputs
           end
            
        end     
        function out=script_path(~)
            if ispc            
                out='Z:\Users\Sina\Rule Representation Project\Server Analysis forlders\Cluster Analsis General Functions\';
            else
                out='/Volumes/buschman/Users/Sina/Rule Representation Project/Server Analysis forlders/Cluster Analsis General Functions/';
            end
        end
        function out=function_path_PC(~)
            if ispc
                out='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Analysis Pipeline\Code\Pipeline\Primary Analysis\';
            else
                out='/Volumes/buschman/Projects/Rule_Representation/ElecPhys_Analysis/Rule Representation Project/Analysis Pipeline/Code/Pipeline/Primary Analysis/';
            end
        end
        function out=clust_status_path_PC(~)
            if ispc
                out='Z:\Projects\Rule_Representation\ElecPhys_Analysis\ClusterStatus\';
            else
                out='/Volumes/buschman/Projects/Rule_Representation/ElecPhys_Analysis/ClusterStatus/';
            end
        end
        function out=base_path(~)
             if ispc
                 out='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Core functions\Cluster functions\';
             else
                 out='/Volumes/buschman/Projects/Rule_Representation/ElecPhys_Analysis/Core functions/Cluster functions/';
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
        end               
        function script_name = WriteBashScript(obj,uniqID,func_name,input_val,input_type,varargin)
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            %Writes a spock bash script that will run func_name with variables defined
            %by input_val.
            % Array is the SLURM_ARRAY_TASK_ID value range (default is 1) can be
            % [1-10] and is a char
            % Memory is the memory in GBs(Default is 25) 
            % Time is how much we want this to run in mins
            %input_val and input_type are equal length cell arrays.
            %example: input_type = {'"%s"'}; sprintf(input_type{1},'$SLURM_ARRAY_TASK_ID');
            fprintf('\nSending %s to %s ...\n',uniqID,obj.CompResource)
            %Get the text from the base
            text = fileread([obj.base_path 'spock_base.sh']);
            %keep this local for right now;
            script_name = sprintf('run_%s.sh',uniqID);
            file_path = [obj.script_path];
            % if ~exist(file_path)
            %    mkdir(file_path);
            % end
            FullFileName=[file_path script_name];
%             if exist(FullFileName,'file')
%                 fprintf('\nDeleted previous sh file...Writing a new sh file...\n')
%                 delete(FullFileName);  
%                 pause(5)
%             end % if it already exists delete it 
            fid = fopen(FullFileName, 'w');
            fprintf(fid,'%s','#!/usr/bin/env bash'); %add the first line of the header

            %give the script a funny name
            temp = 'buckersjunk';
            fprintf(fid,"\n%s'%s'",'#SBATCH -J ',temp);
            fprintf(fid,'\n#SBATCH --array=%s',obj.Array);            % add the arrays 
         %   fprintf(fid,'\n#SBATCH --mem-per-cpu=%iG',obj.Memory);            % add memory size 
            fprintf(fid,'\n#SBATCH -t %i',obj.Time);            % add time
            fprintf(fid,'\n#SBATCH --cpus-per-task=%i',obj.Ncores); % how many CPUs per task
            if obj.ExcludeRedshirt
                 fprintf(fid,'\n#SBATCH --exclude=redshirt-n[12-49]'); % exclude redshirt
            end

            %Add the base script
            fprintf(fid,'%s',text);
            %Add the script to go to the path
            fprintf(fid,'\ncd ''%s''\n',obj.function_path);

            %Add the specific function call
            try %make sure to close the fid even if crash
               %Create the variable lengthed inputs
               temp = {};
               for i = 1:numel(input_val)
                   if i~=numel(input_val)
                       temp{i} = [sprintf([input_type{i}],input_val{i}),','];
                   else %final round, no comma
                       temp{i} = sprintf([input_type{i}],input_val{i});
                   end
               end
               ErrMessage='sprintf(''Error in function %s() at line %d.\n\nError Message:\n%s'', ME.stack(1).name, ME.stack(1).line, ME.message);';
               fprintf(fid, ['xvfb-run -d matlab -nosplash -nodisplay -nodesktop -r ',...
                   sprintf('"addpath(''/jukebox/buschman/Projects/Rule_Representation/ElecPhys_Analysis/Core functions/Cluster functions/'');CF=ClusterFuncsTemp;'),... % add cluster functions
                   sprintf('CF.UpdateClusterStatus(''%s'',$SLURM_ARRAY_JOB_ID,$SLURM_ARRAY_TASK_ID,%i,''STR'',''%s'',''%s'');',uniqID,obj.Ncores,AnalysisOpts.Spock_TaskName,AnalysisOpts.Spock_ExtraInfo),...
                   sprintf('try;%s(',func_name),... %the function
                   sprintf('%s',[temp{:}]),...%');catch ME;\nsprintf(''Error in function %s() at line %d.\n\nError Message:\n%s'', ME.stack(1).name, ME.stack(1).line, ME.message);\nend;exit;"']); %the inputs
                   sprintf(');CF.UpdateClusterStatus(''%s'',$SLURM_ARRAY_JOB_ID,$SLURM_ARRAY_TASK_ID,%i,''END'',''%s'',''%s'');',uniqID,obj.Ncores,AnalysisOpts.Spock_TaskName,AnalysisOpts.Spock_ExtraInfo),...
                   'catch me;Message=[me.message, ''Line:'' num2str(me.stack(1).line) '' Func:'' me.stack(1).name];disp(Message);',...
                     sprintf('CF.UpdateClusterStatus(''%s'',$SLURM_ARRAY_JOB_ID,$SLURM_ARRAY_TASK_ID,%i,Message,''%s'',''%s'');k=0;end;exit"',...
                     uniqID,obj.Ncores,AnalysisOpts.Spock_TaskName,AnalysisOpts.Spock_ExtraInfo)]); % send start value); %the inputs              

               fclose(fid);
               %convert to unix
               unix2dos([file_path script_name],1)
            catch
              % fclose(fid);
               error('Failed generating bash script');
            end 
        end
        function [response,job_id,script_name]=RunFunctiononCluster(obj,job_name,func_name,input_val,input_type,varargin)
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            %% runs an aribitrary function on cluster
            if ~obj.RunThisFunc
                response=[];job_id=[];script_name=[];
                return
            end
            [obj.RunWithDepenency]=obj.RunWithDepFunc(obj.job_id_dep);
            cd(obj.script_path);

            %% Spock Preprocessing
            s_conn=obj.Connect2Cluster;% conncet to cluster first
            job_id = cell(1,numel(func_name));
            script_name = obj.WriteBashScript(sprintf('%s',job_name),func_name,input_val,input_type);
            if obj.RunWithDepenency %run thejob with dependency 
                if ~strcmp(obj.DepenencyType,'singleton')
                    response = ssh2_command(s_conn,...
                        ['cd ' obj.local_bucket  ';',... %cd to directory
                        sprintf('sbatch --dependency=%s:%s %s',obj.DepenencyType,obj.job_id_dep,script_name)]);
                else
                    response = ssh2_command(s_conn,...
                        ['cd ' obj.local_bucket  ';',... %cd to directory
                        sprintf('sbatch --dependency=%s',obj.DepenencyType)]);
                end
            else
            %Run job
                response = ssh2_command(s_conn,...
                 ['cd ' obj.local_bucket  ';',... %cd to directory
                 sprintf('sbatch %s',script_name)]);
            end
            
            %get job id
            job_id = erase(response.command_result{1},'Submitted batch job ');
                         
            % Close the connection
            ssh2_close(s_conn);
            % update submitted jobs 
            obj.UpdateSubmittedJobs(job_name,str2double(job_id),obj.Array,obj.Ncores,'SUBMITTED',AnalysisOpts.Spock_TaskName,AnalysisOpts.Spock_ExtraInfo)
            %clear personal information
            clear s_conn password username
            %%
        end
        % aux funcs
        function out=ConcatinateJobIdDep(obj,JobIdArray) % concatinates array job ids for dependancy 
            %@JobIdArray is a cell array of job ids in char
            out=[];
            for i=1:length(JobIdArray)
                if ~isempty(JobIdArray{i})
                    out=[out JobIdArray{i} ':'];
                end
            end
            out=out(1:end-1);
        end
        function s_conn=Connect2Cluster(obj)
            % Open ssh connection
            username = obj.username;%input(' Spock Username: ', 's');
            if isempty(obj.password)
                password = passwordEntryDialog('WindowName','Enter your spock password');
            else
                password=obj.password;
            end

            s_conn = ssh2_config([obj.CompResource '.princeton.edu'],username,password);
        end
        function CancelJobs(obj,JobID) %cancels  jobs             
            s_conn=obj.Connect2Cluster;
            if exist('JobID','var')
                if ischar(JobID);JobID=str2double(JobID);end
                txt=sprintf('scancel %i',JobID);
                obj.UpdateSubmittedJobs([],JobID,[],[],'CANCELLED',[],[])
            else
                txt=sprintf('scancel -u tafazoli');
            end
            response = ssh2_command(s_conn,txt);
            ssh2_close(s_conn);
            fprintf(2,response.command_result{1});
        end
         function RunningStatusCluster(obj) %gets running status on cluster 
            s_conn=obj.Connect2Cluster;
            
            txt=sprintf('squeue -u tafazoli');
            response = ssh2_command(s_conn,txt);
            ssh2_close(s_conn);
          %  fprintf(2,response.command_result{1});
        end
        function out=ConvMat2Char(~,Mat) % converts mat into char with ',' (to be used for cluster)
            if ~iscell(Mat)
                if size(Mat,1)>1;Mat=Mat';end
                if length(Mat)>200
                    out=[num2str(Mat(1)) '-' num2str(Mat(end))];
                else
                    out=cell2mat(arrayfun(@(x) [num2str(x) ','],Mat,'UniformOutput',0));
                    out=out(1:end-1);
                end
                
            else 
                % remove empty cells 
                EmptCell=cellfun(@(x) ~isempty(x),Mat);
                Mat=Mat(EmptCell);
                 if length(Mat)>200
                    out=[num2str(Mat(1)) '-' num2str(Mat(end))];
                 else
                     out=cell2mat(cellfun(@(x) [x ','],Mat,'UniformOutput',0));
                     out=out(1:end-1);
                 end
            end
        end
        function [RunWithDepenency,job_id_dep]=RunWithDepFunc(~,job_id_dep) % checks if this dependency var exists if not returns [] and turns off the dependency
               if  ~isempty(job_id_dep)
                   RunWithDepenency=1;
               else
                   RunWithDepenency=0;
               end
        end
        function out=EmptyOutErrFolders(obj) % deletes the contents of out and err folders 
            cd([obj.script_path filesep 'out'])
            delete *.*
            cd([obj.script_path filesep 'err'])
            delete *.*
        end
        function obj=UpdateClusterPath(obj) % changes function path to PC version if we are runnig this on a PC
            if ispc
                obj.clust_status_path=obj.clust_status_path_PC;
            end
        end
        
        function UpdateClusterStatusOld(obj,JobName,SLURM_ARRAY_JOB_ID,SLURM_ARRAY_TASK_ID,NUM_CORES,Message,TaskName,ExtraInfo) % updates cluster status in a local file so we keep track of what is running and errors in there 
            obj=obj.UpdateClusterPath;

            fprintf('Updating Cluster Status file...\n')
            % generate a random delay so we don't overlap writing the tasks
            if isempty(SLURM_ARRAY_TASK_ID);SLURM_ARRAY_TASK_ID=1;end
            rng(SLURM_ARRAY_JOB_ID);            
            if SLURM_ARRAY_TASK_ID>1000; MaxTaskID=2500; else; MaxTaskID=1000;end
            PauseTimeOrder=randsample(MaxTaskID,MaxTaskID);
            pause(PauseTimeOrder(SLURM_ARRAY_TASK_ID)/20+randsample([0.1:0.1:3],1)); % pause so that eveything is not written at the same time 
           
            %% this is the text version of the file 
            ClusterStatusFileName=[obj.clust_status_path obj.CompResource 'ClusterStatus.txt'];
            fid = fopen(ClusterStatusFileName, 'a');
            %Text=fread(fid);
            %fprintf(fid,'%s',Text);
            fprintf(fid,'\n JobName=%s,ARRAY_JOB_ID=%i,ARRAY_TASK_ID=%i,NUM_CORES=%i,Task=%s,ExtraInfo=%s,MSG=%s',JobName,SLURM_ARRAY_JOB_ID,SLURM_ARRAY_TASK_ID,...
                NUM_CORES,TaskName,ExtraInfo,Message);          
            %fclose(fid);
            
            %% this is the mat file version where I get errors such as /spockClusterStatus.mat as a valid MAT-file.
            %   ClusterStatus=matfile(ClusterStatusFileName,'Writable',true);
            %% generate cluster status mat file 
            ClusterStatusFileName=[obj.clust_status_path obj.CompResource 'ClusterStatus.mat'];  
            try % try loading the cluster status file 
                load(ClusterStatusFileName,'ClusterStatus')
            catch % if you catch error reset the file 
                ClusterStatus=struct('JobName',[],'TaskName',[],'ExtraInfo',[],'SLURM_ARRAY_JOB_ID',nan,...
                    'SLURM_ARRAY_TASK_ID',nan,'NUM_CORES',nan,'Message',['Starting Status File'],'TimeStamp',char(datetime('now','format','HH:mm:ss')),...
                    'DateStamp',char(datetime('now','format','MMMM d, yyyy')),'ClusterTime',[],'ClusterMemory',[]);
                save(ClusterStatusFileName,'ClusterStatus')
            end
            % search if we have an identical entry here for this job
            ID=find([ClusterStatus.SLURM_ARRAY_JOB_ID]==SLURM_ARRAY_JOB_ID & [ClusterStatus.SLURM_ARRAY_TASK_ID]==SLURM_ARRAY_TASK_ID,1,'last');
            if isempty(ID)  % if it is empty just add it to the end of it 
                ID=length(ClusterStatus)+1;
            end
            
            % update parameters based on this input
            ClusterStatus(ID).JobName=JobName;
            ClusterStatus(ID).TaskName=TaskName;
            ClusterStatus(ID).ExtraInfo=ExtraInfo;
            ClusterStatus(ID).SLURM_ARRAY_JOB_ID=SLURM_ARRAY_JOB_ID;
            ClusterStatus(ID).SLURM_ARRAY_TASK_ID=SLURM_ARRAY_TASK_ID;
            ClusterStatus(ID).NUM_CORES=NUM_CORES;
            ClusterStatus(ID).Message=Message;
            ClusterStatus(ID).TimeStamp=char(datetime('now','format','HH:mm:ss'));
            ClusterStatus(ID).DateStamp=char(datetime('now','format','MMMM d, yyyy'));
            ClusterStatus(ID).ClusterTime=obj.Time;
            ClusterStatus(ID).ClusterMemory=obj.Memory;
            
            pause(randsample([0.1:0.1:1],1));
            try
                save(ClusterStatusFileName,'ClusterStatus')
                fprintf('Updated Cluster Status Mat file...\n')
            catch
                fprintf('Could not update Cluster Status Mat file...\n')
            end           
        end
        function UpdateClusterStatus(obj,JobName,SLURM_ARRAY_JOB_ID,SLURM_ARRAY_TASK_ID,NUM_CORES,Message,TaskName,ExtraInfo) % updates cluster status in a local file so we keep track of what is running and errors in there 
            obj=obj.UpdateClusterPath;

            fprintf('Updating Cluster Status file...\n')
            if isempty(SLURM_ARRAY_TASK_ID);SLURM_ARRAY_TASK_ID=1;end
            
            % generate a random delay so we don't overlap writing the tasks
%            rng(SLURM_ARRAY_JOB_ID);            
%             if SLURM_ARRAY_TASK_ID>1000; MaxTaskID=2500; else; MaxTaskID=1000;end
%             PauseTimeOrder=randsample(MaxTaskID,MaxTaskID);
%             pause(PauseTimeOrder(SLURM_ARRAY_TASK_ID)/20+randsample([0.1:0.1:3],1)); % pause so that eveything is not written at the same time 
           
            %% this is the text version of the file 
            ClusterStatusFileName=[obj.clust_status_path obj.CompResource 'ClusterStatus.txt'];
            fid = fopen(ClusterStatusFileName, 'a');
            %Text=fread(fid);
            %fprintf(fid,'%s',Text);
            fprintf(fid,'\n JobName=%s,ARRAY_JOB_ID=%i,ARRAY_TASK_ID=%i,NUM_CORES=%i,Task=%s,ExtraInfo=%s,MSG=%s',JobName,SLURM_ARRAY_JOB_ID,SLURM_ARRAY_TASK_ID,...
                NUM_CORES,TaskName,ExtraInfo,Message);          
            %fclose(fid);
            
            %% this is the mat file version where I get errors such as /spockClusterStatus.mat as a valid MAT-file.
            %   ClusterStatus=matfile(ClusterStatusFileName,'Writable',true);
            %% generate cluster status mat file 
            ClusterStatusFileName=[obj.clust_status_path obj.CompResource 'ClusterStatus' num2str(SLURM_ARRAY_JOB_ID) '_' num2str(SLURM_ARRAY_TASK_ID) ];  
            try % try loading the cluster status file 
                load(ClusterStatusFileName,'ClusterStatus')
            catch % if you catch error reset the file 
                ClusterStatus=struct('JobName',[],'TaskName',[],'ExtraInfo',[],'SLURM_ARRAY_JOB_ID',nan,...
                    'SLURM_ARRAY_TASK_ID',nan,'NUM_CORES',nan,'Message',['Starting Status File'],'TimeStamp',char(datetime('now','format','HH:mm:ss')),...
                    'DateStamp',char(datetime('now','format','MMMM d, yyyy')),'ClusterTime',[],'ClusterMemory',[]);
            end
            % search if we have an identical entry here for this job
            ID=find([ClusterStatus.SLURM_ARRAY_JOB_ID]==SLURM_ARRAY_JOB_ID & [ClusterStatus.SLURM_ARRAY_TASK_ID]==SLURM_ARRAY_TASK_ID,1,'last');
            if isempty(ID) & ~isempty(ClusterStatus(1).JobName)  % if it is empty just add it to the end of it 
                ID=length(ClusterStatus)+1;
            elseif isempty(ID) & isempty(ClusterStatus(1).JobName)
                ID=1;
            end
                
            % update parameters based on this input
            ClusterStatus(ID).JobName=JobName;
            ClusterStatus(ID).TaskName=TaskName;
            ClusterStatus(ID).ExtraInfo=ExtraInfo;
            ClusterStatus(ID).SLURM_ARRAY_JOB_ID=SLURM_ARRAY_JOB_ID;
            ClusterStatus(ID).SLURM_ARRAY_TASK_ID=SLURM_ARRAY_TASK_ID;
            ClusterStatus(ID).NUM_CORES=NUM_CORES;
            ClusterStatus(ID).Message=Message;
            ClusterStatus(ID).TimeStamp=char(datetime('now','format','HH:mm:ss'));
            ClusterStatus(ID).DateStamp=char(datetime('now','format','MMMM d, yyyy'));
            ClusterStatus(ID).ClusterTime=obj.Time;
            ClusterStatus(ID).ClusterMemory=obj.Memory;
            
            try
                save(ClusterStatusFileName,'ClusterStatus')
                fprintf('Updated Cluster Status Mat file...\n')
            catch
                fprintf('Could not update Cluster Status Mat file...\n')
            end           
        end
        function UpdateSubmittedJobs(obj,JobName,SLURM_ARRAY_JOB_ID,SLURM_ARRAY_TASK_ID,NUM_CORES,Message,TaskName,ExtraInfo) % updates what jobs we have submitted in a local file so we keep track of what is running and errors in there 
            obj=obj.UpdateClusterPath;

            % generate a random delay so we don't overlap writing the tasks
            if isempty(SLURM_ARRAY_TASK_ID);SLURM_ARRAY_TASK_ID=1;end
                       
            %% generate cluster status mat file 
            ClusterStatusFileName=[obj.clust_status_path obj.CompResource 'ClusterStatus.mat'];  
            try % try loading the cluster status file 
                load(ClusterStatusFileName,'SubmittedJobs')
            catch % if you catch error reset the file 
                SubmittedJobs=struct('JobName',[],'TaskName',[],'ExtraInfo',[],'SLURM_ARRAY_JOB_ID',nan,...
                    'SLURM_ARRAY_TASK_ID',nan,'NUM_CORES',nan,'Message',['Starting Status File'],'TimeStamp',char(datetime('now','format','HH:mm:ss')),...
                    'DateStamp',char(datetime('now','format','MMMM d, yyyy')),'ClusterTime',[],'ClusterMemory',[]);
            end
            % search if we have an identical entry here for this job
            ID=find([SubmittedJobs.SLURM_ARRAY_JOB_ID]==SLURM_ARRAY_JOB_ID,1,'last');
            if isempty(ID) & ~isempty(SubmittedJobs(1).JobName)  % if it is empty just add it to the end of it 
                ID=length(SubmittedJobs)+1;
            elseif isempty(ID) & isempty(SubmittedJobs(1).JobName)
                ID=1;
            end
            if isempty(JobName) % then we are just updating the message                
                SubmittedJobs(ID).Message=Message;
                SubmittedJobs(ID).TimeStamp=char(datetime('now','format','HH:mm:ss'));
                SubmittedJobs(ID).DateStamp=char(datetime('now','format','MMMM d, yyyy'));                
            else
                % update parameters based on this input
                SubmittedJobs(ID).JobName=JobName;
                SubmittedJobs(ID).TaskName=TaskName;
                SubmittedJobs(ID).ExtraInfo=ExtraInfo;
                SubmittedJobs(ID).SLURM_ARRAY_JOB_ID=SLURM_ARRAY_JOB_ID;
                SubmittedJobs(ID).SLURM_ARRAY_TASK_ID=SLURM_ARRAY_TASK_ID;
                SubmittedJobs(ID).NUM_CORES=NUM_CORES;
                SubmittedJobs(ID).Message=Message;
                SubmittedJobs(ID).TimeStamp=char(datetime('now','format','HH:mm:ss'));
                SubmittedJobs(ID).DateStamp=char(datetime('now','format','MMMM d, yyyy'));
                SubmittedJobs(ID).ClusterTime=obj.Time;
                SubmittedJobs(ID).ClusterMemory=obj.Memory;
            end
            try
                save(ClusterStatusFileName,'SubmittedJobs','-append')
                fprintf('Updated SubmittedJobs records...\n')
            catch
                fprintf('Could not update SubmittedJobs records...\n')
            end           
        end
       
        function GenClusterStatusFile(obj)% generate a cluster status file in the current folder
            obj=obj.UpdateClusterPath;
            ClusterStatusFileName=[obj.clust_status_path obj.CompResource 'ClusterStatus.txt'];
            s=input('Type YES if you want to clear the cluster status txt and mat file:','s');
            if strcmp(s,'YES');else;return;end
            fid = fopen(ClusterStatusFileName, 'w');
            fprintf(fid,'Starting this file  %s',datetime);
            fclose(fid);
            % reset the mat file as well
            ClusterStatusMatFileName=[obj.clust_status_path obj.CompResource 'ClusterStatus.mat'];
            ClusterStatus=struct('JobName',[],'TaskName',[],'ExtraInfo',[],'SLURM_ARRAY_JOB_ID',nan,...
                'SLURM_ARRAY_TASK_ID',nan,'NUM_CORES',nan,'Message',['Starting Status File'],'TimeStamp',char(datetime('now','format','HH:mm:ss')),...
                'DateStamp',char(datetime('now','format','MMMM d, yyyy')),'ClusterTime',[],'ClusterMemory',[]);
           % do we want to clear Submitted jobs as well?
            s=input('Type YES if you want to clear the Submitted Jobs  and mat file:','s');
            if strcmp(s,'YES')
                SubmittedJobs=ClusterStatus; % this is to track what jobs we have submitted
                save(ClusterStatusMatFileName,'ClusterStatus','SubmittedJobs')
            else
                try 
                    save(ClusterStatusMatFileName,'ClusterStatus','-append');
                catch % generate a completely new file
                    SubmittedJobs=ClusterStatus;
                    save(ClusterStatusMatFileName,'ClusterStatus','SubmittedJobs');
                end
            end
        end       
        function OutputFileClust(obj,SLURM_ARRAY_JOB_ID,SLURM_ARRAY_TASK_ID,OpenFile) % opens cluster output file 
            if ischar(SLURM_ARRAY_JOB_ID) %then we are feeding everything in a text 
                DashInd=strfind(SLURM_ARRAY_JOB_ID,'_');
                SLURM_ARRAY_TASK_ID=str2double(SLURM_ARRAY_JOB_ID(DashInd+1:end));
                SLURM_ARRAY_JOB_ID=str2double(SLURM_ARRAY_JOB_ID(1:DashInd-1));
            end
            outfilename=sprintf('%sout%sbuckersjunk_%i_%i.out' ,obj.script_path,filesep,SLURM_ARRAY_JOB_ID,SLURM_ARRAY_TASK_ID);
            if exist('OpenFile','var')
                if OpenFile==1
                    edit(outfilename)
                end
            else
                edit(outfilename)
            end
        end
        function out=ErrFileClust(obj,SLURM_ARRAY_JOB_ID,SLURM_ARRAY_TASK_ID) % opens cluster output file 
            if ischar(SLURM_ARRAY_JOB_ID) %then we are feeding everything in a text 
                DashInd=strfind(SLURM_ARRAY_JOB_ID,'_');
                SLURM_ARRAY_TASK_ID=str2double(SLURM_ARRAY_JOB_ID(DashInd+1:end));
                SLURM_ARRAY_JOB_ID=str2double(SLURM_ARRAY_JOB_ID(1:DashInd-1));
            end
            outfilename=sprintf('%serr%sbuckersjunk_%i_%i.err' ,obj.script_path,filesep,SLURM_ARRAY_JOB_ID,SLURM_ARRAY_TASK_ID);
            % now check if this file is empty or not if it is empty then just say it is empty
            % if not copy the content
            ErrorFile=sprintf('%i_%i',SLURM_ARRAY_JOB_ID,SLURM_ARRAY_TASK_ID);
            if ~exist(outfilename,'file'),fprintf(2,'\nErr File:%s does not exist',ErrorFile);out=0;return;end
            fileID=fopen(outfilename,'r');
            if fseek(fileID, 1, 'bof') == -1
               % fprintf(2,'\nErr File:%s is empty',ErrorFile);
                out=-1;
            else                
                frewind(fileID)
                out=fgetl(fileID);
                %edit(outfilename);
                fprintf(2,'\nErr File %s:%s',ErrorFile,out);
            end
            fclose(fileID);
        end
        function OpenClusterStatusFile(obj)
            obj=obj.UpdateClusterPath;
            ClusterStatusFileName=[obj.clust_status_path obj.CompResource 'ClusterStatus.txt'];
            edit(ClusterStatusFileName);
        end
        function [ClusterStatus]=LoadClusterStatusMatFile(obj)
            obj=obj.UpdateClusterPath;
            ClusterStatusMatFileName=[obj.clust_status_path obj.CompResource 'ClusterStatus.mat'];
            load(ClusterStatusMatFileName,'ClusterStatus');
        end
        function ShowErrorRuns(obj,DateTh,OpenFile) % opens all of the runs that have resulted in error
            %@ DateTh the date threshold we care about default(30 days) 
            %@OpenFile are we opennign the output files
            obj=obj.UpdateClusterPath;
            ClusterStatusMatFileName=[obj.clust_status_path obj.CompResource 'ClusterStatus.mat'];

            [ClusterStatus,SubmittedJobs]=obj.ClusterStatusReport;
            nID=length(ClusterStatus);
            ErrorID=arrayfun(@(x) ~(strcmp(ClusterStatus(x).Message,'STR') | strcmp(ClusterStatus(x).Message,'END') | ...
                isempty(ClusterStatus(x).Message) | strcmp(ClusterStatus(x).Message,'Starting Status File')),1:nID);
            if ~exist('DateTh','var');DateTh=datetime-30;end
            if isempty(DateTh);DateTh=datetime-30;end
            DateID=arrayfun(@(x) ClusterStatus(x).DateStamp>DateTh,1:nID);
            % open all output files (that have code error in them)
            arrayfun(@(x) obj.OutputFileClust(ClusterStatus(x).SLURM_ARRAY_JOB_ID,ClusterStatus(x).SLURM_ARRAY_TASK_ID,OpenFile),...
                find(ErrorID & DateID))
            % now go through cluster status report and see if there are
            % files that have any error (TIME OUT etc) in them 
            TimeOutErrID=find(arrayfun(@(x) strcmp(ClusterStatus(x).Message,'STR'),1:nID) & DateID) ;
            TimeOutErrTxt=arrayfun(@(x) obj.ErrFileClust(ClusterStatus(x).SLURM_ARRAY_JOB_ID,...
                ClusterStatus(x).SLURM_ARRAY_TASK_ID),TimeOutErrID,'UniformOutput',0);
            % now update the cluster status file for files that have TIME out Error
            TimeOutErrIDind=cellfun(@(x) ~(strcmp(num2str(x),'0') | strcmp(num2str(x),'-1')),TimeOutErrTxt);          
            for ID= find(TimeOutErrIDind)% find the index that has problem and replace it
                ClusterStatus(TimeOutErrID(ID)).Message=TimeOutErrTxt{ID};
            end
            % now look for the same ErrorIDs in the SubmittedJobs and update if there has been any error
            AllErrIDs=[TimeOutErrID(TimeOutErrIDind) find(ErrorID & DateID)];
            if ~isempty(AllErrIDs)
                SJIds=arrayfun(@(x) [SubmittedJobs.SLURM_ARRAY_JOB_ID]==ClusterStatus(x).SLURM_ARRAY_JOB_ID,AllErrIDs,'UniformOutput',0);
                for i=1:length(AllErrIDs)
                    for k=find(SJIds{i})
                        SubmittedJobs(k).Message=ClusterStatus(AllErrIDs(i)).Message;
                    end
                end
            end
            % save the cluster status file again but first delete the current file 
            save(ClusterStatusMatFileName,'ClusterStatus','SubmittedJobs','-append')
            pause(5)
        end
        function [ClusterStatus,SubmittedJobs]=ClusterStatusReport(obj)% creates a report from the cluster statuses
            obj=obj.UpdateClusterPath;
            ClusterStatusMatFolder=[obj.clust_status_path  ];% obj.CompResource 'ClusterStatus.mat'
            if ~strcmp(ClusterStatusMatFolder,'Z:\Projects\Rule_Representation\ElecPhys_Analysis\ClusterStatus\')
                error('Wrong clusterstatus folder')
            end
            a=dir(ClusterStatusMatFolder);na=length(a);
            UseableFiles=find(arrayfun(@(x) contains(a(x).name,[obj.CompResource 'ClusterStatus']) & ...
                contains(a(x).name,['.mat']) & ~contains(a(x).name,'old','IgnoreCase',1) & ...
                ~strcmp(a(x).name,[obj.CompResource 'ClusterStatus.mat']),1:na));
            if ~isempty(UseableFiles)
                ClusterStatusMats=arrayfun(@(x) load([a(1).folder filesep a(x).name],'ClusterStatus'),UseableFiles);
                for i=1:length(ClusterStatusMats)
                    ThisClusterStatus(i)=ClusterStatusMats(i).ClusterStatus;
                end
            else
                ThisClusterStatus=[];
            end
            ClusterStatusMatFileName=[obj.clust_status_path obj.CompResource 'ClusterStatus.mat'];
            % load any cluster status we alrady have and add stuff to it 
            if exist(ClusterStatusMatFileName,'file')
                load(ClusterStatusMatFileName,'ClusterStatus','SubmittedJobs')
            else
                ClusterStatus=[];SubmittedJobs=[];
            end            
            ClusterStatus=[ClusterStatus ThisClusterStatus];
            if isempty(ClusterStatus) & isempty(ThisClusterStatus);obj.GenClusterStatusFile;return;end
            save(ClusterStatusMatFileName,'ClusterStatus','-append')
            if ~isempty(UseableFiles)
                fprintf(2,'\nDeleting all status files...\n')
                arrayfun(@(x) delete([a(1).folder filesep a(x).name]),UseableFiles);
            end
            %openvar('ClusterStatus')
        end
        function ClusterStatus=LimitDateTimeClusterStatus(~,ClusterStatus,RefDate,RefTime) % Limits date and time for cluster status structure 
            % @RefDate and @RefTime are relative date/time with respect to now
            % e.g RefDate=1 RefTime=1 means yesterday 1 hours ago from now
            nID=length(ClusterStatus);      
            if isempty(RefDate);RefDate=0;end;if isempty(RefTime);RefTime=0;end

          %  DateID=arrayfun(@(x) ClusterStatus(x).DateStamp>=(datetime('now'-RefDate)),1:nID);
          %  TimeID=arrayfun(@(x) ClusterStatus(x).TimeStamp>=timeofday(datetime('now')-hours(RefTime)),1:nID);
            RefID=arrayfun(@(x)  [ClusterStatus(x).DateStamp ' ' ClusterStatus(x).TimeStamp]>=(datetime('now')-RefDate-hours(RefTime)),1:nID);
            ClusterStatus=ClusterStatus(RefID);
        end
        function Classifier_TaskSpockTime=GetClassifierSpockTime(~) % gets requested time for running classifer 
            global AnalysisOpts
            
            % define timing based on the processing steps                                 
            Classifier_TaskSpockTime=AnalysisOpts.PopulationAna.Classifier_TaskSpockTime;
             
            if AnalysisOpts.DividSpockClassifier==1 % if we are running all reps in one then calculate based on reps
                FullFoldRunTime=floor(AnalysisOpts.CurrentClassifierOpts.Nrep*5*25/60)*(AnalysisOpts.ExchangeableCalShuffClassifier+1); % this 5 folds 25sec for fold
                Classifier_TaskSpockTime=FullFoldRunTime*ones(1,length(Classifier_TaskSpockTime));
            elseif AnalysisOpts.DividSpockClassifier==3
                if AnalysisOpts.CurrentClassifierOpts.UseCV;nCV=10;else;nCV=1;end
                if ~AnalysisOpts.CalShuffleClassifier
                     if AnalysisOpts.AreaNum==1 % for PFC adjust timing
                         FullFoldRunTime=floor(AnalysisOpts.CurrentClassifierOpts.Nrep*nCV*12/60)+20; % this 10 folds 12 sec for fold 60 mins to start
                     else
                         FullFoldRunTime=floor(AnalysisOpts.CurrentClassifierOpts.Nrep*nCV*8/60)+20; % this 10 folds 6 sec for fold 60 mins to start
                     end
                else
                   % % trying to keep this under 2880 mins, otherwise it gets into very long task catgeory( this timing didn't work with redshirt
                   % FullFoldRunTime=floor(AnalysisOpts.CurrentClassifierOpts.NrepShufperFold*10*17/60)+45; % this 10 folds 17 sec for fold 60 mins to start
                   if AnalysisOpts.UseRep4Cluster & contains(AnalysisOpts.CurrentClassifierOpts.Name,'Learning3D') % we are running each shuffle repetition seperately
                       FullFoldRunTime=floor(AnalysisOpts.CurrentClassifierOpts.NrepShuf*200*16/60); % this 10 folds 20 sec for fold 60 mins to start
                   elseif AnalysisOpts.UseRep4Cluster & contains(AnalysisOpts.CurrentClassifierOpts.Name,'BalInCongV')
                       if AnalysisOpts.AreaNum==1 % for PFC adjust timing 
                           FullFoldRunTime=floor(AnalysisOpts.CurrentClassifierOpts.NrepShuf*160*7/60); % this 10 folds 20 sec for fold 60 mins to start
                       else% for other areas adjust timing for shorter period
                           FullFoldRunTime=floor(AnalysisOpts.CurrentClassifierOpts.NrepShuf*140*7/60); % this 10 folds 20 sec for fold 60 mins to start                           
                       end
                   else
                       if AnalysisOpts.CalShuffTrlOrderClassifier
                           FullFoldRunTime=floor(AnalysisOpts.CurrentClassifierOpts.NrepShuf*nCV*28/60); % this 10 folds 20 sec for fold 60 mins to start
                       else
                           FullFoldRunTime=floor(AnalysisOpts.CurrentClassifierOpts.NrepShuf*nCV*28/60); % this 10 folds 20 sec for fold 60 mins to start
                       end
                   end
                end
                if AnalysisOpts.GetOnlyShuffLabelsClassifier;FullFoldRunTime=60;end % if are only saving the shuffle labels
                Classifier_TaskSpockTime=FullFoldRunTime*ones(1,length(Classifier_TaskSpockTime));
            end
            % if we are running cross temporal classifier 
            if AnalysisOpts.RunCrossTemporalClassifer
                NConds=AnalysisOpts.CurrentClassifierOpts.NConds;
                if AnalysisOpts.DividSpockClassifier==3
                    error('we recommend using DividSpockClassifier=2 for RunCrossTemporalClassifer=1 ')
                    Classifier_TaskSpockTime=5000*ones(1,length(Classifier_TaskSpockTime));
                elseif AnalysisOpts.DividSpockClassifier==2  % assuming that we are not running shuffle on this
                    % assume NConds conditions and 10 CV 100 sec per run
                    Classifier_TaskSpockTime=floor((NConds*10*100)/60+60)*ones(1,length(Classifier_TaskSpockTime));
                end
            end
            % if we are sweeping classifier conditions
            if AnalysisOpts.SweepClassifierConds;Classifier_TaskSpockTime=3*AnalysisOpts.PopulationAna.Classifier_TaskSpockTime;end
            % if we are running dummy file then put everything to 30 mins 
          %  if AnalysisOpts.RunDummyFile;Classifier_TaskSpockTime=30*ones(1,length(Classifier_TaskSpockTime));end
            % if we are plotting files 
            if AnalysisOpts.ProcessingStep==4
                if AnalysisOpts.RunCrossTemporalClassifer
                    Classifier_TaskSpockTime=600*ones(1,length(Classifier_TaskSpockTime));
                else
                    Classifier_TaskSpockTime=240*ones(1,length(Classifier_TaskSpockTime));
                end               
            end
            % if we are plotting comparisions
            if AnalysisOpts.ProcessingStep==5
                Classifier_TaskSpockTime=AnalysisOpts.PopulationAna.ClassifierComparision_SpockTime;
            end
            % if we are concatinating conditions 
            if AnalysisOpts.ProcessingStep==8
                Classifier_TaskSpockTime=AnalysisOpts.PopulationAna.Classifier_TaskConcateSpockTime;
                if AnalysisOpts.RunCrossTemporalClassifer 
                    Classifier_TaskSpockTime=2*Classifier_TaskSpockTime;
                end
            end
        end
        function GLM_SpockTime=GetGLMSpockTime(~)% gets requested time for running GLM analysis
            global AnalysisOpts
            
            switch AnalysisOpts.ProcessingStep
                case 11 % fit and run GLM shuffle for less than 100
                    GLM_SpockTime=AnalysisOpts.GLMSpockFitTime*length(AnalysisOpts.SingCellAna.GLMMdlName2Test);
                case 16 % plot results for each area
                    GLM_SpockTime=1000;
                case 15
                    GLM_SpockTime=5;
                otherwise
                    GLM_SpockTime=200;
            end
        end
        
        function password = maskedInput(obj,prompt)
            % Initialize an empty password string
            password = '';

            % Display the prompt
            fprintf(prompt);

            % Disable the command window's echo to mask the input
            oldEchoState = get(0, 'Echo');
            set(0, 'Echo', 'off');

            % Create a loop to handle each character individually
            while true
                % Capture the user's input one character at a time
                c = obj.getChar();

                % Break the loop when the Enter key is pressed
                if c == 13
                    break;
                end

                % Append the character to the password string and display an asterisk
                password = [password, c];
                fprintf('*');
            end

            % Re-enable the command window's echo
            set(0, 'Echo', oldEchoState);

            % Add a new line character to match the behavior of input()
            fprintf('\n');
        end

        function c = getChar(obj)
            % This function captures a single character from the keyboard
            c = obj.getCharJava;
        end

        function c = getCharJava(~)
            % This function captures a single character using Java's System.in
            import java.io.*;
            import java.util.*;

            br = BufferedReader(InputStreamReader(System.in));
            c = char(br.read());
        end


    end
end

