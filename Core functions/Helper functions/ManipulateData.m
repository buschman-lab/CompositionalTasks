classdef ManipulateData

    properties
        Dim=3;
        smoothingFactor=4;
        DimSmoothing; % what dimension are we perforing the smoothing 
        WantedDate=''; % look at this specific date instead
        WantedCh=''; % which channel we are looking at 
        SaveInResults=''; % save this in results folder
        Deg=1; % if we calculate angle in degree
        SelfName=0; % name of the file is defined by user
        SelectivityIndexMethod='norm'; % can be 'norm','mean','morph'
        SmoothingMethod='gauss';
        SavedinSameNameFolder=0; % if we have saved a file in a folder with the same name of the file (usually for images)
        SaveAnalysisOpts=1; % are we savign AnalysisOpts variable whenever we save any file?
        CheckifFileExistsOnly=0; % when we load a file do we only check if the file exists and don't load or save it?
        ShowClustCorrectionPlot=0; % are we shoing the plot for details of cluster correction
    end
    properties (Access=private)
    end
    methods
        function  obj=ManipulateData(varargin)
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
        end
        function out=Power(~,data)
            out=abs(data).^2;
        end  % Calcultes power
        function out=NormPower(~,data,f) % calculates the normlaized power (Power normalized by 1/f)
            % data is Frequency by Time
            power=abs(data).^2;
            if size(f,1)==1;f=f';end
            out=power.*repmat(f,size(data,1)/length(f),size(data,2));
        end
        function out=NormMedianPower(~,data)
            power=abs(data).^2;
            MedianPower=median(power,2);
            out=power./repmat(MedianPower,1,size(data,2));
        end
        function out=NormMeanPower(~,data,f)% takes percent change of power around the mean
            if ~isreal(data) % then calculate power
                data=data.*conj(data);
            end
            if ~isempty(f) % then normalize by f as well
                if size(f,1)==1;f=f';end
                data=data.*f;
            end
            dataRS=reshape(data,[size(data,1) size(data,2) size(data,3)*size(data,4)]);

            AvgNormPower=nanmean(nanmean(dataRS,3),2);
            out=(data-repmat(AvgNormPower,[1 size(data,2) size(data,3) size(data,4)]))./(repmat(AvgNormPower,[1 size(data,2) size(data,3) size(data,4)]));
        end
        function out=RealImag(obj,data)
            RealData=real(data);RealData=obj.NormalizeMinMax(RealData);
            ImagData=imag(data);ImagData=obj.NormalizeMinMax(ImagData);
            out=[RealData;ImagData];
        end
        function out=Cell2Mat(obj,data)
            if iscell(data)
                out=[];
                if size(data,1)==1 %% we are working on one channel only
                    for i=1:length(data)
                        X=arrayfun(@(x) data{i}(:,:,x),1:size(data{i},obj.Dim),'UniformOutput',0);
                        out=[out cell2mat(X)];
                    end
                elseif size(data,1)>1 %% we have multiple channels so we concatinate them vertically

                    for j=1:size(data,1)
                        temp=[];
                        for i=1:size(data,2)
                            X=arrayfun(@(x) data{j,i}(:,:,x),1:size(data{j,i},obj.Dim),'UniformOutput',0);
                            temp=[temp cell2mat(X)];
                        end
                        out=[out;temp];
                    end
                end
            else
                X=arrayfun(@(x) data(:,:,x),1:size(data,obj.Dim),'UniformOutput',0);
                out=cell2mat(X);
            end
        end % concatinates the cell trial data into a matrix in time
        function [out,Time]=SmoothData(obj,data,width,varargin) % smooth data with different methods
            obj=obj.ParseParams(varargin) ; %%Process optional inputs

            switch obj.SmoothingMethod
                case 'binom'
                    h = [1/2 1/2];
                    binomialCoeff = conv(h,h);
                    for n = 1:width % smoothingFactor
                        binomialCoeff = conv(binomialCoeff,h);
                    end
                    out = filter(binomialCoeff, 1, data);
                    Time=1:length(data);
                case 'gauss'
                 %   out=obj.GaussianSmooth(data',width);
                 %   out=out';
                    out=smoothdata(data,2,'gaussian',width);
                    Time=1:length(data);
                case 'movmean'
                   % [out,TimeInds]=obj.MoveMean(data,width,'valid');
                   % Time=TimeInds(:,2);
                   if isempty(obj.DimSmoothing);Dim=2;else;Dim=obj.DimSmoothing;end
                    out=smoothdata(data,Dim,'movmean',width);
                    Time=1:length(data);                    
            end
        end  % smooths data
        function out=GaussianSmooth(~,data,width) % gaussian smooth the data
            w=gausswin(width);
            out=filtfilt(w/sum(w),1,data);
        end
        function out=ReshapeTrials(~,data,Ntrials)
            % Ntrials is the number trials for each rule for example
            %[ 51 100 40]
            Ntrials=[0 Ntrials];
            NTim=size(data,2);
            NtrialsTot=sum(Ntrials);
            NTrialTime=NTim/NtrialsTot;
            TrialNumInd=arrayfun(@(x) sum(Ntrials(1:x)), 1:length(Ntrials));
            for l=1:size(data,1)
                ThisL= reshape(data(l,:),NtrialsTot,NTrialTime);
                for i=1:length(Ntrials)-1
                    out{l,i}(:,:)=ThisL(TrialNumInd(i)+1:TrialNumInd(i+1),:);
                end
            end

        end  % reshapes the H matrix into trials for each rule
        function out=Turn2Cell(~,data)  % turns the data into cell
            out=arrayfun(@(x) data(x),1:length(data),'UniformOutput',0);
        end
        function out=Reshape2DMat2Cell(~,data,Dim) % reshpaes a 2D mat into 2D matrixes and puts each raw into a cell
            out=arrayfun(@(x) reshape(data(x,:),Dim),1:size(data,1),'UniformOutput',0);
        end
        function out=Reshape4DMat(~,data,Dim,Resize) % reshapes a 2D image mat into 4d for VAE training
            % resize if we need to
            if isempty(Resize)
                Resize=Dim;
            end
            for i=1:size(data,1)
                out(:,:,1,i)= imresize(reshape(data(i,:),Dim),Resize);
            end

        end
        function out=ReshapeCell2Mat(obj,data,Dim,varargin) % reshapes contents of cell into matrix in the specified dimension
            out=[];
            switch Dim
                case 1
                    for i=1:length(data)
                        f=size(data{i},1);
                        d=size(out,2);
                        out(:,d+1:d+f)=data{i}';
                    end
                case 2
                    for i=1:length(data)
                        f=size(data{i},1);
                        d=size(out,1);
                        out(d+1:d+f,:)=data{i}';
                    end
                case 3
                    for i=1:length(data)
                        out(:,:,i)=data{i};
                    end
                case 4
                    for i=1:length(data)
                        out(:,:,:,i)=data{i};
                    end
                case 5
                    for i=1:length(data)
                        out(:,:,:,:,i)=data{i};
                    end
                case 32 % copy dimensions of cell into second dimension
                    for i=1:length(data)
                        out(:,i,:)=data{i};
                    end
                case 42 % copy raw dimension of cell into 4th dimension
                    [ii,jj]=size(data);
                    for j=1:jj
                        for i=1:ii
                            out{j}(:,:,:,i)=data{i,j};
                        end
                    end
                case 52 % concatinates data from all of the conditions from a specified field of a structure
                    % as an input from varargin
                    % structs them so that each cell has its data into 4
                    % dimensions (Motif*Time*Freq*Block)
                    FieldName=varargin{1};

                    NRecs=length(data); % number of recordings
                    NchsTot=sum(arrayfun(@(x) size(data(x).(FieldName),1),1:NRecs)); % total number of channels
                    out=cell(1,NchsTot);ChCnt=1;
                    % iterate and put the data for each channel into one
                    % cell(Motif*Time*Freq*Block)
                    for Rec=1:NRecs
                        NChs=size(data(Rec).(FieldName),1);
                        Nblks=size(data(Rec).(FieldName),2);
                        for ch=1:NChs
                            for Blk=1:Nblks
                                out{ChCnt}=cat(4,out{ChCnt},data(Rec).(FieldName){ch,Blk});
                            end
                            ChCnt=ChCnt+1; % increase the number of channels
                        end
                    end
                case 62 % concatinate the content of the cells into a two dimensions along rows (mostly used for concatinating PSTHs)
                    out=[];
                    for i=1:length(data)
                        out=cat(1,out,data{i});
                    end
                case 63 % concatinate the content of the cells into a three dimensions to indicate blocks (mostly used for concatinating PSTHs)
                    out=[];
                    for i=1:length(data)
                        out=cat(3,out,data{i});
                    end  
                 case 64 % concatinate the content of the cells into a two dimensions along columns(mostly used for concatinating PSTHs)
                    out=[];
                    for i=1:length(data)
                        out=cat(2,out,data{i});
                    end   
                otherwise
                    out=data;
            end
        end
        function out=ReshapeStruct2Mat(~,data,fieldname,Dim,varargin) % reshapes contents of struct into matrix in the specified dimension
            switch Dim
                case 1
                    for i=1:length(data)
                        out(:,i)=data(i).(fieldname);
                    end
                case 2
                    for i=1:length(data)
                        out(i,:)=data(i).(fieldname);
                    end
                case 3
                    for i=1:length(data)
                        out(:,:,i)=data(i).(fieldname);
                    end
                case 4
                    for i=1:length(data)
                        out(:,:,:,i)=data(i).(fieldname);
                    end
                case 5
                    for i=1:length(data)
                        out(:,:,:,:,i)=data(i).(fieldname);
                    end
                case 32 % copy dimensions of struct into second dimension
                    for i=1:length(data)
                        out(:,i,:)=data(i).(fieldname);
                    end
                case 42 % copy raw dimension of struct into 4th dimension
                    [ii,jj]=size(data);
                    for j=1:jj
                        for i=1:ii
                            out{j}(:,:,:,i)=data(i,j).(fieldname);
                        end
                    end
                case 52 % concatinates data from all of the conditions from a specified field of a structure
                    % as an input from varargin
                    % structs them so that each cell has its data into 4
                    % dimensions (Motif*Time*Freq*Block)
                    FieldName=varargin{1};

                    NRecs=length(data); % number of recordings
                    NchsTot=sum(arrayfun(@(x) size(data(x).(FieldName),1),1:NRecs)); % total number of channels
                    out=cell(1,NchsTot);ChCnt=1;
                    % iterate and put the data for each channel into one
                    % cell(Motif*Time*Freq*Block)
                    for Rec=1:NRecs
                        NChs=size(data(Rec).(FieldName),1);
                        Nblks=size(data(Rec).(FieldName),2);
                        for ch=1:NChs
                            for Blk=1:Nblks
                                out{ChCnt}=cat(4,out{ChCnt},data(Rec).(FieldName){ch,Blk});
                            end
                            ChCnt=ChCnt+1; % increase the number of channels
                        end
                    end
                case 62 % concatinate the content of the cells into a two dimensions (mostly used for concatinating PSTHs)
                    out=[];
                    for i=1:length(data)
                        out=cat(1,out,data(i).(fieldname));
                    end
                case 63 % concatinate the content of the cells into a three dimensions to indicate blocks (mostly used for concatinating PSTHs)
                    out=[];
                    for i=1:length(data)
                        out=cat(3,out,data(i).(fieldname));
                    end
                otherwise
                    out=data;
            end
        end

        function out=ReshapeCellStruct2Mat(obj,data,Dim,varargin) % reshape each field of cell struct to a matrix
            if isempty(varargin)
                FieldNames=fieldnames(data{1});
            else
                FieldNames= varargin;
            end
            for f=1:length(FieldNames)
                Field=FieldNames{f};
                out.(Field)=[];
                for i=1:length(data)
                    if size(data{i}.(Field),3)==1 && size(data{i}.(Field),4)==1 && Dim==3                        
                        out.(Field)(:,:,i)=data{i}.(Field);
                    elseif size(data{i}.(Field),3)==1 && size(data{i}.(Field),4)==1 && Dim==1
                        out.(Field)=cat(1,out.(Field),data{i}.(Field));
                    elseif size(data{i}.(Field),3)>1 && size(data{i}.(Field),4)==1
                        out.(Field)(:,:,:,i)=data{i}.(Field);
                    elseif   size(data{i}.(Field),4)>1
                        out.(Field)(:,:,:,:,i)=data{i}.(Field);
                    end
                end
            end
        end
        function out=Normalize_0_1(~,data)  %normalizes data to 0 and 1
            if size(data,3)~=1  % if we are 3dim apply on first 2
                for i=1:size(data,3)
                    temp=data(:,:,i);
                    X=mapminmax(temp(:)',0,1);
                    out(:,:,i)=reshape(X,size(X));
                end
            else   %if 2 dim appy on each raw
                out=mapminmax(data,0,1);
            end

        end
        function out=NormalizeMinMax(~,data) % nromalizes data using man and max of data without using mapminmax func
            MinData=min(data(:));MaxData=max(data(:));
            out=(data-MinData)/(MaxData-MinData);
        end
        function out=NormalizeMotifs(~,data) % normalizes motifs that are in cell
            SizeW=size(data{1});
            out=cellfun(@(x) reshape(mapminmax(x(:),0,1),SizeW(1),SizeW(2)),data,'UniformOutput',0);

        end
        function out = gauss2d(obj,mat, sigma, center)
            gsize = size(mat);
            [R,C] = ndgrid(1:gsize(1), 1:gsize(2));
            out = obj.gaussC(R,C, sigma, center);

        end
        function val = gaussC(obj,x, y, sigma, center)
            xc = center(1);
            yc = center(2);
            exponent = ((x-xc).^2 + (y-yc).^2)./(2*sigma);
            val = (exp(-exponent));
        end
        function out = ReshapeSquareMatrix(~,Pairs,Vals) % arranges pais of masurmeants into a simularity matrix
            Npairs=size(Pairs,1);
            Maxi=max(Pairs(:));
            out=nan(Maxi);
            for i=1:Npairs
                out(Pairs(i,1),Pairs(i,2))=Vals(i);
            end
        end
        function out = MakeASymmetricMatrix(~,Mat) % replaces matrix with a symmetiric one
            if iscell(Mat)
                LMat=size(Mat{1},1);
                for i=1:length(Mat)
                    X=zeros(size(LMat));
                    out{i}=X+triu(Mat{i},1)+triu(Mat{i},1)';
                end
            else
                LMat=size(Mat,1);
                X=zeros(size(LMat));
                out=X+triu(Mat,1)+triu(Mat,1)';
            end
        end
        function Struct = rmfieldExept(~,Struct,WantedFields) % removes all fields of a struct exept the wanted
            % WantedFields is cell array containg fields we want to keep
            StrctFeilds=fieldnames(Struct);
            for i=1:length(StrctFeilds)
                if ~sum(strcmp(WantedFields,StrctFeilds{i}))
                    Struct=rmfield(Struct,StrctFeilds{i});
                end
            end
        end
        function out=TrimMatrix(~,data,TrimSiz) % trims a matrix in a given dimension (used for Motifs)
            if ~iscell(data)
                if size(data,3)==3*TrimSiz
                    out=data(:,:,TrimSiz+1:2*TrimSiz);
                else
                    out=data;
                end
            else
                if size(data{1},2)==3*TrimSiz
                    out=cellfun(@(x) x(:,TrimSiz+1:2*TrimSiz),data,'UniformOutput',0);
                else
                    out=data;
                end
            end
        end
        function out=IsVarExistinFile(~,FileName,VarName) % checks if a var exists in a file without openning it
            if exist(FileName,'file')
                variableInfo = who('-file', FileName);
                out=ismember(VarName, variableInfo); % returns true
            else
                out=0;
            end
        end
        function out=AggregateCellVals(~,data) % this is function to aggregate acorss blocks of Xcorrs
            Ndata=length(data);
            Ncells=length(data{1});
            N2cells=length(data{1}{1});
            out=cell(1,Ncells);

            for i= 1:Ndata
                for j=1:Ncells
                    for k=1:N2cells
                        out{1,j}{1,k}(:,:,i)=data{1, i}{1, j}{1, k};
                    end
                end
            end
        end
        function out=ConcatinateCellValls(~,data)
            Ndata=length(data);
            out=[];
            for i=1:Ndata
                if iscell(data{i});temp=data{i};else;temp=data(i);end
                out=cat(2,out,temp);
            end
        end
        function out=GenMovAvgInds(~,NPoints,NAvg,Type) % generates indices of moving average

            switch Type
                case 'valid' % only take the valid portions of average 'Trailing'
                    NSteps=NPoints-NAvg+1;
                    for i=1:NSteps
                        out(i,:)=[i i+NAvg-1];
                    end
                case 'same'  % shrink the window in the last parts to match the window size
                    NSteps=NPoints-NAvg+1;
                    for i=1:NSteps
                        out(i,:)=[i i+NAvg-1];
                    end
                    for i=NSteps+1:NPoints
                        out(i,:)=[i NPoints];
                    end
            end
        end
        function out=GenMovingInds(~,Str,End,Length,Step) % generates indices of moving indices
            if Str==0 & End==0;out={NaN};return;end %if Str and End are 0 then return NaN 
            if Str<0 || End<0;Str=-Str;End=-End;Sign=-1;else;Sign=1;end
            Npoints=(Length:Step:End)-Length;if isempty(Npoints);Npoints=0;end
            out=arrayfun(@(x) Sign*((Str+x):(Str+Length-1+x)),Npoints,'UniformOutput',0);
        end
        function [out,AvgInds]=MoveMean(obj,data,NAvg,Type)

            if size(data,2)==1 & size(data,1)>1
                data=data';
            end

            NPoints=size(data,2);
            AvgInds=obj.GenMovAvgInds(NPoints,NAvg,Type);
            out=cell2mat(arrayfun(@(x) nanmean(data(:,AvgInds(x,1):AvgInds(x,2)),2),1:size(AvgInds,1),'UniformOutput',false));

        end
        function SaveFileFullPath=SaveVar(obj,AnalysisPathName,VarVal,VarName,ExtraTxt,varargin) % saves a specific variable into file and folders
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs

            if ~obj.SelfName
               
                if ~isempty(obj.WantedDate)
                    AnalysisOpts.CurrentRecDate=AnalysisOpts.RecDate;
                    AnalysisOpts.RecDate=obj.WantedDate;
                    CurrentCh_RecDate=AnalysisOpts.CurrentCh_RecDate;
                    AnalysisOpts.CurrentCh_RecDate=obj.WantedDate;
                end
                if ~isempty(obj.WantedCh)
                    CurrentCh=AnalysisOpts.CurrentCh;CurrentChClust=AnalysisOpts.CurrentChClust;
                    AnalysisOpts.CurrentCh=obj.WantedCh;
                    AnalysisOpts.CurrentChClust=[]; 
                end
                if isempty(AnalysisPathName);AnalysisPathName='Analysis_Data';end % just save it as analysis data
                [~,~,SaveFileFullPath]=GenerateFileName(AnalysisOpts.FS,AnalysisOpts.DataSavePath,AnalysisPathName,...
                    AnalysisOpts.CurrentCh_Animal,AnalysisOpts.CurrentCh_RecDate,...
                    [AnalysisOpts.CurrentCh AnalysisOpts.CurrentChClust],ExtraTxt,'ext','.mat','MakeFolder',1,'ExtFolderName',AnalysisOpts.RecDate);
            else
                SaveFileFullPath=ExtraTxt;
            end

            eval([ VarName '=VarVal;']);
            if ~exist(SaveFileFullPath,'file');GenNewFile=1;save(SaveFileFullPath,'GenNewFile','-v7.3');end

            DateTime=datetime; % save the datetime this was saved
            save(SaveFileFullPath,VarName,'DateTime','-append')
            if obj.SaveAnalysisOpts % are we saving analysis opts as well
                save(SaveFileFullPath,'AnalysisOpts','-append')
            end
            fprintf('\nSaving variable:%s in file:%s',VarName,SaveFileFullPath);
            if ~isempty(obj.WantedDate)
                AnalysisOpts.RecDate=AnalysisOpts.CurrentRecDate;
                AnalysisOpts.CurrentCh_RecDate=CurrentCh_RecDate;
            end
            if ~isempty(obj.WantedCh)
                AnalysisOpts.CurrentCh=CurrentCh;
                AnalysisOpts.CurrentChClust=CurrentChClust;
            end
        end
        function SaveFileFullPath=SaveVar2Struct(obj,AnalysisPathName,StructName,VarVal,VarName,ExtraTxt,varargin) % add a variable to an existing structure in the file
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs

            if ~obj.SelfName
                if ~isempty(obj.WantedDate)
                    AnalysisOpts.CurrentRecDate=AnalysisOpts.RecDate;
                    AnalysisOpts.RecDate=obj.WantedDate;
                    CurrentCh_RecDate=AnalysisOpts.CurrentCh_RecDate;
                    AnalysisOpts.CurrentCh_RecDate=obj.WantedDate;
                end
                if ~isempty(obj.WantedCh)
                    CurrentCh=AnalysisOpts.CurrentCh;CurrentChClust=AnalysisOpts.CurrentChClust;
                    AnalysisOpts.CurrentCh=obj.WantedCh;
                    AnalysisOpts.CurrentChClust=[]; 
                end
                if isempty(AnalysisPathName);AnalysisPathName='Analysis_Data';end % just save it as analysis data
                [~,~,SaveFileFullPath]=GenerateFileName(AnalysisOpts.FS,AnalysisOpts.DataSavePath,AnalysisPathName,...
                    AnalysisOpts.CurrentCh_Animal,AnalysisOpts.CurrentCh_RecDate,...
                    [AnalysisOpts.CurrentCh AnalysisOpts.CurrentChClust],ExtraTxt,'ext','.mat','MakeFolder',1,'ExtFolderName',AnalysisOpts.RecDate);
            else
                SaveFileFullPath=ExtraTxt;
            end

            if ~exist(SaveFileFullPath,'file');warning(sprintf('\n%s does not exist',SaveFileFullPath));return;end

            load(SaveFileFullPath,StructName);
            eval([StructName '.' VarName '=VarVal;']);
            DateTime=datetime; % save the datetime this was saved
            save(SaveFileFullPath,StructName,'DateTime','AnalysisOpts','-append')
            fprintf('\nSaving variable:%s in structure:%s in file:%s',VarName,StructName,SaveFileFullPath);
            if ~isempty(obj.WantedDate)
                AnalysisOpts.RecDate=AnalysisOpts.CurrentRecDate;
                AnalysisOpts.CurrentCh_RecDate=CurrentCh_RecDate;
            end
            if ~isempty(obj.WantedCh)
                AnalysisOpts.CurrentCh=CurrentCh;
                AnalysisOpts.CurrentChClust=CurrentChClust;
            end
        end

        function [FileName,Path,FullPath,FullPathALL]=GetFileName(obj,AnalysisPathName,ExtraTxt,varargin) % gets the file name
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs

            if isempty(AnalysisPathName);AnalysisPathName='Analysis_Data';end % just save it as analysis data
            if ~isempty(obj.WantedDate)
                    AnalysisOpts.CurrentRecDate=AnalysisOpts.RecDate;
                    AnalysisOpts.RecDate=obj.WantedDate;
                    CurrentCh_RecDate=AnalysisOpts.CurrentCh_RecDate;
                    AnalysisOpts.CurrentCh_RecDate=obj.WantedDate;
            end
            if ~isempty(obj.WantedCh)
                    CurrentCh=AnalysisOpts.CurrentCh;
                    CurrentChClust=AnalysisOpts.CurrentChClust;
                    AnalysisOpts.CurrentCh=obj.WantedCh;
                    AnalysisOpts.CurrentChClust=[];
            end
            if obj.SaveInResults
                SavePath=AnalysisOpts.ResultsSavePath;
            else
                SavePath=AnalysisOpts.DataSavePath;
            end
            if strcmpi(AnalysisOpts.RecDate,'ALL')
                d=1;
                for ThisDate=AnalysisOpts.DateSet_2look
                    if ~ischar(ThisDate)
                        AnalysisOpts.RecDate=[AnalysisOpts.PreDateTxt AnalysisOpts.DateSet{ThisDate} ];
                    end
                    [~,~,FullPathALL{d}] =GenerateFileName(AnalysisOpts.FS,AnalysisOpts.DataSavePath,AnalysisPathName,...
                        AnalysisOpts.CurrentCh_Animal,AnalysisOpts.CurrentCh_RecDate,[AnalysisOpts.CurrentCh AnalysisOpts.CurrentChClust],ExtraTxt,'ext','.mat');
                    d=d+1;
                end
                AnalysisOpts.RecDate='ALL';
                [FileName,Path,FullPath] =GenerateFileName(AnalysisOpts.FS,SavePath,AnalysisPathName,...
                    AnalysisOpts.CurrentCh_Animal,AnalysisOpts.CurrentCh_RecDate,[AnalysisOpts.CurrentCh AnalysisOpts.CurrentChClust],ExtraTxt,'ext','.mat');
            else

                [FileName,Path,FullPath] =GenerateFileName(AnalysisOpts.FS,SavePath,AnalysisPathName,...
                    AnalysisOpts.CurrentCh_Animal,AnalysisOpts.CurrentCh_RecDate,[AnalysisOpts.CurrentCh AnalysisOpts.CurrentChClust],ExtraTxt,'ext','.mat');
                FullPathALL=FullPath;
            end
            if obj.SavedinSameNameFolder % if this file is saved in the a folder with the name of the file then add a folder with the same name.
                FullPath=[Path FileName(1:end-4) filesep FileName];
            end
            if ~isempty(obj.WantedDate)
                AnalysisOpts.RecDate=AnalysisOpts.CurrentRecDate;
                AnalysisOpts.CurrentCh_RecDate=CurrentCh_RecDate;
            end
            if ~isempty(obj.WantedCh)
                AnalysisOpts.CurrentCh=CurrentCh;
                AnalysisOpts.CurrentChClust=CurrentChClust;
            end
        end
        function out=DeleteFile(obj,AnalysisPathName,ExtraTxt,DoDelete,varargin) % deletes this file
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs

            if ~obj.SelfName
                if ~isempty(obj.WantedDate)
                    AnalysisOpts.CurrentRecDate=AnalysisOpts.RecDate;
                    AnalysisOpts.RecDate=obj.WantedDate;
                    CurrentCh_RecDate=AnalysisOpts.CurrentCh_RecDate;
                    AnalysisOpts.CurrentCh_RecDate=obj.WantedDate;
                end
                if ~isempty(obj.WantedCh)
                    CurrentCh=AnalysisOpts.CurrentCh;CurrentChClust=AnalysisOpts.CurrentChClust;
                    AnalysisOpts.CurrentCh=obj.WantedCh;
                    AnalysisOpts.CurrentChClust=[]; 
                end
                if isempty(AnalysisPathName);AnalysisPathName='Analysis_Data';end % just save it as analysis data
                [~,~,SaveFileFullPath]=GenerateFileName(AnalysisOpts.FS,AnalysisOpts.DataSavePath,AnalysisPathName,...
                    AnalysisOpts.CurrentCh_Animal,AnalysisOpts.CurrentCh_RecDate,...
                    [AnalysisOpts.CurrentCh AnalysisOpts.CurrentChClust],ExtraTxt,'ext','.mat','MakeFolder',1,'ExtFolderName',AnalysisOpts.RecDate);
            else
                SaveFileFullPath=ExtraTxt;
            end

            if ~exist(SaveFileFullPath,'file');out=0;return;end
            if DoDelete
                delete(SaveFileFullPath);
                fprintf('\nFile %s is being deleted ...',SaveFileFullPath);
            end
            if ~isempty(obj.WantedDate)
                AnalysisOpts.RecDate=AnalysisOpts.CurrentRecDate;
                AnalysisOpts.CurrentCh_RecDate=CurrentCh_RecDate;
            end
            if ~isempty(obj.WantedCh)
                AnalysisOpts.CurrentCh=CurrentCh;
                AnalysisOpts.CurrentChClust=CurrentChClust;
            end
        end
        function [out,FileExist,FileFullPath,filesize]=LoadVar(obj,AnalysisPathName,VarName,ExtraTxt,UseDataPointer,varargin) % load a specific variable from file and folders
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs

            if ~obj.SelfName
                if ~isempty(obj.WantedDate)
                    AnalysisOpts.CurrentRecDate=AnalysisOpts.RecDate;
                    AnalysisOpts.RecDate=obj.WantedDate;
                    CurrentCh_RecDate=AnalysisOpts.CurrentCh_RecDate;
                    AnalysisOpts.CurrentCh_RecDate=obj.WantedDate;
                end
                if ~isempty(obj.WantedCh)
                    CurrentCh=AnalysisOpts.CurrentCh;CurrentChClust=AnalysisOpts.CurrentChClust;
                    AnalysisOpts.CurrentCh=obj.WantedCh;
                    AnalysisOpts.CurrentChClust=[]; 
                end
                if isempty(AnalysisPathName);AnalysisPathName='Analysis_Data';end % just save it as analysis data
                [~,~,FileFullPath]=GenerateFileName(AnalysisOpts.FS,AnalysisOpts.DataSavePath,AnalysisPathName,...
                    AnalysisOpts.CurrentCh_Animal,AnalysisOpts.CurrentCh_RecDate,[AnalysisOpts.CurrentCh AnalysisOpts.CurrentChClust],ExtraTxt,'ext','.mat','ExtFolderName',AnalysisOpts.RecDate);
            else
                FileFullPath=ExtraTxt;
            end
        %    fprintf('\nLoading variable:%s in file:%s',VarName,FileFullPath);
            if exist(FileFullPath,'file')
                FileExist=true;
                % get size of the file 
                s = dir(FileFullPath);
                filesize=s.bytes;
                % if we are only checking if this file exists then return
                if obj.CheckifFileExistsOnly;out='';return;end

                if UseDataPointer
                    out=matfile(FileFullPath,'Writable',true);
                else
                    try
                        out=load(FileFullPath,VarName);
                        eval(['out=out.' VarName ';']);
                    catch ME
                        warning(ME.message)
                        out='';
                        FileExist=false;
                    end
                end
            else
                FileExist=false;
                warning(['File ' FileFullPath ' does not exist']);
                out='';
                filesize=[];
            end            
            if ~isempty(obj.WantedDate)
                AnalysisOpts.RecDate=AnalysisOpts.CurrentRecDate;
                AnalysisOpts.CurrentCh_RecDate=CurrentCh_RecDate;
            end
            if ~isempty(obj.WantedCh)
                AnalysisOpts.CurrentCh=CurrentCh;
                AnalysisOpts.CurrentChClust=CurrentChClust;
            end
        end

        function out=LoadVarSeries(~,AnalysisPathName,VarName,VarName2,ExtraTxt,Range) % load a series of vars and puts them in a matrix
            % range is the range of values at the end of var
            % VarName2 is the variable inside the VarName structure
            global AnalysisOpts
            fprintf('\nLoading variable series ...')
            if isempty(AnalysisPathName);AnalysisPathName='Analysis_Data';end % just save it as analysis data
            error('This function is not maintained; use LoadVarSeries2 instead')
            for i=Range
                fprintf('\nLoading variable series %i ...',i)
                ind=find(i==Range);
                [~,~,SaveFileFullPath]=GenerateFileName(AnalysisOpts.FS,AnalysisOpts.DataSavePath,AnalysisPathName,AnalysisOpts.CurrentCh_Animal,AnalysisOpts.CurrentCh_RecDate,[AnalysisOpts.CurrentCh AnalysisOpts.CurrentChClust],[ExtraTxt '_' num2str(i)],'ext','.mat');                
                eval(['temp=load(SaveFileFullPath,''' VarName '_' num2str(i) ''');']);
                if isempty(VarName2)
                    eval(['out{ind}=temp.' VarName '_' num2str(i) ';']);
                else
                    eval(['out{ind}=temp.' VarName '_' num2str(i) '.' VarName2 ';']);
                end
            end
        end
        function [out,FileExist,SaveFileFullPath]=LoadVarSeries2(obj,AnalysisPathName,VarName,VarName2,ExtraTxt,Range,varargin) % load a series of vars each saved in a seperate file  and puts them in a matrix
            % range is the range of values at the end of var
            % VarName2 is the variable inside the VarName structure
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs

            fprintf('\nLoading variable series from %i to %i',Range(1),Range(end))
            if isempty(AnalysisPathName);AnalysisPathName='Analysis_Data';end % just save it as analysis data
            
            if ~isempty(obj.WantedDate)
                AnalysisOpts.CurrentRecDate=AnalysisOpts.RecDate;
                AnalysisOpts.RecDate=obj.WantedDate;
                CurrentCh_RecDate=AnalysisOpts.CurrentCh_RecDate;
                AnalysisOpts.CurrentCh_RecDate=obj.WantedDate;
            end
            if ~isempty(obj.WantedCh)
                CurrentCh=AnalysisOpts.CurrentCh;CurrentChClust=AnalysisOpts.CurrentChClust;
                AnalysisOpts.CurrentCh=obj.WantedCh;
                AnalysisOpts.CurrentChClust=[]; 
            end            
            fprintf('\nLoading variable series %s ...',VarName)
            for i=Range
                ind=find(i==Range);
                [~,~,SaveFileFullPath]=GenerateFileName(AnalysisOpts.FS,AnalysisOpts.DataSavePath,AnalysisPathName,...
                    AnalysisOpts.CurrentCh_Animal,AnalysisOpts.CurrentCh_RecDate,[AnalysisOpts.CurrentCh AnalysisOpts.CurrentChClust],[ExtraTxt  num2str(i)],'ext','.mat','ExtFolderName',AnalysisOpts.RecDate);
                if exist(SaveFileFullPath,'file')
                    FileExist(ind)=true;
                    if obj.CheckifFileExistsOnly % if we only want to know if the file exists
                        out(ind)=nan;
                    else
                        eval(['temp=load(SaveFileFullPath,''' VarName ''');']);
                        if isempty(VarName2)
                            eval(['out(ind)=temp.' VarName  ';']);
                        else
                            eval(['out(ind)=temp.' VarName '.' VarName2 ';']);
                        end
                    end
                else
                    FileExist(ind)=false;
                end
            end
            if ~exist('out','var');out=[];FileExist=false;end %we return empty matrix
        end
        function [FileExist,SaveFileFullPath]=DeleteFileSeries(obj,AnalysisPathName,ExtraTxt,Range) % deletes a series of file
            % range is the range of values at the end of var
            % VarName2 is the variable inside the VarName structure
            global AnalysisOpts
            fprintf('\nLoading variable series from %i to %i',Range(1),Range(end))
            if isempty(AnalysisPathName);AnalysisPathName='Analysis_Data';end % just save it as analysis data
            
            if ~isempty(obj.WantedDate)
                AnalysisOpts.CurrentRecDate=AnalysisOpts.RecDate;
                AnalysisOpts.RecDate=obj.WantedDate;
                CurrentCh_RecDate=AnalysisOpts.CurrentCh_RecDate;
                AnalysisOpts.CurrentCh_RecDate=obj.WantedDate;
            end
            if ~isempty(obj.WantedCh)
                CurrentCh=AnalysisOpts.CurrentCh;CurrentChClust=AnalysisOpts.CurrentChClust;
                AnalysisOpts.CurrentCh=obj.WantedCh;
                AnalysisOpts.CurrentChClust=[]; 
            end            
            fprintf('\nLoading variable series %s ...')
            for i=Range
                ind=find(i==Range);
                [~,~,SaveFileFullPath]=GenerateFileName(AnalysisOpts.FS,AnalysisOpts.DataSavePath,AnalysisPathName,...
                    AnalysisOpts.CurrentCh_Animal,AnalysisOpts.CurrentCh_RecDate,[AnalysisOpts.CurrentCh AnalysisOpts.CurrentChClust],[ExtraTxt  num2str(i)],'ext','.mat','ExtFolderName',AnalysisOpts.RecDate);
                if exist(SaveFileFullPath,'file')
                    FileExist(ind)=true;
                    fprintf('\nFile %s is being deleted ...',SaveFileFullPath);
                    delete(SaveFileFullPath);
                else
                    fprintf(2,'\nFile %s does not exist to be deleted ...',SaveFileFullPath);
                    FileExist(ind)=false;
                end
            end
        end

        function [RecData,RecSpecs]=LoadDataFromRec(obj,DateNum,AnalysisPathName,VarName,VarName2,ExtraTxt,Range,UseDataPointer)
            global AnalysisOpts
            TrialFunc=TrialFuncs;
            % save the current Animal and Recdate so we don't mess with
            % that
            AnalysisOpts.CurrentDate=AnalysisOpts.RecDate;
            AnalysisOpts.CurrentAnimal=AnalysisOpts.Animal;
            AnalysisOpts.CurrentChSet=AnalysisOpts.Ch;
            AnalysisOpts.CurrentChArea=AnalysisOpts.ChArea;

            % first get information about this recording
            for Rec=1:length(DateNum)
                % replace current values
                AnalysisOpts.RecDate=['18' AnalysisOpts.DateSet{DateNum(Rec)} ];
                AnalysisOpts.Animal=AnalysisOpts.AnimalSet{DateNum(Rec)};
                GenerateClusterPath(AnalysisOpts.AnalysisType) % genereate the path we need for this analysis
                fprintf('\nLoading data from %s %s',AnalysisOpts.RecDate,AnalysisOpts.Animal)

                [RecSpecs(Rec).TrialTimes,RecSpecs(Rec).RuleBlockTrials,RecSpecs(Rec).ChannelInfo,...
                    RecSpecs(Rec).ChannelArea,RecSpecs(Rec).ChsSet,RecSpecs(Rec).BlockSpec]=TrialFunc.InitializeTrialFuncs; % run common trial functions
                RecSpecs(Rec).RecDate=AnalysisOpts.RecDate;RecSpecs(Rec).Animal=AnalysisOpts.Animal;

                % if we are loadign coherence data then find the range first
                if (contains(VarName,'coh','IgnoreCase',true) || contains(VarName,'corr','IgnoreCase',true)) && strcmpi(Range,'ALL') % if we are thinking about coherence
                    NChsSet=size(nchoosek(RecSpecs(Rec).ChsSet,2),1);
                    ThisRange=1:NChsSet;
                else
                    ThisRange=Range;
                end
                if contains(VarName,'coh','IgnoreCase',true)
                    AnalysisOpts.CurrentCh=RecSpecs(Rec).ChsSet(end); % should be removed in later versions
                end
                if isempty(Range) % means we are only looking at one variable
                    RecData(Rec).out=obj.LoadVar(AnalysisPathName,VarName,ExtraTxt,UseDataPointer);
                else
                    RecData(Rec).out=obj.LoadVarSeries(AnalysisPathName,VarName,VarName2,ExtraTxt,ThisRange);
                end
            end
            % change everything back to default
            AnalysisOpts.RecDate=AnalysisOpts.CurrentDate;
            AnalysisOpts.Animal =AnalysisOpts.CurrentAnimal;
            AnalysisOpts.Ch=AnalysisOpts.CurrentChSet;
            AnalysisOpts.ChArea=AnalysisOpts.CurrentChArea;
            if ~isempty(AnalysisOpts.Ch)
                AnalysisOpts.CurrentCh=AnalysisOpts.Ch(end);
            else
                AnalysisOpts.CurrentCh=AnalysisOpts.Ch;
            end
            GenerateClusterPath(AnalysisOpts.AnalysisType) % genereate the path we need for this analysis
        end
        function [FuncOut,RecSpecs]=RunFuncOnRec(obj,ClassName,FuncName,DateNum,Inputs,Nargout,varargin) % runs a general function from any class for specific recording
            % Inputs: is organzie as cell array for each recording and shoul
            % dmatch the number of recordings. If the number of inputs is
            % one then all of the recordings are going to have the same
            % input
            % Nargout : number of outputs we expect from this function
            % leave empty if it is one

            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            %% load necessary classes
            TrialFunc=TrialFuncs;
            if isempty(ClassName)
                eval(['func=@' FuncName]);
            else
                eval(['func=@ClassName.' FuncName ';']);
            end
            %% sort out inputs
            if size(Inputs,1)==1;Inputs=repmat(Inputs,[length(DateNum),1]);end
            if isempty(Nargout);Nargout=1;end

            %% save the current Animal and Recdate so we don't mess with
            %% that
            AnalysisOpts.CurrentDate=AnalysisOpts.RecDate;
            AnalysisOpts.CurrentAnimal=AnalysisOpts.Animal;
            AnalysisOpts.CurrentChSet=AnalysisOpts.Ch;
            AnalysisOpts.CurrentChArea=AnalysisOpts.ChArea;

            %% first get information about this recording and then run the function of each recording
            FuncOut=cell(length(DateNum),Nargout);
            for Rec=1:length(DateNum)
                % replace current values
                AnalysisOpts.RecDate=[AnalysisOpts.PreDateTxt AnalysisOpts.DateSet{DateNum(Rec)} ];
                AnalysisOpts.Animal=AnalysisOpts.AnimalSet{DateNum(Rec)};
                GenerateClusterPath(AnalysisOpts.AnalysisType) % genereate the path we need for this analysis
                fprintf('\nRunning %s on data from %s %s',FuncName,AnalysisOpts.RecDate,AnalysisOpts.Animal)
                % get the specs for this recording
                [RecSpecs(Rec).TrialTimes,RecSpecs(Rec).RuleBlockTrials,RecSpecs(Rec).ChannelInfo,...
                    RecSpecs(Rec).ChannelArea,RecSpecs(Rec).ChsSet,RecSpecs(Rec).BlockSpec]=TrialFunc.InitializeTrialFuncs; % run common trial functions
                RecSpecs(Rec).RecDate=AnalysisOpts.RecDate;RecSpecs(Rec).Animal=AnalysisOpts.Animal;

                %% run this function
                eval(sprintf('[FuncOut{%i,:}]=func(Inputs{%i,:});',Rec,Rec));
            end
            % change everything back to default
            AnalysisOpts.RecDate=AnalysisOpts.CurrentDate;
            AnalysisOpts.Animal =AnalysisOpts.CurrentAnimal;
            AnalysisOpts.Ch=AnalysisOpts.CurrentChSet;
            AnalysisOpts.ChArea=AnalysisOpts.CurrentChArea;
            %             if ~isempty(AnalysisOpts.Ch)
            %                 AnalysisOpts.CurrentCh=AnalysisOpts.Ch(end);
            %             else
            %                 AnalysisOpts.CurrentCh=AnalysisOpts.Ch;
            %             end
            GenerateClusterPath(AnalysisOpts.AnalysisType) % genereate the path we need for this analysis
        end

        function out=IsVarExistinFileAnalysis(~,AnalysisPathName,VarName,ExtraTxt) % checks if a var exists in a file without openning it- for analysis purposes
            global AnalysisOpts
            if isempty(AnalysisPathName);AnalysisPathName='Analysis_Data';end % just save it as analysis data
            [~,~,SaveFileFullPath]=GenerateFileName(AnalysisOpts.FS,AnalysisOpts.DataSavePath,[AnalysisPathName],AnalysisOpts.CurrentCh_Animal,AnalysisOpts.CurrentCh_RecDate,[AnalysisOpts.CurrentCh AnalysisOpts.CurrentChClust],ExtraTxt,'ext','.mat');
            if exist(SaveFileFullPath,'file')
                variableInfo = who('-file', SaveFileFullPath);
                out=ismember(VarName, variableInfo); % returns true
            else
                out=0; % if the file doesn't exist then put it to zero
            end
        end
        function out=CheckSizeofVarInFile(~,AnalysisPathName,VarName1,VarName2,ExtraTxt) % gets the size of A variable in file (VarName2 is nested in VarName1)
            global AnalysisOpts
            if isempty(AnalysisPathName);AnalysisPathName='Analysis_Data';end % just save it as analysis data
            [~,~,SaveFileFullPath]=GenerateFileName(AnalysisOpts.FS,AnalysisOpts.DataSavePath,[AnalysisPathName],AnalysisOpts.CurrentCh_Animal,AnalysisOpts.CurrentCh_RecDate,[AnalysisOpts.CurrentCh AnalysisOpts.CurrentChClust],ExtraTxt,'ext','.mat');
            if exist(SaveFileFullPath,'file')
                MatFile=matfile(SaveFileFullPath);
                if isempty(VarName2)
                    out=size(MatFile.(VarName1));
                else
                    % copy this into a temp var and then read it
                    temp=MatFile.(VarName1);
                    out=size(temp.((VarName2)));
                end
            else
                out=0; % if the file doesn't exist then put it to zero
            end
        end
        function [SigTest,p_SigTest]=StatTeston3Dmat(~,data,mu,Tail) % performs statistical test on 3D matrix; x and y as dimensions Z as samples
            [subi,subj]=arrayfun(@(x) ind2sub([size(data,1) size(data,2)],x),1:size(data,1)*size(data,2));
            [SigTest,p_SigTest]=arrayfun(@(x) ttest(squeeze(data(subi(x),subj(x),:)),mu,'tail',Tail),1:length(subi));
            SigTest=reshape(SigTest,[size(data,1) size(data,2)]);
            p_SigTest=reshape(p_SigTest,[size(data,1) size(data,2)]);
        end
        function CopyVars2AnalysisData(~,varargin)
            global AnalysisData

            if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs for options.'); end
            for i = 1:2:length(varargin)
                try
                    AnalysisData.(varargin{i}) = varargin{i+1};
                catch
                    error('Couldn''t set option ''%s''.', varargin{2*i-1});
                end
            end
        end
        function CopyVars2AnalysisOpts(~,varargin)
            global AnalysisOpts

            if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs for options.'); end
            for i = 1:2:length(varargin)
                try
                    AnalysisOpts.(varargin{i}) = varargin{i+1};
                catch
                    error('Couldn''t set option ''%s''.', varargin{2*i-1});
                end
            end
        end
        function Struct=CopyVars2Struct(~,Struct,varargin)

            if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs for options.'); end
            for i = 1:2:length(varargin)
                try
                    Struct.(varargin{i}) = varargin{i+1};
                catch
                    error('Couldn''t set option ''%s''.', varargin{2*i-1});
                end
            end
        end
        function StructTarg=CopyStructFields(~,StructSource,StructTarg,Fields)% copys variables from one struct to the other
            for Field=Fields'
                StructTarg.(Field{1})=StructSource.(Field{1});
            end
        end
        function out=CopyCell2Struct(~,CellData,StructFieldNames)% copies data from cell to each strcuture with it's field name
            for i=1:length(StructFieldNames)
                out.(StructFieldNames{i})=CellData{i};
            end
        end
        function out=ConcatinateStructs(obj,data,dim,varargin) % concatinates the strutures across recordings
            % data is the structure where each raw is set of structure for
            % a recording .
            switch dim
                case 1 % this is when we concatinate across recordings
                    out=cell(1,size(data,2));
                    for col=1:size(data,2)
                        for Rec=1:size(data,1)
                            out{col}=[out{col} data{Rec,col}];
                        end
                    end

                case 2 % this is when we concatinate across neurons
                    fields=fieldnames(data)';
                    for f=fields
                        % get the dimensions of the data first
                        if isempty(varargin)
                            fsize=size(data(1).(f{1}));
                            ReshapeDim=find(fsize>1,1,'last')+1;
                        else
                            ReshapeDim=varargin{1};
                        end
                        temp=arrayfun(@(x) data(x).(f{1}),1:length(data),'UniformOutput','false');
                        out.(f{1})=obj.ReshapeCell2Mat(temp,ReshapeDim);
                    end
            end
        end
        function [I,Ish,Iz]=CalculateInformation(~,Stim,Resp) % calculates information between stimulus and response
            %             infopt.method='gs';
            %             infopt.bias='gsb';
            %         %    infopt.method='dr';
            %         %    infopt.bias='naive';
            %             infopt.btsp = 100;
            %             infopt.verbose=0;
            %             nb=3;
            %
            %        %     Resp=Resp+0.01*rand(size(Resp)); % add some noise to make sure this is correct
            %             Resp_bin=bin_eqpop(Resp,nb);
            %             [BuildRespMat,BuildNt]=buildr(Stim,Resp_bin);
            %             infopt.nt=BuildNt;
            %             I=information(BuildRespMat,infopt,'I');
            %%%%%%%%%%%%%%%%%%%% another way of calculating information with a
            %%%%%%%%%%%%%%%%%%%% different toolbox
            nb=3;
            infopt.btsp = 100;
            %   Resp=Resp+0.01*rand(size(Resp));
            Resp_bin=bin_equidist(Resp,nb);
            N=length(Stim);
            I=mutInfo(Resp_bin,Stim);
            Ish=arrayfun(@(x) mutInfo(Resp_bin(randperm(N,N)),Stim(randperm(N,N))),1:infopt.btsp);
            I(isinf(I))=nan;
            Iz=zscore([I Ish]);Iz=Iz(1); % get the zscore value of information against shuffle
        end
        function Selctivity=SelectivityIndex(obj,Stim,Resp,varargin) % calculates selectivity to stimulus and response
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            %   Resp=abs(Resp); % make sure this is correct
            switch obj.SelectivityIndexMethod
                case 'norm'  % take the norm of response
                    Selctivity=norm(Resp);
                case 'mean'  % take the mean of response
                    MeanResp=arrayfun(@(x) mean(Resp(Stim==x)),unique(Stim));
                    Selctivity=mean(MeanResp);
                case 'morph' % calculate morph selectivity
                    MeanResp=arrayfun(@(x) mean(Resp(Stim==x)),unique(Stim));
                    NB=length(MeanResp);
                    Rmax=max(MeanResp(:));
                    if NB==1;Selctivity=MeanResp;else
                        Selctivity=(NB -(sum(MeanResp(:))/Rmax))/(NB-1);
                    end
            end
        end
        function [SI,zSI,SIshuffle]=CalSelectivityIndex(obj,Stim,Resp,Nshuffle)
            n=length(Resp);
            SI=obj.SelectivityIndex(Stim,Resp);
            SIshuffle=arrayfun(@(x) obj.SelectivityIndex(Stim(randperm(n,n)),Resp),1:Nshuffle);
            SIzscore=zscore([SI SIshuffle]);
            zSI=SIzscore(1);
        end
        function MI=ModulationIndex(~,Stim,Resp) % calculates the modulation of FR based on the stimulus(color or shape) (cit Paninchenlo Buschman 2020)

            UniqStim=unique(Stim);
            N=length(UniqStim);
            rc=arrayfun(@(x) mean(Resp(Stim==x)),UniqStim);
            Zc=rc./sum(rc);
            MI=nansum(arrayfun(@(x) Zc(x)*log(N*Zc(x)),1:N))/log(N);
            %   if isnan(MI) || isinf(MI);MI=0;end
        end
        function [MI,zMI,MIshuffle]=CalModulationIndex(obj,Stim,Resp,Nshuffle)% calculate modulation index
            n=length(Resp);
            MI=obj.ModulationIndex(Stim,Resp);
            MIshuffle=arrayfun(@(x) obj.ModulationIndex(Stim(randperm(n,n)),Resp(randperm(n,n))),1:Nshuffle);
            MIzscore=zscore([MI MIshuffle]);
            zMI=MIzscore(1);
        end
        function [TargColBin,BinColWheel]=BinColorSpace(~,nb,TargCol) % devides colorspace in nb bins and finds the bin of target color
            ColWheel=0:0.001:2*pi+0.001;
            BinColWheel=bin_equidist(ColWheel,nb);
            for TC=TargCol
                TargColBin(TC==TargCol)=BinColWheel(find(TC<=ColWheel,1,'first'));
            end
        end
        function out=BinData(~,binsiz,nbin,str,stp)
            if binsiz & nbin;warning('you can only specify  binsize or number of bins ... ignoring bin size value');end
            if ~isempty(nbin)
                binsiz=(stp-str)/nbin;
            end
            out=str:binsiz:stp;
        end
        function [out] = CalAngle(obj,data,varargin)
            %CALANGLE Calculates angle in radian and degree corrects for jumps
            % Deg=1 reports in degree
            obj=obj.ParseParams(varargin) ;      %%Process optional inputs
            out=angle(data);
            out(out<0)=out(out<0)+2*pi;
            if obj.Deg
                out=180*out/pi;
            end
            % test
            %plot(obj.CalAngle((cos(0:0.01:2*pi)+1i*sin(0:0.01:2*pi)),1));
        end
        function out=DetermineDate(~,DateNum)
            global AnalysisOpts
            if ~ischar(DateNum)
                out=[AnalysisOpts.DateSet{DateNum}]; % revert to this date
            else
                out=DateNum; % revert to this date
            end
        end
        function Animal=DetermineDateAnimal(obj,Date)
            global AnalysisOpts
            DateNum=obj.DetermineDateNum(Date);
            Animal=AnalysisOpts.AnimalSet{DateNum};
        end
        function [out,uniqdata]=CategorizeData(~,data) % bins unqiue data points in the data and assigns them numbers (to use with ANOVA)
            uniqdata=sort(unique(data));
            if size(uniqdata,1)>1;uniqdata=uniqdata';end
            out=zeros(size(data));
            cnt=1;
            for d=uniqdata
                out(data==d)=cnt;
                cnt=cnt+1;
            end
        end
        function [out,levels,levelVal]=GenCategoricalVars(obj,data) % takes categorical data and generates seperats category levels(to use with costum GLM)
            % data has to be a vector
            if size(data,1)==1;data=data';end

            [levels,levelVal]=obj.CategorizeData(data);
            out=cell2mat(arrayfun(@(x) levels==x,unique(levels),'UniformOutput',false)');

        end
        function R2=CalR2(~,X,Xhat) % calculate coefficient of determiantion R2
            R2=1-(sum((X(:)-Xhat(:)).^2)/sum((X(:)-mean(X(:))).^2));
            %     residuals=X-Xhat;
            %    R2 = 1-nanvar(residuals(:))./nanvar(X(:));
            %  R2 = sqrt(sum((X(:)-Xhat(:)).^2));
            if isinf(R2);R2=nan;end
        end
        function CatTuneIndex=CalGLMCategoryTuningIndex(obj,BetaMLSing,ML) % calculates category tuning index GLM morph levels
            % BetaMLSing: is the beta value for each morph level MLxTime
            % ML Morphlevel values 
            nTim=size(BetaMLSing,2);
            CatML=obj.CategorizeMorphlevel(ML);
            CatML1=CatML==1;
            CatML2=CatML==2;
            
            for t=1:nTim % loop on time points
                WithinDistCat1=mean(pdist(BetaMLSing(CatML1,t),'euclidean'));
                WithinDistCat2=mean(pdist(BetaMLSing(CatML2,t),'euclidean'));
                WithinDist=mean([WithinDistCat1 WithinDistCat2]); %WCD
                BetweenDist=pdist2(BetaMLSing(CatML1,t),BetaMLSing(CatML2,t),'euclidean'); %BCD  
                BetweenDist=mean(BetweenDist(:));
                CatTuneIndex(t)=(BetweenDist-WithinDist)/(BetweenDist+WithinDist);
            end
        end
        function [R2time,R2trl,R2TimePoint]=CalR2time(~,X,Xhat) % calculate R2 when there is time in the trial
            % defined by 1-SSres/SStot=1-var(X-Xhat)/var(X) sumed over
            % trials

            % R2=1-(sum(((X-Xhat,2)).^2)/sum(mean(X-repmat(mean(X,1),[size(X,1) 1]),2).^2));
            R2time=1-(sum(var(X-Xhat,0,1))/sum(var(X,0,1)));
            R2trl =1-(sum(var(X-Xhat,0,2))/sum(var(X,0,2)));
            % calculate R2 per each time point 
            %R2TimePoint=max(arrayfun(@(x) 1-((var(X(:,x)-Xhat(:,x),0,1))/(var(X(:,x),0,1))),1:size(X,2)));
            R2TimePoint = 1 - sum((X - Xhat).^2, 'all') / sum((X - mean(X)).^2,"all");

        end

        function [out,out2]=FindFieldinStruct(~,S,Name) % finds specific fields in a structure
            FN=fieldnames(S)';
            temp=contains(FN,Name,'IgnoreCase',true);
            out=FN(temp);
            out2=FN(~temp);
        end
        function err = MeanSquaredError(~,resp, target) %MSE
            err = nanmean((resp(:) - target(:)).^2);
        end
        function err = MeanSquaredErrorTrial(obj,resp, target) %MSE for set of trials
            % err=mean(arrayfun(@(y) obj.MeanSquaredError(resp(y,:),target(y,:)),1:size(target,1),'UniformOutput',1));
            err=mean(arrayfun(@(y) obj.MeanSquaredError(resp(:,y),target(:,y)),1:size(target,2),'UniformOutput',1));

        end
        function err = SumSquaredError(~,resp, target) % SSE
            err = sum((resp(:) - target(:)).^2);
        end
        function err = EuclideanDistance(~,resp, targ)
            err = sqrt(sum((resp(:) - targ(:)).^2));
        end

        function MLrad=ConvertML2rad(~,ML)
            MLrad=ML*2*pi./200;
        end

        function CPD=CallCPD(obj,Y,Yfull,Yred) % calculates percent of coefficient of partial determination
            % Y is true FR matrix
            %Yfull firing rate estimated with full model
            %Yred firing rate estimated using reduced model with a specific
            %factor
            % CPDi= SSE(x−i) − SSE(x)/SSE(x−i); SSE(x-i): sum squared error
            % for full model that lacks factor i
            % from Chiang, Feng-Kuei, and Joni D Wallis. “Neuronal Encoding in Prefrontal Cortex during Hierarchical Reinforcement Learning�? 30, no. 8 (n.d.): 12.

            NTim=size(Y,2);
            for t=1:NTim % calculate CPD for each time point of the sliding window
                SSExi=obj.SumSquaredError(Yred(:,t),Y(:,t));
                SSEx=obj.SumSquaredError(Yfull(:,t),Y(:,t));
                CPD(t)=(SSExi-SSEx)*100/SSExi;
            end
        end


        function out=GetTrlsFromData(~,data,Trls) % gets trials from a struct using the dimension that is maximum
            sz=size(data);
            [~,MaxDim]=max(sz);% find the dimension of number of trials
            inds=arrayfun(@(x) 1:x,sz,'UniformOutput',0);
            inds{MaxDim} = Trls;
            out=data(inds{:});
        end
        function [TimeAxis,StrTimeAxis]=GenerateTimeAxis(obj,PSTHTimRef,varargin) % generates time axis based on parameters set in AnalysisOpts
            global AnalysisOpts

            if nargin==2 % use timing from AnalysisOpts
                Bin=AnalysisOpts.PopulationAna.PSTHbin*0.001;
                BinShift=AnalysisOpts.SpkParams.PSTH_BinShift;
                PeriodLength=AnalysisOpts.SpkParams.PeriodLength;
                BaselineDelay=AnalysisOpts.SpkParams.BaselineDelay;
            else % use timing in varargin
                Bin=varargin{1}.PopulationAna.PSTHbin*0.001;
                BinShift=varargin{1}.SpkParams.PSTH_BinShift;
                PeriodLength=varargin{1}.SpkParams.PeriodLength;
                BaselineDelay=varargin{1}.SpkParams.BaselineDelay;
            end
            numbins=(PeriodLength-Bin)/BinShift+1;
            %% calculate PSTH
            TimeAxis=arrayfun(@(x) Bin+BinShift*(x-1),1:numbins)-BaselineDelay;
            if ~exist('PSTHTimRef','var')
                StrTimeAxis=TimeAxis-Bin;
            else
                switch PSTHTimRef
                    case 'leading' % start of spike count
                        StrTimeAxis=TimeAxis-Bin;
                    case 'trailing' % this is what we are comuting
                        StrTimeAxis=TimeAxis;
                    case 'centered'  % centered
                        StrTimeAxis=TimeAxis-Bin/2;
                end
            end
        end
        function [duplicate_indices,nonduplicate_indices]=FindIndDuplicateValues(~,data) % finds indices of duplicate values in a vector

            [~, w] = unique( data, 'stable' );
            duplicate_indices = setdiff( 1:numel(data), w );
            nonduplicate_indices=setdiff( 1:numel(data),duplicate_indices);
        end

        function data=ReplaceEmptywithNaN(~,data)% replaces empty values of a matrix with NaN
            if isempty(data);data=NaN;end
        end
        function [DateNum,Animal]=GetDateNum(~,Date)  % gets the number of the date for the specifed text of the date
            global AnalysisOpts

            if iscell(Date)
                DateNum=cellfun(@(x) find(strcmp(AnalysisOpts.DateSet,x)),Date);
                Animal=AnalysisOpts.AnimalSet(DateNum);
            else
                DateNum= find(strcmp(AnalysisOpts.DateSet,Date)) ;
                Animal=AnalysisOpts.AnimalSet(DateNum);
            end
        end
        function ReportAnalysisDetails(obj,varargin)
            global AnalysisOpts
            if ~isempty(varargin)
                OptTxt=varargin{1};
            else
                OptTxt=[];
            end
            fprintf(2,'\n*************************************************\n')
            fprintf(2,'Anlyzing data:Project:%s, %s,\nAnimal:%s \nTrlSpkTimeFieldName:%s,\nSpkCntStartFieldName:%s,\nPeriodLength=%.2f,\nBaselineDelay=%.2f,\nPSTH Bin=%i,\nProcessing Step=%i,%s\nPairNum=%i,\nChannel:%s,%s',AnalysisOpts.Project,AnalysisOpts.AnalysisType,AnalysisOpts.Animal,...
                AnalysisOpts.TrlSpkTimeFieldName,AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.SpkParams.PeriodLength,AnalysisOpts.SpkParams.BaselineDelay,AnalysisOpts.PopulationAna.PSTHbin,AnalysisOpts.ProcessingStep,...
                AnalysisOpts.PopulationAna.ProcessingStepNames{AnalysisOpts.ProcessingStep+1},AnalysisOpts.PairNum,obj.ConvMat2Char(AnalysisOpts.Ch_2look),OptTxt);
            fprintf(2,'\n*************************************************\n')

        end
        function out=CovertMorphLvl2Angle(~,data) % convert the morph level of shapes and colors that change from (0 to 200) to (0 to 2pi)
            out=data*2*pi./200;

        end
        function pval=CalpValShuffle(obj,ShuffleData,Data) % outdated 
            for i=1:size(ShuffleData,2)
                [n,xout]=hist([ShuffleData(:,i) ; Data(i)]);
                pval(i)=sum(n(xout>=Data(i)))/size([ShuffleData(:,i) ; Data(i)],1);
            end
        end
        function [pval,Zval]=CalPval(obj,ShuffleData,Observed,tail) % calculates p-value and Z value
            if strcmp(tail,'one')
                NullDist=[ShuffleData Observed];nNull=length(NullDist);
                pval=sum(NullDist>=Observed)/nNull;
            elseif  strcmp(tail,'two')
                MeanObsv=abs(Observed-mean(ShuffleData));
                NullDist=[ShuffleData-mean(ShuffleData) MeanObsv];nNull=length(NullDist);
                pval=sum(abs(NullDist)>=MeanObsv)/(2*nNull);
            end
                % calculate Z value too 
                Z=@(x,y) (x-mean(y))/std(y);
                Zval=Z(Observed,ShuffleData);
        end
        function [pvals, clustMass, clustIdx,nullMass] = ClusterMassCorrection(obj,x,nIter,thresh,dependent_samples,varargin)
            
            nT = size(x,1);
            
            %thresh = 0.05;
            if dependent_samples==0; TTestFunc='ttest2';else;TTestFunc='ttest';end
            
            if ~isempty(varargin) && strcmpi(varargin{1},'tail')
                eval(['[h, p, ~, stat] =' TTestFunc '(x,0,''tail'',varargin{2});']);
            else
                eval(['[h, p, ~, stat] =' TTestFunc '(x);']);
            end
            
            h(isnan(h)) = 0;
            p(isnan(p)) = 1;
            p = p <= thresh;
            
            [clustIdx, ~, clustMass] = obj.getClust(p,stat.tstat);
            
            nullMass = nan(nIter,1);
            for iIter = 1:nIter
                tmpX = x;
                idx = rand(nT,1)>.5;
                tmpX(idx,:,:) = -1.*tmpX(idx,:,:);
                
                if ~isempty(varargin) && strcmpi(varargin{1},'tail')
                    eval(['[tmpH, tmpP, ~, tmpStat] =' TTestFunc '(tmpX,0,''tail'',varargin{2});']);
                else
                    eval(['[tmpH, tmpP, ~, tmpStat] =' TTestFunc '(tmpX);']);
                end
                
                tmpH(isnan(tmpH)) = 0;
                tmpP(isnan(tmpP)) = 1;
                tmpP = tmpP <= thresh;
                                
                [~, ~, tmpMass] = obj.getClust(tmpP,tmpStat.tstat);
                if isempty(tmpMass)
                    nullMass(iIter) = 0;
                else
                    nullMass(iIter) = max(abs(tmpMass));
                end
            end
            
            pvals = sum(nullMass>=abs(clustMass)')./nIter;
            
        end
        function [clustIdx, nClusts, clustMass] = getClust(obj,h,statistic)
            
            if size(h,3)>1;h=squeeze(h);statistic=squeeze(statistic);end
            [clustIdx, nClusts] = bwlabel(h);
            clustMass = nan(nClusts,1);
            for iClust = 1:nClusts
                clustMass(iClust) = sum(statistic(clustIdx==iClust));
            end
        end
        
        function [pvals, clustMass, clustIdx,nullMass] = ClusterMassCorrection_permutation2D(obj,Shuffle,observed,thresh,two_sided,varargin)
            error('This function is not correct');
            if two_sided;tail='two';else;tail='one';end
            nT = size(Shuffle,1);
            nTim=size(Shuffle,3);
            Ntrl=size(Shuffle,2);
            NItems=nTim*Ntrl;
            
            [row,col] = ind2sub([Ntrl nTim],1:NItems);
            [pval,zval]=arrayfun(@(x) obj.CalPval(Shuffle(:,row(x),col(x))',observed(1,row(x),col(x)),tail),1:NItems);
            pval=reshape(pval,[Ntrl nTim]);zval=reshape(zval,[Ntrl nTim]);
            pval(isnan(pval)) = 1;
            p = pval <= thresh;
            
            % find related clusters
            [clustIdx, ~, clustMass] = obj.getClust(p,zval);
            
            nullMass = nan(nT,1);
            for iIter = 1:nT % iterate on shuffle trials and include the observed 
                ObservedTemp=Shuffle(iIter,:,:); % take shuffle
                ShuffleTemp=[Shuffle(setdiff(1:nT,iIter),:,:)];%;observed'];
                %ShuffleTemp=cat(1,Shuffle(setdiff(1:nT,iIter),:,:),observed);
                [pval_Temp,zval_Temp]=arrayfun(@(x) obj.CalPval(ShuffleTemp(:,row(x),col(x))',ObservedTemp(1,row(x),col(x)),tail),1:NItems);
                pval_Temp=reshape(pval_Temp,[Ntrl nTim]);zval_Temp=reshape(zval_Temp,[Ntrl nTim]);
                
                pval_Temp(isnan(pval_Temp)) = 1;
                tmpP = pval_Temp <= thresh;
                                
                [~, ~, tmpMass] = obj.getClust(tmpP,zval_Temp);
                if isempty(tmpMass)
                    nullMass(iIter) = 0;
                else
                    nullMass(iIter) = mean(abs(tmpMass));
                end
            end           
            pvals = sum(nullMass>=abs(clustMass)')./size(ShuffleTemp,1);            
        end
        
        function [obs_clust_p, clustMass, clustIdx,max_clustMass_dist,clusters,statsummery] = ClusterMassCorrection_permutation(obj,Shuffle,observed,p_thresh,two_sided,varargin)
            
             obj=obj.ParseParams(varargin) ;      %%Process optional inputs

            if two_sided;tail='two';else;tail='one';end
            %% Given your observed and shuffled values over time -- lets call it NullDist which is a matrix of (Nrand + 1) x T 
            %where the first row is your observed and the remaining are shuffled.
            %You first need to find the threshold for defining a moment in time as significant.
            %This can either be done non-parametrically (i.e., take the (1-p)th percentile at each moment in time) or
            %parametrically (i.e., mean + X*standard deviation where X is taken from norminv(1-p, 0, 1)). 
            %Either way you end up with a 1xT vector that is your threshold.
            Nrand = size(Shuffle,1);
            T=size(Shuffle,3);
            Ntrl=size(Shuffle,2);
            NItems=T*Ntrl;
            % obeserved and shiffle are organized as Nrep*NTrl*NTim
            % get the percentile of the null dustribution.
            NullDist=cat(1,observed,Shuffle);
            ZNullDist=zscore(NullDist,0,1);
          %  p_threshold=prctile(NullDist,(1-p_thresh)*100,1); 
            p_threshold=mean(NullDist,1)+norminv(1-p_thresh,0,1)*std(NullDist,1,1);

            %% Now, use this to find all moments in time that are greater than threshold (threshold = D > repmat(p_threshold, (Nrand+1), 1)). 
            threshold = NullDist > repmat(p_threshold, (Nrand+1), 1,1);

            %% Sum up the value of contiguous clusters within each time series. This can either be on the z-transformed data (zD) or the raw data (D).
            % And it can either be done on the values (zD or D) or on just the part above threshold (e.g., zD(D > repmat(threshold, (Nrand+1), 1))).     
            
            [clustIdx, ~, clustMass]=arrayfun(@(x)  obj.getClust(threshold(x,:,:),ZNullDist(x,:,:)),1:(Nrand+1),'UniformOutput',0);
            %% This will now give you a list of the cluster sizes for each row in your data. Usually I store this as a cell array. Let's call it clustMass
            % Now you want to take the maximum cluster size across time for each of the shuffles. This will give you a max_cluster_size=@cellfun(@max, clustMass(2:end))that is Nrandx1 (edited)
            
            max_clustMass=cellfun(@max, clustMass(2:end),'UniformOutput',0);
            emptyClustmass=cellfun(@isempty,max_clustMass);
            max_clustMass_dist=[cell2mat(max_clustMass(~emptyClustmass)) zeros(1,sum(emptyClustmass))];

            %% Now calculate the p-value of each of the observed clusters (index by i) by comparing to the
            %distribution of max_cluster_size: obs_clust_p(i) = 1-sum(cluster_size{1}(i)>max_cluster_size)/Nrand
            obs_clust_p=arrayfun(@(x) 1-sum(clustMass{1}(x)>max_clustMass_dist)/Nrand,1:length(clustMass{1}));
            clustIdx=clustIdx{1};clustMass=clustMass{1};

            %  generate a summery for stat test
            clusters=arrayfun(@(x) find(clustIdx==x),unique(clustIdx(clustIdx~=0)),'UniformOutput',0);

            % create a summery matrix where for each time point in the clusters we assign a p-value
            if ~isempty(clusters)
                statsummery=arrayfun(@(x) obs_clust_p(x)*ones(size(clusters{x})),1:length(clusters),'UniformOutput',0);
            else
                statsummery=[];
            end
            
            if obj.ShowClustCorrectionPlot
                ClusterMassCorrectionCheck % if you want to plot values to double check
            end

        end
       
        function [obs_clust_p, clustMass, clustIdx,max_clustMass_dist,clusters,statsummery] = ClusterMassCorrection_permutationTwoTail(obj,Shuffle,observed,p_thresh,two_sided,varargin)
            
             obj=obj.ParseParams(varargin) ;      %%Process optional inputs

            if two_sided;tail='two';else;tail='one';end
            %% Given your observed and shuffled values over time -- lets call it NullDist which is a matrix of (Nrand + 1) x T 
            %where the first row is your observed and the remaining are shuffled.
            %You first need to find the threshold for defining a moment in time as significant.
            %This can either be done non-parametrically (i.e., take the (1-p)th percentile at each moment in time) or
            %parametrically (i.e., mean + X*standard deviation where X is taken from norminv(1-p, 0, 1)). 
            %Either way you end up with a 1xT vector that is your threshold.
            Nrand = size(Shuffle,1);
            T=size(Shuffle,3);
            Ntrl=size(Shuffle,2);
            NItems=T*Ntrl;
            % obeserved and shiffle are organized as Nrep*NTrl*NTim
            % get the percentile of the null dustribution.
            NullDist=cat(1,observed,Shuffle);
            ZNullDist=zscore(NullDist,0,1);
            p_thresholdUpper=prctile(NullDist,(1-p_thresh)*100,1); 
            p_thresholdLower=prctile(NullDist,(p_thresh)*100,1); 
            
         %   p_thresholdUpper=mean(NullDist,1)+norminv(1-p_thresh,0,1)*std(NullDist,1,1);
         %   p_thresholdLower=mean(NullDist,1)-norminv(p_thresh,0,1)*std(NullDist,1,1);
            %% Now, use this to find all moments in time that are greater than threshold (threshold = D > repmat(p_threshold, (Nrand+1), 1)). 
            threshold = NullDist > repmat(p_thresholdUpper, (Nrand+1), 1,1) | NullDist < repmat(p_thresholdLower, (Nrand+1), 1,1); 

            %% Sum up the value of contiguous clusters within each time series. This can either be on the z-transformed data (zD) or the raw data (D).
            % And it can either be done on the values (zD or D) or on just the part above threshold (e.g., zD(D > repmat(threshold, (Nrand+1), 1))).     
            
            [clustIdx, ~, clustMass]=arrayfun(@(x)  obj.getClust(threshold(x,:,:),ZNullDist(x,:,:)),1:(Nrand+1),'UniformOutput',0);
            %% This will now give you a list of the cluster sizes for each row in your data. Usually I store this as a cell array. Let's call it clustMass
            % Now you want to take the absolute value of  cluster size across time for each of the shuffles. This will give you a max_cluster_size=@cellfun(@max, clustMass(2:end))that is Nrandx1 (edited)
            clustMass=cellfun(@abs,clustMass,'UniformOutput',0);
            max_clustMass=cellfun(@max, clustMass(2:end),'UniformOutput',0);
            emptyClustmass=cellfun(@isempty,max_clustMass);
            max_clustMass_dist=[cell2mat(max_clustMass(~emptyClustmass)) zeros(1,sum(emptyClustmass))];
 
            %% Now calculate the p-value of each of the observed clusters (index by i) by comparing to the
            %distribution of max_cluster_size: obs_clust_p(i) = 1-sum(cluster_size{1}(i)>max_cluster_size)/Nrand
          %  obs_clust_p=arrayfun(@(x) 1-(sum(clustMass{1}(x)>=max_clustMass_dist)+1)/(Nrand+1),1:length(clustMass{1}));
            obs_clust_p=arrayfun(@(x) (sum(max_clustMass_dist>=clustMass{1}(x))+1)/(Nrand+1),1:length(clustMass{1}));
           
            clustIdx=clustIdx{1};clustMass=clustMass{1};

            %  generate a summery for stat test
            clusters=arrayfun(@(x) find(clustIdx==x),unique(clustIdx(clustIdx~=0)),'UniformOutput',0);

            % create a summery matrix where for each time point in the clusters we assign a p-value
            if ~isempty(clusters)
                statsummery=arrayfun(@(x) obs_clust_p(x)*ones(size(clusters{x})),1:length(clusters),'UniformOutput',0);
            else
                statsummery=[];
            end
            
            if obj.ShowClustCorrectionPlot
          %      ClusterMassCorrectionCheckTwoTail % if you want to plot values to double check
            end

        end
       
        function out=ComparePval(~,pval,Alpha,n) % compares significance level
            if isempty(n);n=1;end
            out=pval(:)<(Alpha/n);
            out=reshape(out,size(pval));
        end
        function nSig=GetPvalClusterTime(obj,Time,Cluster,Pval,Alpha,StrTim,EndTim) %connects the number of significant point based on cluster and time
             % Time is the time axis
             % Cluster is the cluster correction 
             % Pval is from the cluster correction 
             % Alpha is the Alpha value from the cluster correction 
             % TimeAfter is the threshold of time that we want to look at the significant points after that

             TimeX=cellfun(@(x) Time(x),Cluster,'UniformOutput',0);
             [~,~,XInc]=cellfun(@(x) obj.RemoveEntryFromVec(x,x(x<=StrTim | x>=EndTim)),TimeX,'UniformOutput',0);
             X=arrayfun(@(x) Cluster{x}(XInc{x}),1:length(XInc),'UniformOutput',0);
             P=arrayfun(@(x) Pval{x}(XInc{x})',1:length(XInc),'UniformOutput',0);
             nSig=sum(cell2mat(P)<Alpha);


        end
        function PhenoClusterPlot(obj,data,MarkerInd) % clusters the data based on phenocluster and plots them
            global AnalysisOpts
            FigParams=fig_params;
            opts.Marker='*';
            opts.MarkerSize=5;
            opts.MarkerSet={'*','d','o','s','.'};
            data=double(data); % if the data is not scalar
            % cluster the data
            [idx, idx_knn, ovr_q]= PhenoCluster(data);
            if isempty(idx);fprintf('\n Phenoclust found no cluster in the data');return;end
            ClustInds=idx(:,1);
            Clusts=unique(ClustInds);
            nClusts=length(Clusts);
            Col=FigParams.getColorPalet(nClusts);
            npoints=size(idx,1);

            % do dimensionality reduction with TSNE
            fprintf('\nUsing TSNE to project the data ....')
            tnsedata=tsne(data,'Algorithm','barneshut','NumPCAComponents',10,'NumDimensions',2);
            % go through each community and plot them
            subplot(2,2,2)
            if isempty(MarkerInd)
                hold on
                arrayfun(@(x) plot(tnsedata(ClustInds==x,1),tnsedata(ClustInds==x,2),'Marker',....
                    opts.Marker,'Color',Col(x,:),'MarkerSize',opts.MarkerSize),Clusts,'UniformOutput',0);
            else % plot each point by its own marker
                hold on
                arrayfun(@(x) plot(tnsedata(x,1),tnsedata(x,2),'Marker',....
                    opts.MarkerSet{MarkerInd(x)},'Color',Col(ClustInds(x),:),'MarkerSize',opts.MarkerSize),1:npoints,'UniformOutput',0);

%                 arrayfun(@(x) plot3(tnsedata(x,1),tnsedata(x,2),tnsedata(x,3),'Marker',....
%                     opts.MarkerSet{MarkerInd(x)},'Color',Col(ClustInds(x),:),'MarkerSize',opts.MarkerSize),1:npoints,'UniformOutput',0);
%                 xlabel('TSNE1');ylabel('TSNE2');zlabel('TSNE3');

                subplot(2,2,3)
                for i=1:5;ThisArea=MarkerInd==i;Sum(i,:)=arrayfun(@(x) sum(ClustInds(ThisArea)==x)/sum(ThisArea),Clusts);end
                plot(Sum');
                legend(AnalysisOpts.AreaNames)
                xlabel('ClusterNumber');ylabel('Proportion of Neurons');
            end
            FigParams.FormatAllAxesFig(gcf)

            MeanClust=arrayfun(@(x) mean(data(ClustInds==x,:),1)',Clusts,'UniformOutput',0)';


        end
        function OrthoVector=GetOrthogonalVector2Plane(~,PlaneVec) % calculates orthogonal vector to a plane % should be equivalent to cross command
            % @PlaneVec is two vectors specifying a plane
            OrthoVector=[PlaneVec(1,2)*PlaneVec(2,3)-PlaneVec(1,3)*PlaneVec(2,2);...
                         PlaneVec(1,3)*PlaneVec(2,1)-PlaneVec(1,1)*PlaneVec(2,3);
                         PlaneVec(1,1)*PlaneVec(2,2)-PlaneVec(1,2)*PlaneVec(2,1)];
        end
        function [angd,angr,CosTheta,DotVects]=GetAngleBetVectors(obj,v1,v2) % calculates angle between two vectors 
            
            DotVects =dot(v1,v2);                
          %  CosTheta=DotVects/(norm(v1)*norm(v2));    
          CosTheta = max(min(DotVects/(norm(v1)*norm(v2)),1),-1);
% 
%           if length(v1)==3
%               angd = atan2d(norm(cross(v1,v2)),dot(v1,v2));
%               angr = atan2(norm(cross(v1,v2)),dot(v1,v2));
%           else
              angd=(acosd(CosTheta));
              angr=(acos(CosTheta));
%           end

        end

        function [Vec,Inds,IndInc]=RemoveEntryFromVec(~,Vec,Entry) % removes entries from a Vector
            if iscell(Vec)
                for i=1:length(Vec)
                    if size(Vec{i},2)>1;Vec{i}=Vec{i}';end
                    if isnan(Entry)
                        Inds{i}=sum(cell2mat(arrayfun(@(x) isnan(Vec{i}),Entry,'UniformOutput',0)),2);
                    else
                        Inds{i}=sum(cell2mat(arrayfun(@(x) Vec{i}==x,Entry,'UniformOutput',0)),2);
                    end
                    if ~isempty(Inds{i})
                        Vec{i}=Vec{i}(~Inds{i})';
                    end
                    IndInc{i}=find(~Inds{i});
                end
            else                
                if size(Vec,2)>1;Vec=Vec';end
                if isnan(Entry)
                    Inds=sum(cell2mat(arrayfun(@(x) isnan(Vec),Entry,'UniformOutput',0)),2);
                else
                    Inds=sum(cell2mat(arrayfun(@(x) Vec==x,Entry,'UniformOutput',0)),2);
                end
                if ~isempty(Inds)
                    Vec=Vec(~Inds)';
                    IndInc=(~Inds);
                else
                    IndInc=logical(ones(size(Vec)));
                end
             end
        end
        function out=CategorizeMorphlevel(~,ML)
            if sum(ismember(ML,[1 2 -1]));warning('Skipping Categorization alread category data');out=ML;return;end
            out=zeros(size(ML));
            out(ML>150 | ML<50)=1;
            out(ML<150 & ML>50)=2;
            out(ML==150 | ML==50)=-1;
        end
        function out=ConvMat2Char(~,Mat) % converts mat into char with ',' (to be used for cluster)
            if ~iscell(Mat)
                if size(Mat,1)>1;Mat=Mat';end
                if length(Mat)>100
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
        function cellArray=CovertDouble2CellStr(~,doubleVector) % converts a double vector into  a cell with str values 
            % Example vector of double values
         %   doubleVector = [1.23, 4.56, 7.89, 10.11];

            % Initialize a cell array to store character arrays
            cellArray = cell(size(doubleVector));

            % Convert each double value to a character array and store in the cell array
            for i = 1:numel(doubleVector)
                cellArray{i} = num2str(doubleVector(i));
            end        

        end
        function Congruency=DetermineStimCongruency(obj,Dim1,Dim2) % determines if the stimulus is congrunet or incongruent 
           Congruency=NaN*ones(1,length(Dim1));
            for i=1:length(Dim1)
                %@Dim input can be shape or color category or morphlevel
                if sum(Dim1(i)==0:10:180)  % we dealing with ML
                    % convert ML to category
                    Dim1(i)=obj.CategorizeMorphlevel(Dim1(i));
                    Dim2(i)=obj.CategorizeMorphlevel(Dim2(i));
                end

                % find out congruency level
                Congruency(i)=~(Dim1(i)==Dim2(i));
                if Dim1(i)==-1 || Dim2(i)==-1
                    Congruency(i)=-1;
                end
            end
%             a=histc(Dim2,[-1 1 2]);
%             fprintf(2,'\n-1:%0.4f 0:%0.4f 1:%0.4f',a(1)/sum(a),a(2)/sum(a),a(3)/sum(a))
        end
        function FactorData=UpdateFactorDataCongruency(obj,FactorData,ClassifierOpts) % updates the congruency on the factor data
            for Neu=1:length(FactorData)
                for i=1:length(FactorData(Neu).AllFactors)
                    FactorData(Neu).AllFactors{i}(:,ClassifierOpts.Congruency.Ind)=...
                        obj.DetermineStimCongruency(FactorData(Neu).AllFactors{i}(:,ClassifierOpts.ShapeCat.Ind),FactorData(Neu).AllFactors{i}(:,ClassifierOpts.ColorCat.Ind));
                end
            end
        end
        function y=RandSamplePopulation(obj,population,k) % randomly samples from a population . k<=length(population)
           opts.IgnoreErrors=0; % don't throw error if population is empty or  k>length(population)
            if k>length(population) | isempty(population)
                if opts.IgnoreErrors
                    y=population;
                else
                error('RandSample can not be used when k>length(population)')
                end
            elseif length(population)==1 & k==1
                y=population;
            else
                y=randsample(population,k);
            end
        end
        function [sampData,sampInds]=RandSampleWithCondition(obj,N,data,cond) % randomly samples N data points from data matrix that match cond condition
             sampInds=cell2mat(arrayfun(@(x) obj.RandSamplePopulation(find(data==cond(x))',N(x)),1:length(N),'uniformoutput',0));
             sampData=data(sampInds);
        end
        
        function [sampData,sampInds]=RandSampleWithCondition_MaxMatch(obj,N,data,cond,MatchFactors,CurrentFactors) % randomly samples N data points from data matrix that match cond condition
             NSamps=size(data,1);
             % loop on MatchFactors and see if you can find any stimlus
             % that exactly matches these conditions 
             CondVals=cell2mat(arrayfun(@(x) cond(x)*ones(1,N(x)),1:length(N),'UniformOutput',0));
             CurrentSamps=1:NSamps;
             for i=1:sum(N)
                 if isempty(MatchFactors)
                     MatchStim=[];
                 else
                     ThisStim=MatchFactors(i,:);
                     MatchStimIndAll=cell2mat(arrayfun(@(x) CurrentFactors(CurrentSamps,x)==ThisStim(x)',1:length(ThisStim),'UniformOutput',0));
                     MatchStim=CurrentSamps(sum(MatchStimIndAll,2)==length(ThisStim));
                 end
                 if ~isempty(MatchStim) % then we have a matching stimulus sample from these stimulus conditions
                     sampInds(i)=obj.RandSamplePopulation(MatchStim,1);
                 else
                     sampInds(i)=obj.RandSamplePopulation(CurrentSamps(find(data(CurrentSamps)==CondVals(i))'),1);
                 end
                 % remove this sample from the population now 
                 CurrentSamps=setdiff(1:NSamps,sampInds);
             end
             sampData=data(sampInds);
             if ~isempty(setdiff(sampData,cond))  
                 error('Conditions are not fully matching here');
             end
        end

        function [sampData,sampInds]=RandSampleWithDoubleCondition(obj,N,data,cond) % randomly samples N data points from data matrix that match two conditions
            %@N number of samples for  each condition(should match the number of rows in cond 
            %@data is a two column matrix with each column representing condition values for samples 
            %@cond is a two column matrix with corresponding columns to data 
                
            sampInds=cell2mat(arrayfun(@(x) obj.RandSamplePopulation(find(data(:,1)==cond(x,1) & data(:,2)==cond(x,2))',...
                N(x)),1:length(N),'uniformoutput',0));
            
            sampData=data(sampInds,:);
        end
        function [sampData,sampInds]=RandSampleWithDoubleCondition_MaxMatch(obj,N,data,cond,MatchFactors,CurrentFactors) % randomly samples N data points from data matrix that match two conditions
            %@N number of samples for  each condition(should match the number of rows in cond 
            %@data is a two column matrix with each column representing condition values for samples 
            %@cond is a two column matrix with corresponding columns to data 
             NSamps=size(data,1);
             % loop on MatchFactors and see if you can find any stimlus
             % that exactly matches these conditions 
             CondVals=cell2mat((arrayfun(@(x) repmat(cond(x,:),[N(x),1])',1:length(N),'UniformOutput',0)))';
             CurrentSamps=1:NSamps;
             for i=1:sum(N)
                 if isempty(MatchFactors)
                     MatchStim=[];
                 else
                     ThisStim=MatchFactors(i,:);
                     MatchStimIndAll=cell2mat(arrayfun(@(x) CurrentFactors(CurrentSamps,x)==ThisStim(x)',1:length(ThisStim),'UniformOutput',0));
                     MatchStim=CurrentSamps(sum(MatchStimIndAll,2)==length(ThisStim));
                 end
                 if ~isempty(MatchStim) % then we have a matching stimulus sample from these stimulus conditions
                     sampInds(i)=obj.RandSamplePopulation(MatchStim,1);
                 else
                     try
                         sampInds(i)=obj.RandSamplePopulation(CurrentSamps(find(data(CurrentSamps,1)==CondVals(i,1) & data(CurrentSamps,2)==CondVals(i,2)))',1);
                     catch
                         error('check this code')
                     end
                    % sampInds(i)=obj.RandSamplePopulation(CurrentSamps(find(data(CurrentSamps)==CondVals(i))'),1);
                 end
                 % remove this sample from the population now 
                 CurrentSamps=setdiff(1:NSamps,sampInds);
             end
             sampData=data(sampInds,:);
             if ~isempty(setdiff(sampData(:,1),cond(:,1))) | ~isempty(setdiff(sampData(:,2),cond(:,2))) 
                 error('Conditions are not fully matching here');
             end
        end
        function FileNameSyntax=GetFileNameSyntax(obj,AnalysisName)% naming conventions of differnt analysis 
           global AnalysisOpts
            %% classifier results 
            switch AnalysisName
                case 'ClassifierData'
         %   obj.SaveVar('Classifier',ClassifierResults,'ClassifierResults',...
                FileNameSyntax=['_' ClassifierOpts.Name '_' AnalysisOpts.CurrentCh_AreaName '_' AnalysisOpts.SpkCntStartFieldName '_' AnalysisOpts.TrlSpkTimeFieldName];
                case 'ImportantData'
                %% ImportantData
                 FileNameSyntax=[AnalysisOpts.DataSavePath 'Core Data' filesep 'SpikingData' filesep AnalysisOpts.Animal '_' AnalysisOpts.Area2look{1} '_' AnalysisOpts.TrlSpkTimeFieldName '_' AnalysisOpts.SpkCntStartFieldName '.mat']; 
            end
        end
        function String=RepDashSpace(~,String)% replaces dash with space
            if isempty(String);return;end
            if iscell(String)
                String= cellfun(@(x) strrep(x,'_',' '),String,'UniformOutput',0);
            else
                String=strrep(String,'_',' ');
            end
        end
        function [MatInclude,MatWhold,IncludeInd,WholdInd]=RandSampleMatrixwithProportion(~,Mat,Prop_Withhold,Dim,IncludeInd,WholdInd) 
            %randomly samples data from a matrix in a specific dimension with the Prop_Withhold of number of samples
          if isempty(Dim); Dim=1;end
          if isempty(IncludeInd) % if we provide the inds then don't recalculate it.
              N=size(Mat,Dim);
              nWhold=ceil(N*Prop_Withhold);
              WholdInd=randsample(1:N,nWhold);
              IncludeInd=setdiff(1:N,WholdInd);
          end
          if Dim==1
              MatWhold=Mat(WholdInd,:);
              MatInclude=Mat(IncludeInd,:);
          elseif Dim==2
              MatWhold=Mat(:,WholdInd);
              MatInclude=Mat(:,IncludeInd);
          end        
        end
        function ind=GetInd(~,Mat,Val) %finds a value in a cell or matrix
            if ischar(Val)
                ind=find(strcmp(Mat,Val));
            else
                ind=find(Mat==Val);
            end
        end
        function out=FlipRepDirection(~,data)
            out=zeros(size(data));
            out(data==1)=2;
            out(data==2)=1;
            out(data==3)=4;
            out(data==4)=3;
        end
        function yCI95=Calculate95thPercentile(obj,y,BootstrapdataFlag) % calculates 95th confidance inteval 
            %@Y is the matrix (NxT) N is the number of observations and T
            %is independant variables
            % taken from https://www.mathworks.com/matlabcentral/answers/414039-plot-confidence-interval-of-a-signal?s_tid=answers_rc1-2_p2_MLT
            % BootstrapdataFlag are the repetitions boot strap?

            N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
            yMean = mean(y,1);                                  % Mean Of All Experiments At Each Value Of ‘x’
            if BootstrapdataFlag
                ySEM = nanstd(y,1,1);
            else
                ySEM = nanstd(y,1,1)/sqrt(N);                       % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
            end
            CI95 = tinv([0.975], N-1);                          % Calculate 95% Probability Intervals Of t-Distribution
            yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
        end
        function [X,P]=DifferentiateSigClusters(obj,x,p) % finds seperate clusters in significance level 
            %x matrix of index of significant points 
            %p significant level 
            if isempty(x);X=[];P=[];return;end
            DiffInd=diff(x);
            ClustInds=unique([x(1) x(DiffInd>1)+1 x(end)]); % find discontiuties
            X=arrayfun(@(y)  x(x>=ClustInds(y) & x<=ClustInds(y+1)),1:(size(ClustInds,2)-1),'uniformoutput',0);
            P=arrayfun(@(y)  p(x>=ClustInds(y) & x<=ClustInds(y+1)),1:(size(ClustInds,2)-1),'uniformoutput',0);
 
        end
        function [X,P]=DifferentiateSigClusters2(obj,p) % finds seperate clusters in significance level (we are only given a pval per point)
            %p significant level per time point 
            if isempty(p);X=[];P=[];return;end
            % create a time matrix where we set different values based on the p value threshold 
            L=length(p);
            x=1:L;
            [~,~,pdesc]=arrayfun(@(x) obj.pvalueStar(x),p,'UniformOutput',0);
            pdesc=cell2mat(pdesc);

            DiffInd=find(diff(pdesc));
            ClustInds=unique([x(1) x(DiffInd+1) x(end)]); % find discontiuties
            for y=1:(size(ClustInds,2)-1)
                if y<(size(ClustInds,2)-1)
                    X{y}=x(x>=ClustInds(y) & x<ClustInds(y+1));
                    P{y}=(pdesc(x>=ClustInds(y) & x<ClustInds(y+1)));
                else
                    X{y}=x(x>=ClustInds(y) & x<=ClustInds(y+1));
                    P{y}=(pdesc(x>=ClustInds(y) & x<=ClustInds(y+1)));
                end
            end
        end
        function [symb,col,fp]=pvalueStar(~,pval,varargin)
            if ~isempty(varargin);Alpha=varargin{1};else;Alpha=0.05;end

            if pval>Alpha
                symb='';
                col=1;
                fp=1;
            elseif pval<Alpha & pval>=Alpha*0.2
                symb='*';
                col=2;
                fp=Alpha-0.0001;
            elseif pval<Alpha*0.2 & pval>=Alpha*0.02
                symb='**';
                col=2;
                fp=Alpha*0.2-0.0001;
            elseif pval<Alpha*0.02
                symb='***';
                col=2;
                fp=Alpha*0.02-0.0001;
            elseif isnan(pval)
                symb='';
                col=1;
                fp=nan;
            end
        end
        function [azimuth,elevation,r]=ConvertCart2Spherical(obj,x,y,z) % converts cartisian to spherical coordinates
            Angles360 = @(a) rem(360+a, 360); % turn every angle to the range of 0 to 360

            azimuth = wrapTo360(Angles360(atan2d(y,x)));
            elevation = wrapTo360(Angles360(atan2d(z,sqrt(x.^2 + y.^2))));
            r = sqrt(x.^2 + y.^2 + z.^2);
        end
        function [xhat,P,p]=ProjectIntoSubspace(obj,b,a1,a2) % project into vector b onto closest point p in a plane that is defined by two basis vectors a1 and a2
            % this solution is taken from https://ocw.mit.edu/courses/18-06sc-linear-algebra-fall-2011/00e9c8f0eafedeab21a3d079a17ed3d8_MIT18_06SCF11_Ses2.2sum.pdf
            
            % if a1 and a2 form a basis for a plane then that plane is the
            % clumn space of the matrix A=[a1 a2]
            A=[a1 a2];
            % we know that p=xhat1*a1+xhat2*a2=Axhat and we want to find
            % xhat
            xhat=inv(A'*A)*A'*b;
            P=A*inv(A'*A)*A';
            p=P*b;               
        end
        function PC=CalPercentChange(~,V1,V2) % calculates percente change
            PC=(V2-V1)./abs(V1);
        end
        function [Data1,Data2]=SplitDatainTwo(~,Data) % splits data in two sets of trials
            % Data is cell array where each array is a Ntrl*NNeu
           
            Inds1=cellfun(@(x) randperm(size(x,1),floor(size(x,1)/2)),Data,'UniformOutput',0);
            Inds2=arrayfun(@(x) setdiff(1:size(Data{x},1),Inds1{x}),1:length(Data),'UniformOutput',0);
            Data1=arrayfun(@(x) Data{x}(Inds1{x},:),1:length(Data),'UniformOutput',0);
            Data2=arrayfun(@(x) Data{x}(Inds2{x},:),1:length(Data),'UniformOutput',0);
        end
        function out=Diag3D(~,A)% retuns diagonal values of a 3D and 2D matrix 
            if size(A,1)~=size(A,2)
                out=A'; 
            else
                out=reshape(A(logical(repmat(eye(size(A(:,:,1))),[1 1 size(A,3)]))),[size(A,1) size(A,3)])';
            end

        end
        function [a,p,Pairs]=DoTtestonCell(~,Mat,PairdttestFlag) % performs ttest on all of the pairs of cells in Mat
            Pairs=nchoosek(1:length(Mat),2);
            if PairdttestFlag
                arrayfun(@(x) ttest2(Mat{Pairs(x,1)},Mat{Pairs(x,2)}),1:size(Pairs,1));
            else
                [a,p]=arrayfun(@(x) ttest(Mat{Pairs(x,1)},Mat{Pairs(x,2)}),1:size(Pairs,1));
            end
            Pairs=mat2cell(Pairs,[ones(1,size(Pairs,1))],[2]);
        end
        function nUniq=CalNumUniqueCombsofConds(obj,Conds) % find the number of uniqe coniditons for this matrix
            % Conds is a N by B matrix where N is the number of
            % observations
            nUniq=size(unique(Conds,'rows'),1);
        end
        function DateNum=DetermineDateNum(obj,DateNum)
            global AnalysisOpts
            
            if DateNum==0 | strcmp(DateNum,'ALL') | strcmp(DateNum,'A');DateNum=0;AnalysisOpts.Animal='ALL';
            elseif strcmp(DateNum,'Chico') | strcmp(DateNum,'C');DateNum=AnalysisOpts.ChicoRecNums;
            elseif strcmp(DateNum,'Silas') | strcmp(DateNum,'S');DateNum=AnalysisOpts.SilasRecNums;
            elseif ischar(DateNum); DateNum=find(strcmp(AnalysisOpts.DateSet,DateNum));
            elseif isnumeric(DateNum) % datenum is a number
                IsChico=sum(arrayfun(@(x) sum(x==AnalysisOpts.ChicoRecNums),DateNum));
                IsSilas=sum(arrayfun(@(x) sum(x==AnalysisOpts.SilasRecNums),DateNum));
                if IsChico & IsSilas;AnalysisOpts.Animal='ALL';end
                if ~IsChico & IsSilas;AnalysisOpts.Animal='Silas';end
                if IsChico & ~IsSilas;AnalysisOpts.Animal='Chico';end
            end
            if sum(DateNum)==sum(AnalysisOpts.ChicoRecNums);AnalysisOpts.Animal='Chico';DateNum=AnalysisOpts.ChicoRecNums;end
            if sum(DateNum)==sum(AnalysisOpts.SilasRecNums);AnalysisOpts.Animal='Silas';DateNum=AnalysisOpts.SilasRecNums;end           
        end
        function Inds=FindValsinMat(obj,Mat,Vals) % find indexes for values in a matrix             
            Inds=obj.ReshapeCell2Mat(arrayfun(@(x) Mat==x,Vals,'UniformOutput',0),3);
            Inds=sum(Inds,3)>0;
        end
        function Inds=FindValsinMat2D(obj,Mat,Vals) % find indexes for values in a matrix but in two dimensions            
            Inds=obj.ReshapeCell2Mat(arrayfun(@(x) Mat(:,1)==Vals(x,1) & Mat(:,2)==Vals(x,2),1:size(Vals,1),'UniformOutput',0),3);
            Inds=sum(Inds,3)>0;
        end
        function [Dup,Count,UV,NDup,DupLoc,nCell]=FindDuplicateVals(~,MatCell) % Mat is a cell array
            Mat=cell2mat(MatCell);
            UV=unique(Mat(:));
            Count=arrayfun(@(x) sum(Mat==x),UV);
            Dup=UV(Count>1);
            NDup=length(Dup);
            % now find out which array had the dupicataes 
            DupLoc=arrayfun(@(y) cellfun(@(x) sum(x==y),MatCell,'UniformOutput',1),Dup,'UniformOutput',0);
            nCell=cellfun(@length,MatCell);
        end
        function [a,p]=CorrelateMats(~,M1,M2)
            if size(M1)==size(M2)
                [a,p]=corr(M1,M2);
            else
                a=nan*ones(size(M1));
            end
            
        end
        function ThrowErrorTxt(~,me) % try catch error text. me is the error strcuture from catch
            fprintf(2,'\nERROR **********%s***********',me.message)
        end
        function [OnSetLatencySec,PeakLatencySec]=CalPeakBaselineLatency(~,Data,Time) % calculates latency to peak and baseline
            global AnalysisOpts
            % Data needs to be normalizaed already
            % TH is the treshold where we consider the rise time
            baseline_percentage = AnalysisOpts.baseline_percentage_RiseTime; % You can adjust this value
            peak_percentage     = AnalysisOpts.peak_percentage_RiseTime;     % You can adjust this value
           
            % Interpolation points
            interpolation_points = linspace(1, numel(Data), 1000); % Increase the number for higher resolution
            
            % Interpolate the normalized PSTH
            interpolated_psth = interp1(1:numel(Data), Data, interpolation_points, 'linear');
            interpolated_Time=  interp1(1:numel(Data), Time, interpolation_points, 'linear');
            
            SampleOnset=find(interpolated_Time>=0,1,'first');           
            OnSetLatencyInd=find(interpolated_psth(SampleOnset:end)>=baseline_percentage,1,'first')+SampleOnset-1;
            OnSetLatencySec=interpolated_Time(OnSetLatencyInd);
            PeakLatencyInd=find(interpolated_psth(SampleOnset:end)>=peak_percentage,1,'first')+SampleOnset-1;
             PeakLatencySec=interpolated_Time(PeakLatencyInd);                
            if isempty(OnSetLatencySec);OnSetLatencySec=nan;end
            if isempty(PeakLatencySec);PeakLatencySec=nan;end
        end
        function [outTxt,p]=Perform_Jttred_StatTest(obj,groups,Data) % performs JONCKHEERE-TERPSTRSA statistical test using jtternd
          % @ Data is organized as rep*Stage

            % example data organization is 
            % Data=[0 0 1 1 2 2 4 9 0 0 5 7 8 11 13 23 25 97 2 3 6 9 10 11 11 12 21 ...
            %       0 3 5 6 10 19 56 100 132 2 4 6 6 6 7 18 39 60];
            %  Ind=[ones(1,8) 2.*ones(1,10) 3.*ones(1,9) 4.*ones(1,9) 5.*ones(1,9)];
            %    x=[Data' Ind'];
            L=size(Data,1);
            groups=transpose(cell2mat(arrayfun(@(x) groups(x)*ones(1,L),1:length(groups),'UniformOutput',0)));
            Data=Data(:);
            p=jttrend([Data groups]);
            outTxt=sprintf('p_Jttred=%0.4f',p);
        end

        function [outTxt,p,tau,z_score,H]=Perform_Modified_MannKendall_StatTest(obj,Data,varargin) % performs Modified_MannKendall statistical test using jtternd
          % @ Data is organized as rep*Stage
           global AnalysisOpts
           
            if isempty(varargin)
                num_stages=1:size(Data,2);
            else
                num_stages=varargin{1};
            end
           % if size(Data,1)>1;Data=mean(Data,1);end
           try
            if AnalysisOpts.Classifier_TrlShuff_TrendCorrMethodCode % which code we are using
                [tau, z_score, p, H] =  Modified_MannKendall_test(num_stages, Data, ...
                    AnalysisOpts.Modified_MannKendall_significance_value_tau, ...
                   AnalysisOpts.Modified_MannKendall_significance_value_ac);%,'uniformoutput',1);
            else
                [H,p]=Mann_Kendall_Modified(Data,AnalysisOpts.Modified_MannKendall_significance_value_tau);
            end
           catch me
               obj.ThrowErrorTxt(me);
               tau=nan; z_score=nan; p=nan; H=nan;
           end
           % p=cell2mat(p);
           % tau=tau{1};
            outTxt=sprintf('p_MMK=%0.4f',p);
        end
        function [outTxt,F_statistic,critical_f_value,p_anova,df,dferr,tbl,stats,alpha]=Perform_OneWayANOVA_StatTest(obj,groups,Data,alpha) 
            % @ Data is organized as rep*Stage

            L=size(Data,1);
            groups=(cell2mat(arrayfun(@(x) groups(x)*ones(1,L),1:length(groups),'UniformOutput',0)));
            
            Data=[Data(:)]';
            
            % Perform one-way ANOVA
            [p_anova, tbl, stats] = anova1(Data, groups, 'off');
            
            % Get the test statistic from ANOVA results
            F_statistic = tbl{2, 5};
            df=tbl{2,3};
            dferr=tbl{3,3};
            % Set significance level (e.g., 0.05)
           % alpha = 0.05;
            
            % Compute critical F-value
            critical_f_value = finv(1 - alpha, numel(unique(groups)) - 1, numel(Data) - numel(unique(groups)));
            
            % Make a decision based on the test statistic and critical F-value
            outTxt=sprintf('F(%i,%i)=%0.4f,p=%0.4f',df,dferr,F_statistic,p_anova);            
        end
        function [outTxt,p_values,significant_stage]=Perform_PairTtestTrend_StatTest(obj,Data,alpha) % performs pairs of ttest compare when there is a difference between the first group in a trend
            % @ Data is organized as rep*Stage
            
            % Number of stages
            num_stages = size(Data,2);
            
            % Calculate pairwise t-test p-values
            p_values = zeros(num_stages - 1, 1);
            
            for i = 2:num_stages
                [~, p_values(i - 1)] = ttest2(Data(:,1), Data(:,i)); % using a two sampel ttest for these data 
            end
            ncomparisions=num_stages-1;
            % Set significance level (e.g., 0.05)
            % alpha = 0.05;
            
            % Find the first stage where the performance is significantly different
            significant_stage = find(p_values < (alpha/ncomparisions))+1;      
            outTxt=sprintf('Paird %s', obj.ConvMat2Char(significant_stage));
        end
        function outTxt=PerformAllTrend_StatTest(obj,groups,Data,alpha) % performs series of statistical tests on ordinal data to show the trend of significance 
           % performs JONCKHEERE-TERPSTRSA, One-way ANOVA and ttest 
           % [outTxt{1,1}]=obj.Perform_Jttred_StatTest(groups,Data);
           % performs Modified_MannKendall, One-way ANOVA and ttest
           try
            [outTxt{1,1}]=obj.Perform_Modified_MannKendall_StatTest(Data);
         %   [outTxt{2,1}]=obj.Perform_OneWayANOVA_StatTest(groups,Data,alpha);
         %   [outTxt{3,1}]=obj.Perform_PairTtestTrend_StatTest(Data,alpha);
           catch me
               obj.ThrowErrorTxt(me)
               outTxt=[];
           end
        end
        function [x_no_outliers,y_no_outliers]=removeoutliers(~,x,y,z_score_threshold)
            
            % Calculate z-scores for each data point
            z_scores_x = (x - mean(x)) / std(x);
            z_scores_y = (y - mean(y)) / std(y);
            
            % Define z-score threshold for outlier detection
            if isempty(z_score_threshold)
            z_score_threshold = 3;  % Adjust as needed
            end
            
            % Identify outliers based on z-scores
            outliers = abs(z_scores_x) > z_score_threshold | abs(z_scores_y) > z_score_threshold;
            
            % Remove outliers from data
            x_no_outliers = x(~outliers);
            y_no_outliers = y(~outliers);
        end
        function [a,p]=Correlation(obj,X,Y,RemoveOutLierFlag,varargin)
           global AnalysisOpts
            % the observations should be columns
            if RemoveOutLierFlag
                [X,Y]=obj.removeoutliers(X,Y,[]);
            end
            if size(Y,2)>size(X,2) 
                Ncomp=size(Y,2);
            else
                Ncomp=1;
            end
            if isempty(varargin)
                [a,p]=corr(X,Y);
            elseif strcmp(varargin{1},'Modified_MannKendall')
                [outTxt,p,a,z_score,H]=arrayfun(@(x) obj.Perform_Modified_MannKendall_StatTest(Y(:,x),X),1:Ncomp,'UniformOutput',0);
                a=cell2mat(a);p=cell2mat(p);
%                 % perform correction as well
%                 if AnalysisOpts.Classifier_TrlShuff_UseBonferroni==1 & Ncomp>1
%                     p=p*Ncomp; % do the bonferroni correction
%                 elseif AnalysisOpts.Classifier_TrlShuff_UseBonferroni==2 & Ncomp>1% then do benjamini-Hotchbery correction
%                     [~, ~, ~, p]=fdr_bh(p,AnalysisOpts.Classifier_TrlShuff_BenjaminiHochberg_FalseDiscoveryRate,'pdep','yes');
%                 end
            elseif strcmp(varargin{1},'Spearman')
                [a,p]=corr(X,Y,'type','spearman');
            end
        end
        function pval=PerformStattestMultipleComparisionCorrection(obj,pval)
            global AnalysisOpts

            Ncomp=length(pval);
            if AnalysisOpts.Classifier_TrlShuff_UseBonferroni==1  & Ncomp>1
                pval=pval*Ncomp; % do the bonferroni correction
            elseif AnalysisOpts.Classifier_TrlShuff_UseBonferroni==2 & Ncomp>1% then do benjamini-Hotchbery correction
                [~, ~, ~, pval]=fdr_bh(pval,AnalysisOpts.Classifier_TrlShuff_BenjaminiHochberg_FalseDiscoveryRate,'pdep','yes');
            elseif AnalysisOpts.Classifier_TrlShuff_UseBonferroni==3 & Ncomp>1% then do bonf-holm correction
                [pval]=bonf_holm(pval,0.05);
            end
        end
        function [Y,opts]=NormalizeData(obj,X,Y,opts)% normalizes data with subtracting the baseline and then by max or by mean
            global AnalysisOpts
            
            % adjust point base line subtraction 
            if strcmp(opts.NPnts_SubtractBaseLine,'auto') % then adjuct it based on the time limit we have put
               opts.NPnts_SubtractBaseLine=X>AnalysisOpts.PaperSpec.NPnts_SubtractBaseLine(1) & X<AnalysisOpts.PaperSpec.NPnts_SubtractBaseLine(2);
            elseif ~islogical(opts.NPnts_SubtractBaseLine) 
               opts.NPnts_SubtractBaseLine=1:opts.NPnts_SubtractBaseLine;
            end
            % do we need to normalize the the mean data by their max?
            if opts.SubtractBaseLine
                meanY=nanmean(Y,1);
                BaseLineMean=mean(meanY(opts.NPnts_SubtractBaseLine));
                Y=Y-BaseLineMean;           
            end
            if opts.NormalizebyMax 
                meanY=nanmean(Y,1);
                MaxMeanY=max(meanY);
                Y=Y./MaxMeanY;             
            end
            if opts.NormalizebyMean % divide by mean  
                meanY=nanmean(Y,1);
                BaseLineMean=mean(meanY(opts.NPnts_SubtractBaseLine));           
                Y=Y./BaseLineMean;             
            end     
        end
        function [Yorg,Y,YSTD]=CalDataSTD(~,Y,YSTD,opts) % calculates standard diviation of the data 
            Yorg=Y;
            if isempty(YSTD)                
                if opts.CircularData  % are these circular data
                    YSTD=rad2deg(circ_std(Y,[],[],1));
                    Y=rad2deg(circ_mean(Y,[],1));
                else
                    YSTD=nanstd(Y,1,1);
                    Y(isinf(Y))=nan;
                    Y=nanmean(Y,1);
                end
                if opts.UseCos % are we using cosine
                    if sum(Y>180);error('This can only be used with anlges 0-180');end
                    YSTD=acosd(YSTD);
                    Y=acosd(Y);
                end
                switch opts.STD_method
                    case 'sem' % standard error of the mean
                        YSTD=YSTD/sqrt(size(Yorg,1));
                        %    warning('STD plot method is not bootstrap')
                    case '95p' % 95 percentile
                        YSTD=obj.Calculate95thPercentile(Yorg,0);
                        %   warning('STD plot method is not bootstrap')
                    case 'bootstrap' % in this case sem=std
                        YSTD=YSTD;
                        %  YSTD=obj.ManData.Calculate95thPercentile(Yorg,1);
                end
            end
        end
        function [OnSetLatencySec1,OnSetLatencySec2,pval]=CompareRiseTimesPaired(obj,Y1,Y2,Time,NormalizebyMax,NPnts_SubtractBaseLine,SubtractBaseLine) % compares rise times of two distibutions 
            % Y1 nd Y2 is matrix of Trl by Time
            global AnalysisOpts
                      
            % prepare vars
            opts.NormalizebyMax=NormalizebyMax;
            opts.SubtractBaseLine=SubtractBaseLine;
            opts.NPnts_SubtractBaseLine=NPnts_SubtractBaseLine;            
            opts.NormalizebyMean=0; % we don't use normalize by mean here
            num_bootstrap_iterations = 1000;
            
            NTrls1=size(Y1,1);NTrls2=size(Y2,1);
            
           NBootStrapSamp=100;
            % create a bootstrap resampleing to better estimate the mean
            for i = 1:num_bootstrap_iterations
                % Randomly draw a subsample with replacement Y1
                subsample1 = randsample(NTrls1, NBootStrapSamp, true);
                Y1BS{i}=mean(Y1(subsample1,:),1);
                % Randomly draw a subsample with replacement  Y2              
                Y2BS{i}=mean(Y2(subsample1,:),1);                               
            end

              % normalize the data first
            [NormY1]=(arrayfun(@(x) obj.NormalizeData(Time,Y1BS{x},opts)',1:num_bootstrap_iterations,'UniformOutput',0));
            [NormY2]=(arrayfun(@(x) obj.NormalizeData(Time,Y2BS{x},opts)',1:num_bootstrap_iterations,'UniformOutput',0));
            
            % calculate latency
            [OnSetLatencySec1]=arrayfun(@(x) obj.CalPeakBaselineLatency(NormY1{x},Time),1:num_bootstrap_iterations,'UniformOutput',1);
            [OnSetLatencySec2]=arrayfun(@(x) obj.CalPeakBaselineLatency(NormY2{x},Time),1:num_bootstrap_iterations,'UniformOutput',1);
         
            
%             % normalize the data first
%             [NormY1]=(arrayfun(@(x) obj.NormalizeData(Time,Y1(x,:),opts)',1:NTrls1,'UniformOutput',0));
%             [NormY2]=(arrayfun(@(x) obj.NormalizeData(Time,Y2(x,:),opts)',1:NTrls2,'UniformOutput',0));
%             
%             % calculate latency
%             [OnSetLatencySec1]=arrayfun(@(x) obj.CalPeakBaselineLatency(NormY1{x},Time),1:NTrls1,'UniformOutput',1);
%             [OnSetLatencySec2]=arrayfun(@(x) obj.CalPeakBaselineLatency(NormY2{x},Time),1:NTrls2,'UniformOutput',1);
%            
            valid_indices = ~isnan(OnSetLatencySec1) & ~isnan(OnSetLatencySec2);
            OnSetLatencySec1 = OnSetLatencySec1(valid_indices);
            OnSetLatencySec2 = OnSetLatencySec2(valid_indices);

            %  Perform a statistical test (Mann-Whitney U-test) to compare rise times
            Diff=(OnSetLatencySec1-OnSetLatencySec2)';
            pval=1-obj.CalpValShuffle(Diff,0);
            
            % Display the results
           % fprintf('Mann-Whitney U-test p-value: %.4f\n', pval);
        end
       
        function [OnSetLatencySec1,OnSetLatencySec2,pval]=CompareRiseTimes(obj,Y1,Y2,Time,NormalizebyMax,NPnts_SubtractBaseLine,SubtractBaseLine) % compares rise times of two distibutions 
            % Y1 nd Y2 is matrix of Trl by Time
            global AnalysisOpts
                      
            % prepare vars
            opts.NormalizebyMax=NormalizebyMax;
            opts.SubtractBaseLine=SubtractBaseLine;
            opts.NPnts_SubtractBaseLine=NPnts_SubtractBaseLine;            
            opts.NormalizebyMean=0; % we don't use normalize by mean here
            num_bootstrap_iterations = 1000;
            
            NTrls1=size(Y1,1);NTrls2=size(Y2,1);
            
           
            % create a bootstrap resampleing to better estimate the mean
            for i = 1:num_bootstrap_iterations
                % Randomly draw a subsample with replacement Y1
                subsample1 = randsample(NTrls1, NTrls1, true);
                Y1BS{i}=mean(Y1(subsample1,:),1);
                % Randomly draw a subsample with replacement  Y2              
                subsample2= randsample(NTrls2, NTrls2,true);
                Y2BS{i}=mean(Y2(subsample2,:),1);                               
            end
            
            % normalize the data first
            [NormY1]=(arrayfun(@(x) obj.NormalizeData(Time,Y1BS{x},opts),1:num_bootstrap_iterations,'UniformOutput',0));
            [NormY2]=(arrayfun(@(x) obj.NormalizeData(Time,Y2BS{x},opts),1:num_bootstrap_iterations,'UniformOutput',0));
            
            % calculate latency
            [OnSetLatencySec1]=arrayfun(@(x) obj.CalPeakBaselineLatency(NormY1{x},Time),1:num_bootstrap_iterations,'UniformOutput',1);
            [OnSetLatencySec2]=arrayfun(@(x) obj.CalPeakBaselineLatency(NormY2{x},Time),1:num_bootstrap_iterations,'UniformOutput',1);
           
            valid_indices = ~isnan(OnSetLatencySec1) & ~isnan(OnSetLatencySec2);
            OnSetLatencySec1 = OnSetLatencySec1(valid_indices);
            OnSetLatencySec2 = OnSetLatencySec2(valid_indices);

            %  Perform a statistical test (Mann-Whitney U-test) to compare rise times
            [pval, hval] = ranksum(OnSetLatencySec1, OnSetLatencySec2);
            
            % Display the results
           % fprintf('Mann-Whitney U-test p-value: %.4f\n', pval);
        end
       
        
        
        function AreaNum=getAreaNum(~,AreaName) % conevrts areaname to area num
            global AnalysisOpts
            
            if ischar(AreaName);AreaName={AreaName};end
            AreaNum=cellfun(@(x) find(strcmp(AnalysisOpts.AreaNames,x)),AreaName);
        end
        function population_significance=CalPvalPopulationUsingZscore(obj,ShuffleDist,z_scores_individual,Pval) % pvalue based on the zscore
            nShuffle = size(ShuffleDist, 1);
            numTimePoints = size(ShuffleDist, 2);

            % Preallocate for efficiency
            population_significance = zeros(1, numTimePoints);

            % Generate a population shuffle distribution for each time point
            numResamples = 1000;

            for t = 1:numTimePoints
                population_shuffle = zeros(numResamples, 1);

                for i = 1:numResamples
                    sampled_indices = randsample(nShuffle, 2, true);
                    population_shuffle(i) = mean(ShuffleDist(sampled_indices, t));
                end

                % Calculate the population z-score for this time point
                aggregated_z = mean(z_scores_individual(:, t));
                population_mean = mean(population_shuffle);
                population_std = std(population_shuffle);
                population_z = (aggregated_z - population_mean) / population_std;

                % Calculate the p-value (two-tailed test)
                p_value = 2 * min(normcdf(population_z), 1 - normcdf(population_z));

                % Determine if the time point's effect is significant
                population_significance(t) = p_value;
            end


        end
        function [clusters, p_values, t_sums, permutation_distribution,statsummery ]=CalPvalPopulation(obj,ShuffleDist,Observed,Pval)
            nShuffle = size(ShuffleDist, 1);
            numTimePoints = size(ShuffleDist, 2);
            NNeu=size(Observed,1);
            NSampPerShuff=nShuffle/NNeu; % how many repetitions of shuffle per neuron


            % Generate a population shuffle distribution for each time point
            numResamples = 1000;
            population_shuffle = zeros(numResamples, numTimePoints);

            for i = 1:numResamples
                sampled_indices = obj.SampleOnNumfromEachInterval(nShuffle,NSampPerShuff); % sample one shuffle value from each neuron
                population_shuffle(i,:) = mean(ShuffleDist(sampled_indices, :),1);
            end
            population_observed=mean(Observed,1);
            % make sure it matches the test structure
            population_shuffle=permute(population_shuffle,[1 3 2]);
            population_observed=permute(population_observed,[1 3 2]);

            % use cluster correction now for significance 
            [p_values,t_sums, clustIdx,permutation_distribution,clusters,statsummery] = obj.ClusterMassCorrection_permutation(population_shuffle,population_observed,Pval,0);
                       
        end
        function sampledValues=SampleOnNumfromEachInterval(obj,N,Int) % samples one number from each interval of values
            % Define the range and interval
           % N = 100;  % Maximum value
           % Int = 5;  % Interval

            % Calculate the number of intervals
            numIntervals = floor(N / Int);

            % Initialize an array to store sampled values
            sampledValues = zeros(1, numIntervals);

            % Sample a random number from each interval
            for i = 1:numIntervals
                lowerBound = 1 + (i - 1) * Int;
                upperBound = i * Int;
                sampledValues(i) = randi([lowerBound, upperBound]);
            end

        end
        function p=LearningTimeTrendStatTest(obj,MetricValsOrg,varargin)
            obj=obj.ParseParams(varargin) ; %%Process optional inputs

            NonOverlappingBins=[1 size(MetricValsOrg,3)]; % find what are our non overlapping bins

            % loop on time points
            NTim=size(MetricValsOrg,1);
            for t=1:NTim
                [~,p(t)]=obj.Perform_Jttred_StatTest(1:length(NonOverlappingBins),squeeze(MetricValsOrg(t,:,NonOverlappingBins)));
            end
            % correct p value for bon ferroni
            p=p*NTim;

        end
        function NewFactorData = SwapShuffleTrialsFactorData(obj, ClassifierOpts, FactorData, ShuffRep, DimNum, FactorLevelComb, TargetFactorInd, DimMethod, TargetFactorInd1ndD,Cond,TrlRng)
            % SwapShuffleTrialsFactorData - This function swaps trials in shuffled data for factor data.
            global AnalysisOpts
            tic
            NumLevelComb = size(FactorLevelComb, 1);  % Calculate the number of factor level combinations
            NNeu = length(FactorData);  % Get the number of neurons
            BlkOrderInd=(strcmp(AnalysisOpts.factornames,'BlkOrder'));
            RevTrlNumInd=(strcmp(AnalysisOpts.factornames,'TrialNumReverse'));
            
            if find(TargetFactorInd)==find(strcmp(AnalysisOpts.factornames,'TrialNum'))
                ThisIsTrlShuff=1;
            else 
                ThisIsTrlShuff=0;
            end % if we are doing trial shuffle
            
            for Neu = 1:NNeu
                % Initialize new matrices for the current neuron
                NewFactorData(Neu).data = cell(1, NumLevelComb);
                NewFactorData(Neu).AllFactors = cell(1, NumLevelComb);

                % Initialize fields for special cases (spike count periods)
                for nSkpCnt = 1:size(ClassifierOpts.SpkCountPeriod, 1)
                    NewFactorData(Neu).FactorData_SpkCnt{nSkpCnt} = cell(1, NumLevelComb);
                end

                % Get relevant data and labels
%                 TrialIndex = ClassifierOpts.ClassifierShuffleTrialIndex{Cond}{TrlRng}{DimNum}{Neu};
%                 Label = ClassifierOpts.ClassifierShuffleLabel{Cond}{TrlRng}{DimNum}{Neu}(:, ShuffRep);
%                 TrainCondIndex = ClassifierOpts.ClassifierShuffleTrainCondIndex{Cond}{TrlRng}{DimNum}{Neu};
                TrialIndex = ClassifierOpts.ClassifierShuffleTrialIndex{DimNum}{Neu};
                Label = ClassifierOpts.ClassifierShuffleLabel{DimNum}{Neu}(:, ShuffRep); % ShuffRep is 1 
                TrainCondIndex = ClassifierOpts.ClassifierShuffleTrainCondIndex{DimNum}{Neu};
                BlkOrder = ClassifierOpts.ClassifierShuffleTrainCondBlkOrder{DimNum}{Neu}(:, ShuffRep);
                RevTrlNum=ClassifierOpts.ClassifierShuffleRevTrlNum{DimNum}{Neu}(:, ShuffRep);
              
                % Get unique labels and train condition indices
                UniqLabel = unique(Label)';
                UniqCondInd = unique(TrainCondIndex);
                UsedIndFLC=[]; % gather what FLC we have used here
                if DimMethod == 1
                    % Swap trials based on label and condition
                    for L = UniqLabel
                        i = find(Label == L)';
                        % Determine where the current trial belongs
                        IndFLC = arrayfun(@(x) find(FactorLevelComb(:, 1) == L & FactorLevelComb(:, 2) == FactorLevelComb(TrainCondIndex(x), 2)),i);
                        UniqIndFLC=unique(IndFLC);

                        % Copy data and factors across the matrix
                        for flc=UniqIndFLC
                            NewFactorData(Neu).data{flc}      =cell2mat(arrayfun(@(x) FactorData(Neu).data{TrainCondIndex(x)}(TrialIndex(x), :)',i(IndFLC==flc),'UniformOutput',0))';
                            NewFactorData(Neu).AllFactors{flc}=cell2mat(arrayfun(@(x) FactorData(Neu).AllFactors{TrainCondIndex(x)}(TrialIndex(x), :)',i(IndFLC==flc),'UniformOutput',0))';

                            % Update for special cases (spike count periods)
                            for nSkpCnt = 1:size(ClassifierOpts.SpkCountPeriod, 1)
                                NewFactorData(Neu).FactorData_SpkCnt{nSkpCnt}{flc} = cell2mat(arrayfun(@(x) FactorData(Neu).FactorData_SpkCnt{nSkpCnt}{TrainCondIndex(x)}(TrialIndex(x)),i(IndFLC==flc),'UniformOutput',0))';
                            end

                            % Change the label for all trials in the same condition
                            NewFactorData(Neu).AllFactors{flc}(:, TargetFactorInd) = L;
                        end
                        UsedIndFLC=[UsedIndFLC UniqIndFLC];
                    end
                elseif DimMethod == 2
                    % Copy labels into the target factor
                    NewFactorData(Neu) = FactorData(Neu);
                    for IndFLC = UniqCondInd
                        NewFactorData(Neu).AllFactors{IndFLC}(:, TargetFactorInd) = Label(TrainCondIndex == IndFLC);
                    
                        if ThisIsTrlShuff % if this is trial order shuffle condition then keep the information about blockorder
                            NewFactorData(Neu).AllFactors{IndFLC}(:, BlkOrderInd) = BlkOrder(TrainCondIndex == IndFLC);
                            NewFactorData(Neu).AllFactors{IndFLC}(:, RevTrlNumInd) = RevTrlNum(TrainCondIndex == IndFLC);
                        end
                    end
                    UsedIndFLC=UniqCondInd;
                elseif DimMethod == 3

                    % Swap trials based on label and condition
                    for L = UniqLabel
                        i = find(Label == L)';
                        FirstDimFactorVal = arrayfun(@(x) FactorData(Neu).AllFactors{TrainCondIndex(x)}(TrialIndex(x), TargetFactorInd1ndD),i,'UniformOutput',1);

                        % Determine where the current trial belongs
                        IndFLC = arrayfun(@(x) find(FactorLevelComb(:, 1) == FirstDimFactorVal(x) & FactorLevelComb(:, 2) == FactorLevelComb(TrainCondIndex(i(x)), 2)),1:length(i));
                        UniqIndFLC=unique(IndFLC);

                        % Copy data and factors across the matrix
                        for flc=UniqIndFLC
                            NewFactorData(Neu).data{flc}=cell2mat(arrayfun(@(x) FactorData(Neu).data{TrainCondIndex(x)}(TrialIndex(x), :)',i(IndFLC==flc),'UniformOutput',0))';
                            NewFactorData(Neu).AllFactors{flc} = cell2mat(arrayfun(@(x) FactorData(Neu).AllFactors{TrainCondIndex(x)}(TrialIndex(x), :)',i(IndFLC==flc),'UniformOutput',0))';

                            % Update for special cases (spike count periods)
                            for nSkpCnt = 1:size(ClassifierOpts.SpkCountPeriod, 1)
                                NewFactorData(Neu).FactorData_SpkCnt{nSkpCnt}{flc} = cell2mat(arrayfun(@(x) FactorData(Neu).FactorData_SpkCnt{nSkpCnt}{TrainCondIndex(x)}(TrialIndex(x)),i(IndFLC==flc),'UniformOutput',0))';
                            end

                            % Change the label for all trials in the same condition
                            NewFactorData(Neu).AllFactors{flc}(:, TargetFactorInd) = L;
                        end
                        UsedIndFLC=[UsedIndFLC UniqIndFLC];
                    end
                end
               
                % now copy the data that were not touched back into the matrix
                NotusedIndFLC=setdiff(1:size(FactorLevelComb,1),unique(UsedIndFLC));
                for NUFLCind=NotusedIndFLC
                    NewFactorData(Neu).data(NUFLCind)=FactorData(Neu).data(NUFLCind);
                    NewFactorData(Neu).AllFactors(NUFLCind)=FactorData(Neu).AllFactors(NUFLCind);

                    for nSkpCnt = 1:size(ClassifierOpts.SpkCountPeriod, 1)
                        NewFactorData(Neu).FactorData_SpkCnt{nSkpCnt}(NUFLCind) =FactorData(Neu).FactorData_SpkCnt{nSkpCnt}(NUFLCind);
                    end
                end
               

            end
%             % if one of the target fields is congruency then<<<THIS IS NOT CORRECT>> refer to slack chat with Tim on 2/23/2024
%             % recalculate congruency for this
%             TargetFields=ClassifierOpts.StimCongruencyFactorName{ClassifierOpts.StimCongruency(DimNum)};
%             if sum(strcmp(TargetFields,'Congruency'))
%                 NewFactorData=obj.UpdateFactorDataCongruency(NewFactorData,ClassifierOpts);
%             end
            fprintf('Elapsed Time to swap Dim %i=%0.2f',DimNum,toc)
        end
       
        function NewFactorData = SwapShuffleTrialsFactorDataOld(~, ClassifierOpts, FactorData, ShuffRep, DimNum, FactorLevelComb, TargetFactorInd, DimMethod, TargetFactorInd1ndD)
            % SwapShuffleTrialsFactorData - This function swaps trials in shuffled data for factor data.

            NumLevelComb = size(FactorLevelComb, 1);  % Calculate the number of factor level combinations
            NNeu = length(FactorData);  % Get the number of neurons

            for Neu = 1:NNeu
                % Initialize new matrices for the current neuron
                NewFactorData(Neu).data = cell(1, NumLevelComb);
                NewFactorData(Neu).AllFactors = cell(1, NumLevelComb);

                % Initialize fields for special cases (spike count periods)
                for nSkpCnt = 1:size(ClassifierOpts.SpkCountPeriod, 1)
                    NewFactorData(Neu).FactorData_SpkCnt{nSkpCnt} = cell(1, NumLevelComb);
                end

                % Get relevant data and labels
                TrialIndex = ClassifierOpts.ClassifierShuffleTrialIndex{DimNum}{Neu};
                Label = ClassifierOpts.ClassifierShuffleLabel{DimNum}{Neu}(:, ShuffRep);
                TrainCondIndex = ClassifierOpts.ClassifierShuffleTrainCondIndex{DimNum}{Neu};

                % Get unique labels and train condition indices
                UniqLabel = unique(Label)';
                UniqCondInd = unique(TrainCondIndex);

                if DimMethod == 1
                    % Swap trials based on label and condition
                    for L = UniqLabel
                        for i = find(Label == L)'
                            % Determine where the current trial belongs
                            IndFLC = find(FactorLevelComb(:, 1) == L & FactorLevelComb(:, 2) == FactorLevelComb(TrainCondIndex(i), 2));

                            % Copy data and factors across the matrix
                            NewFactorData(Neu).data{IndFLC} = cat(1, NewFactorData(Neu).data{IndFLC}, FactorData(Neu).data{TrainCondIndex(i)}(TrialIndex(i), :));
                            NewFactorData(Neu).AllFactors{IndFLC} = cat(1, NewFactorData(Neu).AllFactors{IndFLC}, FactorData(Neu).AllFactors{TrainCondIndex(i)}(TrialIndex(i), :));

                            % Update for special cases (spike count periods)
                            for nSkpCnt = 1:size(ClassifierOpts.SpkCountPeriod, 1)
                                NewFactorData(Neu).FactorData_SpkCnt{nSkpCnt}{IndFLC} = cat(1, NewFactorData(Neu).FactorData_SpkCnt{nSkpCnt}{IndFLC}, FactorData(Neu).FactorData_SpkCnt{nSkpCnt}{TrainCondIndex(i)}(TrialIndex(i)));
                            end
                        end
                        % Change the label for all trials in the same condition
                        NewFactorData(Neu).AllFactors{IndFLC}(:, TargetFactorInd) = L;
                    end
                elseif DimMethod == 2
                    % Copy labels into the target factor
                    NewFactorData(Neu) = FactorData(Neu);
                    for IndFLC = UniqCondInd
                        NewFactorData(Neu).AllFactors{IndFLC}(:, TargetFactorInd) = Label(TrainCondIndex == IndFLC);
                    end
                elseif DimMethod == 3
                    % Swap trials based on labels and rules
                    for L = UniqLabel
                        for i = find(Label == L)'
                            % Determine the first dimension factor value
                            FirstDimFactorVal = FactorData(Neu).AllFactors{TrainCondIndex(i)}(TrialIndex(i), TargetFactorInd1ndD);

                            % Determine where the current trial belongs
                            IndFLC = find(FactorLevelComb(:, 1) == FirstDimFactorVal & FactorLevelComb(:, 2) == FactorLevelComb(TrainCondIndex(i), 2));

                            % Copy data and factors across the matrix
                            NewFactorData(Neu).data{IndFLC} = cat(1, NewFactorData(Neu).data{IndFLC}, FactorData(Neu).data{TrainCondIndex(i)}(TrialIndex(i), :));
                            NewFactorData(Neu).AllFactors{IndFLC} = cat(1, NewFactorData(Neu).AllFactors{IndFLC}, FactorData(Neu).AllFactors{TrainCondIndex(i)}(TrialIndex(i), :));

                            % Update for special cases (spike count periods)
                            for nSkpCnt = 1:size(ClassifierOpts.SpkCountPeriod, 1)
                                NewFactorData(Neu).FactorData_SpkCnt{nSkpCnt}{IndFLC} = cat(1, NewFactorData(Neu).FactorData_SpkCnt{nSkpCnt}{IndFLC}, FactorData(Neu).FactorData_SpkCnt{nSkpCnt}{TrainCondIndex(i)}(TrialIndex(i)));
                            end
                        end
                    end
                end
            end
        end
       
        function combinations=CreatPairCombsCell(~,Mat1,Mat2) % create paired combination of values for two matrices
            % Sample cell arrays with numeric values
            %             cell1 = {1, 2, 3};
            %             cell2 = {4, 5};

            % Initialize an empty cell array to store the combinations
            combinations = [];

            % Loop through the elements in cell1 and cell2
            for i = 1:length(Mat1)
                for j = 1:length(Mat2)
                    % Create a pair by combining one element from cell1 and one from cell2
                    pair = [Mat1(i), Mat2(j)];

                    % Add the pair to the combinations cell array
                    combinations = [combinations; pair];
                end
            end
        end
        function combinations=CreatTripleCombsCell(~,Mat1,Mat2,Mat3) % create Triple combination of values for three matrices
            % Sample cell arrays with numeric values
            %             cell1 = {1, 2, 3};
            %             cell2 = {4, 5};cell3 = {4, 5 6};

            % Initialize an empty cell array to store the combinations
            combinations = [];

            % Loop through the elements in cell1 and cell2
            for i = 1:length(Mat1)
                for j = 1:length(Mat2)
                    for k = 1:length(Mat3)
                        % Create a triple by combining one element from
                        % cell1 and one from cell2 and cell3
                        pair = [Mat1(i), Mat2(j) Mat3(k)];

                        % Add the pair to the combinations cell array
                        combinations = [combinations; pair];
                    end
                end
            end
        end
        function chunks=ChunkMat(~,yourMatrix,maxChunkSize)
            % Your matrix (replace this with your actual matrix)
 
            % Maximum chunk size
         %   maxChunkSize = 2000;

            % Get the size of the matrix
            [m] = length(yourMatrix);

            % Calculate the number of chunks
            numChunks = ceil(m / maxChunkSize);

            % Create a cell array to store the chunks
            chunks = cell(1, numChunks);

            % Split the matrix into chunks
            for i = 1:numChunks
                startIdx = (i - 1) * maxChunkSize + 1;
                endIdx = min(i * maxChunkSize, m);
                chunks{i} = yourMatrix(startIdx:endIdx);
            end


        end
        
        function [chunks,UniqFactor,Reconstrcut]=ChuckMatValue(~,yourMatrix,maxChunkSize)

            SortedValues=sort(yourMatrix);
            Remaining=mod(SortedValues,maxChunkSize);
            Factor=floor(SortedValues/maxChunkSize);
            UniqFactor=unique(Factor);
            chunks=arrayfun(@(x) Remaining(Factor==x),UniqFactor,'UniformOutput',0);
            % reconstruct the orignial data 
            Reconstrcut=sort(cell2mat(arrayfun(@(x) chunks{x}+maxChunkSize*UniqFactor(x),1:length(UniqFactor),'UniformOutput',0)));
            if sum(sort(yourMatrix)==Reconstrcut)~=length(yourMatrix);error('check this code');end
        end
        function GetWorkspaceVars(~)% shows all of the variables in the workspace
            workspace_vars = whos;
            for i = 1:length(workspace_vars)
                var_name = workspace_vars(i).name;
                var_size_bytes = workspace_vars(i).bytes;
                var_size_MB = var_size_bytes / 1e6; % Convert bytes to megabytes
                fprintf('Variable: %s, Size: %.2f MB\n', var_name, var_size_MB);
            end
        end
        function [MeanProjDiag,SumProjDiag,MeanProjDiagPos,MeanProjDiagNeg]=ProjectontoMatDiag(obj,Mat,AntiDiagFlag) % projects the values onto the diagonal of a matrix
            % if you want the anti diagonal elements of the matrix then use the AntiDiagFlag
           
            [m,n]=size(Mat);
            MeanProjDiag=arrayfun(@(k) mean(obj.getDiagonalElementsofMat(Mat,k,AntiDiagFlag)),(-m+1):(n-1));
            SumProjDiag =arrayfun(@(k) mean(obj.getDiagonalElementsofMat(Mat,k,AntiDiagFlag)),(-m+1):(n-1));
            % get mean of the positive and negative values 
            k=1;
            for i=(-m+1):(n-1)
                DiagVals=diag(Mat,i);
                MeanProjDiagPos(k)=mean(DiagVals(DiagVals>0));
                MeanProjDiagNeg(k)=mean(DiagVals(DiagVals<=0));
                k=k+1;
            end
            MeanProjDiagPos(isnan(MeanProjDiagPos))=0;
            MeanProjDiagNeg(isnan(MeanProjDiagNeg))=0;
        end
        function p_value=BinomialTest(~,numTrials,numSuccess,alpha) %runs binomial test on the data based on total and number of success
            % Example data
            %             numTrials = 50;  % Total number of trials
            %             numSuccess = 10; % Number of successes (selective cells)
            %             alpha = 0.05;    % Significance level

            % Binomial test
            p_value =1- binocdf(numSuccess, numTrials, alpha);
        end
        function DiagonalElements = getDiagonalElementsofMat(~,matrix, k,AntiDiagFlag)
            % if you want the anti diagonal elements of the matrix then use the AntiDiagFlag

            if AntiDiagFlag
                DiagonalElements = diag(fliplr(matrix), k);
            else
                DiagonalElements = diag((matrix), k);
            end
        end

    end
end