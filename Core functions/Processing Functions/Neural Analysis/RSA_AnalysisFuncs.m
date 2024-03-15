classdef RSA_AnalysisFuncs
    %RSA_ANALYSISFUNCS functions to do RSA analysis on Neural Data
    % this is getting a lot of help from Flesch et al 2020 code
    % Flesch T, Juechems K, Dumbalska T, Saxe A, Summerfield C: Orthogonal representations for robust context-dependent task performance in brains and neural networks. Neuron 2022, 110:1258-1270.e11.
    
    
    properties
        RDMdistancemetric='euclidean' %'correlation'  refer to pdist help
        fmincon_N_ITERS=500;
        % parametric model options
        ParamMdlRDM_type % type of ParamMdlRDM
        ParamMdlRDM_OrthogonalContexts; % 0 if the two contexts are about same features 1 if the two contexts are about orthogonal features
        ParamMdlRDM_FreezRot % are we freezing rotation?
        ParamMdlRDM_FreezCompression % are we freezing compression
        ParamMdlRDM_FreezBias % are we freezing bias
        ParamMdlRDM_TestFreezMdl % are we tesing all of the freeze conditions
        ParamMdlRDM_FreezMdlList={'FreezRot'} %,'FreezCompression','FreezBias'}; % list of freeze models
        ParamMdlRDM_RotBounds=[-180 180]; % limits of rotation       
        % 
        Correction_RotStep=10; % degrees for dividing rotations
        % plotting options
        MDSview=[0 30] % view angle of the MDS 3D plot
        CorrectRSAresults4Loss=1; % are we correcting for Loss
    end
    properties (Access=private)
        ManData=ManipulateData;
        TrialFunc=TrialFuncs;
        FigParams=fig_params;
    end
    methods
        function obj = RSA_AnalysisFuncs(varargin)
            %RSA_ANALYSISFUNCS Construct an instance of this class
            %   Detailed explanation goes here
            if nargin~=0 % initialize vars
                obj=obj.ParseParams(varargin) ; %%Process optional inputs
            end
        end
        function obj=ParseParams(obj,InputArgs)
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
        function results=RunRSAanalysisNeuralData(obj,FactorLvLInds,FactorLvLData,FactorLevels,Opts,varargin) % runs RSA analysis on Neural data
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            FactorLvLData=obj.OrganizeQuadrants(FactorLvLInds,FactorLvLData,Opts);
            % take mean of the data and concatinate them in a matrix
            MeanFactorLvLData=cell2mat(cellfun(@mean ,FactorLvLData,'UniformOutput',0)');
            
            % calculate RDM of Neural Data
            BrainRDM=obj.CalculateRDM(MeanFactorLvLData);
            
            % do RSA analysis
            results=obj.rsa_parammod_fmincon(BrainRDM);
            
            % now look if need to test the freeze conditions as well
            if obj.ParamMdlRDM_TestFreezMdl
                for i=1:length(obj.ParamMdlRDM_FreezMdlList) % run for each freeze condition
                    FreezMdlTxt=sprintf('ParamMdlRDM_%s',obj.ParamMdlRDM_FreezMdlList{i});
                    result_freez{i}=obj.rsa_parammod_fmincon(BrainRDM,FreezMdlTxt,1);
                    obj.(FreezMdlTxt)=0;
                end
                % results.ParamMdlRDM_FreezMdlList=obj.ParamMdlRDM_FreezMdlList;
                results.freezeConds=result_freez;
            end
        end
        
        function RDM=CalculateRDM(obj,X,varargin) % calculates Representation Dissimilarity metrix
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            %X Input data, specified as a numeric matrix of size m-by-n.Rows correspond to individual observations, and columns correspond to individual neurons.
            RDM=squareform(pdist(X,obj.RDMdistancemetric));
            
        end
        
        function  results=rsa_parammod_fmincon(obj,brainRDM,varargin)
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            %% rsa_roi_parammod_fmincon()
            %
            % fits fully parametrised model at single subject level
            % using fmincon (mse between constructed rdm and subject rdm)
            % saves the best fitting parameters (and their indices) for each subject
            
            % NOTE: don't forget to use -v7.3
            %
            % Timo Flesch, 2021,
            % Human Information Processing Lab,
            % Experimental Psychology Department
            % University of Oxford
            params.names.modelset = 'parametrised';
            params.num.RecSessions=1; % for now we are looking at all of the recording sessions together
            maskName='';
            N_ITERS=obj.fmincon_N_ITERS;
            
            % determine if contexts are orthogonal
            switch obj.ParamMdlRDM_type
                case 'ParallelColor'
                    obj.ParamMdlRDM_OrthogonalContexts=0;
                case 'OrthogonalShapeColor'
                    obj.ParamMdlRDM_OrthogonalContexts=1;
            end
            
            
            % %context
            ctx_min = 0;
            ctx_max = 2;
            
            % compression:
            comp_min = 0.0;
            comp_max = 1.0;
            % north can be shape and south can be color
            
            % rotation:
            bounds_ctx1_rel = [comp_min, comp_max];
            bounds_ctx2_rel = [comp_min, comp_max];
            bounds_ctx1_irrel = [comp_min, comp_max];
            bounds_ctx2_irrel = [comp_min, comp_max];
            bounds_rot = [obj.ParamMdlRDM_RotBounds(1), obj.ParamMdlRDM_RotBounds(2)];
            bounds_ctx = [ctx_min, ctx_max];
            constraints = [bounds_ctx1_rel; bounds_ctx2_rel; bounds_ctx1_irrel; bounds_ctx2_irrel; bounds_rot; bounds_ctx];
            
            betas_hat = nan(N_ITERS, params.num.RecSessions, 6);
            
            for it = 1:N_ITERS
                % initvals = [.0,.0,1,1,0,2]; % optim [.0,.0,1,1,0,2]
                
                %if contains(maskName,'EVC')
                %    initvals = [rand(1),rand(1),rand(1),rand(1),randsample(-90:90,1),2*rand(1)];
                %else
                initvals = [.6*rand(1),.6*rand(1),.5+(1-.5)*rand(1),.5+(1-.5)*rand(1),randsample(bounds_rot(1):bounds_rot(2),1),ctx_max*rand(1)];
                %end
                
                %%  reorganize  RDM from ePhys/trueMdl data
                y_sub = scale01(vectorizeRDM(brainRDM)');
                
                % eval(sprintf('loss = @(initvals) sum((y_sub - obj.ParamModelRDM_%s(initvals)).^2);',obj.ParamMdlRDM_type));
                loss = @(x) sum((y_sub - obj.ParamModelRDM_General(x)).^2);
                %  tic
                [betas, minloss] = fmincon(loss, initvals, [], [], [], [], constraints(:, 1), constraints(:, 2), [], optimoptions('fmincon', 'Display', 'off'));
                %  ElapsedTime(it)=toc;
                betas_hat(it, :) = betas;
                minlossAll(it,:) = minloss;
                
                % use general model for this
                % eval(sprintf('rdms(it, :, :) = squareform(obj.ParamModelRDM_%s(betas));',obj.ParamMdlRDM_type));
                rdms(it, :, :) = squareform(obj.ParamModelRDM_General(betas));
            end
            %  fprintf('\nMean Elapsed time:%0.4f',mean(ElapsedTime))
            
            % average over independent iterations:
            ModelRDM = squeeze(nanmean(rdms,1));
            Mean_betas_hat = [squeeze(nanmean(betas_hat,1))];
            
            % store results
            results=struct('BrainRDM',nan,'minlossAll',nan,'betas_hat',nan,'rdms',nan,'TrueBeta',nan);
            
            %   results.Mean_betas_hat = Mean_betas_hat;
            %   results.ModelRDM = ModelRDM;
            results.BrainRDM=brainRDM;
            % save raw values as well
            results.minlossAll=minlossAll;
            results.betas_hat=squeeze(betas_hat);
            results.rdms=rdms;
           % if obj.CorrectRSAresults4Loss % get only the correct loss solutions
           %     results=obj.GetCorrectRSAresults4Loss(results);
           % end
            % if we testign freeze conditions then remove some of the field
            if  obj.ParamMdlRDM_FreezRot | obj.ParamMdlRDM_FreezCompression | obj.ParamMdlRDM_FreezBias
                results=obj.ManData.rmfieldExept(results,{'betas_hat','minlossAll'});
            end
            
            
        end
        function [RSAresults,Angle,MinLoss,AngleEdges]=GetCorrectRSAresults4Loss(obj,RSAresults,varargin) % looks at the results and selects the meaningful estimates based on the loss function
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            Angle=RSAresults.betas_hat(:,5);
            MinLoss=RSAresults.minlossAll;
            %            EpsAng=1; % 1deg of epsilon
            % MeanLoss=[ mean(RSAresults.minlossAll(Angle>=-EpsAng & Angle<=EpsAng)) mean(RSAresults.minlossAll(Angle<-EpsAng)) mean(RSAresults.minlossAll(Angle>EpsAng))];
            %[Loss,RegAngleInd]=min(MeanLoss);
            %             if RegAngleInd==1
            %                 CorrectEstimateInd=Angle>=-EpsAng & Angle<=EpsAng;
            %             elseif RegAngleInd==2
            %                 CorrectEstimateInd=Angle<-EpsAng;
            %             elseif RegAngleInd==3
            %                 CorrectEstimateInd=Angle>EpsAng;
            %             end
            opts.Method='minLoss'; % can be "minLoss" too
            if obj.CorrectRSAresults4Loss
                % get a histogram of angle values with resolution of 10 degrees
                EpsAngRes=obj.Correction_RotStep;
                [a,b,bin]=histcounts(Angle,obj.ParamMdlRDM_RotBounds(1):EpsAngRes:obj.ParamMdlRDM_RotBounds(2));
                AngleEdges=movmean(b,2,"Endpoints","discard");
                Bins=unique(bin);
                BinLoss=arrayfun(@(x) mean(RSAresults.minlossAll(bin==x)),Bins);
                [~,MinLossInd]=min(BinLoss);
                MinLoss=nan*ones(1,length(AngleEdges));
                MinLoss(Bins)=BinLoss; 
                switch opts.Method
                    case 'minLoss' % take minimum loss as metric
                          CorrectEstimateInd=bin==Bins(MinLossInd);
                    case  'MaxNumber' % take maximum number of repetitions
                        [~,MaxNumInd]=max(a);
                        CorrectEstimateInd=bin==MaxNumInd;
                end
            else
                CorrectEstimateInd=logical(ones(1,length(Angle)));
                BinLoss=mean(RSAresults.minlossAll(CorrectEstimateInd));
            end
            
            if isfield(RSAresults,'rdms')
                RSAresults.rdms=squeeze(mean(RSAresults.rdms(CorrectEstimateInd,:,:),1));
            end
            RSAresults.betas_hat=squeeze(mean(RSAresults.betas_hat(CorrectEstimateInd,:),1));
            RSAresults.minlossAll=mean(RSAresults.minlossAll(CorrectEstimateInd),1);            
          %  obj.ShowAllMDS(RSAresults,[],[],Angle,MinLoss,AngleEdges);
          %  pause
        end
        %% Different model RDMs
        function [x_sub,rdm] = ParamModelRDM_OrthogonalShapeColor(obj,theta,varargin) % Parametrical Model RDM based on our task for two rules only that are orthogonal in there stimulus representation
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % based on the paper ‘north garden’’/blue rectangle//Context B
            % respond to leafiness ;Shape
            % versus ‘‘south garden’’/orange rectangle//Context A),
            % branchiness ;Color
            % Task 1 respVect(1:4,:) the relevant is shape the irrelevant is
            % color; Task 2 the relevant is color and irrelavant is shape
            c_rel_north = theta(1);
            c_rel_south = theta(2);
            c_irrel_north = theta(3);
            c_irrel_south = theta(4);
            a1 = 0;
            a2 = theta(5);
            ctx = theta(6);
            % note: north=90 and south=0 are the optimal, i.e. ground-truth
            % values % north=color  south=shape
            [b, l] = meshgrid([-1 1], [-1 1]); % -1 and 1 are the categories of shape and color
            b = b(:); % color
            l = l(:); % shape
            % compress irrelevant dimension:
            respVect = [[(1 - c_irrel_north) .* b, (1 - c_rel_north) .* l]; [(1 - c_rel_south) .* b, (1 - c_irrel_south) .* l]]; %l=north,b=south
            % rotate vector
            respVect(1:4, :) = respVect(1:4, :) * [cos(deg2rad(a1)), -sin(deg2rad(a1)); sin(deg2rad(a1)), cos(deg2rad(a1))];
            respVect(5:end, :) = respVect(5:end, :) * [cos(deg2rad(a2)), -sin(deg2rad(a2)); sin(deg2rad(a2)), cos(deg2rad(a2))];
            respVect = [respVect [zeros(length(b), 1); ctx + zeros(length(b), 1)]];
            rdm = squareform(pdist(respVect));
            %   rdm = helper_expandRDM(rdm,n_runs);
            x_sub = scale01(vectorizeRDM(rdm)');
        end
        function [x_sub,rdm] = ParamModelRDM_ParallelColor(obj,theta,varargin) % Parametrical Model RDM based on our task for two rules only that are Parallel in there stimulus representation
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % based on the paper ‘north garden’’/blue rectangle//Context B
            % respond to leafiness ;Color
            % versus ‘‘south garden’’/orange rectangle//Context A),
            % branchiness ;Color
            % both tasks respond to color
            % Task 1 respVect(1:4,:) the relevant is color the irrelevant is
            % shape; Task 2 the relevant is color and irrelavant is shape
            c_rel_north = theta(1);
            c_rel_south = theta(2);
            c_irrel_north = theta(3);
            c_irrel_south = theta(4);
            a1 = 0;
            a2 = theta(5);
            ctx = theta(6);
            % values % north=color  south=shape
            [b, l] = meshgrid([-1 1], [-1 1]); % -1 and 1 are the categories of shape and color
            b = b(:); % color
            l = l(:); % shape
            % compress irrelevant dimension:
            respVect = [[(1 - c_rel_north) .* b, (1 - c_irrel_north) .* l]; [(1 - c_rel_south) .* b, (1 - c_irrel_south) .* l]]; %l=north,b=south
            % rotate vector
            respVect(1:4, :) = respVect(1:4, :) * [cos(deg2rad(a1)), -sin(deg2rad(a1)); sin(deg2rad(a1)), cos(deg2rad(a1))];
            respVect(5:end, :) = respVect(5:end, :) * [cos(deg2rad(a2)), -sin(deg2rad(a2)); sin(deg2rad(a2)), cos(deg2rad(a2))];
            respVect = [respVect [zeros(length(b), 1); ctx + zeros(length(b), 1)]];
            rdm = squareform(pdist(respVect));
            %   rdm = helper_expandRDM(rdm,n_runs);
            x_sub = scale01(vectorizeRDM(rdm)');
        end
        function [x_sub,rdm,respVect] = ParamModelRDM_General(obj,theta,varargin) % a general parametric RDM that has rotation in 3 dimensions,compression and context bais
            % theta is a vector [c_rel_ctx1,c_rel_ctx2,c_irrel_ctx1,c_irrel_ctx2,a2,ctx]
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            opts.StimDim1=[-1 1]; % stimulus values for feature 1
            opts.StimDim2=[-1 1]; % stimulus values for feature 2
            opts.OrthogonalContexts=obj.ParamMdlRDM_OrthogonalContexts; % 0 if the two contexts are about same features 1 if the two contexts are about orthogonal features
            opts.nStim=2*length(opts.StimDim1);
            opts.nTotStim=2*opts.nStim;
            
            c_rel_ctx1 = theta(1);    % compression for relevant dimension of context 1
            c_rel_ctx2 = theta(2);    % compression for relevant dimension of context 2
            c_irrel_ctx1 = theta(3);  % compression for irrelevant dimension of context 1
            c_irrel_ctx2 = theta(4);  % compression for irrelevant dimension of context 2
            a1 = 0;                   % angle of rotation of context 1 to context 2
            a2 = theta(5);            % angle of rotation of context 2 to context 1
            ctx = theta(6);           % context bias between the two contexts
            
            % freeze components
            if obj.ParamMdlRDM_FreezRot
                a2=0;
            end
            if obj.ParamMdlRDM_FreezCompression
                c_rel_ctx1=0;  c_rel_ctx2=0; c_irrel_ctx1=0; c_irrel_ctx2=0;
            end
            if obj.ParamMdlRDM_FreezBias
                ctx=0;
            end
            
            % values
            [b, l] = meshgrid(opts.StimDim1, opts.StimDim2);
            b = b(:); % feature 1
            l = l(:); % feature 2
            
            % compress irrelevant dimension:
            if opts.OrthogonalContexts
                % feature 1  is irrelavant for context 1 and relavant for context 2; feature 2 is irrelavant for context 2 and relavant for context 1
                respVect = [[(1 - c_irrel_ctx1 ) .* b, (1 - c_rel_ctx1) .* l]; [(1 - c_rel_ctx2) .* b, (1 - c_irrel_ctx2) .* l]];
            else
                % feature 1 and 2 are is relavant for context 1 and context 2
                respVect = [[(1 - c_rel_ctx1   ) .* b, (1 - c_irrel_ctx1) .* l]; [(1 - c_rel_ctx2) .* b, (1 - c_irrel_ctx2) .* l]];
            end
            
            % rotate vector
            respVect(1:opts.nStim, :) = respVect(1:opts.nStim, :) * [cos(deg2rad(a1)), -sin(deg2rad(a1)); sin(deg2rad(a1)), cos(deg2rad(a1))];
            respVect((opts.nStim+1):end, :) = respVect((opts.nStim+1):end, :) * [cos(deg2rad(a2)), -sin(deg2rad(a2)); sin(deg2rad(a2)), cos(deg2rad(a2))];
            % add context bias
            respVect = [respVect [zeros(length(b), 1); ctx + zeros(length(b), 1)]];
            % calculate RDM
            rdm = squareform(pdist(respVect,obj.RDMdistancemetric));
            %   rdm = helper_expandRDM(rdm,n_runs);
            x_sub = scale01(vectorizeRDM(rdm)');
        end
        function [x_sub,rdm,respVect] = ParamModelRDM3DRot_General(obj,theta,varargin) % a general parametric RDM that has rotation,compression and context bais
            % theta is a vector [c_rel_ctx1,c_rel_ctx2,c_irrel_ctx1,c_irrel_ctx2,a2,ctx]
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            opts.StimDim1=[-1 1]; % stimulus values for feature 1
            opts.StimDim2=[-1 1]; % stimulus values for feature 2
            opts.OrthogonalContexts=obj.ParamMdlRDM_OrthogonalContexts; % 0 if the two contexts are about same features 1 if the two contexts are about orthogonal features
            opts.nStim=2*length(opts.StimDim1);
            opts.nTotStim=2*opts.nStim;
            
            c_rel_ctx1 = theta(1);    % compression for relevant dimension of context 1
            c_rel_ctx2 = theta(2);    % compression for relevant dimension of context 2
            c_irrel_ctx1 = theta(3);  % compression for irrelevant dimension of context 1
            c_irrel_ctx2 = theta(4);  % compression for irrelevant dimension of context 2
            a1 = 0;                   % angle of rotation of context 1 to context 2
            rotx = theta(5);            % angle of rotation of context 2 to context 1 around x axis
            roty = theta(6);            % angle of rotation of context 2 to context 1 around y axis
            rotz = theta(7);            % angle of rotation of context 2 to context 1 around z axis            
            ctxbias_x = theta(8);            % context bias between the two contexts with regards to x           
            ctxbias_y = theta(9);            % context bias between the two contexts            
            ctxbias_z = theta(10);           % context bias between the two contexts            
            
            % freeze components
            if obj.ParamMdlRDM_FreezRot
                b1=0;b2=0;b3=0;
            end
            if obj.ParamMdlRDM_FreezCompression
                c_rel_ctx1=0;  c_rel_ctx2=0; c_irrel_ctx1=0; c_irrel_ctx2=0;
            end
            if obj.ParamMdlRDM_FreezBias
                ctx=0;
            end
            
            % values
            [b, l] = meshgrid(opts.StimDim1, opts.StimDim2);
            b = b(:); % feature 1
            l = l(:); % feature 2
            
            % compress irrelevant dimension:
            if opts.OrthogonalContexts
                % feature 1  is irrelavant for context 1 and relavant for context 2; feature 2 is irrelavant for context 2 and relavant for context 1
                respVect = [[(1 - c_irrel_ctx1 ) .* b, (1 - c_rel_ctx1) .* l]; [(1 - c_rel_ctx2) .* b, (1 - c_irrel_ctx2) .* l]];
            elseif opts.OrthogonalContexts==0
               % feature 1 and 2 are is relavant for context 1 and context 2
               respVect = [[(1 - c_rel_ctx1   ) .* b, (1 - c_irrel_ctx1) .* l]; [(1 - c_rel_ctx2) .* b, (1 - c_irrel_ctx2) .* l]];
            elseif opts.OrthogonalContexts==-1
                % feature 1 and 2 are replcated in both rules 
                 respVect = [[(1 - c_rel_ctx1   ) .* b, (1 - c_irrel_ctx1) .* l];[(1 - c_rel_ctx1   ) .* b, (1 - c_irrel_ctx1) .* l]] ;
                 
            end
            
            % add the third dimension
            respVect = [respVect [zeros(length(b), 1); zeros(length(b), 1)]];
            if opts.OrthogonalContexts==-1 || opts.OrthogonalContexts==0 || opts.OrthogonalContexts==1
                respVect([opts.nStim opts.nStim*2],3)=respVect([opts.nStim opts.nStim*2],3)+0.0001;
            end
             
            % add context bias to each dimension
            respVect=respVect+[zeros(length(b), 3) ;[ctxbias_x+zeros(length(b), 1) ctxbias_y+zeros(length(b), 1) ctxbias_z+zeros(length(b), 1)]];
           
            %% rotate both vectors around y axis remove any 0 0 0 componenent 
%             ThetaRotBoth=30;
%             RotyBoth=[cos(deg2rad(ThetaRotBoth)) 0 -sin(deg2rad(ThetaRotBoth));0 1, 0;sin(deg2rad(ThetaRotBoth)),0, cos(deg2rad(ThetaRotBoth))];
%             respVect(1:opts.nStim, :) =[RotyBoth*respVect(1:opts.nStim, :)']';
%             respVect((opts.nStim+1):end, :) =[RotyBoth*respVect((opts.nStim+1):end, :)']';
%             
            % rotate vector
            % respVect(1:opts.nStim, :) = respVect(1:opts.nStim, :);
            Rotx=[1 0 0;0 cos(deg2rad(rotx)), sin(deg2rad(rotx));0 -sin(deg2rad(rotx)), cos(deg2rad(rotx))];
            Roty=[cos(deg2rad(roty)) 0 -sin(deg2rad(roty));0 1, 0;sin(deg2rad(roty)),0, cos(deg2rad(roty))];
            Rotz=[cos(deg2rad(rotz)) sin(deg2rad(rotz)) 0;-sin(deg2rad(rotz)) cos(deg2rad(rotz)) 0;0 ,0, 1];                     
            
            respVect((opts.nStim+1):end, :) =[Rotz*Roty*Rotx*respVect((opts.nStim+1):end, :)']';
            
                               
            % calculate RDM
            rdm = squareform(pdist(respVect,obj.RDMdistancemetric));
            %   rdm = helper_expandRDM(rdm,n_runs);
            x_sub = scale01(vectorizeRDM(rdm)');
        end

        %% Demo functions
        function ParametrizedRSAanaDemo(obj,varargin) % Demo to fit simulated data to parametrized RDMs and test differnt hypothetis
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            obj.CorrectRSAresults4Loss=1; % don't correct for loss
            %% Simulate a model that rotates a representation with respect to another representation
            Demo1.Rot=-180:10:180;h=[];Sp=[];
            RSAresults=struct('BrainRDM',nan,'minlossAll',nan,'betas_hat',nan,'rdms',nan,'TrueBeta',nan);
            PointOrder=[1 2 3 4 1];%[1 2 4 3 1]
            for r=1:length(Demo1.Rot)
                % construct our ground truth RDM and then fit it
                TrueBeta=[0.2 0.2 0.8 0.8 Demo1.Rot(r) 0.5];
                [~,GT_ModelRDM,respVect]=obj.ParamModelRDM_General(TrueBeta,'ParamMdlRDM_OrthogonalContexts',1);
                RSAresults(r)=obj.rsa_parammod_fmincon(GT_ModelRDM,'ParamMdlRDM_type','General',...
                    'ParamMdlRDM_OrthogonalContexts',1,'fmincon_N_ITERS',100);
                RSAresults(r).TrueBeta=TrueBeta;
                if 0 % let's check the analysis using PCA now
                    respVect=[0 0 0;1,0,0;1,1,0;0,1,0;0,0,0;cosd(Demo1.Rot(r)),0,sind(Demo1.Rot(r));...
                        cosd(Demo1.Rot(r)),1,sind(Demo1.Rot(r));0,1,0];
                    
                    PCAresults(r)=obj.DoPCAanalysisonMdlData(respVect);
                    cla
                    hold on
                    plot3(respVect(PointOrder,1),respVect(PointOrder,2),respVect(PointOrder,3))
                    plot3(respVect(PointOrder+4,1),respVect(PointOrder+4,2),respVect(PointOrder+4,3))
                    
                    plot3(respVect([1],1),respVect([1],2),respVect([1],3),'r*')
                    plot3(respVect([5],1),respVect([5],2),respVect([5],3),'r*')
                    title(sprintf('Inferred Angle:%0.1f, True:%i',PCAresults(r).PairAngled),Demo1.Rot(r))
                    view(30,30)
                    pause
                end
            end
            [h,Sp,mvFrame]=obj.ShowAllMDS(RSAresults,h,Sp);
            [~,~,MDSFigFileName]=obj.ManData.GetFileName(['Subspace'],['_Demo_RSA_AngRot' num2str(obj.ParamMdlRDM_RotBounds(2))],'SaveInResults',1,'WantedDate','ALL');
            obj.FigParams.MakeMovieFromFrames(mvFrame,1,MDSFigFileName)
            
            %% Simulate a model that rotates a representation, changes bias and changes compression
            Demo1.Rot=-180:30:180;h=[];Sp=[];nRots=length(Demo1.Rot);
            RSAresults=struct('BrainRDM',nan,'minlossAll',nan,'betas_hat',nan,'rdms',nan,'TrueBeta',nan);
            for r=1:length(Demo1.Rot)
                % construct our ground truth RDM and then fit it
                TrueBeta=[r/(nRots+1) 0.5 0.9 0.6 Demo1.Rot(r) 1-(r/(nRots+1))];
                GT_ModelRDM=obj.ParamModelRDM_General(TrueBeta,'ParamMdlRDM_OrthogonalContexts',1);
                RSAresults(r)=obj.rsa_parammod_fmincon(GT_ModelRDM','ParamMdlRDM_type','General',...
                    'ParamMdlRDM_OrthogonalContexts',1,'fmincon_N_ITERS',50);
                RSAresults(r).TrueBeta=TrueBeta;
            end
            [h,Sp,mvFrame]=obj.ShowAllMDS(RSAresults,h,Sp);
            [~,~,MDSFigFileName]=obj.ManData.GetFileName(['Subspace'],['_Demo_RSA_AngRot' num2str(obj.ParamMdlRDM_RotBounds(2)) '_Comp_Bias'],'SaveInResults',1,'WantedDate','ALL');
            obj.FigParams.MakeMovieFromFrames(mvFrame,1,MDSFigFileName)
            
            
        end
        %% preprocessing functions
        function FactorLvLData=OrganizeQuadrants(obj,FactorLvLInds,FactorLvLData,Opts,varargin) % organzied Quadrants data to match the RDM
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % if our main target factors are quadrants then we need to fix
            % the order of objects to match the RDM model
            % the Quadrants are organized as [Shape Color]
            %[1 1][2 2][1 2][2 1]->[RB,GT,GB,RT]
            % the model RDM is organized as
            %[1 1][2 1][1 2][2 2]->[RB,RT,GB,GT]
            % so the transformation we need is [1,4,3,2] for each rule
            if ~strcmp(Opts.MainTargetFactor,'Quadrants');return;end
            nRule=1:length(unique(FactorLvLInds));
            TransInds=cell2mat(arrayfun(@(x) [1,4,3,2]+(x-1)*4,nRule,'UniformOutput',0));
            FactorLvLData=FactorLvLData(TransInds);
        end
        %% PCA functions 
        function [PairAngled,PairAngler,CosTheta,Pairs]=CalculateAnglesBetSubspaces(obj,PCAscores,Factors,nD,varargin) % calculates angles between the planes that are fit to each condition
            %@PCAscores is a cell with scores of different conditions we want to compare  
            %@Factors cellarray of factors if we are using quadrants then we collapse based on
            %Shape or Color and calculate angle between lines
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
           
            NConds=1:length(PCAscores);
                % fit a plane to data per second condition now
                [PCA,PCAcoeff,OriginPlane,FittedPlane,PlaneEq]=arrayfun(@(x) obj.FitPlane2Data(PCAscores{x}(:,1:3)), NConds,'UniformOutput',0);
                
                % calculate perpenducular vector
                PerpendVec=cellfun(@(x) cross(x(1,:),x(2,:)),PCAcoeff,'UniformOutput',0);
                
                % calculate angle between pairs of conditions
                PairsInd=nchoosek(NConds,2);
                [PairAngled,PairAngler,CosTheta]=arrayfun(@(x) obj.ManData.GetAngleBetVectors(PerpendVec{PairsInd(x,1)},PerpendVec{PairsInd(x,2)}),1:size(PairsInd,1));
                Pairs=NConds(PairsInd);
            
        end
        function [PCA,PCAcoeff,OriginPlane,FittedPlane,PlaneEq,OrthoVector,OriginPlaneEq]=FitPlane2Data(obj,data,varargin) % finds the best 2D plane that fits the data using PCA
            %@ data points where we want to find the best plane fitting them
            %@ PCAcoeff vectors specifying the plane
            %@ OriginPlane: plane that passes through origin
            %@ FittedPlane: plane that is shifted to the mean of the points
            %@ PlaneEq: equation of the plane using
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs

            rMean=squeeze(mean(data,1)); % take the mean across factor values

            % take PCA per each time point now
            [PCA.coeff,PCA.score,PCA.latent,PCA.tsquared,PCA.explained,PCA.mu]=...
                pca(data-repmat(rMean,[size(data,1) 1]));
            % now the plane is the two principal component vectors , orogin and their sum
            PCAcoeff=PCA.coeff(:,1:2)';
            NDims=size(PCAcoeff,2);
            if NDims==4
                % find orthogonal complement plane of these vectors
                % this will be the null space of Ax=0 where the two vectors
                % PCAcoeff are the row spaces 
                NullSpace=null(PCAcoeff);
                % PCAcoeff*NullSpace=0 should be
                % if you want to solve it with linear algebra way 
                % RREF=rref([PCAcoeff;zeros(2,4)]);
                % AlgebraNull=[-RREF(1:2,3)' 1 0;-RREF(1:2,4)' 0 1];
                OrthoVector=NullSpace';
                OriginPlane=[];FittedPlane=[];PlaneEq=[];OriginPlaneEq=[];
            elseif NDims==3
                % now if the vectors are 3D
                sumPCAvectors=sum(PCAcoeff,1);
                OriginPlane=[0 0 0;PCAcoeff(1,:);sumPCAvectors;PCAcoeff(2,:)];
                % add mean to origin palce to get the fitted plane
                FittedPlane=rMean+OriginPlane;
                % get the equation for this plane as well
                % calculate orthogonal vector
                OrthoVector=cross(PCAcoeff(1,:),PCAcoeff(2,:));
                PlaneEq=@(x,y) (OrthoVector(1)*(rMean(1)-x)+OrthoVector(2)*(rMean(2)-y))./OrthoVector(3) + rMean(3);
                OriginPlaneEq=@(x,y) (OrthoVector(1)*(0-x)+OrthoVector(2)*(0-y))./OrthoVector(3) + 0;% plane that passes through origin
            else
                error('This function only works with 3 and 4 dimensions')
            end
        end
        function PCAresults=DoPCAanalysisonMdlData(obj,respVect) % Cross checks resutls with PCA analysis
           
            
            [PCAresults.PairAngled,~,PCAresults.CosTheta]=obj.CalculateAnglesBetSubspaces({respVect(1:4,:) respVect(5:8,:)},[]); % calculates angles between the planes that are fit to each condition
           % [PCAresults.SubspaceCompression1]=obj.CalculateCompressionBetSubspaces(SubspaceSpkCnt{IndResp}.CondScore{1}(1));  % Calculates compression for each subspace
           % [PCAresults.SubspaceCompression2] =obj.CalculateCompressionBetSubspaces(SubspaceSpkCnt{Tim}.CondScore{2}(1));      % Calculates compression for each subspace
            
        end
        %% plot functions
        function ShowMDS(obj,rdm,ndims,Sp,CtxLegTxt,varargin)% % projects rdm into ndims-D space and visualises results as scatter plot
            %@rdm is the rdm of the data
            %@ndims number for dimensions should be 3
            %@CtxLegTxt is the txt for legend
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            ctxMarkerEdgeCol = [0,11,229; 229,143,0]./255;
            nObjsPerRule=2;nRule=2;
            nStims = 8;
            if isempty(CtxLegTxt);CtxLegTxt={'Ctx1','Ctx2'};end
            
            try
                xyz = mdscale(rdm,ndims,'Criterion','strain');
            catch
                xyz = mdscale(rdm,ndims,'Criterion','strain','Start','random');
            end
            ctxMarkerCol     = 'w';
            ctxMarkerSize = 20;
            
            scat_leafiness = [20:-(10/(nObjsPerRule-1)):10];
            scat_branchiness= [1 0 0;0 1 0];%(RDMcolormap_hack(nObjsPerRule));
            
            
            [b,l] = meshgrid(1:nObjsPerRule,1:nObjsPerRule);
            b = b(:);l = l(:);
            x = xyz(:,1);
            y = xyz(:,2);
            z = xyz(:,3);
            if ~isempty(Sp);subplot(Sp);end
            cla
            % plot grid connecting adjacent points
            h1=obj.DisplayMDSgrid(xyz(1:nStims/2,:),ctxMarkerEdgeCol(1,:));
            h2=obj.DisplayMDSgrid(xyz(nStims/2+1:end,:),ctxMarkerEdgeCol(2,:));
            
            % scatterplot (size= dim1, col= dim2)
            [b,l] = meshgrid(1:nObjsPerRule,1:nObjsPerRule);
            b = [b(:);b(:)];l = [l(:);l(:)];
            for ii = 1:nStims/2
                %                 plot3(x(ii),y(ii),z(ii),'MarkerFaceColor',ctxMarkerCol,'MarkerSize',ctxMarkerSize, ...
                %                     'Marker','square','MarkerEdgeColor',ctxMarkerEdgeCol(1,:),'LineWidth',2);
                hold on;
                plot3(x(ii),y(ii),z(ii),'MarkerFaceColor',scat_branchiness(b(ii),:),'MarkerSize',scat_leafiness(l(ii)), ...
                    'Marker','diamond','MarkerEdgeColor','None');
            end
            for ii = nStims/2+1:nStims
                %                 plot3(x(ii),y(ii),z(ii),'MarkerFaceColor',ctxMarkerCol,'MarkerSize',ctxMarkerSize, ...
                %                     'Marker','square','MarkerEdgeColor',ctxMarkerEdgeCol(2,:),'LineWidth',2);
                hold on;
                plot3(x(ii),y(ii),z(ii),'MarkerFaceColor',scat_branchiness(b(ii),:),'MarkerSize',scat_leafiness(l(ii)), ...
                    'Marker','diamond','MarkerEdgeColor','None');
            end
            xlabel('MDS1');ylabel('MDS2');zlabel('MDS3');
            legend([h1 h2],CtxLegTxt,'Location','best')
            view(obj.MDSview(1),obj.MDSview(2))
            obj.FigParams.FormatAxes(gca);
        end
        function h=DisplayMDSgrid(obj,xyz,edgeCol,edgeWidth,varargin)% draws lines between adjacent points on a plane
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            nObjsPerRule=2;
            
            if ~exist('edgeCol','var')
                edgeCol = [1,1,1].*.8;
            end
            
            if ~exist('edgeWidth','var')
                edgeWidth = 1.5;
            end
            
            [l,b] = meshgrid(1:nObjsPerRule,1:nObjsPerRule);
            b = b(:);
            l = l(:);
            bl = [b,l];
            if ~exist('xyz','var')
                xyz = [b,l];
            end
            
            for i = 1:(nObjsPerRule-1)
                for j = 1:(nObjsPerRule-1)
                    p1 = xyz(bl(:,1)==i & bl(:,2)==j,:);
                    p2 = xyz(bl(:,1)==i+1 & bl(:,2)==j,:);
                    plot3([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],'LineWidth',edgeWidth,'Color',edgeCol);
                    hold on;
                    p2 = xyz(bl(:,1)==i & bl(:,2)==j+1,:);
                    plot3([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],'LineWidth',edgeWidth,'Color',edgeCol);
                end
            end
            
            i = (nObjsPerRule);
            for j = 1:(nObjsPerRule-1)
                p1 = xyz(bl(:,1)==i & bl(:,2)==j,:);
                p2 = xyz(bl(:,1)==i & bl(:,2)==j+1,:);
                plot3([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],'LineWidth',edgeWidth,'Color',edgeCol);
            end
            
            j = (nObjsPerRule);
            for i = 1:(nObjsPerRule-1)
                p1 = xyz(bl(:,1)==i & bl(:,2)==j,:);
                p2 = xyz(bl(:,1)==i+1 & bl(:,2)==j,:);
                h=plot3([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],'LineWidth',edgeWidth,'Color',edgeCol);
            end
        end
        function [h,Sp,mvFrame]=ShowAllMDS(obj,RSAresults,h,Sp,Angle,MinLoss,AngleEdges,varargin) % shows the results for all of the mds
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            Ntrlrng=length(RSAresults);
            
            if isempty(h)
                h=obj.FigParams.RenderFigure(1,[]);
                [h,Sp]=obj.FigParams.RenderSubplots(2,5,h{1},[]);
            end
            
            
            CtxLegTxt={'Ctx1','Ctx2'};
            for ntrlrng=1:Ntrlrng
                if ~exist('Angle','var')
                    % correct results for this trial range now
                    [RSAresults(ntrlrng),Angle,MinLoss,AngleEdges]=obj.GetCorrectRSAresults4Loss(RSAresults(ntrlrng));
                end
                if isnan(RSAresults(ntrlrng).TrueBeta);RSAresults(ntrlrng).TrueBeta=RSAresults(ntrlrng).betas_hat;end
                %obj.ShowMDS(RSAresults(ntrlrng).BrainRDM,3,Sp(1),CtxLegTxt);               
                % show Brain/TrueMdl RDM
                obj.ShowMDS(RSAresults(ntrlrng).BrainRDM,3,Sp(1),CtxLegTxt,'MDSview',[0 90]);
                obj.ShowMDS(RSAresults(ntrlrng).BrainRDM,3,Sp(2),CtxLegTxt,'MDSview',[90 0]);
                title('Brain/TrueMdl RDM')

                %    obj.ShowMDS(RSAresults(ntrlrng).BrainRDM,3,Sp(4),CtxLegTxt,'MDSview',[0 0]);
                % show model RDM
                % obj.ShowMDS(RSAresults(ntrlrng).rdms,3,Sp(3),CtxLegTxt);
                obj.ShowMDS(RSAresults(ntrlrng).rdms,3,Sp(3),CtxLegTxt,'MDSview',[0 90]);
                obj.ShowMDS(RSAresults(ntrlrng).rdms,3,Sp(4),CtxLegTxt,'MDSview',[90 0]);
                title('Model RDM')

                % obj.ShowMDS(RSAresults(ntrlrng).rdms,3,Sp(8),CtxLegTxt,'MDSview',[0 0]);
                                
                subplot(Sp(5));cla
                hold on
                bar(AngleEdges,MinLoss);
                plot(AngleEdges,MinLoss,'r*');
                title('Mean Loss')
                axis square
                ylabel('Loss')
%                 xticklabels({'Loss'})
%                 xtickangle(-45)
                
                subplot(Sp(6));cla
                hist(Angle,-180:10:180);hold on;axis square
                title('Dist Angle')
                xlabel('Angle')
                ylabel('Count')
                v=axis;
                plot([RSAresults(ntrlrng).TrueBeta(5) RSAresults(ntrlrng).TrueBeta(5)],[v(3) v(4)],'r','linewidth',1);
                legend('Ang Hat','True Ang','Location','best')
                
                subplot(Sp(7))
                hold on
                plot(ntrlrng,RSAresults(ntrlrng).betas_hat(5),'r*')
                plot(ntrlrng,RSAresults(ntrlrng).TrueBeta(5),'k*')
                ylabel('Angle')
                xlabel('Step')
                legend('Ang Hat','True Ang','Location','best')
                title([{'Comp estimated'}; {'vs true rot'}])
                axis square
                
                subplot(Sp(8));
                bar([RSAresults(ntrlrng).betas_hat([1:4 6])',RSAresults(ntrlrng).TrueBeta([1:4 6])']);axis square
                legend('Beta Hat','True Beta','Location','best')
                xticklabels({'Cmp Rel Ctx1','Cmp Rel Ctx2','Cmp irRel Ctx1','Cmp irRel Ctx2','Ctx Bias'})
                xtickangle(-45)
                title([{'Comp estimated'};{'vs true betas'}])
                
                subplot(Sp(9));hold on
                plot(ntrlrng,RSAresults(ntrlrng).betas_hat(1),'r*')
                plot(ntrlrng,RSAresults(ntrlrng).TrueBeta(1),'rd')
                axis square
                legend('Beta CmpRelCtx1','True CmpRelCtx1','Location','best')
                title([{'Comp of estimated'}; {'vs true Compressions'}])
                xlabel('Step')
                ylabel('Compression')
                
                subplot(Sp(10));hold on
                plot(ntrlrng,RSAresults(ntrlrng).betas_hat(6),'k*')
                plot(ntrlrng,RSAresults(ntrlrng).TrueBeta(6),'kd')
                axis square
                legend('Beta BiasCtx','True BiasCtx','Location','best')
                title([{'Comp of estimated'}; {'vs true Bias'}])
                ylabel('Bias')
                xlabel('Step')
                % get frames for a movie
                mvFrame(ntrlrng) = getframe(gcf);
                clear Angle
            end
        end
        
    end
end

