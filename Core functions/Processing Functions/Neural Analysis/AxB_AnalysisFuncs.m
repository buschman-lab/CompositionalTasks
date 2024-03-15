classdef AxB_AnalysisFuncs
    %AXB_ANALYSISFUNCS functions to perfrom Ax=B analysis
    %   Detailed explanation goes here
    
    properties
        PlotAxBDetails=1; % are we plotting AxB transformation details
        PlotAxBCircSubspace=1; % are we plotting circles that span the subspaces
        AxBGeneralizeTest=0 % are we testing the generalization of AxB
        ProjectXY2Plane=0; % are we projecting X and Y into a plane before proceeding with analysis
        PointOrder=[1 2 4 3 1]; % order of objects 
        UsePCA4AxB=0; % are we using PCA before doing the AxB analysis 
        SubtractFactorMean=1; % subtract factor mean before doing AxB
    end
    properties (Access=private)
        ManData=ManipulateData;
        TrialFunc=TrialFuncs;
        FigParams=fig_params;
        RSAana=RSA_AnalysisFuncs;
    end
    
    methods
        function obj = AxB_AnalysisFuncs(varargin)
            %AxB_ANALYSISFUNCS Construct an instance of this class
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
        
        %% Demo functions
        function ParametrizedAxBanaDemoOld(obj,varargin) % Demo to fit simulated data with AxB and test differnt hypothetis
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            %% Simulate a model that rotates a representation with respect to another representation
            Demo1.Rot=-180:10:180;
            PointOrder=[1 2 4 3 1];%[1 2 3 4 1];%[1 2 4 3 1]
            for r=1:length(Demo1.Rot)
                % construct our ground truth RDM and then fit it
                %x y z             Bx By Bz
                [~,~,RespVect]=obj.RSAana.ParamModelRDM3DRot_General([0.2 0.2 0.8 0.8 0 0 Demo1.Rot(r) 0.0 0.0 0],...
                    'ParamMdlRDM_OrthogonalContexts',1);
                % add a tiny bit of noise to this so we avoid low rank matrixes
                %RespVect=RespVect+0.0000005*rand(size(RespVect));
                TrueComp=obj.CalFakeDataCompression(RespVect,PointOrder); % calculate ture compression
                
                % now decode what is in the space
                [AxBresults(r)] =obj.RunAxBanalysisNeuralData([],RespVect,[],[]);
                
                % use PCA to double check the results as well
                PCAresults(r)=obj.RSAana.DoPCAanalysisonMdlData(RespVect);
                
                RotAxis=AxBresults(r).RotAxis;
                RespVectExmp=RespVect(PointOrder(1)+4,:);
                PRotAxis=[RespVectExmp(1) RespVectExmp(2) -(RespVectExmp(1)*RotAxis(1)+RespVectExmp(2)*RotAxis(2))/RotAxis(3)];
                dot(RotAxis,PRotAxis);
                TransformPRotAxis=[PRotAxis*AxBresults(r).RotMat'];
                [azimuth_PRotAxis,elevation_PRotAxis,r_PRotAxis]=obj.ManData.ConvertCart2Spherical(PRotAxis(1),PRotAxis(2),PRotAxis(3));
                [azimuth_TransformPRotAxis,elevation_TransformPRotAxis,r_TransformPRotAxis]=obj.ManData.ConvertCart2Spherical(TransformPRotAxis(1),TransformPRotAxis(2),TransformPRotAxis(3));
                azimuth_PRot_Diff=azimuth_TransformPRotAxis-azimuth_PRotAxis;
                elevation_PRot_Diff=elevation_TransformPRotAxis-elevation_PRotAxis;
                
                N=360;
                RandVectsXY=[cosd(1:N)' sind(1:N)'];
                PRotAxisRand=[RandVectsXY [-(RandVectsXY(:,1)*RotAxis(1)+RandVectsXY(:,2)*RotAxis(2))/RotAxis(3)]];
                DotRand=arrayfun(@(x) dot(PRotAxisRand(x,:),RotAxis),1:N);
                
                TransformPRotRand=cell2mat(arrayfun(@(x) [PRotAxisRand(x,:)*AxBresults(r).RotMat']',1:N,'UniformOutput',0))';
                
                [azimuth_PRotAxisRand,elevation_PRotAxisRand]=arrayfun(@(x) obj.ManData.ConvertCart2Spherical(PRotAxisRand(x,1),PRotAxisRand(x,2),PRotAxisRand(x,3)),1:N,'uniformoutput',1);
                [azimuth_TransformPRotAxisRand,elevation_TransformPRotAxisRand]=arrayfun(@(x) obj.ManData.ConvertCart2Spherical(TransformPRotRand(x,1),TransformPRotRand(x,2),TransformPRotRand(x,3)),1:N,'uniformoutput',1);
                
                AzimuthRandDiff=azimuth_TransformPRotAxisRand-azimuth_PRotAxisRand;
                ElevationRandDiff=elevation_TransformPRotAxisRand-elevation_PRotAxisRand;
                
                RotRandVects=arrayfun(@(x) obj.ManData.GetAngleBetVectors(PRotAxisRand(x,:),TransformPRotRand(x,:)),1:N);
                
                %% now fit a plane to the X and project the vectors onto that plane
                [PCA,PCAcoeff,OriginPlane,FittedPlane,~,OrthoVector,OrgPlaneEq]=obj.RSAana.FitPlane2Data(RespVect(5:8,:));
                
                
                Ns=360;
                Base1=PCAcoeff(1,:);Base2=PCAcoeff(2,:);
                %   SubspaceXvects=[cosd(1:Ns)' sind(1:Ns)' OrgPlaneEq(cosd(1:Ns),sind(1:Ns))']; % generate vectors on this subspace
                SubspaceXvects=cell2mat(arrayfun(@(x) [cosd(x)*Base1+(sind(x))*Base2;cosd(x)*Base1-(sind(x))*Base2 ]',1:Ns/2,'uniformoutput',0))';
                DotSubVectsN=arrayfun(@(x) dot(SubspaceXvects(x,:),OrthoVector),1:Ns);
                if sum(DotSubVectsN>0.01);error('Plane equation is not correct');end
                % now transform the SubspaceXvects with rotation matrix and
                % project them back into this plane
                TransformSubspaceXvects=cell2mat(arrayfun(@(x) [SubspaceXvects(x,:)*AxBresults(r).RotMat']',1:Ns,'UniformOutput',0))';
                % project the back into the plane
                [~,ProjMat,projBackSubspaceXvects]=arrayfun(@(x) obj.ManData.ProjectInto3DSubspace(TransformSubspaceXvects(x,:)',PCAcoeff(1,:)',PCAcoeff(2,:)'),1:Ns,'UniformOutput',0);
                projBackSubspaceXvects=cell2mat(projBackSubspaceXvects)';
                Dot_projBackSubspaceXvects=arrayfun(@(x) dot(projBackSubspaceXvects(x,:),OrthoVector),1:Ns);
                if sum(Dot_projBackSubspaceXvects>0.01);error('Error in projecting back into plane');end
                % now calculate in-plane and out of plane rotations
                RotInplane_SubspaceVects=arrayfun(@(x) obj.ManData.GetAngleBetVectors(SubspaceXvects(x,:),projBackSubspaceXvects(x,:)),1:Ns);
                
                % out of plane rotation would be the angle between the
                % projected vector and normal vector of the plane
                [RotOutplane_SubspaceVects,~,DotTransOrtho]=arrayfun(@(x) obj.ManData.GetAngleBetVectors(OrthoVector,TransformSubspaceXvects(x,:)),1:Ns);
                RotOutplane_SubspaceVects(RotOutplane_SubspaceVects<=90)=90-RotOutplane_SubspaceVects(RotOutplane_SubspaceVects<=90);
                RotOutplane_SubspaceVects(RotOutplane_SubspaceVects>90)=RotOutplane_SubspaceVects(RotOutplane_SubspaceVects>90)-90;
                
                % calculate the angle between plane's normal vector and
                % rotation axis: that is the out of plane rotation
                RotOutplane=obj.ManData.GetAngleBetVectors(OrthoVector,RotAxis);
                
                %   %% do the same process for rotation axis
                %   [~,ProjMat_RotAxis,projBackRotAxis]=obj.ManData.ProjectInto3DSubspace(RotAxis,PCAcoeff(1,:)',PCAcoeff(2,:)');
                %   RotInplane_Rot=arrayfun(@(x) obj.ManData.GetAngleBetVectors(SubspaceXvects(x,:),projBackSubspaceXvects(x,:)),1:Ns);
                
                subplot(1,4,[1 2])
                cla
                hold on
                plot3(RespVect(PointOrder,1),RespVect(PointOrder,2),RespVect(PointOrder,3))
                plot3(RespVect(PointOrder+4,1),RespVect(PointOrder+4,2),RespVect(PointOrder+4,3))
                
                plot3(RespVect([1],1),RespVect([1],2),RespVect([1],3),'r*')
                plot3(RespVect([5],1),RespVect([5],2),RespVect([5],3),'r*')
                % plot Rotation vector from origin as well
                xlabel('X Axis'); ylabel('Y Axis'); zlabel('Z axis')
                view(30,30)
                
                if 1
                    % plot a sample vector and its transformation
                    RotAxisLeg=quiver3(0,0,0,AxBresults(r).RotAxis(1),AxBresults(r).RotAxis(2),AxBresults(r).RotAxis(3),'g');% rotation axis
                    NormAxisLeg=quiver3(0,0,0,OrthoVector(1),OrthoVector(2),OrthoVector(3),'r');% normal vector to plane
                    %  quiver3(0,0,0,PRotAxis(1),PRotAxis(2),PRotAxis(3),'r');
                    %  quiver3(0,0,0,TransformPRotAxis(1),TransformPRotAxis(2),TransformPRotAxis(3),'k');
                    %  plot3(FittedPlane(:,1),FittedPlane(:,2),FittedPlane(:,3),'g*')
                    % PointMarker=['*' cellfun(@(x) '.',1:Ns-1)];
                    OrgSubSp=plot3(SubspaceXvects(:,1),SubspaceXvects(:,2),SubspaceXvects(:,3),'r.');
                    TransFor= plot3(TransformSubspaceXvects(:,1),TransformSubspaceXvects(:,2),TransformSubspaceXvects(:,3),'g.');
                    ProjBck=plot3(projBackSubspaceXvects(:,1),projBackSubspaceXvects(:,2),projBackSubspaceXvects(:,3),'k.');
                    NTest=30; % give some examples
                    quiver3(0,0,0,SubspaceXvects(NTest,1),SubspaceXvects(NTest,2),SubspaceXvects(NTest,3),'r','LineStyle',':','LineWidth',4);
                    quiver3(0,0,0,TransformSubspaceXvects(NTest,1),TransformSubspaceXvects(NTest,2),TransformSubspaceXvects(NTest,3),'g','LineStyle',':','LineWidth',4);
                    quiver3(0,0,0,projBackSubspaceXvects(NTest,1),projBackSubspaceXvects(NTest,2),projBackSubspaceXvects(NTest,3),'k','LineStyle',':','LineWidth',4);
                    NTest=130; % give some examples
                    quiver3(0,0,0,SubspaceXvects(NTest,1),SubspaceXvects(NTest,2),SubspaceXvects(NTest,3),'r','LineStyle','--','LineWidth',4);
                    quiver3(0,0,0,TransformSubspaceXvects(NTest,1),TransformSubspaceXvects(NTest,2),TransformSubspaceXvects(NTest,3),'g','LineStyle','--','LineWidth',4);
                    quiver3(0,0,0,projBackSubspaceXvects(NTest,1),projBackSubspaceXvects(NTest,2),projBackSubspaceXvects(NTest,3),'k','LineStyle','--','LineWidth',4);
                    
                    legend([RotAxisLeg NormAxisLeg OrgSubSp TransFor ProjBck],{'RotAxis','NormPlane','inSub','Trans','Proj'},'Location','best')
                    
                end
                
                sgtitle({sprintf('Inferred Angle PCA:%0.2f, True:%i',PCAresults(r).PairAngled,Demo1.Rot(r));...
                    sprintf('Inferred Rot AxB:%0.2f/%0.2f/%0.2f',AxBresults(r).Rot(1),AxBresults(r).Rot(2),AxBresults(r).Rot(3));...
                    sprintf('Inferred Rot Tot AxB:%0.2f',AxBresults(r).RotTot);
                    sprintf('Inferred Rot Axis AxB:%0.2f/%0.2f/%0.2f',AxBresults(r).RotAxis(1),AxBresults(r).RotAxis(2),AxBresults(r).RotAxis(3));
                    sprintf(' Comp:%0.2f/%0.2f/%0.2f',AxBresults(r).CompTot(1),AxBresults(r).CompTot(2),AxBresults(r).CompTot(3));...
                    sprintf(' True Comp:%0.2f/%0.2f',TrueComp(1),TrueComp(2));...
                    sprintf('Angle test vect:AZ%0.2f/EL%0.2f',azimuth_PRotAxis,elevation_PRotAxis);...
                    sprintf('Angle trans test vect:AZ%0.2f/EL%0.2f',azimuth_TransformPRotAxis,elevation_TransformPRotAxis);...
                    sprintf('Rot test vect:AZ%0.2f/EL%0.2f',wrapTo180(azimuth_TransformPRotAxis-azimuth_PRotAxis),wrapTo180(elevation_TransformPRotAxis-elevation_PRotAxis));...
                    sprintf('Rot ang NP:%0.2f',AxBresults(r).RotAng_NP)});
                axis equal
                subplot(143)
                %polarhistogram(deg2rad(AzimuthRandDiff),0:2*pi/99:2*pi);title('Azimuth');
                polarhistogram(deg2rad((RotInplane_SubspaceVects)),0:2*pi/99:2*pi);title('InPlane');
                subplot(144)
                % polarhistogram(deg2rad(ElevationRandDiff),0:2*pi/99:2*pi);title('Elevation');
                polarhistogram(deg2rad((RotOutplane_SubspaceVects)),0:2*pi/99:2*pi);title('OutofPlane');
                pause
                mvFrame(r) = getframe(gcf);
            end
            %  [h,Sp,mvFrame]=obj.ShowAllMDS(RSAresults,h,Sp);
            [~,~,MDSFigFileName]=obj.ManData.GetFileName(['Subspace'],['_Demo_AxB_RotX' ],'SaveInResults',1,'WantedDate','ALL');
            obj.FigParams.MakeMovieFromFrames(mvFrame,1,MDSFigFileName)
            
        end
        function ParametrizedAxBanaDemo(obj,varargin) % Demo to fit simulated data with AxB and test differnt hypothetis
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            %% Simulate a model that rotates a representation with respect to another representation
            Demo1.Rot=0:10:180;
            PointOrder=obj.PointOrder;%[1 2 4 3 1];
            for r=1:length(Demo1.Rot)
                % construct our ground truth RDM and then fit it
                %                                                    CR1  CR2 CIR1 CIR2  x y z            Bx  By  Bz
                [~,~,RespVect]=obj.RSAana.ParamModelRDM3DRot_General([0.2 0.2 0.8  0.8   30 Demo1.Rot(r) 70 0   0   0],...
                    'ParamMdlRDM_OrthogonalContexts',-1);
                % add a tiny bit of noise to this so we avoid low rank matrixes
                %RespVect=RespVect+0.0000005*rand(size(RespVect));
                TrueComp=obj.CalFakeDataCompression(RespVect,PointOrder); % calculate ture compression
                
                % now decode what is in the space
                [AxBresults(r)] =obj.RunAxBanalysisNeuralData([],RespVect,[],[]);
                
                % use PCA to double check the results as well
                PCAresults(r)=obj.RSAana.DoPCAanalysisonMdlData(RespVect);
                RotAxis=AxBresults(r).RotAxis;
                %% now fit a plane to the X and project the vectors onto that plane
                [PCA,PCAcoeff,OriginPlane,FittedPlane,~,OrthoVector,OrgPlaneEq]=obj.RSAana.FitPlane2Data(RespVect(5:8,:));
                
                Ns=360;
                Base1=PCAcoeff(1,:);Base2=PCAcoeff(2,:);
                SubspaceXvects=cell2mat(arrayfun(@(x) [cosd(x)*Base1+(sind(x))*Base2;cosd(x)*Base1-(sind(x))*Base2 ]',1:Ns/2,'uniformoutput',0))';
                DotSubVectsN=arrayfun(@(x) dot(SubspaceXvects(x,:),OrthoVector),1:Ns);
                if sum(DotSubVectsN>0.01);error('Plane equation is not correct');end
                % now transform the SubspaceXvects with rotation matrix and
                % project them back into this plane
               % TransformSubspaceXvects=cell2mat(arrayfun(@(x) [SubspaceXvects(x,:)*AxBresults(r).RotMat']',1:Ns,'UniformOutput',0))';
                TransformSubspaceXvects=cell2mat(arrayfun(@(x) [AxBresults(r).RotMat*SubspaceXvects(x,:)'],1:Ns,'UniformOutput',0))';
               
                % project the back into the plane
                [~,ProjMat,projBackSubspaceXvects]=arrayfun(@(x) obj.ManData.ProjectInto3DSubspace(TransformSubspaceXvects(x,:)',PCAcoeff(1,:)',PCAcoeff(2,:)'),1:Ns,'UniformOutput',0);
                projBackSubspaceXvects=cell2mat(projBackSubspaceXvects)';
                Dot_projBackSubspaceXvects=arrayfun(@(x) dot(projBackSubspaceXvects(x,:),OrthoVector),1:Ns);
                if sum(Dot_projBackSubspaceXvects>0.01);error('Error in projecting back into plane');end
                % now calculate in-plane and out of plane rotations
                RotInplane_SubspaceVects=arrayfun(@(x) obj.ManData.GetAngleBetVectors(SubspaceXvects(x,:),projBackSubspaceXvects(x,:)),1:Ns);
                
                % out of plane rotation would be the angle between the
                % projected vector and normal vector of the plane
                [RotOutplane_SubspaceVects,~,DotTransOrtho]=arrayfun(@(x) obj.ManData.GetAngleBetVectors(OrthoVector,TransformSubspaceXvects(x,:)),1:Ns);
                RotOutplane_SubspaceVects(RotOutplane_SubspaceVects<=90)=90-RotOutplane_SubspaceVects(RotOutplane_SubspaceVects<=90);
                RotOutplane_SubspaceVects(RotOutplane_SubspaceVects>90)=RotOutplane_SubspaceVects(RotOutplane_SubspaceVects>90)-90;
                
                % calculate the angle between plane's normal vector and
                % rotation axis: that is the out of plane rotation
                RotOutplane=obj.ManData.GetAngleBetVectors(OrthoVector,RotAxis);
                
                %   %% do the same process for rotation axis
                %   [~,ProjMat_RotAxis,projBackRotAxis]=obj.ManData.ProjectInto3DSubspace(RotAxis,PCAcoeff(1,:)',PCAcoeff(2,:)');
                %   RotInplane_Rot=arrayfun(@(x) obj.ManData.GetAngleBetVectors(SubspaceXvects(x,:),projBackSubspaceXvects(x,:)),1:Ns);
                
                subplot(1,4,[1 2])
                cla
                hold on
                plot3(RespVect(PointOrder,1),RespVect(PointOrder,2),RespVect(PointOrder,3))
                plot3(RespVect(PointOrder+4,1),RespVect(PointOrder+4,2),RespVect(PointOrder+4,3))
                
                plot3(RespVect([1],1),RespVect([1],2),RespVect([1],3),'r*')
                plot3(RespVect([5],1),RespVect([5],2),RespVect([5],3),'r*')
                % plot Rotation vector from origin as well
                xlabel('X Axis'); ylabel('Y Axis'); zlabel('Z axis')
                view(30,30)
                
                if 1
                    % plot a sample vector and its transformation
                    RotAxisLeg=quiver3(0,0,0,AxBresults(r).RotAxis(1),AxBresults(r).RotAxis(2),AxBresults(r).RotAxis(3),'g');% rotation axis
                    NormAxisLeg=quiver3(0,0,0,OrthoVector(1),OrthoVector(2),OrthoVector(3),'r');% normal vector to plane
                    OrgSubSp=plot3(SubspaceXvects(:,1),SubspaceXvects(:,2),SubspaceXvects(:,3),'r.');
                    TransFor= plot3(TransformSubspaceXvects(:,1),TransformSubspaceXvects(:,2),TransformSubspaceXvects(:,3),'g.');
                    ProjBck=plot3(projBackSubspaceXvects(:,1),projBackSubspaceXvects(:,2),projBackSubspaceXvects(:,3),'k.');
                    NTest=30; % give some examples
                    quiver3(0,0,0,SubspaceXvects(NTest,1),SubspaceXvects(NTest,2),SubspaceXvects(NTest,3),'r','LineStyle',':','LineWidth',4);
                    quiver3(0,0,0,TransformSubspaceXvects(NTest,1),TransformSubspaceXvects(NTest,2),TransformSubspaceXvects(NTest,3),'g','LineStyle',':','LineWidth',4);
                    quiver3(0,0,0,projBackSubspaceXvects(NTest,1),projBackSubspaceXvects(NTest,2),projBackSubspaceXvects(NTest,3),'k','LineStyle',':','LineWidth',4);
                    NTest=130; % give some examples
                    quiver3(0,0,0,SubspaceXvects(NTest,1),SubspaceXvects(NTest,2),SubspaceXvects(NTest,3),'r','LineStyle','--','LineWidth',4);
                    quiver3(0,0,0,TransformSubspaceXvects(NTest,1),TransformSubspaceXvects(NTest,2),TransformSubspaceXvects(NTest,3),'g','LineStyle','--','LineWidth',4);
                    quiver3(0,0,0,projBackSubspaceXvects(NTest,1),projBackSubspaceXvects(NTest,2),projBackSubspaceXvects(NTest,3),'k','LineStyle','--','LineWidth',4);
                    legend([RotAxisLeg NormAxisLeg OrgSubSp TransFor ProjBck],{'RotAxis','NormPlane','inSub','Trans','Proj'},'Location','best')
                end
                
                sgtitle({sprintf('Inferred Angle PCA:%0.2f, True:%i',PCAresults(r).PairAngled,Demo1.Rot(r));...
                    sprintf('Inferred Rot AxB:%0.2f/%0.2f/%0.2f',AxBresults(r).Rot(1),AxBresults(r).Rot(2),AxBresults(r).Rot(3));...
                    sprintf('Inferred Rot Tot AxB:%0.2f',AxBresults(r).RotTot);
                    sprintf('Inferred Rot Axis AxB:%0.2f/%0.2f/%0.2f',AxBresults(r).RotAxis(1),AxBresults(r).RotAxis(2),AxBresults(r).RotAxis(3));
                    sprintf(' Comp:%0.2f/%0.2f/%0.2f',AxBresults(r).CompTot(1),AxBresults(r).CompTot(2),AxBresults(r).CompTot(3));...
                    sprintf(' True Comp:%0.2f/%0.2f',TrueComp(1),TrueComp(2));...
                    sprintf('Rot ang NP:%0.2f',AxBresults(r).RotAng_NP)});
                
                axis equal
                subplot(143)
                %polarhistogram(deg2rad(AzimuthRandDiff),0:2*pi/99:2*pi);title('Azimuth');
                polarhistogram(deg2rad((RotInplane_SubspaceVects)),0:2*pi/180:2*pi);title('InPlane');
                subplot(144)
                % polarhistogram(deg2rad(ElevationRandDiff),0:2*pi/99:2*pi);title('Elevation');
                polarhistogram(deg2rad((RotOutplane_SubspaceVects)),0:2*pi/180:2*pi);title('OutofPlane');
                pause
                mvFrame(r) = getframe(gcf);
            end
            %  [h,Sp,mvFrame]=obj.ShowAllMDS(RSAresults,h,Sp);
            [~,~,MDSFigFileName]=obj.ManData.GetFileName(['Subspace'],['_Demo_AxB_RotX' ],'SaveInResults',1,'WantedDate','ALL');
            obj.FigParams.MakeMovieFromFrames(mvFrame,1,MDSFigFileName)
            
        end
        
        function AxBresults=NeuralDataAxBanaDemo(obj,FactorLvLInds,FactorLvLData,FactorLevels,Opts,TimePoint,AnalysisTitle,varargin) % Demo to fit Neural data with AxB and test differnt hypothetis
           %@TimePoint is the time point we are looking at right now
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            PointOrder=[1 2 3 4 1]; 
            NDims='Both';
            % organizied by [RB,RT,GB,GT]
            MarkerSizeSet=[15 5 15 5 15 5 15 5];MarkerColorSet={'r','r','g','g','r','r','g','g'};
            NNeu=size(FactorLvLData{1},2);
                       
            % calculate in-plane and out-of-plane angles
            if NDims==3
                % now decode what is in the space
                [AxBresults,X,Y] =obj.RunAxBanalysisNeuralData(FactorLvLInds,FactorLvLData,FactorLevels,Opts,NDims);
                RespVect=[Y';X'];
                [AxBresults,PCAresults,SubspaceXvects,TransformSubspaceXvects,projBackSubspaceXvects,FullTransformSubspaceXvects,OrthoVector]=...
                    obj.CalAnglesin3D(RespVect,AxBresults,X,Y);
                 RespVect=[Y';X'];
            elseif NDims==4
                % now decode what is in the space
                [AxBresults,X,Y] =obj.RunAxBanalysisNeuralData(FactorLvLInds,FactorLvLData,FactorLevels,Opts,NDims);
                RespVect=[Y';X'];
                [AxBresults]=obj.CalAnglesin4D(RespVect,AxBresults,X,Y); 
                 
            else % do both
                % now decode what is in the space
                [AxBresults,X,Y] =obj.RunAxBanalysisNeuralData(FactorLvLInds,FactorLvLData,FactorLevels,Opts,4);
                RespVect=[Y';X'];
                [AxBresults4D]=obj.CalAnglesin4D(RespVect,AxBresults,X,Y);
                % now decode what is in the space
                [AxBresults,X,Y] =obj.RunAxBanalysisNeuralData(FactorLvLInds,FactorLvLData,FactorLevels,Opts,3);
                RespVect=[Y';X'];
                [AxBresults,PCAresults,SubspaceXvects,TransformSubspaceXvects,projBackSubspaceXvects,FullTransformSubspaceXvects,OrthoVector,X,Y]=...
                    obj.CalAnglesin3D(RespVect,AxBresults,X,Y);
                RespVect=[Y';X']; % redefine X RespVect
                
                % run everything with full dimensions now to check it later on 
                [AxBresultsFull] =obj.RunAxBanalysisNeuralData(FactorLvLInds,FactorLvLData,FactorLevels,Opts,NNeu,'UsePCA4AxB',1);
                 
                
                AxBresults=obj.ManData.CopyVars2Struct(AxBresults,'RotInplane_SubspaceVects4D',AxBresults4D.RotInplane_SubspaceVects,...
                    'RotOutplane_SubspaceVects4D',AxBresults4D.RotOutplane_SubspaceVects,...
                    'Compression_X4D',AxBresults4D.Compression_X,'Compression_Y4D',...
                    AxBresults4D.Compression_Y,'Compression_XCompressed4D',AxBresults4D.Compression_XCompressed,...
                    'RotInplane_X4D',AxBresults4D.RotInplane_X,'RotOutplane_X4D',AxBresults4D.RotOutplane_X,...
                    'GeneralizationErr4D',AxBresults4D.GeneralizationErr,'RotTotEig4D',AxBresults4D.RotTotEig,...
                    'RotTotFull',AxBresultsFull.RotTot,'RotAng_NPFull', AxBresultsFull.RotAng_NP,...
                    'Compression_XFull', AxBresultsFull.Compression_X,...
                    'Compression_XCompressedFull',AxBresultsFull.Compression_XCompressed,'RotTotEigFull',AxBresultsFull.RotTotEig);
            end
            TrueComp=obj.CalFakeDataCompression(RespVect,PointOrder); % calculate ture compression
            if obj.PlotAxBDetails
                obj.FigParams.RenderFigure(1,[]);
                subplot(2,4,[1 2 5 6])
                cla
                obj.PlotAxBXYplanes(RespVect,[],AnalysisTitle); % plot X and Y planes 

                if obj.PlotAxBCircSubspace
                    % multiply compression to transformationed to see if
                    % they match with what it should be
                    CompressTransomed=[AxBresults.CompMatprime* TransformSubspaceXvects']';
                    % multiply only the compression to the data and see
                    % where it taks 
                    OnlyCompress=[AxBresults.CompMat*SubspaceXvects']';
                    % plot a sample vector and its transformation
                    RotAxisLeg=quiver3(0,0,0,AxBresults.RotAxis(1),AxBresults.RotAxis(2),AxBresults.RotAxis(3),'g');% rotation axis
                    NormAxisLeg=quiver3(0,0,0,OrthoVector(1),OrthoVector(2),OrthoVector(3),'r');% normal vector to plane
                    OrgSubSp=plot3(SubspaceXvects(:,1),SubspaceXvects(:,2),SubspaceXvects(:,3),'r.');
                    TransFor= plot3(TransformSubspaceXvects(:,1),TransformSubspaceXvects(:,2),TransformSubspaceXvects(:,3),'g.');
                    ProjBck=plot3(projBackSubspaceXvects(:,1),projBackSubspaceXvects(:,2),projBackSubspaceXvects(:,3),'k.');
                    FullTransLeg=plot3(FullTransformSubspaceXvects(:,1),FullTransformSubspaceXvects(:,2),FullTransformSubspaceXvects(:,3),'c.');
                    CompTransLeg=plot3(CompressTransomed(:,1),CompressTransomed(:,2),CompressTransomed(:,3),'y.');
                    OnlyCompressLeg=plot3(OnlyCompress(:,1),OnlyCompress(:,2),OnlyCompress(:,3),'m.');

                    
                    NTest=30; % give some examples
                    quiver3(0,0,0,SubspaceXvects(NTest,1),SubspaceXvects(NTest,2),SubspaceXvects(NTest,3),'r','LineStyle',':','LineWidth',4);
                    quiver3(0,0,0,TransformSubspaceXvects(NTest,1),TransformSubspaceXvects(NTest,2),TransformSubspaceXvects(NTest,3),'g','LineStyle',':','LineWidth',4);
                    quiver3(0,0,0,projBackSubspaceXvects(NTest,1),projBackSubspaceXvects(NTest,2),projBackSubspaceXvects(NTest,3),'k','LineStyle',':','LineWidth',4);
                    NTest=130; % give some examples
                    quiver3(0,0,0,SubspaceXvects(NTest,1),SubspaceXvects(NTest,2),SubspaceXvects(NTest,3),'r','LineStyle','--','LineWidth',4);
                    quiver3(0,0,0,TransformSubspaceXvects(NTest,1),TransformSubspaceXvects(NTest,2),TransformSubspaceXvects(NTest,3),'g','LineStyle','--','LineWidth',4);
                    quiver3(0,0,0,projBackSubspaceXvects(NTest,1),projBackSubspaceXvects(NTest,2),projBackSubspaceXvects(NTest,3),'k','LineStyle','--','LineWidth',4);
                    legend([XsubspacePlot YsubspacePlot RotAxisLeg NormAxisLeg OrgSubSp TransFor ProjBck FullTransLeg CompTransLeg OnlyCompressLeg],...
                        {'X','Y','RotAxis','NormPlane','inSub','Trans','Proj','Full','CompTran','Comp'},'Location','best')
                end
                
                sgtitle({AnalysisTitle;sprintf('Inferred Angle PCA:%0.2f',PCAresults.PairAngled);...
                    sprintf('Inferred Rot AxB:%0.2f/%0.2f/%0.2f',AxBresults.Rot(1),AxBresults.Rot(2),AxBresults.Rot(3));...
                    sprintf('Inferred Rot Tot AxB:%0.2f',AxBresults.RotTot);
                    sprintf('Inferred Rot Axis AxB:%0.2f/%0.2f/%0.2f',AxBresults.RotAxis(1),AxBresults.RotAxis(2),AxBresults.RotAxis(3));
                    sprintf(' Comp:%0.2f/%0.2f/%0.2f',AxBresults.CompTot(1),AxBresults.CompTot(2),AxBresults.CompTot(3));...
                    sprintf(' True Comp:%0.2f/%0.2f',TrueComp(1),TrueComp(2));...
                    sprintf('Rot ang NP:%0.2f RconErr:%0.4f',AxBresults.RotAng_NP,AxBresults.ReconstErr);...
                    sprintf('Time: %0.2f %s',TimePoint,AnalysisOpts.SpkCntStartFieldName)});
                
                axis equal
                
                % plot distribution of in-plane and out-of-plane angles  in
                % 3D
                subplot(243);cla;
                polarhistogram(deg2rad((AxBresults.RotInplane_SubspaceVects)),0:2*pi/180:2*pi);
                v=axis;hold on
                arrayfun(@(x) polarplot(deg2rad(AxBresults.RotInplane_X(x)),v(4),'d','markersize',MarkerSizeSet(x),'MarkerFaceColor',MarkerColorSet{x}),1:4);title('InPlane');
                
                subplot(244);cla;
                polarhistogram(deg2rad((AxBresults.RotOutplane_SubspaceVects)),0:2*pi/180:2*pi);v=axis;hold on
                arrayfun(@(x) polarplot(deg2rad(AxBresults.RotOutplane_X(x)),v(4),'d','markersize',MarkerSizeSet(x),'MarkerFaceColor',MarkerColorSet{x}),1:4);title('OutofPlane');
               
                % plot distribution of in-plane and out-of-plane angles  in
                % 4D
                subplot(2,4,7);cla
                polarhistogram(deg2rad((AxBresults4D.RotInplane_SubspaceVects)),0:2*pi/180:2*pi);
                v=axis;hold on
                arrayfun(@(x) polarplot(deg2rad(AxBresults4D.RotInplane_X(x)),v(4),'d','markersize',MarkerSizeSet(x),'MarkerFaceColor',MarkerColorSet{x}),1:4);title('InPlane4D');
                
                subplot(2,4,8);cla;
                polarhistogram(deg2rad((AxBresults4D.RotOutplane_SubspaceVects)),0:2*pi/180:2*pi);v=axis;hold on
                arrayfun(@(x) polarplot(deg2rad(AxBresults4D.RotOutplane_X(x)),v(4),'d','markersize',MarkerSizeSet(x),'MarkerFaceColor',MarkerColorSet{x}),1:4);title('OutofPlane4D');
            end
            % only keep important variables
            if ~AnalysisOpts.CalShuffleAxB
                %                 AxBresults=obj.ManData.rmfieldExept(AxBresults,{'RotMat','ReconstErr','CompMat','A','RotInplane_X','RotOutplane_X','RotInplane_X4D',...
                %                     'Compression_X4D','Compression_Y4D','Compression_XCompressed4D','Compression_X','Compression_Y','Compression_XCompressed',...
                %                     'RotOutplane_X4D','RotTot','RotAngPCA'});
                if AnalysisOpts.GetFullAxBdata
                    AxBresults=obj.ManData.rmfieldExept(AxBresults,{'RotTotEig','RotTotEig4D','RotTotEigFull','RotX','Xproj2plane','Yproj2plane','RotInplane_X','RotOutplane_X','RotInplane_X4D',...
                        'Compression_X4D','Compression_Y4D','Compression_XCompressed4D','Compression_X','Compression_Y','Compression_XCompressed',...
                        'RotOutplane_X4D','GeneralizationErr4D','GeneralizationErr','RotTot','RotAngPCA','RotTotFull','RotAng_NPFull','Compression_XFull','Compression_XCompressedFull'});
                    
                else
                    AxBresults=obj.ManData.rmfieldExept(AxBresults,{'RotTotEig','RotTotEig4D','RotTotEigFull','RotInplane_X','RotOutplane_X','RotInplane_X4D',...
                        'Compression_X4D','Compression_Y4D','Compression_XCompressed4D','Compression_X','Compression_Y','Compression_XCompressed',...
                        'RotOutplane_X4D','GeneralizationErr4D','GeneralizationErr','RotTot','RotAngPCA','RotTotFull','RotAng_NPFull','Compression_XFull','Compression_XCompressedFull'});
                end
            else
                AxBresults=obj.ManData.rmfieldExept(AxBresults,{'RotTotEig','RotTotEig4D','RotTotEigFull','RotInplane_X','RotOutplane_X','RotInplane_X4D','RotOutplane_X4D','RotTot',...
                    'Compression_X4D','Compression_Y4D','Compression_XCompressed4D','RotAngPCA','Compression_X','Compression_Y','Compression_XCompressed'...
                    'GeneralizationErr4D','GeneralizationErr','RotTotFull','RotAng_NPFull','Compression_XFull','Compression_XCompressedFull'});
            end
        end
        function [AxBresults,PCAresults,SubspaceXvects,TransformSubspaceXvects,projBackSubspaceXvects,FullTransformSubspaceXvects,OrthoVector,X,Y]=CalAnglesin3D(obj,RespVect,AxBresults,X,Y,varargin) % calculates all angles in 3 dimensions
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
           
            % use PCA to double check the results as well
            PCAresults=obj.RSAana.DoPCAanalysisonMdlData(RespVect);
            RotAxis=AxBresults.RotAxis;
           
            %% now fit a plane to the X and project the vectors onto that plane
            [PCA,PCAcoeff,OriginPlane,FittedPlane,~,OrthoVector,OrgPlaneEq]=obj.RSAana.FitPlane2Data(X');
            Base1=PCAcoeff(1,:);Base2=PCAcoeff(2,:);
            
            %% project X and Y into planes that fit them (PAY ATTENTION TO THIS )
            [~,~,Xproj2plane]=arrayfun(@(x) obj.ManData.ProjectIntoSubspace(X(:,x),Base1',Base2'),1:4,'UniformOutput',0);
            Xproj2plane=cell2mat(Xproj2plane);
            % project Y into its best fitting plane 
            [~,PCAcoeff_Y]=obj.RSAana.FitPlane2Data(Y');
            Base1Y=PCAcoeff_Y(1,:);Base2Y=PCAcoeff_Y(2,:);
            [~,~,Yproj2plane]=arrayfun(@(x) obj.ManData.ProjectIntoSubspace(Y(:,x),Base1Y',Base2Y'),1:4,'UniformOutput',0);
            Yproj2plane=cell2mat(Yproj2plane);
            if obj.ProjectXY2Plane % Are we using the data from projections 
                X=Xproj2plane;
                Y=Yproj2plane;
            end                        
            AxBresults.Yproj2plane=Yproj2plane;
            AxBresults.Xproj2plane=Xproj2plane;
            
            % generate a circle using the basis functions on that subspace
            Ns=360;
            Base1=PCAcoeff(1,:);Base2=PCAcoeff(2,:);
            SubspaceXvects=cell2mat(arrayfun(@(x) [cosd(x)*Base1+(sind(x))*Base2;cosd(x)*Base1-(sind(x))*Base2 ]',1:Ns/2,'uniformoutput',0))';
            DotSubVectsN=arrayfun(@(x) dot(SubspaceXvects(x,:),OrthoVector),1:Ns);
            if sum(DotSubVectsN>0.01);error('Plane equation is not correct');end
           
            % now transform the SubspaceXvects with rotation matrix and
            % project them back into this plane
            [RotInplane_SubspaceVects,RotOutplane_SubspaceVects,TransformSubspaceXvects,projBackSubspaceXvects,FullTransformSubspaceXvects]=...
                obj.CalInPlaneOutPlaneSubspaceRotation(AxBresults.A,AxBresults.RotMat,SubspaceXvects,Base1',Base2',OrthoVector); % calculates in-plane and out-of-plane rotation with respect to a subspace
            
            % get the in-plane and out-of-plane rotation for acutual data
            % points of X
           [RotInplane_X,RotOutplane_X,TransformX,projBackX,FullTransformX]=...
                obj.CalInPlaneOutPlaneSubspaceRotation(AxBresults.A,AxBresults.RotMat,X',Base1',Base2',OrthoVector); % calculates in-plane and out-of-plane rotation with respect to a subspace
           
            % copy results into results matrix     
            AxBresults=obj.ManData.CopyVars2Struct(AxBresults,'RotInplane_SubspaceVects',RotInplane_SubspaceVects,'RotOutplane_SubspaceVects',RotOutplane_SubspaceVects,...
                'RotInplane_X',RotInplane_X,'RotOutplane_X',RotOutplane_X,'RotAngPCA',PCAresults.PairAngled);

            % calculate the angle between plane's normal vector and
            % rotation axis: that is the out of plane rotation
            RotOutplane=obj.ManData.GetAngleBetVectors(OrthoVector,RotAxis);
            
            %   %% do the same process for rotation axis
            %   [~,ProjMat_RotAxis,projBackRotAxis]=obj.ManData.ProjectInto3DSubspace(RotAxis,PCAcoeff(1,:)',PCAcoeff(2,:)');
            %   RotInplane_Rot=arrayfun(@(x) obj.ManData.GetAngleBetVectors(SubspaceXvects(x,:),projBackSubspaceXvects(x,:)),1:Ns);
        end       
        function [AxBresults,PCAresults,SubspaceXvects,TransformSubspaceXvects,projBackSubspaceXvects,FullTransformSubspaceXvects]=CalAnglesin4D(obj,RespVect,AxBresults,X,Y,varargin) % calculates all angles in 4 dimensions
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
           
            % use PCA to double check the results as well
            PCAresults=obj.RSAana.DoPCAanalysisonMdlData(RespVect);
            RotAxis=AxBresults.RotAxis;
           
            %% now fit a plane to the X and project the vectors onto that plane
            [PCA,PCAcoeff,~,~,~,OrthoPlane,~]=obj.RSAana.FitPlane2Data(X');
            Base1=PCAcoeff(1,:);Base2=PCAcoeff(2,:);
            %% project X and Y into planes that fit them (PAY ATTENTION TO THIS )
            [~,~,Xproj2plane]=arrayfun(@(x) obj.ManData.ProjectIntoSubspace(X(:,x),Base1',Base2'),1:4,'UniformOutput',0);
            Xproj2plane=cell2mat(Xproj2plane);
            % project Y into its best fitting plane 
            [~,PCAcoeff_Y]=obj.RSAana.FitPlane2Data(Y');
            Base1Y=PCAcoeff_Y(1,:);Base2Y=PCAcoeff_Y(2,:);
            [~,~,Yproj2plane]=arrayfun(@(x) obj.ManData.ProjectIntoSubspace(Y(:,x),Base1Y',Base2Y'),1:4,'UniformOutput',0);
            Yproj2plane=cell2mat(Yproj2plane);
            if obj.ProjectXY2Plane % Are we using the data from projections 
                X=Xproj2plane;
                Y=Yproj2plane;
            end
                
            %%
            % generate a circle using the basis functions on original subspace
            Ns=360;
            SubspaceXvects=cell2mat(arrayfun(@(x) [cosd(x)*Base1+(sind(x))*Base2;cosd(x)*Base1-(sind(x))*Base2 ]',1:Ns/2,'uniformoutput',0))';
            DotSubVectsN=SubspaceXvects*OrthoPlane';
            if sum(DotSubVectsN(:)>0.01);error('Plane equation is not correct');end
           
            % now transform the SubspaceXvects with rotation matrix and
            % project them back into this plane
            [RotInplane_SubspaceVects,RotOutplane_SubspaceVects,TransformSubspaceXvects,projBackSubspaceXvects,FullTransformSubspaceXvects]=...
                obj.CalInPlaneOutPlaneSubspaceRotation(AxBresults.A,AxBresults.RotMat,SubspaceXvects,Base1',Base2',OrthoPlane); % calculates in-plane and out-of-plane rotation with respect to a subspace
            
            % get the in-plane and out-of-plane rotation for acutual data
            % points of X
           [RotInplane_X,RotOutplane_X,TransformX,projBackX,FullTransformX]=...
                obj.CalInPlaneOutPlaneSubspaceRotation(AxBresults.A,AxBresults.RotMat,X',Base1',Base2',OrthoPlane); % calculates in-plane and out-of-plane rotation with respect to a subspace
           
            % copy results into results matrix     
            AxBresults=obj.ManData.CopyVars2Struct(AxBresults,'RotInplane_SubspaceVects',RotInplane_SubspaceVects,'RotOutplane_SubspaceVects',RotOutplane_SubspaceVects,...
                'RotInplane_X',RotInplane_X,'RotOutplane_X',RotOutplane_X,'RotAngPCA',PCAresults.PairAngled);

            % calculate the angle between plane's normal vector and
            % rotation axis: that is the out of plane rotation
          %  RotOutplane=obj.ManData.GetAngleBetVectors(OrthoVector,RotAxis);
            
            %   %% do the same process for rotation axis
            %   [~,ProjMat_RotAxis,projBackRotAxis]=obj.ManData.ProjectInto3DSubspace(RotAxis,PCAcoeff(1,:)',PCAcoeff(2,:)');
            %   RotInplane_Rot=arrayfun(@(x) obj.ManData.GetAngleBetVectors(SubspaceXvects(x,:),projBackSubspaceXvects(x,:)),1:Ns);
        end
        function Compression=CalNeuralDataCompression(obj,X,TransOrder,varargin)            
            % @X PC principal components of the quadrants
            % calculates compression if the input is PCs for Quadrants the
            % distance is euclidean distance between mean([Dist(GB-GT)
            % Dist(RB-RT)]/mean([GB-RB]-[RT-GT)
            % if we are organizing quadrants then 
            %@ Transorder is the order of tranformations in the case of neural data it is            
            %[1 1][2 2][1 2][2 1]->[RB,GT,GB,RT] to
            % the neural data is organized as
            %[1 1][2 1][1 2][2 2]->[RB,RT,GB,GT]
            
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            Dist= @(x,y) sqrt(sum((x(:) - y(:)).^2)); %Euclidean Distance
            A1=X(TransOrder(1),:);
            A2=X(TransOrder(2),:);
            A3=X(TransOrder(3),:);
            A4=X(TransOrder(4),:);
            
            Compression=mean([Dist(A2,A3) Dist(A1,A4)])/mean([Dist(A3,A1) Dist(A4,A2)]);
        end
        
        function Comp=CalFakeDataCompression(obj,RespVect,PointOrder) % calculate compression for fake data
            
            % calculte Euclidean distance from first point to second
            Dist11=obj.ManData.EuclideanDistance(RespVect(PointOrder(1),:),RespVect(PointOrder(2),:));
            Dist12=obj.ManData.EuclideanDistance(RespVect(PointOrder(2),:),RespVect(PointOrder(3),:));
            
            Dist21=obj.ManData.EuclideanDistance(RespVect(PointOrder(1)+4,:),RespVect(PointOrder(2)+4,:));
            Dist22=obj.ManData.EuclideanDistance(RespVect(PointOrder(2)+4,:),RespVect(PointOrder(3)+4,:));
            
            Comp(1)=Dist11/Dist21;
            Comp(2)=Dist12/Dist22;
        end        
        function [AxBresults,X,Y,A] = RunAxBanalysisNeuralData(obj,FactorLvLInds,FactorLvLData,FactorLevels,Opts,NDims,varargin) % runs AxB analysis on Neural data
            %FactorLvLInds index of factor levels
            %FactorLvLData Data from factors
            %Opts SubspaceAnalysisOpts
            %NDims Number of dimensions to anlayze neural data. Fake data's
            %dimension is always 3
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
          
            [A,X,Y,TransOrder]=obj.SolveAxBEq(FactorLvLInds,FactorLvLData,Opts,NDims);
            
            %% compute rotational and compression of A
            [R,P,Pprime]=obj.PolarDecomposeA(A);

            % compute rotation angle and axis using eigen value and eigen
            % vectors
            [EigVect,EigVal]=eig(R);EigVal=diag(EigVal);
            AxBresults.Rot=rad2deg(angle(EigVal));
            %AxBresults.RotEig=sum(abs(rad2deg(angle(EigVal))))/2; % get the rotation angle with the eigen values
            % determine the rotation axis
            % AxBresults.RotAxis=EigVect(:,imag(EigVal)==0);
            AxBresults.RotTot=acosd((trace(R)-1)/2); 
            AxBresults.Rot(floor(AxBresults.Rot)==180)=0;
            AxBresults.RotTotEig=wrapTo180(sum(abs(rad2deg(angle(EigVal))))/2);%                     
            AxBresults.DetRotMat=det(R);
            AxBresults.RRT=R*R';
            
            % comput rotation axis and rotation angle that is not proncipal
            % angle
            [AxBresults.RotAxis,AxBresults.RotAng_NP]=obj.ComputRotAxisandAngle(R);
            AxBresults.RotLabRef=acosd(R);
            % compression
            %AxBresults.Comp=det(V);
            AxBresults.ReconstErr= (norm(A*X-Y));
            if AxBresults.ReconstErr>0.01; error('High reconstrcution error');end
                        
            % [U1,S1,V1] = svd(U); %U1 and V1 should be equal
            AxBresults.CompTot=diag(P);
            AxBresults.CompMat=P;
            AxBresults.CompMatprime=Pprime;
            AxBresults.RotMat=R;
            AxBresults.A=A;
            AxBresults.RotX=R*X; % get the rotation of X
            % tranform X into Y using only compression P then calculate compression value 
            XCompressed=[AxBresults.CompMat*X];
            AxBresults.Compression_X=obj.CalNeuralDataCompression(X',TransOrder);
            AxBresults.Compression_Y=obj.CalNeuralDataCompression(Y',TransOrder);
            AxBresults.Compression_XCompressed=obj.CalNeuralDataCompression(XCompressed',TransOrder);
            
            %% test generalization ability of A
            AxBresults.GeneralizationErr=obj.TestAxBgeneralization(FactorLvLInds,FactorLvLData,Opts,NDims);
                                  
            %% get a test vector and see what are the Theta and Phi angles ( not using this currently)
%             RotAxis=AxBresults.RotAxis;
%             RespVectExmp=X(:,1);
%             PRotAxis=[RespVectExmp(1) RespVectExmp(2) -(RespVectExmp(1)*RotAxis(1)+RespVectExmp(2)*RotAxis(2))/RotAxis(3)];
%             TransformPRotAxis=[PRotAxis*AxBresults.RotMat'];
%             %                 acosd(dot(PRotAxis,TransformPRotAxis)/(norm(PRotAxis)*norm(TransformPRotAxis)))
%             %                 acosd(TransformPRotAxis/norm(TransformPRotAxis))-acosd(PRotAxis/norm(PRotAxis))
%             %
%             [azimuth_PRotAxis,elevation_PRotAxis,r_PRotAxis]=obj.ManData.ConvertCart2Spherical(PRotAxis(1),PRotAxis(2),PRotAxis(3));
%             [azimuth_TransformPRotAxis,elevation_TransformPRotAxis,r_TransformPRotAxis]=obj.ManData.ConvertCart2Spherical(TransformPRotAxis(1),TransformPRotAxis(2),TransformPRotAxis(3));
%             AxBresults.RotAximuth=azimuth_TransformPRotAxis-azimuth_PRotAxis;
%             AxBresults.RotElevation=elevation_TransformPRotAxis-elevation_PRotAxis;
            
        end
        function [R,P,Pprime]=PolarDecomposeA(obj,A) % applies polar decomposition to square matrix a

            % [R P V] = poldecomp(A); % we are not using this because it
            % can't handle singular matrices 
            [U,S,V] = svd(A);
            R=U*V';
            P=V*S*V';      %A=R*P*X
            Pprime=U*S*U'; %A=Pprime*R*X
        end
        function [GeneralizationErr,A,Xtest,Ytest,GeneralizationCorr]=TestAxBgeneralization(obj,FactorLvLInds,FactorLvLData,Opts,NDims,varargin) % splits the data into half and runs AxB analysis then tests on other half
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
           [FactorLvLData1,FactorLvLData2]= obj.ManData.SplitDatainTwo(FactorLvLData);
          
           % run AxB analysis on one data set
           [A]=obj.SolveAxBEq(FactorLvLInds,FactorLvLData1,Opts,NDims);
           
           % now predict the test data 
           [~,Xtest,Ytest]=obj.SolveAxBEq(FactorLvLInds,FactorLvLData2,Opts,NDims,1);   
            
           Yhat=A*Xtest;
           % take average Euclidean distance 
           GeneralizationErr=mean(arrayfun(@(x) obj.ManData.EuclideanDistance(Ytest(:,x),Yhat(:,x)),1:size(Yhat,2)));
           GeneralizationCorr=mean(arrayfun(@(x) corr(Ytest(:,x),Yhat(:,x)),1:size(Yhat,2)));
            
        end
        function [GeneralizationErr,CompositionErr,GeneralizationErr_Sh,GeneralizationCorr,CompositionCorr,GeneralizationCorr_Sh]=TestAxBCompositionality(obj,FactorLvLInds,FactorLvLData,Opts,varargin)% tests the compositionality in dynamics 
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs

            % get factor level data for all times for rule 1 2 3
            NDims=size(FactorLvLData{1}{1},2); % use maximum dimensions
            nShuffle=10;
            nRun=10;
            AnalysisOpts.CalShuffleAxB=0;
            ComAnalysisTitles={'Ycomp3_R1C2','Ycomp3_TR1C2','Ycomp3_R2C1',...
                'Ycomp3_R1','Ycomp3_R2','Ycomp3_R3','Ycomp3_TR1',...
                'Ycomp3_C1','Ycomp3_C2','Ycomp3_C3',...
                'Ycomp3_R1C1','Ycomp3_R2C2',...
                'Ycomp3_S2A1','Ycomp3_S2A2','Ycomp3_S2A3','Ycomp3_S1A1'};

            for Run=1:nRun
                for Rule=1:3 % loop on data for each rule
                    [GeneralizationErr(Run,Rule),A{Rule},Xtest{Rule},Ytest{Rule},GeneralizationCorr(Run,Rule)]=obj.TestAxBgeneralization(FactorLvLInds,FactorLvLData{Rule},Opts,NDims);
                    % [A{Rule},X{Rule},Y{Rule}]=obj.SolveAxBEq(FactorLvLInds,FactorLvLData,Opts,NDims);
                    [R{Rule},C{Rule},Pprime{Rule}]=obj.PolarDecomposeA(A{Rule});%% compute rotational and compression of A
                end

                % now construct rule 3 using rotation in Rule1 and Compression P in Rule 2
                Ycomp3_R1C2=R{1}*C{2}*Xtest{3}; % build rule 3 from rotation of rule 1 and compression of rule 2
                Ycomp3_TR1C2=transpose(R{1})*C{2}*Xtest{3}; % build rule 3 from transpose of rotation of rule 1 and compression of rule 2
                Ycomp3_R2C1=R{2}*C{1}*Xtest{3}; % build rule 3 from rotation of rule 2 and compression of rule 1
                
                Ycomp3_R1=R{1}*Xtest{3}; % build rule 3 from rotation of rule 1             
                Ycomp3_R2=R{2}*Xtest{3}; % build rule 3 from rotation of rule 2 
                Ycomp3_R3=R{3}*Xtest{3}; % build rule 3 from rotation of rule 3
                
                Ycomp3_TR1=transpose(R{1})*Xtest{3}; % build rule 3 from transpose of rotation of rule 1             
                
                Ycomp3_C1=C{1}*Xtest{3}; % build rule 3 from compression of rule 1 
                Ycomp3_C2=C{2}*Xtest{3}; % build rule 3 from compression of rule 2
                Ycomp3_C3=C{3}*Xtest{3}; % build rule 3 from compression of rule 3

                Ycomp3_R1C1=R{1}*C{1}*Xtest{3}; % build rule 3 from rotation of rule 1 and compression of rule 1
                Ycomp3_R2C2=R{2}*C{2}*Xtest{3}; % build rule 3 from rotation of rule 2 and compression of rule 2

                Ycomp3_S2A1=A{1}*Xtest{2};% build rule 3 from sensory of rule 2 and transformation of rule 1
                Ycomp3_S2A2=A{2}*Xtest{2};% build rule 3 from sensory of rule 2 and transformation of rule 2
                Ycomp3_S2A3=A{3}*Xtest{2};% build rule 3 from sensory of rule 2 and transformation of rule 3
                Ycomp3_S1A1=A{1}*Xtest{1};% build rule 3 from sensory of rule 1 and transformation of rule 1
                
               
                for i=1:length(ComAnalysisTitles)
                    eval(sprintf('CompositionErr(Run,i)=mean(arrayfun(@(x) obj.ManData.EuclideanDistance(Ytest{3}(:,x),%s(:,x)),1:size(%s,2)));',ComAnalysisTitles{i},ComAnalysisTitles{i}));
                    eval(sprintf('CompositionCorr(Run,i)=mean(arrayfun(@(x) corr(Ytest{3}(:,x),%s(:,x)),1:size(%s,2)));',ComAnalysisTitles{i},ComAnalysisTitles{i}));                    
                end
            end

            % compute shuffle error for each rule as well
             AnalysisOpts.CalShuffleAxB=1;
             for Rule=1:3 % loop on data for each rule
                [GeneralizationErr_Sh(:,Rule),~,~,~,GeneralizationCorr_Sh(:,Rule)]=arrayfun(@(x) obj.TestAxBgeneralization(FactorLvLInds,FactorLvLData{Rule},Opts,NDims),1:nShuffle,'uniformoutput',0);                          
             end
             AnalysisOpts.CalShuffleAxB=0;
                
             GeneralizationErr=mean(GeneralizationErr);
             CompositionErr=mean(CompositionErr);
             GeneralizationErr_Sh=mean(cell2mat(GeneralizationErr_Sh));
           
             CompositionCorr=mean(CompositionCorr);
             GeneralizationCorr_Sh=mean(cell2mat(GeneralizationCorr_Sh));
             GeneralizationCorr=mean(GeneralizationCorr);
        end
        function [A,X,Y,TransOrder]=SolveAxBEq(obj,FactorLvLInds,FactorLvLData,Opts,NDims,OnlyXYflag,varargin) % solve euqation of Ax=B
            % the organization of A and B is N*B where N is the number of
            % neurons and B is the number of conditions
            %FactorLvLInds index of factor levels
            %FactorLvLData Data from factors
            %Opts SubspaceAnalysisOpts
            %NDims Number of dimensions to anlayze neural data. Fake data's
            %dimension is always 3
            % OnlyXYflag just give the XY 
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
                       
            if iscell(FactorLvLData) & ~obj.UsePCA4AxB & obj.SubtractFactorMean % then take PCA
                [FactorLvLDataQuad,TransOrder]=obj.OrganizeQuadrants(FactorLvLInds,FactorLvLData,Opts);
                % take mean of the data and concatinate them in a matrix
                MeanFactorLvLData=cell2mat(cellfun(@mean ,FactorLvLDataQuad,'UniformOutput',0)');
                % X=MeanFactorLvLData(1:4,:)';
                % Y=MeanFactorLvLData(5:8,:)';
                ConXY=MeanFactorLvLData;
                [~,Score_PCA_ConXY,~,~,Explained]=pca(ConXY,'Centered',true);
                % mean center these data first
                Score_PCA_ConXY(1:4,:)=Score_PCA_ConXY(1:4,:)-mean(Score_PCA_ConXY(1:4,:),1);
                Score_PCA_ConXY(5:8,:)=Score_PCA_ConXY(5:8,:)-mean(Score_PCA_ConXY(5:8,:),1);                
                X=Score_PCA_ConXY(5:8,1:NDims)'; % is  what is being transformed
                Y=Score_PCA_ConXY(1:4,1:NDims)'; % is the target
            elseif iscell(FactorLvLData) & obj.UsePCA4AxB & obj.SubtractFactorMean% if we have taken the PCA alrady
                [FactorLvLDataQuad,TransOrder]=obj.OrganizeQuadrants(FactorLvLInds,FactorLvLData,Opts);
                MeanFactorLvLData=cell2mat(arrayfun(@(x) mean(FactorLvLDataQuad{x},1)' ,1:length(FactorLvLDataQuad),'UniformOutput',0))';

                Score_PCA_ConXY(1:4,:)=MeanFactorLvLData(1:4,:)-mean(MeanFactorLvLData(1:4,:),1);
                Score_PCA_ConXY(5:8,:)=MeanFactorLvLData(5:8,:)-mean(MeanFactorLvLData(5:8,:),1);
                X=Score_PCA_ConXY(5:8,1:NDims)';  % is what is being transformed
                Y=Score_PCA_ConXY(1:4,1:NDims)';  % is the target
            elseif ~obj.SubtractFactorMean  % don't subtract factor mean 
                [FactorLvLDataQuad,TransOrder]=obj.OrganizeQuadrants(FactorLvLInds,FactorLvLData,Opts);
                MeanFactorLvLData=cell2mat(arrayfun(@(x) mean(FactorLvLDataQuad{x},1)' ,1:length(FactorLvLDataQuad),'UniformOutput',0))';

                Score_PCA_ConXY(1:4,:)=MeanFactorLvLData(1:4,:);
                Score_PCA_ConXY(5:8,:)=MeanFactorLvLData(5:8,:);
                X=Score_PCA_ConXY(5:8,1:NDims)';  % is what is being transformed
                Y=Score_PCA_ConXY(1:4,1:NDims)';  % is the target
            else  % if we don't need to transform the order before
                TransOrder=[1 2 3 4];
                NStimTot=size(FactorLvLData,1)/2;
                X=FactorLvLData(NStimTot+1:end,1:3)';  % is what is being transformed
                Y=FactorLvLData(1:NStimTot,1:3)';  % is the target
            end
            if exist('OnlyXYflag','var')
                if OnlyXYflag==1
                    A=[];
                   return
                end
            end
            % now if we are doing shuffle permute the Y matrix order 
            if AnalysisOpts.CalShuffleAxB==1
                Y=Y(:,randperm(size(Y,2)));
            end
            %   AX=Y->X'A'=Y' % solve the reverse problem first 
            Ap=mldivide(X',Y');            
            A=Ap';
            if norm([A*X-Y])>0.1;warning('mldivide is not finding proper solution');end
        end
        function [RotInplane_SubspaceVects,RotOutplane_SubspaceVects,TransformSubspaceXvects,projBackSubspaceXvects,FullTransformSubspaceXvects]=CalInPlaneOutPlaneSubspaceRotation(obj,A,RotMat,SubspaceXvects,Base1,Base2,OrthoVector) % calculates in-plane and out-of-plane rotation with respect to a subspace
           %@A is the full transformation matrix
           % @RotMat is the rotation matrix
           % @SubspaceXvects vectors in the original subpace that need to
           % be projected into new subspace 
           % @ Base1 and Base2 are the basis vectors for the original
           % subspace
           %@OrthoVector orthogonal vector to the subspace plane
           global AnalysisOpts
           
            Ns=size(SubspaceXvects,1);
            NDim=length(A);
             %   TransformSubspaceXvects=cell2mat(arrayfun(@(x) [SubspaceXvects(x,:)*AxBresults(r).RotMat']',1:Ns,'UniformOutput',0))';
            TransformSubspaceXvects=cell2mat(arrayfun(@(x) [RotMat*SubspaceXvects(x,:)'],1:Ns,'UniformOutput',0))';
            
            % transform these vectors with transformation matrix as well
            FullTransformSubspaceXvects=cell2mat(arrayfun(@(x) [A*SubspaceXvects(x,:)'],1:Ns,'UniformOutput',0))';
            
            % project the back into the plane
            [~,ProjMat,projBackSubspaceXvects]=arrayfun(@(x) obj.ManData.ProjectIntoSubspace(TransformSubspaceXvects(x,:)',Base1,Base2),1:Ns,'UniformOutput',0);
            projBackSubspaceXvects=cell2mat(projBackSubspaceXvects)';
            Dot_projBackSubspaceXvects=projBackSubspaceXvects*OrthoVector';
            if sum(Dot_projBackSubspaceXvects(:)>0.01);error('Error in projecting back into plane');end
            
            
                % if we are in 3 dimensions then the in-plane is the angle of
                % the projection before and after on the original plane and
                % out-of-plane is the angle between transformed and rothognal vector of the plane
                % now calculate in-plane and out of plane rotations
                RotInplane_SubspaceVects=arrayfun(@(x) obj.ManData.GetAngleBetVectors(SubspaceXvects(x,:),projBackSubspaceXvects(x,:)),1:Ns);
               
                if NDim==3
                    % out of plane rotation would be the angle between the
                    % projected vector and normal vector of the plane
                    [RotOutplane_SubspaceVects]=arrayfun(@(x) obj.ManData.GetAngleBetVectors(OrthoVector,TransformSubspaceXvects(x,:)),1:Ns);
                    
                    % we can limit out of plane angles to 0 and 90
                    % RotOutplane_SubspaceVects(RotOutplane_SubspaceVects<=90)=90-RotOutplane_SubspaceVects(RotOutplane_SubspaceVects<=90);
                    % RotOutplane_SubspaceVects(RotOutplane_SubspaceVects>90)=RotOutplane_SubspaceVects(RotOutplane_SubspaceVects>90)-90;
                elseif NDim==4
                    % if we are in 4 dimension the in-plane is similar to 3D
                    % but for out-of-plane the angle the transformed matrix
                    % should be projected into nullspace(Complemtary orhtogonal
                    % subspace) of orginal plane and the angle should be
                    % calculated there
                    BaseNull1=OrthoVector(1,:)';BaseNull2=OrthoVector(2,:)';
                    [~,~,proj2NullSubspaceXvects]=arrayfun(@(x) obj.ManData.ProjectIntoSubspace(TransformSubspaceXvects(x,:)',BaseNull1,BaseNull2),1:Ns,'UniformOutput',0);
                    proj2NullSubspaceXvects=cell2mat(proj2NullSubspaceXvects)';
                    Dot_proj2NullSubspaceXvects=proj2NullSubspaceXvects*[Base1 Base2];
                    if sum(Dot_proj2NullSubspaceXvects(:)>0.01);error('Error in projecting back into null plane');end
                    
                    % out of plane rotation would be the angle between the
                    % projected vectors in the null space and rotated
                    % vector
                    % to add direction to this plane first measure the
                    % angle of the vector with the first basis vector if it
                    % is more than 90 degree then multiply -1 by the
                    % transformationvectors, this is because the angle
                    % between each vector and its projection on a plane is
                    % always less or equal to 90
                     [NullBasis_SubspaceVects]=arrayfun(@(x) obj.ManData.GetAngleBetVectors(TransformSubspaceXvects(x,:),BaseNull1),1:Ns);                                        
                     NullBasis_SubspaceVects(NullBasis_SubspaceVects<=90)=1;NullBasis_SubspaceVects(NullBasis_SubspaceVects>90)=-1;
                     
                    [RotOutplane_SubspaceVects]=arrayfun(@(x) obj.ManData.GetAngleBetVectors(TransformSubspaceXvects(x,:),NullBasis_SubspaceVects(x)*proj2NullSubspaceXvects(x,:)),1:Ns);                                        
                end
        end
        function [FactorLvLData,TransOrder]=OrganizeQuadrants(obj,FactorLvLInds,FactorLvLData,Opts,varargin) % organzied Quadrants data to match the RDM
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % if our main target factors are quadrants then we need to fix
            % the order of objects to match the what we want
            % the Quadrants are organized as [Shape Color]
            %[1 1][2 2][1 2][2 1]->[RB,GT,GB,RT]
            % the model RDM is organized as
            %[1 1][2 1][1 2][2 2]->[RB,RT,GB,GT]
            % so the transformation we need is [1,4,3,2] for each rule
            TransOrder=[1,4,3,2];
            if ~strcmp(Opts.MainTargetFactor,'Quadrants');return;end
            nRule=1:length(unique(FactorLvLInds));
            TransInds=cell2mat(arrayfun(@(x) TransOrder+(x-1)*4,nRule,'UniformOutput',0));
            if iscell(FactorLvLData)
                FactorLvLData=FactorLvLData(TransInds);
            else
                FactorLvLData=FactorLvLData(TransInds,:);
            end
        end
        function [RotAxis,RotAng]=ComputRotAxisandAngle(obj,R,varargin) % algorithm to compute axis and angle of rotation without requiring principal angle
            % taken from chapter 5.4.7 %Brannon, Rebecca M., and H. L. Schreyer. Rotation, Reflection,
            %and Frame Changes: Orthogonal Tensors in Computational Engineering Mechanics.
            %Version: 20180401. IOP Expanding Physics. Bristol, UK: IOP Publishing, 2018.
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % Step 0 compute cosine of the angle of rotation
            % c=0.5*(trace(R)-1)
            c=0.5*(trace(R)-1);
            %Step1 Construct matrix of components for R+R'-2*c*I
            CompMat=R+R'-2*c*eye(size(R));
            %Step 2:let x be the column of Compmat with largest magetude
            %and compute the size of x
            MagX=arrayfun(@(x) norm(CompMat(:,x)),1:size(CompMat,2));
            [MagX,IndX]=max(MagX);
            X=CompMat(:,IndX);
            %Step 3: Compute axis of rotation
            RotAxis=X/MagX;if sum(X)==0;RotAxis=zeros(size(X));end
            %Step 4 compute sine of the angle of rotaion
            s=0.5*(RotAxis(1)*(R(3,2)-R(2,3))+RotAxis(2)*(R(1,3)-R(3,1))+RotAxis(3)*(R(2,1)-R(1,2)));
            
            Angles360 = @(a) rem(360+a, 360); % turn every angle to the range of 0 to 360
            RotAng = Angles360(atan2d(s,c));
        end
        %% plot functions 
        function PlotAxBXYplanes(obj,RespVect,Xproj,AnalysisTitle) % plots planes for X Y and look at it through time 
                % Xproj is the projecttion of X with rotation matrix

                hold on
                YsubspacePlot=plot3(RespVect(obj.PointOrder,1),RespVect(obj.PointOrder,2),RespVect(obj.PointOrder,3),'LineWidth',3);
                XsubspacePlot=plot3(RespVect(obj.PointOrder+4,1),RespVect(obj.PointOrder+4,2),RespVect(obj.PointOrder+4,3),'LineWidth',3);

                % organizied by [RB,RT,GB,GT]
                MarkerSizeSet=[15 5 15 5 15 5 15 5];MarkerColorSet={'r','r','g','g','r','r','g','g'};              
                arrayfun(@(x) plot3(RespVect([x],1),RespVect([x],2),RespVect([x],3),'d','MarkerSize',MarkerSizeSet(x),...
                    'MarkerFaceColor',MarkerColorSet{x}),1:8);
                RX=AnalysisTitle(6);RY=AnalysisTitle(10);
                if ~isempty(Xproj)
                    XprojsubspacePlot=plot3(Xproj(obj.PointOrder,1),Xproj(obj.PointOrder,2),Xproj(obj.PointOrder,3),'g');

                    arrayfun(@(x) plot3(Xproj([x],1),Xproj([x],2),Xproj([x],3),'d','MarkerSize',MarkerSizeSet(x),...
                        'MarkerFaceColor',MarkerColorSet{x}),1:4);
                    legend([XsubspacePlot YsubspacePlot XprojsubspacePlot],{['X:R' RX],['Y:R' RY],'proj'},'Location','best')
                else 
                    legend([XsubspacePlot YsubspacePlot],{['X:R' RX],['Y:R' RY]},'Location','best')
                end
                
                 
                xlabel('X Axis'); ylabel('Y Axis'); zlabel('Z axis')
                view(30,30)
        end
       
        
    end
end

