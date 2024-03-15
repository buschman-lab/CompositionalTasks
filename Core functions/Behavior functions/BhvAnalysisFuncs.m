classdef BhvAnalysisFuncs
    %BHVANALYSISFUNCS put analysis functions for behavior here
    properties 
       NavgPerf=15;   % number of avg for performance
    end
    properties (Access=private)
        fp=fig_params;
        ManData=ManipulateData;
    end
    
    methods
        
        function  obj=BhvAnalysisFuncs(varargin)
            global AnalysisOpts
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
        function [AllPSM,AllTrlPerf,IndSamp,AllTrlCount,AllTrlCountDay,NBlocksDay,RewardPulse,NCorrectTrl,NumRewards,AllSeqHist] = CancatinateInfoDays(obj,Perf) % take average of psychmetric curve and plots it
            warning('off')
            Perf=Perf(arrayfun(@(x) ~isempty(x{1}),Perf));
            NDays=length(Perf);
            AllPSM.Perf=cell(1,3);AllPSM.Count=cell(1,3);InitMat={[],[],[]};
            %  AllTrlPerf=struct('FromSwitch',InitMat,'ToSwitch',InitMat,'FromSwitch_Cong',InitMat,'ToSwitch_Cong',InitMat,'FromSwitch_InCong',InitMat,'ToSwitch_InCong',InitMat);
            AllTrlPerf=struct('FromSwitch',[],'ToSwitch',[],'LogisticFit',[]);
            AllTrlPerf.FromSwitch=struct('All',InitMat,'Congruent',InitMat,'InCongruent',InitMat,'IndTrls',InitMat,'AllTrlCnt',InitMat);
            AllTrlPerf.ToSwitch=struct('All',InitMat,'Congruent',InitMat,'InCongruent',InitMat,'IndTrls',InitMat,'AllTrlCnt',InitMat);
            AllTrlPerf.LogisticFit=struct('All',InitMat,'Congruent',InitMat,'InCongruent',InitMat,'IndTrls',InitMat,'AllTrlCnt',InitMat);
            AllSeqHist=cell(1,3);
            AllTrlCount=cell(1,3);
            IndSamp.Perf=cell(1,3);IndSamp.Count=cell(1,3);IndSamp.RT=cell(1,3);IndSamp.CorrCount=cell(1,3); % individual sample
            AllPSM.Level=Perf{1, 1}.PsycPerfLvl;
            AllTrlCountDay=zeros(1,NDays);
            NBlocksDay=zeros(NDays,3);
            % loop on days and rules
            for i=1:NDays
                for Rule=1:length(Perf{i}.PsycPerf)
                    if ~isempty(Perf{i}.OverAllPerf{Rule})
                        %% concatinate PSM
                        AllPSM.Perf{Rule}=[AllPSM.Perf{Rule};Perf{i}.PsycPerf{Rule}];
                        AllPSM.Count{Rule}=[AllPSM.Count{Rule};Perf{i}.PsycCount{Rule}];
                        %% calculate performance for each block with desired nMovAvg. Tile them forward and backward
                        %TrlPerf=cell2mat(cellfun(@(x) [padarray(x,[0 800-length(x)],NaN,'post')]',Perf{i}.OverAllPerf{Rule},'uniformoutput',0));
                        %AllTrlPerf{Rule}=[AllTrlPerf{Rule} TrlPerf];
                        [TrlPerf]=cellfun(@(x) obj.PerfCalCulateOverAllPerformance(x,obj.NavgPerf),Perf{i}.AllBlkTrls{Rule},'uniformoutput',0);
                        for Field={'All','Congruent','InCongruent','IndTrls','AllTrlCnt'}
                            Field=Field{1};
                            TrlPerf_FromSwitch=cell2mat(cellfun(@(x) [padarray(x.(Field),[0 800-length(x.(Field))],NaN,'post')]',TrlPerf,'uniformoutput',0));
                            TrlPerf_ToSwitch=cell2mat(cellfun(@(x) [padarray(x.(Field),[0 800-length(x.(Field))],NaN,'pre')]',TrlPerf,'uniformoutput',0));
                            AllTrlPerf.FromSwitch(Rule).(Field)= [AllTrlPerf.FromSwitch(Rule).(Field) TrlPerf_FromSwitch];
                            AllTrlPerf.ToSwitch(Rule).(Field)= [AllTrlPerf.ToSwitch(Rule).(Field) TrlPerf_ToSwitch];                           
                            %% fit logistic to the performance from Switch x=(0:length(I)-1);
                            if ~(strcmp(Field,'IndTrl') | strcmp(Field,'AllTrlCnt'))
                                AllTrlPerf.LogisticFit(Rule).(Field)=[AllTrlPerf.LogisticFit(Rule).(Field) cellfun(@(x) FitLogistic((0:length(x.(Field))-1)',x.(Field)'),TrlPerf,'UniformOutput',0)];
                            end
                        end
                        %% concatinate Performance
                        AllTrlCount{Rule}=[AllTrlCount{Rule} sum(~isnan(TrlPerf_FromSwitch))+ obj.NavgPerf]; % 15 is the averaging number make sure it is consistent
                        AllTrlCountDay(i)=Perf{i}.TotAttemptedTrls;%AllTrlCountDay(i)+sum(sum(~isnan(TrlPerf))+15);                     
                        %% determine SeqHist for this day of recording
                         RuleSeq=Perf{i}.ConditionTot;
                         SeqHist=nan(1,size(RuleSeq,2));
                         SeqHist(1:2:end)=[false logical(diff(RuleSeq(1:2:end)))];
                         SeqHist(2:2:end)=[false logical(diff(RuleSeq(2:2:end)))];
                         % change rule sequence for rule 2 based on the
                         % rule before it 0 means rule1 and 1 means rule 3
                         Rule2Inds=find(RuleSeq==2);
                         SeqHist(Rule2Inds==1)=-1;Rule2IndsNotBeg=Rule2Inds(Rule2Inds>1);
                         SeqHist(Rule2IndsNotBeg)=RuleSeq(Rule2IndsNotBeg-1);
                         SeqHist(SeqHist==1 & RuleSeq==2)=0;
                         SeqHist(SeqHist==3 & RuleSeq==2)=1; 
                         AllSeqHist{Rule}=[AllSeqHist{Rule} SeqHist(RuleSeq==Rule)];
                        %% average IndSample Info
                        for j=1:length(Perf{i}.IndSampInfo{Rule})
                            IndSamp.Perf{Rule} =cat(3,IndSamp.Perf{Rule},Perf{i}.IndSampInfo{Rule}{j}.ThisObjPerf);
                            IndSamp.Count{Rule}=cat(3,IndSamp.Count{Rule},Perf{i}.IndSampInfo{Rule}{j}.AllTrialsCount);
                            IndSamp.CorrCount{Rule}=cat(3,IndSamp.CorrCount{Rule},Perf{i}.IndSampInfo{Rule}{j}.CorrTrialsCount);                            
                            IndSamp.RT{Rule}   =cat(3,IndSamp.RT{Rule},Perf{i}.IndSampInfo{Rule}{j}.AvgRT);
                        end
                    end
                end
                temp=cellfun(@(x) size(x,1),Perf{i}.PsycPerf);
                NBlocksDay(i,:)=[temp zeros(1,3-length(temp))];
                RewardPulse(i)=Perf{i}.RewardPulse;
                NCorrectTrl(i)=Perf{i}.TotCorrectTrls;
                NumRewards.Switch(i)=Perf{i}.NumRewards.Switch;
            end
            IndSamp.Pdist=Perf{1}.IndSampInfo{2}{1}.Pdist; % grab one Pdist from here
            IndSamp.Congmap=Perf{1}.IndSampInfo{2}{1}.CongMap;
            %% perform statistical test on each of the morph levels 
            IndSamp.TotalTrials=cellfun(@(x) sum(x,3),IndSamp.Count,'UniformOutput',0); % total number of trials
            IndSamp.TotalCorrTrials=cellfun(@(x) sum(x,3),IndSamp.CorrCount,'UniformOutput',0); % total number of correct trials 
            % perform statistical test on these now 
            [IndSamp.p_BinoSigTest,IndSamp.p_Binostar]=arrayfun(@(x) obj.BinomialStatTeston3Dmat(IndSamp.TotalCorrTrials{x},IndSamp.TotalTrials{x},0.5,''),1:3,'UniformOutput',0);
            %% perform binomial test trial performance during learning 
            AllTrlPerf=obj.PerformBinomialTestAllTrlPerf(AllTrlPerf);
            AllTrlPerf.NavgPerf=obj.NavgPerf; % save this number as well
        end
        function [Perf,CongPerCongruent,CongPerInCongruent]=PerfCalCulateOverAllPerformance(obj,Trials,NMovAvgPerf)% calculates average performance for a block
            sm_kern = ones(1, NMovAvgPerf);
            sm_kern_norm = sm_kern./sum(sm_kern);
            INDOBJ_AllTrials=(([Trials.StopCondition] == 1 | [ Trials.StopCondition] == -1));
            INDOBJ_CorrectTrials=([Trials.StopCondition] == 1);
            CorrectMat=INDOBJ_CorrectTrials(INDOBJ_AllTrials);            
            Perf.All=convn(CorrectMat,sm_kern_norm,'valid');
            Perf.IndTrls=double(CorrectMat);
            Perf.AllTrlCnt=convn(CorrectMat,sm_kern,'valid');
            %% calculate performance for congruent and incongruent trials
            FieldName={'Congruent','InCongruent'};
            CongPerf=struct(FieldName{1},[],FieldName{2},[]);
            ObjectML=[Trials.ObjectMorphLevel];
            ColorML=[Trials.ColorMorphLevel];
            StimCongruency=obj.ManData.DetermineStimCongruency(ObjectML,ColorML);
            k=1;
            AttemptedTrials=[Trials.StopCondition] ==-1 | [Trials.StopCondition] ==1 ;           
            for Cong=[1 0]
                INDOBJ_AllTrials=double(([Trials.StopCondition] == 1 | [Trials.StopCondition] == -1) & StimCongruency == Cong);
                INDOBJ_AllTrials(INDOBJ_AllTrials==0)=NaN;
                INDOBJ_CorrectTrials=([Trials.StopCondition] == 1);
                CorrectMat=INDOBJ_AllTrials.*INDOBJ_CorrectTrials;
                CorrectMat=CorrectMat(AttemptedTrials);
                Perf.(FieldName{k})=movmean(CorrectMat,NMovAvgPerf,2,'omitnan','Endpoints','discard');
                Perf.(FieldName{k})=obj.ManData.RemoveEntryFromVec(Perf.(FieldName{k}),NaN); % remove NaN entries
                k=k+1;
            end
            CongPerCongruent=Perf.Congruent;
            CongPerInCongruent=Perf.InCongruent;
        end
        function Perf=CalculateBhvPerf(obj,CorrectMat,NMovAvgPerf) % calculates behavioral performance with logistic fit givven only the correct mat 
            %@CorrectMat matrix with correct trials 
            %@NMovAvgPerf number of avegrage 
            sm_kern = ones(1, NMovAvgPerf);
            sm_kern = sm_kern./sum(sm_kern);                     
            Perf.All=convn(CorrectMat,sm_kern,'valid');
            x=(0:length(Perf.All)-1);
            cla
            Perf.FitResult=FitLogistic(x',Perf.All');
            %plot(x,Perf.All);hold on
            %plot(x,Perf.FitResult.L./(1+exp(-Perf.FitResult.k*(x-Perf.FitResult.x0))),'--g')
            %title(sprintf('k:%0.2f,L:%0.2f,x0:%0.2f',Perf.FitResult.k,Perf.FitResult.L,Perf.FitResult.x0))
            %pause
        end
        function PerfTot=CalculateFeatureBhvPerf(obj,CorrectMat,NMovAvgPerf) % calculates behavioral performance with logistic fit givven only the correct mat 
            %@CorrectMat matrix with correct trials 
            %@NMovAvgPerf number of avegrage 
            
            % remove nan values first but then put them back
            PerfTot=zeros(1,length(CorrectMat));
            PerfTot(isnan(CorrectMat))=nan;
            CorrectMatNan=CorrectMat(~isnan(CorrectMat));
            sm_kern = ones(1, NMovAvgPerf);
            sm_kern = sm_kern./sum(sm_kern);                     
            Perf=convn(CorrectMatNan,sm_kern,'same');
            PerfTot(~isnan(CorrectMat))=Perf;
        end
        function [PerfShapeTot,PerfColorTot]=CalculateBhvFeaturePerformance(obj,ShapeCat,ColorCat,ResponseLoc,Rule,NMovAvgPerf)% calculates performance with respect to both shape and color features
            % calculate the performance accoring to shape and color
            % contingency matrix is category response correct or incorrect
            % remove off axis responses in calculation of performance here 
            if Rule==1 | Rule ==3;offAxis=(ResponseLoc==3 | ResponseLoc==4);else
                offAxis=(ResponseLoc==1 | ResponseLoc==2);end
            NtrlsTot=length(ResponseLoc);
            onAxis=~offAxis;
            ShapeCat=ShapeCat(onAxis);
            ColorCat=ColorCat(onAxis);
            ResponseLoc=ResponseLoc(onAxis);
                        
            Ntrls=length(ResponseLoc);
            if Rule==2;
                ResponseLoc(ResponseLoc==3)=1;ResponseLoc(ResponseLoc==4)=2;end
            
            ShapeContingency=[1 1 1;...
                1 2 0;...
                2 1 0;...
                2 2 1;...
               -1 1 nan;...
               -1 2 nan];
            ColorContingency=[1 1 0;...
                1 2 1;...
                2 1 1;...
                2 2 0;...
                -1 1 nan;...
                -1 2 nan];
            % go through this and find what is correct and what is incorrect 
            CorrectMat_Shape=arrayfun(@(x) ShapeContingency(ShapeContingency(:,1)==ShapeCat(x) & ShapeContingency(:,2)==ResponseLoc(x),3),1:Ntrls);
            CorrectMat_Color=arrayfun(@(x) ColorContingency(ColorContingency(:,1)==ColorCat(x) & ColorContingency(:,2)==ResponseLoc(x),3),1:Ntrls);
            
            PerfShape=obj.CalculateFeatureBhvPerf(CorrectMat_Shape,NMovAvgPerf); 
            PerfColor=obj.CalculateFeatureBhvPerf(CorrectMat_Color,NMovAvgPerf);    
            PerfShapeTot=nan*ones(NtrlsTot,1);PerfColorTot=PerfShapeTot;
            PerfShapeTot(onAxis)=PerfShape;PerfColorTot(onAxis)=PerfColor;
        end
        function [Trial]=CalIndividualSampleInfo(obj,bhv,Trials,Condition) % calculates essential info for each sample
            
            Pdist=bhv.Pdist;
            if iscell(Pdist)
                Pdist=Pdist{1};
            end
            Pdist.Dim1Objs=[0,20,30,40,50,60,70,80,100,120,130,140,150,160,170,180];
            Dim1Ident= [1 1 1 1 -1 2 2 2 2 2 2 2 -1 1 1 1];
            Dim2Ident= [1 1 1 1 -1 2 2 2 2 2 2 2 -1 1 1 1];
            % build a map of congruency 1: congruent, 2:incongruent,
            % 3:50/150
            CongMap.Trans=[1 1 2;
                           2 2 2;
                           1 2 1;
                           2 1 1;
                          -1 1 3;
                          -1 2 3;
                          1 -1 3;
                          2 -1 3;
                         -1 -1 3];
            for d1=1:length(Pdist.Dim1Objs)
                 for d2=1:length(Pdist.Dim1Objs)
                    ID=find(CongMap.Trans(:,1)==Dim1Ident(d1) & CongMap.Trans(:,2)==Dim1Ident(d2));
                    CongMap.Map(d1,d2)=CongMap.Trans(ID,3);
                end
            end
            
            
            L=1:length(Pdist.Dim1Objs);
            d1=1;% Dimension 1 is shape and 2 is color
            ReactionTime=[Trials.ReactionTime];
            for Dim1ML=Pdist.Dim1Objs
                d2=1;
                for Dim2ML=Pdist.Dim1Objs
                    %%% Calclute performance of the each of the objects and put it
                    %%% in ametrix with all other charachtricstics of that object so
                    %%% that for each block we have a summery matrix
                    ThisObj= [Trials.ObjectMorphLevel]==Dim1ML & [Trials.ColorMorphLevel]==Dim2ML & [Trials.Condition]==Condition;
                    AllTrials=([Trials.StopCondition] == 1 | [Trials.StopCondition] ==-1) & ThisObj;
                    AllTrialsPassed=([Trials.StopCondition] ~= 1 & [Trials.StopCondition] ~=-1) & ThisObj; %trial that animals did fix break or too fast etc
                    CorrTrials=([Trials.StopCondition] == 1) & ThisObj;
                    %%
                    Trial.AllTrialsCount(d1,d2)=sum(AllTrials);
                    Trial.CorrTrialsCount(d1,d2)=sum(CorrTrials);
                    Trial.AllTrialsPassedCount(d1,d2)=sum(AllTrialsPassed);
                    Trial.ThisObjPerf(d1,d2)=sum(CorrTrials)/sum(AllTrials);
                    Trial.RT{d1,d2}=ReactionTime(CorrTrials);
                    Trial.AvgRT(d1,d2)=nanmean(ReactionTime(CorrTrials));
                    d2=d2+1;
                end
                d1=d1+1;
            end
            Trial.Pdist=Pdist;
            Trial.CongMap=CongMap;
        end       
        %% plot functions
        function varargout = PlotAvgPSM(obj,AllPSMPerf) % plots average PSM and fits it
            global AnalysisOpts
            % h is axis
            obj.fp.font_size=AnalysisOpts.Plotting.font_size;
            varargout=obj.fp.RenderFigure(4,[]);
            Col=AnalysisOpts.RuleColors;
            RuleFeatures={'Shape','Color','Color'};
            ErrBarCol='k';
            ErrBarLineWidth=2;
            ErrBarMarkerSize=10;
            PSMMethodTxt={'STD_method','sem'};
            % index levels
            MorphLvl=AllPSMPerf.Level;
            HighPerfInd=1:9;LowPerfInd=10:18;
            HighMorphLvl=MorphLvl(HighPerfInd);
            LowMorphLvl=MorphLvl(LowPerfInd);
            MLs0to1=[0 20 30 40 50 60 70 80 100]/100;
            UsedMLs=~isnan(mean(AllPSMPerf.Perf{1}(:,HighPerfInd),1));
            ML=MLs0to1(UsedMLs);
            HighPerfInd=HighPerfInd(UsedMLs);LowPerfInd=LowPerfInd(UsedMLs);
            LowMorphLvl=LowMorphLvl(UsedMLs);HighMorphLvl=HighMorphLvl(UsedMLs);
            YAxisTicks=[0 0.25 0.5 0.75 1];YAxisTickLabels=arrayfun(@(x) num2str(x),YAxisTicks,'UniformOutput',0);
            % loop on rules and plot each
            for Rule=1:3
                if ~isempty(AllPSMPerf.Perf{Rule})
                    
                    % plot PSM
                    figure(varargout{1})
                    h1{Rule}=subplot(3,2,(Rule-1)*2+1);hold on % first half
                    % fit the PSM and plot the fitted function
                    obj.FitPSM(ML,AllPSMPerf.Perf{Rule}(:,HighPerfInd),Col(Rule,:)); 
                    heb1=obj.fp.PlotMeanStd(ML,AllPSMPerf.Perf{Rule}(:,HighPerfInd),[],[RuleFeatures{Rule} ' Morph Level'],'Performance(%)',Col(Rule,:),0,...
                        ['Rule' num2str(Rule)],'AxesGrid','off',PSMMethodTxt{:},'p_marker','.','p_marker_size',ErrBarMarkerSize,'p_line_width',ErrBarLineWidth,'IsthisAxisTime',0);     
                    %heb1=obj.fp.PlotMeanStd(ML,AllPSMPerf.Perf{Rule}(:,HighPerfInd),[],[RuleFeatures{Rule} ' Morph Level'],'Performance(%)',Col(Rule,:),0,...
                        %['Rule' num2str(Rule)],'p_line_style','none','AxesGrid','off',PSMMethodTxt{:},'p_marker','.','p_marker_size',ErrBarMarkerSize,'p_line_width',ErrBarLineWidth,'IsthisAxisTime',0);                    
                    heb1.CapSize=0;
                    set(gca,'XTick',ML,'XTickLabel',{HighMorphLvl});axis square;
                    set(gca,'YTick',YAxisTicks,'YTickLabel',YAxisTickLabels)
                    ylim([0 1])
                    
                    h2{Rule}=subplot(3,2,(Rule-1)*2+2);hold on % second Half
                    obj.FitPSM(ML,AllPSMPerf.Perf{Rule}(:,LowPerfInd),Col(Rule,:)) ;                
                    heb2=obj.fp.PlotMeanStd(ML,AllPSMPerf.Perf{Rule}(:,LowPerfInd),[],[RuleFeatures{Rule} ' Morph Level'],...
                        'Performance(%)',Col(Rule,:),0,['Rule' num2str(Rule)],'AxesGrid',...
                        'off','p_marker','.','p_marker_size',10,'p_line_style','none','IsthisAxisTime',0,PSMMethodTxt{:});
                    heb2.CapSize=0;
                    set(gca,'XTick',ML,'XTickLabel',{LowMorphLvl});axis square
                    set(gca,'YTick',YAxisTicks,'YTickLabel',YAxisTickLabels)
                    ylim([0 1])
                    
                    % plot combined PSM for each rule 
                    figure(varargout{2})
                    subplot(1,3,Rule);hold on
                    ThisMLPerf=nan*ones([size(AllPSMPerf.Perf{Rule}(:,LowPerfInd)),2]);
                    ThisMLPerf(:,:,1)=AllPSMPerf.Perf{Rule}(:,LowPerfInd);
                    ThisMLPerf(:,:,2)=AllPSMPerf.Perf{Rule}(:,HighPerfInd);
                    ThisMLPerf=mean(ThisMLPerf,3);
                    obj.FitPSM(ML,ThisMLPerf,Col(Rule,:)) ;                     
                    hebc=obj.fp.PlotMeanStd(ML,ThisMLPerf,...
                        [],[RuleFeatures{Rule} 'Morph Level'],'Performance(%)',Col(Rule,:),0,['Rule' num2str(Rule)],...
                        'AxesGrid','off','p_marker','.','p_marker_size',10,'p_line_style','none','IsthisAxisTime',0,PSMMethodTxt{:});
                    hebc.CapSize=0;
                    set(gca,'XTick',ML,'XTickLabel',{HighMorphLvl});axis square
                    set(gca,'YTick',YAxisTicks,'YTickLabel',YAxisTickLabels)
                    ylim([0 1])
                    
                    % Plot superimposed PSM 
                    figure(varargout{3})
                    hold on
                    ThisMLPerf=nan*ones([size(AllPSMPerf.Perf{Rule}(:,LowPerfInd)),2]);
                    ThisMLPerf(:,:,1)=AllPSMPerf.Perf{Rule}(:,LowPerfInd);
                    ThisMLPerf(:,:,2)=AllPSMPerf.Perf{Rule}(:,HighPerfInd);
                    ThisMLPerf=mean(ThisMLPerf,3);
                    obj.FitPSM(ML,ThisMLPerf,Col(Rule,:)) ;                    
                    hebc=obj.fp.PlotMeanStd(ML,ThisMLPerf,...
                        [],[RuleFeatures{Rule} 'Morph Level'],'Performance(%)',Col(Rule,:),0,['Rule' num2str(Rule)],...
                        'AxesGrid','off','p_marker','.','p_marker_size',10,'p_line_style','none','IsthisAxisTime',0,PSMMethodTxt{:});
                    hebc.CapSize=0;
                    set(gca,'XTick',ML,'XTickLabel',{HighMorphLvl});axis square
                    set(gca,'YTick',YAxisTicks,'YTickLabel',YAxisTickLabels)
                    ylim([0 1])                    
                end
            end
             
        end
        function varargout = PlotAvgPSMCongInCong(obj,IndSamp,AllPSMPerf) % plots average PSM and fits it
            CongInCong = obj.CalCongIncongPSM(IndSamp);
            PdistMorphLvls=[0,20,30,40,50,60,70,80,100,120,130,140,150,160,170,180];
            % h is axis            
            varargout=obj.fp.RenderFigure(2,[]);
            Col=distinguishable_colors(3);
            % index levels
            MorphLvl=AllPSMPerf.Level;
            HighPerfInd=1:9;LowPerfInd=10:18;
            HighMorphLvl=MorphLvl(HighPerfInd);
            LowMorphLvl=MorphLvl(LowPerfInd);
            % find the indices corresponding to what we have 
            MorphLvlInd=arrayfun(@(x) find(x==PdistMorphLvls),MorphLvl);
            %% plot for full morphline
            % loop on rules and plot each
            figure(varargout{1});
            for Rule=1:3
                    PerfInCong_ro=CongInCong{Rule}.Block.InCongCorrTrials./CongInCong{Rule}.Block.InCongTrials; % incong block
                    %PerfInCong_ro=CongInCong{Rule}.Block.CongCorrTrials./CongInCong{Rule}.Block.CongTrials;%cong
                    
                    %PerfInCong_ro=CongInCong{Rule}.InCongCorrTrials./CongInCong{Rule}.InCongTrials; % incong all                    
                    PerfInCong=PerfInCong_ro(:,MorphLvlInd);
                    PerfInCong(:,1:4)=1-PerfInCong(:,1:4);
                    PerfInCong(:,10:13)=1-PerfInCong(:,10:13);
                    PvalSigStar=CongInCong{Rule}.BinomSigTest_pvalStar(MorphLvlInd);
                    
                    IncongCountTot=CongInCong{Rule}.InCongTrials(MorphLvlInd);
                    IncongCorrCount=CongInCong{Rule}.InCongCorrTrials(MorphLvlInd);
                    % plot PSM
                    figure(varargout{1})
                    h1{Rule}=subplot(3,2,(Rule-1)*2+1);hold on % first half
                    obj.fp.PlotMeanStd(HighPerfInd,PerfInCong(:,HighPerfInd),[],'Morph Level','Performance',Col(Rule,:),0,['Rule' num2str(Rule)],'IsthisAxisTime',0)
                    arrayfun(@(x) text(x-0.5,nanmean(PerfInCong(:,x),1)+0.05,PvalSigStar{x},'fontsize',15),HighPerfInd)
                    arrayfun(@(x) text(x-0.5,nanmean(PerfInCong(:,x),1)+0.2,[num2str(IncongCorrCount(x)) '/' num2str(IncongCountTot(x))],...
                        'fontsize',7),HighPerfInd)
                   
                    set(gca,'XTick',HighPerfInd,'XTickLabel',{HighMorphLvl});axis square
                    ylim([0 1])
                    h2{Rule}=subplot(3,2,(Rule-1)*2+2);hold on % second Half
                    obj.fp.PlotMeanStd(LowPerfInd,PerfInCong(:,LowPerfInd),[],'Morph Level','Performance',Col(Rule,:),0,['Rule' num2str(Rule)],'IsthisAxisTime',0)
                    arrayfun(@(x) text(x-0.5,nanmean(PerfInCong(:,x),1)+0.05,PvalSigStar(x),'fontsize',15),LowPerfInd)
                    arrayfun(@(x) text(x-0.5,nanmean(PerfInCong(:,x),1)+0.2,[num2str(IncongCorrCount(x)) '/' num2str(IncongCountTot(x))],...
                        'fontsize',7),LowPerfInd)
                    
                    set(gca,'XTick',LowPerfInd,'XTickLabel',{LowMorphLvl});axis square
                    ylim([0 1])
                                                                         
            end
            cellfun(@(x) obj.AddPerfThr2Axis(0.7,'k',x),h1);cellfun(@(x) obj.AddPerfThr2Axis(0.7,'k',x),h2);
            
            
            %% plot combining morphelevls in 30/170 and 70/130
             % loop on rules and plot each
            figure(varargout{2});
            SqMorphLevls=[0 20 30 40 50 60 70 80 100]; % squeezed morphlevels
            nML=length(SqMorphLevls);
            for Rule=1:3
                    % combine morphlevels 
                    
                    InCongCorrTrials=obj.ReorganizetoMinPSM(CongInCong{Rule}.Block.InCongCorrTrials,PdistMorphLvls);
                    InCongTrials=obj.ReorganizetoMinPSM(CongInCong{Rule}.Block.InCongTrials,PdistMorphLvls);
                    
                    SumInCongCorrTrials=sum(InCongCorrTrials,1);SumInCongTrials=sum(InCongTrials,1);
                    CongInCongBinomSigTest=arrayfun(@(x) myBinomTest(SumInCongCorrTrials(x),SumInCongTrials(x),0.5,'two'),1:length(SqMorphLevls));
                    BinomSigTest_pvalStar=arrayfun(@(x) pvalueStar(x),CongInCongBinomSigTest,'UniformOutput',0);
                       
                    
                    PerfInCong_ro=InCongCorrTrials./InCongTrials; % incong block
                    %PerfInCong_ro=CongInCong{Rule}.Block.CongCorrTrials./CongInCong{Rule}.Block.CongTrials;%cong
                    
                    %PerfInCong_ro=CongInCong{Rule}.InCongCorrTrials./CongInCong{Rule}.InCongTrials; % incong all                    
                    PerfInCong=PerfInCong_ro;
                    PerfInCong(:,1:4)=1-PerfInCong(:,1:4);
                    PvalSigStar=BinomSigTest_pvalStar;
                    
                    IncongCountTot=SumInCongTrials;
                    IncongCorrCount=SumInCongCorrTrials;
                    % plot PSM
                    figure(varargout{1})
                    h1{Rule}=subplot(1,3,Rule);hold on % first half
                    obj.fp.PlotMeanStd(1:nML,PerfInCong,[],'Morph Level','Performance',Col(Rule,:),0,['Rule' num2str(Rule)],'IsthisAxisTime',0)
                    arrayfun(@(x) text(x-0.5,nanmean(PerfInCong(:,x),1)+0.05,PvalSigStar{x},'fontsize',15),1:nML)
                    arrayfun(@(x) text(x-0.5,nanmean(PerfInCong(:,x),1)+0.2,[num2str(IncongCorrCount(x)) '/' num2str(IncongCountTot(x))],...
                        'fontsize',7),1:nML)                   
                    set(gca,'XTick',1:nML,'XTickLabel',{SqMorphLevls});axis square
                    ylim([0 1])                                                                                             
            end
            cellfun(@(x) obj.AddPerfThr2Axis(0.7,'k',x),h1);
           
        end
        function out = ReorganizetoMinPSM(obj,PSM,MorphLevels) % reorganize to the minimal PSM
            % target morphlevels would be 0 20 30 40 50 60 70 80 100
            MLs={0,[20 180],[30 170],[40 160],[50 150],[60 140],[70 130],[80 120],100};
            for m=1:length(MLs)
                ML=MLs{m};
                if length(ML)==1
                   out(:,m)=PSM(:,MorphLevels==ML);
                else
                   out(:,m)=sum(PSM(:,(MorphLevels==ML(1) | MorphLevels==ML(2))),2);
                end
            end
            
        end
        function varargout = PlotTrlPerf(obj, AllTrlPerf,AllTrlCount,AllTrlCountDay,NBlocksDay,RewardPulse,NCorrectTrl,NumRewards,AllSeqHist)
            global AnalysisOpts
            
            obj.fp.font_size=AnalysisOpts.Plotting.font_size;
            RuleCol= AnalysisOpts.RuleColors; % color of the rules
            SeqCol=AnalysisOpts.ColorPalettCamden; % sequence colors;
            ColSet={RuleCol SeqCol};
            % pick colors
            PickCol{1}=@(Rule) ColSet{1}(Rule,:);
            PickCol{2}=@(Rule) ColSet{2}(1,:);
            PickCol{3}=@(Rule) ColSet{2}(2,:);
            
            LengthRule=[60 35 60];
            FieldNames={'All'};%fieldnames(AllTrlPerf.FromSwitch)';
            %% titles for different sequences for different rules 
            SeqHistTxt{1}={'',' 121 Sequence',' 321 Sequence'};
            SeqHistTxt{2}={'',' 12 Sequence' ,' 32 Sequence'};
            SeqHistTxt{3}={'',' 323 Sequence',' 123 Sequence'};
            SeqSgTitleTxt={['All Performance ' AnalysisOpts.BhvAna.Kidname],['Sequence Performance ' AnalysisOpts.BhvAna.Kidname]};
            % subplot orders 
            Seq_Sp{1}=@(x) {[2,3,1],[2,3,2],[2,3,3+x]};
            Seq_Sp{2}=@(x) {[2,3,x],[nan],[2,3,3+x]};
            Seq_Sp{3}=@(x) {[2,3,x],[nan],[2,3,3+x]};
            % stattest orders
            ChanceLevel=0.5; % define chance level as 50% 
            Stat_Seq{1}=@(Rule,Perf) ChanceLevel*ones(size(Perf{1}{Rule}'));
            Stat_Seq{2}=@(Rule,Perf) nan;
            Stat_Seq{3}=@(Rule,Perf) Perf{2}{Rule}';
            % specs for point stat test trials(only for from switch) 
            PointStatTestTrials=[1:20]; % this will include the first 30 trials trials 1:15 and trials 16:30
                      
            varargout=obj.fp.RenderFigure(2,[]);            
            for S=[-1 0 1]                
                figure(varargout{double(S>-1)+1})
                SeqHist=S;
                SeqHistSet=[-1,0,1];
                nSeq=[-1 0 1]==S;
                switch SeqHist
                    case {0,1}
                        SeqHistBlks=cellfun(@(x) x==SeqHist,AllSeqHist,'UniformOutput',0);
                      %  SeqHistBlks{2}=logical(ones(1,length(AllSeqHist{2})));
                    case -1
                        SeqHistBlks=cellfun(@(x) logical(ones(1,length(x))),AllSeqHist,'UniformOutput',0);
                end
                for Rule=1:3
               %     sgtitle(SeqSgTitleTxt{double(S>-1)+1})
                    Field=FieldNames;
                    Field=Field{1};n=find(strcmp(FieldNames,Field));
                                      
                    % plot from switch data
                    Sp=Seq_Sp{nSeq}(Rule);Sp=Sp{1};
                    if ~isnan(Sp)
                        subplot(Sp(1),Sp(2),Sp(3))
                        yyaxis left
                        ThisRulePerfFromSwitch{nSeq}{Rule}=AllTrlPerf.FromSwitch(Rule).(Field)(1:LengthRule(Rule),SeqHistBlks{Rule})';
                        ThisRuleTrlCntFromSwitch{nSeq}{Rule}=AllTrlPerf.FromSwitch(Rule).AllTrlCnt(1:LengthRule(Rule),SeqHistBlks{Rule})';
                        ThisRuleIndTrlsFromSwitch{nSeq}{Rule}=AllTrlPerf.FromSwitch(Rule).IndTrls(1:LengthRule(Rule),SeqHistBlks{Rule})';

                        obj.fp.PlotMeanStd([],ThisRulePerfFromSwitch{nSeq}{Rule},[],'Trials from Switch','Perf',PickCol{nSeq}(Rule),6,...
                            ['Rule' num2str(Rule)],'AxesGrid','off','LegendTxt',['Rule' num2str(Rule) SeqHistTxt{Rule}{nSeq}],'p_line_style','-',...
                            'STD_method','sem','IsthisAxisTime',0,'p_marker','none');

                        % add chance level line
                        obj.fp.PlotHorizontalLine(0.25,gca,[0.5 0.5 0.5],'p_line_style','--')    
                        axis tight
                        ax1 = gca;
                        ax1.YAxis(2).Visible = 'off';
                        ax1.YAxis(1).Color = 'k';
                       % AxXtic= ax1.XTickLabel;
                        ax1.XTick=[0 20 40 60 ];
                        ax1.XTickLabel={'0','20','40','60'};%arrayfun(@num2str,ax1.XTick,'UniformOutput',0);
                       
                        % perfrom stat test and plot it 
                        if SeqHist==-1 % we are using all of the trials use the calculated binomial test 
                            StatTest.FromSwitch=AllTrlPerf.FromSwitch(Rule).BinomialStatTest(1:LengthRule(Rule));
                            % compare two points during learning P1=trial 1(15) and P2=trial 60(75)
                            [StatTestPoint.FromSwitch.AllTrls(Rule),StatTestPoint.FromSwitch.AllTrls_P1(Rule),StatTestPoint.FromSwitch.AllTrls_P2(Rule)]=...
                                obj.ChiSquareStatTestAllTrlPerf(ThisRuleTrlCntFromSwitch{1}{Rule}(:,1),ThisRuleTrlCntFromSwitch{1}{Rule}(:,LengthRule(Rule)),AllTrlPerf.NavgPerf);
                        elseif SeqHist==1 % then we are doing a chisquare test on same and different block sequences
                            StatTest.FromSwitch=obj.ChiSquareStatTestAllTrlPerf(ThisRuleTrlCntFromSwitch{2}{Rule}(:,1:LengthRule(Rule)),ThisRuleTrlCntFromSwitch{3}{Rule}(:,1:LengthRule(Rule)),AllTrlPerf.NavgPerf);
                            % compare two points during learning of same and different P1=trial 1(15)(123) and P2=trial 60(75)
                            [StatTestPoint.FromSwitch.BlockSeq(Rule),StatTestPoint.FromSwitch.BlockSeq_P1(Rule),StatTestPoint.FromSwitch.BlockSeq_P2(Rule)]=obj.ChiSquareStatTestAllTrlPerf(sum(ThisRuleIndTrlsFromSwitch{2}{Rule}(:,PointStatTestTrials),2),sum(ThisRuleIndTrlsFromSwitch{3}{Rule}(:,PointStatTestTrials),2),length(PointStatTestTrials));
                        end

                        if sum(SeqHist==[-1 1])
                            % plot statistica test
                            [X,P]=obj.ManData.DifferentiateSigClusters2(StatTest.FromSwitch);
                            obj.fp.plot_significance_level(X,P,1:LengthRule(Rule),[0 1.05 0.05],PickCol{nSeq}(Rule),[],0,...
                                'SmoothingMethod','','WidthSmoothing',obj.NavgPerf,'IsthisAxisTime',0);
                        end                      
                        obj.fp.FormatAxes(gca,'daspect',[2*75 1 1])  
                    end
                    
                    % plot to switch data
                     Sp=Seq_Sp{nSeq}(Rule);Sp=Sp{2};
                     if ~isnan(Sp)
                         subplot(Sp(1),Sp(2),Sp(3))
                         yyaxis right
                         hold on
                         ThisRulePerfToSwitch{nSeq}{Rule}=AllTrlPerf.ToSwitch(Rule).(Field)(end-LengthRule(Rule)+1:end,SeqHistBlks{Rule})';
                         ThisRuleTrlCntToSwitch{nSeq}{Rule}=AllTrlPerf.ToSwitch(Rule).AllTrlCnt(end-LengthRule(Rule)+1:end,SeqHistBlks{Rule})';                       
                         ThisRuleIndTrlsToSwitch{nSeq}{Rule}=AllTrlPerf.ToSwitch(Rule).IndTrls(end-LengthRule(Rule)+1,SeqHistBlks{Rule})';
                         

                         leg(n)=obj.fp.PlotMeanStd((LengthRule(1)-LengthRule(Rule)+1):LengthRule(1),ThisRulePerfToSwitch{nSeq}{Rule},[],'Trials to Switch','Perf',PickCol{nSeq}(Rule),6,['Rule' num2str(Rule)],...
                             'AxesGrid','off','LegendTxt',['Rule' num2str(Rule) SeqHistTxt{Rule}{nSeq}],'p_line_style','-','IsthisAxisTime',0,'STD_method','sem','p_marker','none');
                         % add chance level line
                         obj.fp.PlotHorizontalLine(0.25,gca,[0.5 0.5 0.5],'p_line_style','--')
                         axis tight
                         ax1 = gca;
                         ax1.YAxis(1).Visible = 'off';
                        % axis tight
                         ax1.YAxis(2).Color = 'k';
%                        AxXtic= ax1.XTickLabel;
%                        xticklabels(AxXtic(end:-1:1))
                         ax1.XTick=[0 20 40 60];
                         ax1.XTickLabel={'60','40','20','0'};%arrayfun(@num2str,[75 60:-20:0],'UniformOutput',0);

                         % perfrom stat test and plot it 
                        if SeqHist==-1 % we are using all of the trials use the calculated binomial test 
                            StatTest.ToSwitch=AllTrlPerf.ToSwitch(Rule).BinomialStatTest(end-LengthRule(Rule)+1:end);
                            % compare two points during learning P1=trial -1(15) and P2=trial -60(75)
                            StatTestPoint.ToSwitch.AllTrls(Rule)=obj.ChiSquareStatTestAllTrlPerf(ThisRuleTrlCntToSwitch{1}{Rule}(:,1),ThisRuleTrlCntToSwitch{1}{Rule}(:,LengthRule(Rule)),AllTrlPerf.NavgPerf);
                        elseif SeqHist==1 % then we are doing a chisquare test on same and different block sequences
                            StatTest.ToSwitch=obj.ChiSquareStatTestAllTrlPerf(ThisRuleTrlCntToSwitch{2}{Rule}(:,1:LengthRule(Rule)),ThisRuleTrlCntToSwitch{3}{Rule}(:,1:LengthRule(Rule)),AllTrlPerf.NavgPerf);
                            % compare two points during learning of same and different  first 30 trials
                            StatTestPoint.ToSwitch.BlockSeq(Rule)=obj.ChiSquareStatTestAllTrlPerf(sum(ThisRuleIndTrlsToSwitch{2}{Rule}(:,PointStatTestTrials)),sum(ThisRuleIndTrlsToSwitch{3}{Rule}(:,PointStatTestTrials),2),length(PointStatTestTrials));
                        end
                        
                        if sum(SeqHist==[-1 1])
                            % plot statistica test
                            [X,P]=obj.ManData.DifferentiateSigClusters2(StatTest.ToSwitch);
                            obj.fp.plot_significance_level(X,P,(LengthRule(1)-LengthRule(Rule)+1):LengthRule(1),[0 1.05 0.05],PickCol{nSeq}(Rule),[],0,...
                                'SmoothingMethod','','WidthSmoothing',obj.NavgPerf,'IsthisAxisTime',0);
                        end
                         % perfrom stat test and plot it
%                          [clusters, ~, ~, ~,statsummery]=obj.PerformClusterCorrected_tTest(ThisRulePerfToSwitch{nSeq}{Rule}', Stat_Seq{nSeq}(Rule,ThisRulePerfToSwitch), false,[], [], [], []);
%                          try
%                              obj.fp.plot_significance_level(clusters,statsummery,(LengthRule(1)-LengthRule(Rule)+1):LengthRule(1),[0 1.05 0.05],PickCol{nSeq}(Rule),[],0,...
%                                  'SmoothingMethod','movmean','WidthSmoothing',obj.NavgPerf,'IsthisAxisTime',0);
%                          catch me
%                              disp(me.message)
%                          end
                         obj.fp.FormatAxes(gca,'daspect',[2*75 1 1])
                     end
                     
                    % plot distribution of performances for rules
                     Sp=Seq_Sp{nSeq}(Rule);Sp=Sp{3};
                     if ~isnan(Sp)                         
                         subplot(Sp(1),Sp(2),Sp(3))
                         hold on
                         Begdata{nSeq}{Rule}=ThisRulePerfFromSwitch{nSeq}{Rule}(:,1:25);
                         obj.fp.HistogramPlot(Begdata{nSeq}{Rule}(:),[0:0.1:1],'none','Performance','Proportion(%)',...
                             ['Rule' num2str(Rule) ' Start Block Performance' ],'bar_edgecol',PickCol{nSeq}(Rule),...
                             'bar_width',0.5+find(nSeq)/10,'LegendTxt',['Rule' num2str(Rule) SeqHistTxt{Rule}{nSeq}]);
                     end                     
                end
                grid off
            end
            
            
        end
        function varargout = PlotSampInfo(obj,IndSamp)
            varargout=obj.fp.RenderFigure(1,[]);
            
            % prepare vars
            Col={'b','r','g'};
            Pdist=IndSamp.Pdist;
            L=1:length(Pdist.Dim1Objs);
            figure(varargout{1})
            for Rule=1:3
                ThisRulePerf=nanmean(IndSamp.Perf{Rule},3)';
                ThisRuleTrlCnt=nanmean(IndSamp.Count{Rule},3);
                ThisRuleRT=nanmean(IndSamp.RT{Rule},3)';
                ThisRulepstar=IndSamp.p_Binostar{Rule}';
                TotalTrials=IndSamp.TotalTrials{Rule}';
                TotalCorrTrials=IndSamp.TotalCorrTrials{Rule}';
                %% plot ind samp Perf
                subplot(3,3,Rule)
                hold on
                imagesc(ThisRulePerf);
                set(gca,'YDir','normal')
                set(gca,'Xtick',L,'XtickLabel',{Pdist.Dim1Objs},'Ytick',L,'YtickLabel',{Pdist.Dim1Objs});
                colorbar
                xlabel('Shape');ylabel('Color')
                title(['Rule' num2str(Rule)])
                %                 A=zeros(size(ThisRulePerf',1)+1,size(ThisRulePerf',2)+1);
                %                 A(2:end,2:end)=ThisRulePerf';
                %                 A(2:end,1)=Pdist.Dim1Objs;
                %                 A(1,2:end)=Pdist.Dim1Objs;
                d1=1;
                for Dim1ML=Pdist.Dim1Objs
                    d2=1;
                    for Dim2ML=Pdist.Dim1Objs
                        text(d2-.5,d1,ThisRulepstar{d1,d2});
                        d2=d2+1;
                    end
                    d1=d1+1;
                end
                %% plot counts
                subplot(3,3,Rule+3);hold on
                bar(sum(ThisRuleTrlCnt,2),'facecolor','none','edgecolor','b')
                bar(sum(ThisRuleTrlCnt,1),'facecolor','none','edgecolor','r')
                set(gca,'Xtick',L,'XtickLabel',{Pdist.Dim1Objs},'fontsize',5);
                xlabel('Morph Level');ylabel('Count');axis tight
                legend({'Shape','Color'})
                hold on
                title(['Rule' num2str(Rule)])
                                                
                %% plot RTs
                subplot(3,3,Rule+6)
                hold on
                imagesc(ThisRuleRT);
                set(gca,'YDir','normal')
                set(gca,'Xtick',L,'XtickLabel',{Pdist.Dim1Objs},'Ytick',L,'YtickLabel',{Pdist.Dim1Objs});
                colorbar
                xlabel('Shape');ylabel('Color');axis tight
                title(['Rule' num2str(Rule)])
                d1=1;
                
                for Dim1ML=Pdist.Dim1Objs
                    d2=1;
                    for Dim2ML=Pdist.Dim1Objs
                        text(d2-.5,d1,sprintf('%i/%i',TotalCorrTrials(d1,d2),TotalTrials(d1,d2)),'fontsize',7);
                        d2=d2+1;
                    end
                    d1=d1+1;
                end
                
            end
            
            
            
        end
        function varargout = PlotSampPerfInfo(obj,IndSamp) % plots samples with performance 
            varargout=obj.fp.RenderFigure(1,[]);
            
            % prepare vars
            Col={'b','r','g'};
            Pdist=IndSamp.Pdist;
            L=1:length(Pdist.Dim1Objs);
            figure(varargout{1})
            for Rule=1:3
                ThisRulePerf=nanmean(IndSamp.Perf{Rule},3)';
                ThisRuleTrlCnt=nanmean(IndSamp.Count{Rule},3);
                ThisRuleRT=nanmean(IndSamp.RT{Rule},3)';
                ThisRulepstar=IndSamp.p_Binostar{Rule}';
                TotalTrials=IndSamp.TotalTrials{Rule}';
                TotalCorrTrials=IndSamp.TotalCorrTrials{Rule}';
                %% plot ind samp Perf
                subplot(2,3,Rule)
                hold on
                imagesc(ThisRulePerf);
                set(gca,'YDir','normal')
                set(gca,'Xtick',L,'XtickLabel',{Pdist.Dim1Objs},'Ytick',L,'YtickLabel',{Pdist.Dim1Objs},'fontsize',7);
                colorbar
                xlabel('Shape');ylabel('Color')
                title(['Rule' num2str(Rule)])
                axis tight
                %                 A=zeros(size(ThisRulePerf',1)+1,size(ThisRulePerf',2)+1);
                %                 A(2:end,2:end)=ThisRulePerf';
                %                 A(2:end,1)=Pdist.Dim1Objs;
                %                 A(1,2:end)=Pdist.Dim1Objs;
                d1=1;
                for Dim1ML=Pdist.Dim1Objs
                    d2=1;
                    for Dim2ML=Pdist.Dim1Objs
                        text(d2-.5,d1,ThisRulepstar{d1,d2});
                        d2=d2+1;
                    end
                    d1=d1+1;
                end
                                                              
                %% plot counts
                subplot(2,3,Rule+3)
                hold on
                imagesc(ones(size(ThisRuleRT)));               
                set(gca,'YDir','normal')
                set(gca,'Xtick',L,'XtickLabel',{Pdist.Dim1Objs},'Ytick',L,'YtickLabel',{Pdist.Dim1Objs},'fontsize',7);
                colorbar
                xlabel('Shape');ylabel('Color');axis tight
                title(['Rule' num2str(Rule)])
                d1=1;
                grid on
                
                for Dim1ML=Pdist.Dim1Objs
                    d2=1;
                    for Dim2ML=Pdist.Dim1Objs
             %           text(d2-.5,d1,sprintf('%i/%i',TotalCorrTrials(d1,d2),TotalTrials(d1,d2)),'fontsize',7);
                        text(d2-.5,d1,sprintf('%i',TotalTrials(d1,d2)),'fontsize',7);                        
                        d2=d2+1;
                    end
                    d1=d1+1;
                end
                
            end
            
            
            
        end
        function CongInCong = CalCongIncongPSM(obj,IndSamp) % plots samples with performance 
            
            % prepare vars
            Col={'b','r','g'};
            Pdist=IndSamp.Pdist;
            L=1:length(Pdist.Dim1Objs);
            for Rule=1:3              
                if Rule==1
                    for l=L
                        ThisCong=find(IndSamp.Congmap.Map(l,:)==1);
                        ThisInCong=find(IndSamp.Congmap.Map(l,:)==2);
                        % all trials
                        CongInCong{Rule}.CongTrials(l)=sum(IndSamp.TotalTrials{Rule}(l,ThisCong));
                        CongInCong{Rule}.InCongTrials(l)=sum(IndSamp.TotalTrials{Rule}(l,ThisInCong));
                        % correct trials
                        CongInCong{Rule}.CongCorrTrials(l)=sum(IndSamp.TotalCorrTrials{Rule}(l,ThisCong));
                        CongInCong{Rule}.InCongCorrTrials(l)=sum(IndSamp.TotalCorrTrials{Rule}(l,ThisInCong));
                        CongInCong{Rule}.BinomSigTest(l)=myBinomTest(CongInCong{Rule}.InCongCorrTrials(l),CongInCong{Rule}.InCongTrials(l),0.5,'two');
                        CongInCong{Rule}.BinomSigTest_pvalStar{l}=pvalueStar(CongInCong{Rule}.BinomSigTest(l));
                        
                        CongInCong{Rule}.BinomSigTest_Cong(l)=myBinomTest(CongInCong{Rule}.CongCorrTrials(l),CongInCong{Rule}.CongTrials(l),0.5,'two');
                        CongInCong{Rule}.BinomSigTest_pvalStar_Cong{l}=pvalueStar(CongInCong{Rule}.BinomSigTest_Cong(l));
                       
                        %% do blockwise 
                        %all trials
                        CongInCong{Rule}.Block.CongTrials(:,l)=squeeze(sum(IndSamp.Count{Rule}(l,ThisCong,:),2));
                        CongInCong{Rule}.Block.InCongTrials(:,l)=squeeze(sum(IndSamp.Count{Rule}(l,ThisInCong,:),2));
                        % correct trials
                        CongInCong{Rule}.Block.CongCorrTrials(:,l)  =squeeze(sum(IndSamp.CorrCount{Rule}(l,ThisCong,:),2));
                        CongInCong{Rule}.Block.InCongCorrTrials(:,l)=squeeze(sum(IndSamp.CorrCount{Rule}(l,ThisInCong,:),2));
                    end
                elseif Rule==2 | Rule==3   
                    for l=L
                        ThisCong=find(IndSamp.Congmap.Map(:,l)==1);
                        ThisInCong=find(IndSamp.Congmap.Map(:,l)==2);
                        % all trials
                        CongInCong{Rule}.CongTrials(l)=sum(IndSamp.TotalTrials{Rule}(ThisCong,l));
                        CongInCong{Rule}.InCongTrials(l)=sum(IndSamp.TotalTrials{Rule}(ThisInCong,l));
                        % correct trials
                        CongInCong{Rule}.CongCorrTrials(l)=sum(IndSamp.TotalCorrTrials{Rule}(ThisCong,l));
                        CongInCong{Rule}.InCongCorrTrials(l)=sum(IndSamp.TotalCorrTrials{Rule}(ThisInCong,l));
                        CongInCong{Rule}.BinomSigTest_Cong(l)=myBinomTest(CongInCong{Rule}.CongCorrTrials(l),CongInCong{Rule}.CongTrials(l),0.5,'two');
                        CongInCong{Rule}.BinomSigTest_pvalStar_Cong{l}=pvalueStar(CongInCong{Rule}.BinomSigTest_Cong(l));
                     
                        CongInCong{Rule}.BinomSigTest(l)=myBinomTest(CongInCong{Rule}.InCongCorrTrials(l),CongInCong{Rule}.InCongTrials(l),0.5,'two');
                        CongInCong{Rule}.BinomSigTest_pvalStar{l}=pvalueStar(CongInCong{Rule}.BinomSigTest(l));
                     
                        %% do blockwise 
                        %all trials
                        CongInCong{Rule}.Block.CongTrials(:,l)=squeeze(sum(IndSamp.Count{Rule}(ThisCong,l,:),1));
                        CongInCong{Rule}.Block.InCongTrials(:,l)=squeeze(sum(IndSamp.Count{Rule}(ThisInCong,l,:),1));
                        % correct trials
                        CongInCong{Rule}.Block.CongCorrTrials(:,l)  =squeeze(sum(IndSamp.CorrCount{Rule}(ThisCong,l,:),1));
                        CongInCong{Rule}.Block.InCongCorrTrials(:,l)=squeeze(sum(IndSamp.CorrCount{Rule}(ThisInCong,l,:),1));
                    end
                end
            end
            
            
            
        end
        %% Th Analysis functions 
        function [fh,h1,h2] = PlotAvgPSM_ThAnalysis(obj,ThTrl,AllPSMPerf,Col,fh,h1,h2) % plots average PSM and fits it
            % h is axis
            % index levels
            MorphLvl=AllPSMPerf.Level;
            HighPerfInd=1:9;LowPerfInd=10:18;
            HighMorphLvl=MorphLvl(HighPerfInd);
            LowMorphLvl=MorphLvl(LowPerfInd);
            
            % loop on rules and plot each
            for Rule=1:3
                if ~isempty(AllPSMPerf.Perf{Rule})
                    
                    % plot PSM
                    figure(fh{1});
                    if isempty(h1{Rule})
                        h1{Rule}=subplot(3,2,(Rule-1)*2+1);
                    else
                        axes(h1{Rule})
                        hold(h1{Rule},'on') % first half
                    end
                    obj.fp.PlotMeanStd(1:9,AllPSMPerf.Perf{Rule}(:,HighPerfInd),[],'Morph Level','Performance',Col,3,['Rule' num2str(Rule)])
                    set(gca,'XTick',HighPerfInd,'XTickLabel',{HighMorphLvl});axis square;axis tight
                    ylim([0.2 0.8])
                    if isempty(h2{Rule})
                        h2{Rule}=subplot(3,2,(Rule-1)*2+2); % second Half
                    else
                        axes(h2{Rule})
                        hold(h2{Rule},'on') % first half
                    end
                    obj.fp.PlotMeanStd(LowPerfInd,AllPSMPerf.Perf{Rule}(:,LowPerfInd),[],'Morph Level','Performance',Col,3,['Rule' num2str(Rule)])
                    set(gca,'XTick',LowPerfInd,'XTickLabel',{LowMorphLvl});axis square;axis tight
                   ylim([0.2 0.8])                                      
                end
            end
            cellfun(@(x) obj.AddPerfThr2Axis(0.7,'k',x),h1);cellfun(@(x) obj.AddPerfThr2Axis(0.7,'k',x),h2);
           % legend(num2str(ThTrl),'Location','best')
         end
        function varargout = PlotTrlPerf_ThAnalysis(obj,Col,ThId,fh,AllTrlPerf,AllTrlCount,AllTrlCountDay,NBlocksDay,RewardPulse,NCorrectTrl,NumRewards)
            figure(fh{1})
            Colx=distinguishable_colors(3);
            for Rule=1:3
                subplot(2,3,Rule);hold on
                ThisRulePerf=AllTrlPerf{Rule}';
                obj.fp.PlotMeanStd([],ThisRulePerf,[],'Trial','Performance',Col,3,['Rule' num2str(Rule)])
                ylim([0 1])
            end
            %% plot total number of attempted trails
            figure(fh{2})
            subplot(231) % plot number of trials for rules
            hold on
            for Rule=[1 3]
                plot(Rule*ones(1,length(AllTrlCount{Rule}))+ThId*0.3,AllTrlCount{Rule},'.','color',Col,'MarkerSize',10)
                plot(Rule+ThId*0.3,mean(AllTrlCount{Rule}),'*','color',Col,'MarkerSize',15);
           %     obj.fp.PlotMeanStd(Rule,AllTrlCount{Rule}',[],'Rule','Num Attempted Trls',Col,0,[]);
            end
            xticks([1.7 3.7])
            xticklabels({'1','3'});
            xlim([0 5])
            ylim([50 400])
            xlabel('Rule')
            ylabel('Trials to switch')
            title('Number of trials to switch across days')
            %% plot total number of trials across days
            subplot(232)
            bar(1:length(AllTrlCountDay),AllTrlCountDay,'FaceColor','r')
            text(1:length(AllTrlCountDay),AllTrlCountDay+50,[arrayfun(@num2str,RewardPulse,'UniformOutput',0)]);
         %   text(1:length(AllTrlCountDay),AllTrlCountDay+70,[arrayfun(@num2str,NumRewards.Switch,'UniformOutput',0)] )
            xlabel('Day #')
            ylabel('# attempted trials')
            title('Total Number of trials across days')
           
            %% plot average number of blocks over dats for rules
            subplot(233)
            hold on
            arrayfun(@(x) bar(x,mean(NBlocksDay(:,x),1),'FaceColor','none','EdgeColor',Colx(x,:),'BarWidth',0.5),1:3);
            arrayfun(@(x) errorbar(x,mean(NBlocksDay(:,x),1),std(NBlocksDay(:,x))/sqrt(length(NBlocksDay(:,x))),'color',Colx(x,:)),1:3);
            legend({'R1','R2','R3'})
            ylabel('avg # finished blks')
            title('Avg Number of Finished Blks/Rule')    
            varargout=fh;
        end
        function varargout = PlotAvgPSMCongInCong_ThAnalysis(obj,fh,ThTrl,Col,IndSamp,AllPSMPerf) % plots average PSM and fits it
             CongInCong = obj.CalCongIncongPSM(IndSamp);
             PdistMorphLvls=[0,20,30,40,50,60,70,80,100,120,130,140,150,160,170,180];
            % h is axis
            % index levels
            MorphLvl=AllPSMPerf.Level;
            HighPerfInd=1:9;LowPerfInd=10:18;
            HighMorphLvl=MorphLvl(HighPerfInd);
            LowMorphLvl=MorphLvl(LowPerfInd);
            % find the indices corresponding to what we have 
            MorphLvlInd=arrayfun(@(x) find(x==PdistMorphLvls),MorphLvl);            
            % loop on rules and plot each
            for Rule=1:3
                    PerfInCong_ro=CongInCong{Rule}.Block.InCongCorrTrials./CongInCong{Rule}.Block.InCongTrials; % incong block
                %    PerfInCong_ro=CongInCong{Rule}.Block.CongCorrTrials./CongInCong{Rule}.Block.CongTrials;%cong
                    
                 %   PerfInCong_ro=CongInCong{Rule}.InCongCorrTrials./CongInCong{Rule}.InCongTrials; % incong all
                    
                    PerfInCong=PerfInCong_ro(:,MorphLvlInd);
                    PerfInCong(:,1:4)=1-PerfInCong(:,1:4);
                    PerfInCong(:,10:13)=1-PerfInCong(:,10:13);
                    PvalSigStar=CongInCong{Rule}.BinomSigTest_pvalStar(MorphLvlInd);
                    
                    IncongCountTot=CongInCong{Rule}.InCongTrials(MorphLvlInd);
                    IncongCorrCount=CongInCong{Rule}.InCongCorrTrials(MorphLvlInd);
                    % plot PSM
                    figure(fh{1})
                    h1{Rule}=subplot(3,2,(Rule-1)*2+1);hold on % first half
                    obj.fp.PlotMeanStd(HighPerfInd,PerfInCong(:,HighPerfInd),[],'Morph Level','Performance for Incong',Col,0,['Rule' num2str(Rule) ',' num2str(ThTrl) 'trls'])
                    arrayfun(@(x) text(x-0.5,nanmean(PerfInCong(:,x),1)+0.05,PvalSigStar{x},'fontsize',15),HighPerfInd)
                    arrayfun(@(x) text(x-0.5,nanmean(PerfInCong(:,x),1)+0.2,[num2str(IncongCorrCount(x)) '/' num2str(IncongCountTot(x))],...
                        'fontsize',7),HighPerfInd)
                   
                    set(gca,'XTick',HighPerfInd,'XTickLabel',{HighMorphLvl});axis square
                    ylim([0 1])
                    h2{Rule}=subplot(3,2,(Rule-1)*2+2);hold on % second Half
                    obj.fp.PlotMeanStd(LowPerfInd,PerfInCong(:,LowPerfInd),[],'Morph Level','Performance for Incong',Col,0,['Rule' num2str(Rule) ',' num2str(ThTrl) 'trls'])
                    arrayfun(@(x) text(x-0.5,nanmean(PerfInCong(:,x),1)+0.05,PvalSigStar(x),'fontsize',15),LowPerfInd)
                    arrayfun(@(x) text(x-0.5,nanmean(PerfInCong(:,x),1)+0.2,[num2str(IncongCorrCount(x)) '/' num2str(IncongCountTot(x))],...
                        'fontsize',7),LowPerfInd)
                    
                    set(gca,'XTick',LowPerfInd,'XTickLabel',{LowMorphLvl});axis square
                    ylim([0 1])
                                                                         
            end
            cellfun(@(x) obj.AddPerfThr2Axis(0.7,'k',x),h1);cellfun(@(x) obj.AddPerfThr2Axis(0.7,'k',x),h2);
            varargout=fh;
        end       

        %% AUX functions
        function AddPerfThr2Axis(obj,Th,Col,h) % adds a performance threshold to axis
            
            axes(h);
            v=axis;
            plot([v(1) v(2)],[Th Th],'color',Col);
            plot([v(1) v(2)],[1-Th 1-Th],'color',Col);
            
        end
     %   function TrimPSMbyPerf(obj,Perf)   % remove blocks that the performance        
        %% Pschofit functions (functions used to fit pscychometric curves 
        %% uses Pschofit toolbox made by Matteo Carrandini 
        function pars=FitPSM(obj,xx,Perf,Col) % fit 'erf function' to psychometric data 
             
            if isempty(xx)
               xx=[0 20 30 40 50 60 70 80 100]/100; 
            end
           
            MeanPerf=mean(Perf,1);
            OkInds=~isnan(MeanPerf);
            MeanPerf=MeanPerf(OkInds);
            xx=xx(OkInds);
            nxx = length(xx); 
            % fit to reconstruct the parameters
            pars = mle_fit_psycho([xx; ones(1,nxx);MeanPerf],'erf_psycho');                       
            
            % Plot fitted function  
            obj.fp.PlotMeanStd([0:0.01:1],erf_psycho(pars,[0:0.01:1]),[],'Morph Level','Performance',Col,3,[],'AxesGrid','off')

        end
        
        %% statistics funcs
        function [SigTest,p_SigTest]=TtestStatTeston3Dmat(~,data,mu,Tail) % performs statistical test on 3D matrix; x and y as dimensions Z as samples
            [subi,subj]=arrayfun(@(x) ind2sub([size(data,1) size(data,2)],x),1:size(data,1)*size(data,2));
            [SigTest,p_SigTest]=arrayfun(@(x) ttest(squeeze(data(subi(x),subj(x),:)),mu,'tail',Tail),1:length(subi));
            SigTest=reshape(SigTest,[size(data,1) size(data,2)]);
            p_SigTest=reshape(p_SigTest,[size(data,1) size(data,2)]);
        end
        function [p_SigTest,p_star]=BinomialStatTeston3Dmat(~,correct,total,p,Tail) % performs binomial statistical test on 3D matrix; x and y as dimensions Z as samples
           
            [subi,subj]=arrayfun(@(x) ind2sub([size(total,1) size(total,2)],x),1:size(total,1)*size(total,2));            
            [p_SigTest]=arrayfun(@(x) myBinomTest(correct(subi(x),subj(x)),total(subi(x),subj(x)),p,'two'),1:length(subi));
             p_star=arrayfun(@(x) pvalueStar(p_SigTest(x)),1:length(subi),'UniformOutput',0);
             p_SigTest=reshape(p_SigTest,[size(total,1) size(total,2)]);
             p_star=reshape(p_star,[size(total,1) size(total,2)]);
        end
        function [clusters, p_values, t_sums, permutation_distribution,statsummery ]=PerformClusterCorrected_tTest(obj,trial_group_1, trial_group_2, dependent_samples,p_threshold, num_permutations, two_sided, num_clusters,varargin)
            % Uses cluster-corrected t-test to compare distributions. 
            % Permutation test (non-parametric test for significance) for dependent or independent measures of 1-D or 2-D data.
            % Based on Maris & Oostenveld 2007 for 1-D and 2-D vectors. The test statistic is T-Sum - the total of t-values within a cluster of contingent above-threshold data points. See:
            % Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of EEG-and MEG-data. Journal of Neuroscience Methods, 164(1), 177–190. https://doi.org/10.1016/j.jneumeth.2007.03.024
            % This function is using permutest witten by Edden M.Gerber           
            % refer to permutest.m for description of inputs 
           
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
          
            if isnan(sum(trial_group_1(:))) | isnan(sum(trial_group_2(:))) % if any of these are nan then return
                clusters=nan; p_values=nan; t_sums=nan; permutation_distribution=nan; statsummery=nan;
                return
            end
            % turn empty values into default 
            if isempty(dependent_samples); dependent_samples=false;end%indicate independant samples;
            if isempty(p_threshold);p_threshold=0.05;end
            if isempty(two_sided);two_sided=true;end%do two sided t test
            if isempty(num_clusters);num_clusters=inf;end % detect infinit number of trials
            if isempty(num_permutations);num_permutations=1000;end % number of permutations should be >1000
            
            [clusters, p_values, t_sums, permutation_distribution] =...
                permutest( trial_group_1, trial_group_2, dependent_samples, ...
                p_threshold, num_permutations, two_sided, num_clusters );
            % create a summery matrix where for each time point in the clusters we assign 
            % a p-value 
            if ~isempty(clusters) 
                statsummery=arrayfun(@(x) p_values(x)*ones(size(clusters{x})),1:length(clusters),'UniformOutput',0);
            else 
                statsummery=[];
            end
        end
        function PerformBonferroniCorrected_tTest(obj,trial_group_1, trial_group_2,dependent_samples,p_threshold,two_sided,varargin)
            
        end
        function AllTrlPerf=PerformBinomialTestAllTrlPerf(obj,AllTrlPerf) % performs binomial statistical test on Trl performances 
             Fields={'FromSwitch','ToSwitch'};
                ChanceLevel=0.5;
             for f=Fields
                 for Rule=1:3
                     ThisTrlCnt=AllTrlPerf.(f{1})(Rule).AllTrlCnt;
                     AllTrlsCnt=nan*ones(size(ThisTrlCnt));
                     AllTrlsCnt(~isnan(ThisTrlCnt))=obj.NavgPerf;
                     SumCorrTrls=nansum(ThisTrlCnt,2);
                     SumAllTrl=nansum(AllTrlsCnt,2);
                     % run stat test
                     AllTrlPerf.(f{1})(Rule).BinomialStatTest=obj.BinomialStatTest(SumAllTrl,SumCorrTrls,ChanceLevel);
                 end
             end

        end
        function pval=BinomialStatTest(obj,Ntot,Ncorrect,ChancelLevel)
             pval=arrayfun(@(x) 1-binocdf(Ncorrect(x),Ntot(x),ChancelLevel),1:length(Ntot));
        end
        function [pval,Perf1,Perf2]=ChiSquareStatTestAllTrlPerf(obj,successes_distribution1,successes_distribution2,NumTrsAvg) % compares two distributions duing learning usig chi sqaure test
            % Example data for two binomial distributions
%             successes_distribution1 = 20; % Number of successes for distribution 1
%             trials_distribution1 = 30;    % Total number of trials for distribution 1
% 
%             successes_distribution2 = 25; % Number of successes for distribution 2
%             trials_distribution2 = 30;    % Total number of trials for distribution 2 %        
%            Perform the chi-squared test
%            [h, pval, stats] = chi2gof(observed_data, 'ctrs', [0, 1], 'frequency', 'off');
            for Ntrl=1:size(successes_distribution1,2)
                SumSucess1=sum(successes_distribution1(:,Ntrl),'omitnan');
                SumTot1=sum(~isnan(successes_distribution1(:,Ntrl)))*NumTrsAvg;
                
                SumSucess2=sum(successes_distribution2(:,Ntrl),'omitnan');
                SumTot2=sum(~isnan(successes_distribution2(:,Ntrl)))*NumTrsAvg;
                
                % Combine the data into a matrix
                observed_data = [SumSucess1, SumTot1-SumSucess1;
                    SumSucess2, SumTot2-SumSucess2];               

                % Perform chi-squared test
                [pval(Ntrl)] = chi2test(observed_data);
                Perf1(Ntrl)=SumSucess1/SumTot1;
                Perf2(Ntrl)=SumSucess2/SumTot2;
            end
        end
    end
end

