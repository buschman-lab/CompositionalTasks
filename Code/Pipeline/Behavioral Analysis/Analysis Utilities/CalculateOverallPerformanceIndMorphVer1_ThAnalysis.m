% threshold analaysis where we simulate the monkey behavior using different
% thresold s for number of trials that the performance goes over 0.65
function  [fh,NewTrialInterval]=CalculateOverallPerformanceIndMorphVer1_ThAnalysis(i,NSwitch,Trials,Nback,TrialInterval,ThNumTrl,ThPerf,fh)
global   N_Ind   LastTrialMat_Ind LastTrial_Ind 
                if isempty(fh); fh = figure(20000); end
                figure(fh)
                set(0, 'currentfigure', fh);
                hold on;
                LookMorph={[0],[30],[70],[100],[130],[170],[30,70,130,170]};
                LookMorphLeg={'0','30','70','100','130','170','mm'};
                NMorph=length(LookMorph);
                FirstThTrl=[];
                Col=distinguishable_colors(NMorph);
                subplot(3,1,i);hold on
                
                if i==1
                    Feature='ObjectMorphLevel';
                else
                    Feature='ColorMorphLevel';
                end
                 
                k=1;
                AttemptedTrials=[Trials(TrialInterval).StopCondition] >=-1;
                for ML=1:length(LookMorph)
                    INDOBJ_AllTrials=cell2mat(arrayfun(@(x) double(([Trials(TrialInterval).StopCondition] == 1 | [Trials(TrialInterval).StopCondition] == -1) & [Trials(TrialInterval).(Feature)] == x)',LookMorph{ML},'UniformOutput',0));
                    INDOBJ_AllTrials=sum(INDOBJ_AllTrials,2)';
                 %   INDOBJ_AllTrials=double(([Trials(TrialInterval).StopCondition] == 1 | [Trials(TrialInterval).StopCondition] == -1) & [Trials(TrialInterval).(Feature)] == ML);
                    INDOBJ_AllTrials(INDOBJ_AllTrials==0)=NaN;
                    INDOBJ_CorrectTrials=([Trials(TrialInterval).StopCondition] == 1);
                    CorrectMat=INDOBJ_AllTrials.*INDOBJ_CorrectTrials;
                    CorrectMat=CorrectMat(AttemptedTrials);
                    MeanCorrectMat{k}=movmean(CorrectMat,Nback,2,'omitnan','Endpoints','discard');
                    MeanCorrectMatTh(k,:)=MeanCorrectMat{k}>=ThPerf;
                    k=k+1;
                end
                
                
                    % determine the trial we pass the threshold for al of the
                    % target morphlevels 
                    SumCorrectMatTh=sum(MeanCorrectMatTh([1 4 7],:),1);
                    MinNumTrl=ThNumTrl-Nback;
                    SumCorrectMatThTrl=find(SumCorrectMatTh==3);
                    if sum(AttemptedTrials)>=ThNumTrl && ~isempty(SumCorrectMatThTrl)
                        FirstThTrl=SumCorrectMatThTrl(find(SumCorrectMatThTrl>=MinNumTrl,1,'first'));
                        if ~isempty(FirstThTrl)
                            TargetTrial=FirstThTrl+Nback-1; % trial that we have passed the threshold
                            TargetTrialinBlk=find(cumsum(AttemptedTrials)==TargetTrial,1,'first');
                            NewTrialInterval=TrialInterval(1):TrialInterval(TargetTrialinBlk);
                        else
                            NewTrialInterval=TrialInterval;
                        end                        
                    else
                        NewTrialInterval=TrialInterval;
                    end
                    
                LastTrial_Ind(i)=LastTrial_Ind(i)+length(CorrectMat);
                LastTrialMat_Ind{i}=[LastTrialMat_Ind{i} LastTrial_Ind(i)];
                PlotInterval=LastTrialMat_Ind{i}(N_Ind(i)):(LastTrialMat_Ind{i}(N_Ind(i)+1)-Nback+1);
                if ~isempty(PlotInterval) && length(PlotInterval)>0
                     arrayfun(@(x) plot(PlotInterval(1:length(MeanCorrectMat{x})),MeanCorrectMat{x},'linewidth',2,'color',Col(x,:)),1:NMorph);
                  %  LegendTxt=cellfun(@(x) x,LookMorphLeg,'UniformOutput',0);
                    text(PlotInterval(1),1,['S:' num2str(NSwitch)],'FontSize',15);
                    title(['Rule:' num2str(i)])
                    ylim([0 1.2])    
                    legend(LookMorphLeg,'Location','northeastoutside')           
                end
                v=axis;
                plot([v(1) v(2)],[ThPerf ThPerf],'--k') 
                plot([PlotInterval(end) PlotInterval(end)],[v(3) v(4)],'--k');
                % show target trial 
                if ~isempty(FirstThTrl)
                    plot([PlotInterval(1)+FirstThTrl PlotInterval(1)+FirstThTrl],[v(3) v(4)],'--r');
                end
                N_Ind(i)=N_Ind(i)+1;
end