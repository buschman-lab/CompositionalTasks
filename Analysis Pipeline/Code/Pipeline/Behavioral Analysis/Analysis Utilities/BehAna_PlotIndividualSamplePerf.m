function [Trial]=BehAna_PlotIndividualSamplePerf(bhv,Trials,Condition)
d1=1;% Dimension 1 is shape and 2 is color
Pdist=bhv.Pdist;
Axis=bhv.PlotAxis;L=1:length(Pdist.Dim1Objs);
    figure(1000+Condition)
    subplot(121)
    set(gca,'Xtick',L);
    set(gca,'XtickLabel',{Pdist.Dim1Objs});
    set(gca,'Ytick',L);
    set(gca,'YtickLabel',{Pdist.Dim1Objs});
    colorbar
    xlabel('Shape');ylabel('Color')
    
    subplot(122)
    set(gca,'Xtick',L);
    set(gca,'XtickLabel',{Pdist.Dim1Objs});
    set(gca,'Ytick',L);
    set(gca,'YtickLabel',{Pdist.Dim1Objs});
     xlabel('Shape');ylabel('Color')

  
                            
               
for Dim1ML=Pdist.Dim1Objs
    d2=1;
    for Dim2ML=Pdist.Dim1Objs
            
           %%% Calclute performance of the each of the objects and put it
           %%% in ametrix with all other charachtricstics of that object so
           %%% that for each block we have a summery matrix
           
           ThisObj= [Trials.ObjectMorphLevel]==Dim1ML & [Trials.ColorMorphLevel]==Dim2ML & [Trials.Condition]==Condition ;
           AllTrials=([Trials.StopCondition] == 1 | [Trials.StopCondition] ==-1) & ThisObj; 
           AllTrialsPassed=([Trials.StopCondition] ~= 1 & [Trials.StopCondition] ~=-1) & ThisObj; %trial that animals did fix break or too fast etc
           CorrTrials=([Trials.StopCondition] == 1) & ThisObj;
           %%%%
           Trial.AllTrialsCount(d1,d2)=sum(AllTrials);
           Trial.CorrTrialsCount(d1,d2)=sum(CorrTrials);
           Trial.AllTrialsPassedCount(d1,d2)=sum(AllTrialsPassed);
           Trial.ThisObjPerf(d1,d2)=sum(CorrTrials)/sum(AllTrials);
           
                %%%Ploting objects
               figure(1000+Condition)
               subplot(121)
               hold on
               imagesc(Trial.ThisObjPerf'); 
               set(gca,'YDir','normal')
               subplot(122)
               hold on
               text(d1,d2,[num2str(Trial.AllTrialsCount(d1,d2)) '/' num2str(bhv.TrialSpec{Condition}.Num(d1,d2))])
               xlim([0.5 d1+0.5])
               ylim([0.5 d2+0.5])
       
         
          %%% 
          d2=d2+1;        
    end
    d1=d1+1;
end
 drawnow 

