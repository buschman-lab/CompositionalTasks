

function PlotBlockPerf(BlockNums,ConditionTot,BlockIden,TrialOpts)
%%
global PerfFitTot TrialInfoBlk
optsfileds=fieldnames(TrialOpts);
Noptsfileds=length(optsfileds);
NBlockSet=length(BlockNums);
NBlockRule(1)=sum(BlockIden==1);
NBlockRule(3)=sum(BlockIden==3);
NColmn=5;
NRule=[1 1 1];Page=[1 1 1];
for i=1:NBlockSet
     Thisblocks=BlockNums{i};
     ThisblocksRule=ConditionTot(Thisblocks);
     N=[1 1 1];
     ThsBlkUniqRule=BlockIden(i);
     
    for j=Thisblocks
        
      %  BLKNUMS=Thisblocks(ThisblocksRule~=2);
        
        if ConditionTot(j)~=2
            
            if NRule(ThsBlkUniqRule)>NColmn%3*NBlockRule(ThsBlkUniqRule)
                Page(ThsBlkUniqRule)=Page(ThsBlkUniqRule)+1;
                NRule(ThsBlkUniqRule)=1;
            end
            
            figure(100*ThsBlkUniqRule+ThsBlkUniqRule+Page(ThsBlkUniqRule))
            set(gcf,'Name',['Rule ' num2str(ThsBlkUniqRule) ' Page:' num2str(Page(ThsBlkUniqRule)) ])
            x=0:PerfFitTot.Ntrl(j)-1;
            Col=varycolor(length(Thisblocks));
            
            %%% plot the trajectory 
       %     subplot(3,NBlockRule(ThsBlkUniqRule),NRule(ThsBlkUniqRule));hold on
             subplot(3,NColmn,NRule(ThsBlkUniqRule));hold on
           
            plot(x,PerfFitTot.L(j)./(1+exp(-PerfFitTot.k(j)*(x-PerfFitTot.x0(j)))),'Color',Col(N(ThsBlkUniqRule),:),'LineWidth',3)            
            ylim([0 1])
            title(['Rule' num2str(ConditionTot(j))])
            
            %%% plot the paramters 
        %    subplot(3,NBlockRule(ThsBlkUniqRule),NRule(ThsBlkUniqRule)+NBlockRule(ThsBlkUniqRule));hold on          
             subplot(3,NColmn,NRule(ThsBlkUniqRule)+NColmn);hold on          
          
            bar(N(ThsBlkUniqRule),PerfFitTot.L(j),'FaceColor',Col(N(ThsBlkUniqRule),:))
          %  set(gca,'XTick',(Thisblocks(ThisblocksRule~=2)))
            ylim([0 1])
            %%%%%%%%%%%%%%%%%%%%%%%%%  Now write here what paramters were
            %%%%%%%%%%%%%%%%%%%%%%%%%  used for that sets of blocks
            NRule(ThsBlkUniqRule)
            
            if N(ThsBlkUniqRule)==1

                for y=1:Noptsfileds
                    
   %                subplot(3,NBlockRule(ThsBlkUniqRule),NRule(ThsBlkUniqRule)+2*NBlockRule(ThsBlkUniqRule));hold on   
                    subplot(3,NColmn,NRule(ThsBlkUniqRule)+2*NColmn);hold on   

                    ylim([0 6])
                    temp=unique(TrialOpts.(optsfileds{y})(TrialInfoBlk{i}));
                    text(0,y,[optsfileds{y}(1:4) ':' num2str(temp)]);     
                    
                end
                
            end
            
            N(ThsBlkUniqRule)=N(ThsBlkUniqRule)+1;
            
        end                       
    end
    NRule(ThsBlkUniqRule)=NRule(ThsBlkUniqRule)+1;
    
end
              
end