%% plot Arima fit to compression index results 
ManData=ManipulateData;
Ar=[1 1:9];
Poly=[5 zeros(1,9)];
Path='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Analysis Pipeline\Input Output Data\Code Test Data\';

for ar=1:length(Ar)
    % load file 
    FileName=sprintf('ArimaTestShuffLag%i',Ar(ar));
    if Poly(ar)
        FileName=[FileName 'Poly' num2str(Poly(ar))];
    end
     
    load([Path FileName]);
    AnalysisOpts.Time=-0.15:0.01:0.61; % this has to be modified
    % preprocess the resutls 
    ManData.ClusterMassCorrection_permutationTwoTail(permute(cell2mat(ArShuffle),[3 1 2]),permute(cell2mat(ArObserved),[1 3 2]),0.2,1,'ShowClustCorrectionPlot',1);
    sgtitle(sprintf('ARIMA+Spearman Lag %i, Poly %i',Ar(ar),Poly(ar)));
    
    figure(100)
    subplot(2,5,ar)
    A=squeeze(cell2mat(ArShuffle_p))';
    plot(AnalysisOpts.Time,sum(A<0.05,1)/size(ArShuffle_p,3));
    xlabel('Time');
    ylabel('Probablity of p-shuffle<0.05')
    axis square
    title(FileName)

end
