function [NormInd_Obs,NormInd_Shuf]=ManipulateMeticVal(MetricVal,MetricVal_Shuff,varargin)
global AnalysisOpts

ManData=ManipulateData;

if ~isempty(varargin)
    Method=varargin{1};
else
    Method='Identical';
end

Sq=@(x) squeeze(x);
NrepShuff=size(MetricVal_Shuff,1);
NTim=size(MetricVal,3);
fprintf('\n Manipulating observed and shuffle with %s method',Method)
switch Method
    case 'BaseLine'
        AvgTime=AnalysisOpts.Time<=0.1 & AnalysisOpts.Time>=-0.1;
        NormInd_Obs=MetricVal-repmat(mean(MetricVal(:,:,AvgTime),3),[1 1 NTim]);
        if ~isempty(MetricVal_Shuff)
            NormInd_Shuf=MetricVal_Shuff-repmat(mean(MetricVal_Shuff(:,:,AvgTime),3),[1 1 NTim]);
        else
            NormInd_Shuf=[];
        end
    case 'Identical'
        NormInd_Obs=MetricVal;
        NormInd_Shuf=MetricVal_Shuff;

    case 'NormalizationIndex'
        % subtract 50% perfromance first
       % MetricVal=MetricVal-0.5;MetricVal_Shuff=MetricVal_Shuff-0.5;
                 
        [NormInd_Obs]=ManData.CalNormalizationIndex(Sq(MetricVal(1,2,:)),Sq(MetricVal(1,1,:)));

        NormInd_Shuf=cell2mat(arrayfun(@(x) ManData.CalNormalizationIndex(Sq(MetricVal_Shuff(x,2,:)),Sq(MetricVal_Shuff(x,1,:))),...
            1:NrepShuff,'UniformOutput',0));

        NormInd_Obs=permute(NormInd_Obs,[3 2 1]);
        NormInd_Shuf=permute(NormInd_Shuf,[2 3 1]);

    case 'Subtract50Abs' % subtract 50 and take the absolute value
      %  MetricVal=abs(MetricVal-0.5);MetricVal_Shuff=abs(MetricVal_Shuff-0.5);
        NormInd_Obs=MetricVal(:,2,:)-MetricVal(:,1,:);
        NormInd_Shuf=MetricVal_Shuff(:,2,:)-MetricVal_Shuff(:,1,:);      
end

end
