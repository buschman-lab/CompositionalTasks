function ChangeSpkCntStartFieldName(input)
global AnalysisOpts
if strcmp(input,'p')
AnalysisOpts.SpkCntStartFieldName='SAMPLE_ON';
else
    AnalysisOpts.SpkCntStartFieldName='SACCADE_START';
end

end