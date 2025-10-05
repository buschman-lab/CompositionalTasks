function RenameFiles(NNeu)

global AnalysisOpts

if NNeu==116;return;end

Path=[AnalysisOpts.DataSavePath 'Classifier' filesep 'ALL' filesep];
for i=1:250
    SourceFile=sprintf('%sClassifier_ALL_ALL_NeuRS%i_C%i_3D Color Category Response Xgen BalInCongV5_AreaComp_F_MS_PFC_SAMPLE_ON_AllTrials_100.mat',Path,NNeu,i);
    TargetFile=sprintf('%sClassifier_ALL_ALL_NeuRS%i_C%i_3D Color Category Response Xgen BalInCongV5_AreaComp_MS_PFC_SAMPLE_ON_AllTrials_100.mat',Path,NNeu,i);
    try
        movefile(SourceFile,TargetFile);
        fprintf('\nRenamed %i',i)
    catch
        fprintf('\nErr renaming file %i ---%s',i,SourceFile)
    end
end

end