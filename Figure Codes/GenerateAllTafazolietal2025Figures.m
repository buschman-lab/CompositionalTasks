function GenerateAllTafazolietal2025Figures(DataPath,FigSavePath)
%GENERATEALLTAFAZOLIETAL2025FIGURES Summary of this function goes her

%DataPath is the path to data, download them from 
% https://figshare.com/articles/dataset/Data_for_Tafazoli_et_al_2025_Building_Compositional_Tasks_with_Shared_Neural_Subspaces/30276238
% and copy them into root\Figure Data folder
%FigSavePath is the path you want the resulting figures to be saved
%   Detailed explanation goes here
Figure1(DataPath,FigSavePath)
Figure2(DataPath,FigSavePath)
Figure3(DataPath,FigSavePath)
Figure4(DataPath,FigSavePath)
Figure5(DataPath,FigSavePath)
SuppFigures(DataPath,FigSavePath)
end

