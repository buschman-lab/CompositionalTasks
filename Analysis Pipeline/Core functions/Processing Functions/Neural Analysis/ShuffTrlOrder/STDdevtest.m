
STD_obsv_Shape=squeeze(std(EncodingDist_Obsv{1}.ShapeDistRedReSampAvg,0,1));

zShapeDistAvg_obsv=EncodingDist_Obsv{1}.ShapeDistAvg./STD_obsv_Shape;

for i=1:500
    STD_sh_Shape=squeeze(std(EncodingDist_Shuff{i}.ShapeDistRedReSampAvg,0,1));
    zShapeDistAvg_sh(:,:,i)=EncodingDist_Shuff{i}.ShapeDistAvg./STD_sh_Shape;
end

Observed=ManData.SmoothData(zShapeDistAvg_obsv,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);
Shuffle=ManData.SmoothData(zShapeDistAvg_sh,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);


NStg=16;
Shuffle=permute(Shuffle,[3 1 2]);
TrlRng=[1:NStg]';
col=copper(16);
figure;hold on;arrayfun(@(x) plot(Observed(x,:),'color',col(x,:)),1:16);
MeanShuff=squeeze(mean(Shuffle,1));
figure;hold on;arrayfun(@(x) plot(MeanShuff(x,:),'color',col(x,:)),1:16);
