% finds temporal croscorrelation between two matrixes
function [D,raw,col] = TemporalCrossCorrVer1(ZI,ZJ)
    global AnalysisData
    SizeW=AnalysisData.SizeW;
    ZI_rs=reshape(ZI,SizeW);
    D=arrayfun(@(x) CalCrosCor(ZI_rs,reshape(ZJ(x,:),SizeW),SizeW),1:size(ZJ,1),'UniformOutput',0);
    D=1./transpose(cell2mat(D));
end

function [out,raw,col]=CalCrosCor(X,Y,SizeW)
    bb=arrayfun(@(x) xcorr(X(x,:),Y(x,:)),1:size(Y,1),'UniformOutput',0);
    bb=bb';
    CrossCorr=cell2mat(bb);
    [out,ind]=max(CrossCorr(:)); 
    [raw,col]=ind2sub(SizeW,ind);
end

