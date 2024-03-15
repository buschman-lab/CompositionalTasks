% finds temporal croscorrelation between two matrixes
%using dynamic time warping
%
function [D,Shift] = TemporalCrossCorrDTW(ZI,ZJ)
global AnalysisData

SizeW=AnalysisData.SizeW;
ZI=reshape(ZI,SizeW);

for zj=1:size(ZJ,1)    
    ZJ_ind=reshape(ZJ(zj,:),SizeW);
    % take dtw of Zi and ZJ now
    [dist,ix,iy]=dtw(ZI,ZJ_ind);
    ZI_warp=ZI(:,ix);
    ZJ_warp=ZJ_ind(:,iy);    
    % calculte correlation
    [ssr(zj)] = corr2(ZI_warp,ZJ_warp);
    Shift(zj)=dist; % how many to shift
    
end
D=[ssr]';
end

