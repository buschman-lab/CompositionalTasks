% finds temporal croscorrelation between two matrixes
% uses 1D correlation
% works with PArloop
function [D,Shift] = TemporalCrossCorrDistance1DParLoopDTW(i,j,Motifs,SizeW)

    ZJ=Motifs(j,:);
    ZI=reshape(Motifs(i,:),SizeW);
    
    PlotDetails=0;
    for zj=1:size(ZJ,1)
        
        ZJ_ind=reshape(ZJ(zj,:),SizeW);        
        % take dtw of Zi and ZJ now
       
        [dist,ix,iy]=dtw(ZI,ZJ_ind);
        ZI_warp=ZI(:,ix);
        ZJ_warp=ZJ_ind(:,iy);
        
        % calculte correlation                
        [ssr(zj)] = corr2(ZI_warp,ZJ_warp);          
        Shift(zj)=dist; % how many to shift 
       
        if PlotDetails
            CorrBefore=corr2(ZI,ZJ_ind);
            subplot(221)
            helperCWTTimeFreqPlot(ZI,1:size(ZI,2),1:size(ZI,1),'justplot1',['Motif 1 ' num2str(CorrBefore,2) ],'Time','Frequency(Hz)',0)
            subplot(222)
            helperCWTTimeFreqPlot(ZI_warp,1:size(ZI_warp,2),1:size(ZI_warp,1),'justplot1',['Motif 1 warped ' num2str(ssr(zj),2)],'Time','Frequency(Hz)',0)
            subplot(223)
            helperCWTTimeFreqPlot(ZJ_ind,1:size(ZJ_ind,2),1:size(ZJ_ind,1),'justplot1',['Motif 2'],'Time','Frequency(Hz)',0)
            subplot(224)
            helperCWTTimeFreqPlot(ZJ_warp,1:size(ZJ_warp,2),1:size(ZJ_warp,1),'justplot1',['Motif 2 warped'],'Time','Frequency(Hz)',0)
            pause
        end
    end

    D=[ssr]';

    
end
