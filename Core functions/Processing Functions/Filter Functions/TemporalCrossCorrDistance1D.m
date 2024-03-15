% finds temporal croscorrelation between two matrixes
% uses 1D correlation
function [D,Shift] = TemporalCrossCorrDistance1D(ZI,ZJ,SizeW,PadRef)
    % ZI: Reference motif
    %ZJ: collection of other motifs
    %PadRef: 0 if we dont want to pad and reshape ZI
    if PadRef % prepare and Pad Reference
        ZI=reshape(ZI,SizeW);
        ZI_pad=padarray(ZI,[0 SizeW(2)]);
    else
        ZI_pad=ZI;
    end
   %PlotDetails=0; % plot to make sure everything is correct
    SizePadMat=size(ZI_pad);  % size of final padded matrix
    ShiftRng=-SizeW(2):SizeW(2);
    
    for zj=1:size(ZJ,1)
        
        ZJ_ind=reshape(ZJ(zj,:),SizeW);        
        ZJ_pad=padarray(ZJ_ind,[0 SizeW(2)]);        
        % put all of the shifted values into a matrix( make correlation
        % very fast)         
        Shift_ZIpad=arrayfun(@(x) reshape(circshift(ZI_pad,x,2),[SizePadMat(1)*SizePadMat(2) 1]),ShiftRng,'Uniformoutput',0);
        Shift_ZIpad=cell2mat(Shift_ZIpad);
        % calculte correlation 
        crr=corr(ZJ_pad(:),Shift_ZIpad);       
        [ssr(zj),snd] = max(crr(:));             
        Shift(zj)=ShiftRng(snd); % how many to shift 
       
%         if PlotDetails
%             subplot(131)
%             helperCWTTimeFreqPlot(ZJ_pad,1:size(ZJ_pad,2),AnalysisData.f,'justplot1',['PSD'],'Time','Frequency(Hz)',0)
%             subplot(132)
%             helperCWTTimeFreqPlot(ZI_pad,1:size(ZI_pad,2),AnalysisData.f,'justplot1',['PSD'],'Time','Frequency(Hz)',0)
%             subplot(133)
%             plot(ShiftRng,crr)
%             title('Cross-Correlation')
%             hold on
%             plot(snd,ssr,'or')
%             hold off
%             text(snd*1.05,ssr,'Maximum')
%         end
    end

    D=[ssr]';

    
end

function [out,raw,col]=CalCrosCor(X,Y,SizeW)
    bb=arrayfun(@(x) xcorr(X(x,:),Y(x,:)),1:size(Y,1),'UniformOutput',0);
    bb=bb';
    CrossCorr=cell2mat(bb);
    [out,ind]=max(CrossCorr(:)); 
    [raw,col]=ind2sub(SizeW,ind);
end

