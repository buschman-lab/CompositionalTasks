% finds temporal croscorrelation between two matrixes
% takes the whole matrix elemets and uses corr 
function [D,max_raw,max_col] = TemporalCrossCorrVer2(ZI,ZJ)
    global AnalysisData
    SizeW=AnalysisData.SizeW;
    L=SizeW(2);
    ZI=reshape(ZI,SizeW);
    ZI_pad=padarray(ZI,[0 L]);
    PlotDetails=0; % plot to make sure everything is correct
    
    for zj=1:size(ZJ,1)
        ZJ_ind=reshape(ZJ(zj,:),SizeW);        
        ZJ_pad=padarray(ZJ_ind,[0 L]);
        
        MatRef=ZJ_pad;
        MatComp=ZI_pad;
       
        crr=xcorr2(MatRef,MatComp);
        [ssr,snd] = max(crr(:));        
        [max_raw,max_col] = ind2sub(size(crr),snd);
       
        if PlotDetails
            subplot(141)
            helperCWTTimeFreqPlot(MatRef,1:size(MatRef,2),AnalysisData.f,'justplot1',['PSD'],'Time','Frequency(Hz)',0)
            subplot(142)
            helperCWTTimeFreqPlot(MatComp,1:size(MatComp,2),AnalysisData.f,'justplot1',['PSD'],'Time','Frequency(Hz)',0)
            subplot(143)
            helperCWTTimeFreqPlot(crr,1:size(crr,2),1:size(crr,1),'justplot1',['PSD'],'Time','Frequency(Hz)',0)


            subplot(144)
            plot(crr(:))
            title('Cross-Correlation')
            hold on
            plot(snd,ssr,'or')
            hold off
            text(snd*1.05,ssr,'Maximum')
        end
        D(zj)=ssr;
    end

    D=1./D;

    
end

function [out,raw,col]=CalCrosCor(X,Y,SizeW)
    bb=arrayfun(@(x) xcorr(X(x,:),Y(x,:)),1:size(Y,1),'UniformOutput',0);
    bb=bb';
    CrossCorr=cell2mat(bb);
    [out,ind]=max(CrossCorr(:)); 
    [raw,col]=ind2sub(SizeW,ind);
end

