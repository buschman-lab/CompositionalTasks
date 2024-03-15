% finds temporal croscorrelation between two matrixes
% uses 1D correlation
% works with PArloop
function [D,Shift,NewMotif_corct,crr,ShiftRng] = TemporalCrossCorr(ZI,ZJ,varargin)
global AnalysisData

if ~isempty(varargin)
    SizeW=varargin{1};
else
    SizeW=AnalysisData.SizeW;
end

ZI=reshape(ZI,SizeW);
ZI_pad=padarray(ZI,[0 SizeW(2)]);
PlotDetails=0; % plot to make sure everything is correct
SizePadMat=size(ZI_pad);  % size of final padded matrix
ShiftRng=-SizeW(2):SizeW(2);
% put all of the shifted values into a matrix( make correlation
% very fast)
Shift_ZIpad=arrayfun(@(x) reshape(circshift(ZI_pad,x,2),[SizePadMat(1)*SizePadMat(2) 1]),ShiftRng,'Uniformoutput',0);
Shift_ZIpad=cell2mat(Shift_ZIpad);

for zj=1:size(ZJ,1)
    
    ZJ_ind=reshape(ZJ(zj,:),SizeW);
    ZJ_pad=padarray(ZJ_ind,[0 SizeW(2)]);
    
    % calculte correlation
    crr=corr(ZJ_pad(:),Shift_ZIpad);
    [ssr(zj),snd] = max(crr(:));
    Shift(zj)=-1*ShiftRng(snd); % how many to shift
    
    % shift this motif to match template motif
    NewMotif=circshift(ZJ_pad,Shift(zj),2);
    NewMotif_corct(zj,:)=reshape(NewMotif,[1 size(NewMotif,1)*size(NewMotif,2)]);
    
    
    if PlotDetails
        subplot(141)
        helperCWTTimeFreqPlot(ZJ_pad,1:size(ZJ_pad,2),AnalysisData.cwt_f,'justplot1',['TFR this motifs'],'Time','Frequency(Hz)',0)
        subplot(142)
        helperCWTTimeFreqPlot(ZI_pad,1:size(ZI_pad,2),AnalysisData.cwt_f,'justplot1',['TFR temp motif'],'Time','Frequency(Hz)',0)
        subplot(143)
        plot(ShiftRng,crr)
        title('Cross-Correlation')
        hold on
        plot(ShiftRng(snd),ssr(zj),'or')
        hold off
        text(ShiftRng(snd)*1.05,ssr(zj),'Maximum')
        subplot(144)
        helperCWTTimeFreqPlot(NewMotif,1:size(NewMotif,2),AnalysisData.cwt_f,'justplot1',['PSD'],'Time','Frequency(Hz)',0)
        pause
    end
end

D=[ssr]';


end
