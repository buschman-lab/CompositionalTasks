function CompareCWTRules(CWTRule1,CWTRule2,TrueTime,cwt_f,R2)

   if R2==0
       meanCWTRule1=mean(CWTRule1(:,:,end-99),3); 
       meanCWTRule2=mean(CWTRule2(:,:,end-99),3);
   else
       meanCWTRule1=mean(CWTRule1(:,:,end-49),3); 
       meanCWTRule2=mean(CWTRule2(:,:,end-49),3);
   end
    
    diffmean=abs(meanCWTRule1-meanCWTRule2);
    colormap(jet)
    subplot(2,2,1)
    helperCWTTimeFreqPlot(meanCWTRule1,TrueTime,cwt_f,'power','PSD Rule1','Time from Stim Onset(ms)','Frequency(Hz)')

    subplot(2,2,2)
    helperCWTTimeFreqPlot(meanCWTRule2,TrueTime,cwt_f,'power','PSD Rule3','Time from Stim Onset(ms)','Frequency(Hz)')

    
    subplot(2,2,3)
    helperCWTTimeFreqPlot(diffmean,TrueTime,cwt_f,'power','|PSD Rule1-Rule3|','Time from Stim Onset(ms)','Frequency(Hz)')

    subplot(2,2,4)
    helperCWTTimeFreqPlot(diffmean,TrueTime,cwt_f,'phase','Phase Rule1-Rule3','Time from Stim Onset(ms)','Frequency(Hz)')

    