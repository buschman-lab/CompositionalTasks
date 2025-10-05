ABuff=rand(1,500);
BBuff=rand(1,500);
n=0;
figure
for M=[1 5 20]
A=movmean(ABuff,M,2);
B=movmean(BBuff,M,2);

[pccf, lags,a_ma,p_ma]=ManData.PartialCrossCorr(A,B,10);
[a,p]=corr(A',B');
subplot(3,2,1+(n)*2);
scatter(A,B);
title(['K=' num2str(M) ' Corr before correction:a=' num2str(a,2), ' p=' num2str(p,2)])

subplot(3,2,2+(n)*2);

stem(lags,pccf,"filled");
title(['K=' num2str(M) ' Corr after correction:a=' num2str(a_ma,2), ' p=' num2str(p_ma,2)])


n=n+1;
end



ABuff=rand(1,500);
BBuff=ABuff+0.5*rand(1,500);
n=0;
figure
for M=[1 5 20]
A=movmean(ABuff,M,2);
B=movmean(BBuff,M,2);

[pccf, lags,a_ma,p_ma]=ManData.PartialCrossCorr(A,B,10);
[a,p]=corr(A',B');
subplot(3,2,1+(n)*2);
scatter(A,B);
title(['K=' num2str(M) 'Corr before correction:a=' num2str(a,2), ' p=' num2str(p,2)])

subplot(3,2,2+(n)*2);

stem(lags,pccf,"filled");
title(['K=' num2str(M) ' Corr after correction:a=' num2str(a_ma,2), ' p=' num2str(p_ma,2)])


n=n+1;
end

