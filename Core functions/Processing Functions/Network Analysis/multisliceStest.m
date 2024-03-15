%load playdata MaxXcorr
ReshapeRew=[];
Pairs=nchoosek(1:29,2);
ChArea=[1     1     1     1     1     1     1     1     4     4     4     4     4     4  4     4     5     5     5     5     5     5     5     5     5     5     5     5 5];
Area_2look={'PFC','Striatum','IT','FEF','LIP' };
ManData=ManipulateData;
for i=1:16;ReshapeRew{i}=ManData.ReshapeSquareMatrix(Pairs,AvgTimBoxTrls(5).Reward{1, 1}(:,i));
end
ReshapeRewSym=ManData.MakeASymmetricMatrix(ReshapeRew);
jjj=0:2:5;
 iii=1:2:10;
for jj=1:3
    for ii=1:5
        i=iii(ii);
        j=jjj(jj);
        [S,Q]=Community_OrederedMultiSliceNet(ReshapeRewSym,0.25*i,j*0.2);
        subplot(3,5,ii+5*(jj-1));
        imagesc(S);title(['g:' num2str(0.25*i),'o:' num2str(j*0.2)]);
        set(gca,'YDir','normal')
        colormap('colorcube')
        xticks([1:16])
        XTickLblLL=arrayfun(@(x) num2str(x),-10:5,'UniformOutput',0);
        xticklabels(XTickLblLL);
        YTickLblLL=arrayfun(@(x) Area_2look{x},ChArea,'UniformOutput',0);
        yticks([1:29])
        yticklabels(YTickLblLL);
    end
end
figure
[S,Q]=Community_OrederedMultiSliceNet(ReshapeRewSym,0.25*6,2*0.2);
imagesc(S)
set(gca,'YDir','normal')
colormap('colorcube')
xticks([1:16])
XTickLblLL=arrayfun(@(x) num2str(x),-10:5,'UniformOutput',0);
xticklabels(XTickLblLL);
YTickLblLL=arrayfun(@(x) Area_2look{x},ChArea,'UniformOutput',0);
yticks([1:29])

yticklabels(YTickLblLL);
colorbar