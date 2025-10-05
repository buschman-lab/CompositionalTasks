% simulate linear and non linear encoding 
% define a boundary that is combination of linear and non linear 
Bondary=@(x) x.^2;
X=-100:0.1:100;
nX=length(X);
% generate data for this boundary with some noise
data=Bondary(X)+100*randn(1,nX);
Labels=[data>=Bondary(X)]';
%data(Labels==1)=data(Labels==1)+50;
%data(Labels==0)=data(Labels==0)-50;
% plot the data and the boundary
figure
subplot(221);
hold on
plot(X,Bondary(X),'Color','b','LineWidth',5,'LineStyle','--');
plot(X,data,'.');
plot(X(Labels),data(Labels),'.','color','y');
legend({'boundry','Label 1','Label 2'})
axis square

% do the classification task 
Predictors=[X' data'];
c = cvpartition(nX,'Holdout',0.2);
TrainInds=randsample(1:nX,c.TrainSize);
TestInds=setdiff(1:nX,TrainInds);
PredictorsTrain=Predictors(TrainInds,:);PredictorsTest=Predictors(TestInds,:);
LabelsTrain=Labels(TrainInds,:);LabelsTest=Labels(TestInds,:);
clear Accuracy AccuracyRBF
for polyorder=1:3   
    PolyClassifier=fitcsvm(PredictorsTrain,LabelsTrain,'KernelFunction','polynomial', 'PolynomialOrder', polyorder);
    [PredictLabel]=predict(PolyClassifier,PredictorsTrain);
    Accuracy(polyorder)=sum(PredictLabel==LabelsTrain)/c.TrainSize;
end

RBFClassifier=fitcsvm(PredictorsTrain,LabelsTrain,'KernelFunction','RBF');
[PredictLabel]=predict(RBFClassifier,PredictorsTrain);
AccuracyRBF=sum(PredictLabel==LabelsTrain)/c.TrainSize;

 
% % subtract the residuals
% residuals=double(Labels)-LinearScorePredict(:,1);
% LabelsNonlin=residuals>=median(residuals);
% 
% subplot(222);hold on
% plot(X,residuals,'b.');
% plot(X(LabelsNonlin),residuals(LabelsNonlin),'.','color','y');
% 
% % LinearClassifierResiduals=fitcsvm(Predictors,residuals,'KernelFunction','rbf');
% % [PredictLabelNonLin]=predict(LinearClassifier,Predictors);
% % Accuracy=sum(PredictLabelNonLin==Labels)/nX;
% 
% NonLinFit=fitrgp(Predictors,residuals);
% residualspredit=predict(NonLinFit,Predictors);
% subplot(223);hold on
% plot(X,residualspredit,'b.');
% 
% LinFit=fitrlinear(Predictors,residuals);
% residualsLinpredit=predict(LinFit,Predictors);
% subplot(224);hold on
% plot(X,residualsLinpredit,'b.');
