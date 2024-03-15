for i=1:6
     classificationSVM=fitclinear(trainingPredictors,trainingResponse,'Learner','svm','Lambda',Lambda(i),...
                    'Regularization','ridge','Solver',solver{opts.solver});
                Observed=ProjectData2SVMhyperplane(obj,classificationSVM,validationPredictors,validationResponse,ClassifierOpts);
                Accu(i)=Observed.Accurary;
end
          

CVMdl = fitclinear(trainingPredictors,trainingResponse,'KFold',5,...
    'Learner','svm','Solver','sparsa','Regularization','lasso',...
    'Lambda',Lambda,'GradientTolerance',1e-8)
numCLModels = numel(CVMdl.Trained)
ce = kfoldLoss(CVMdl);
Mdl = fitclinear(trainingPredictors,trainingResponse,...
    'Learner','svm','Solver','sparsa','Regularization','lasso',...
    'Lambda',Lambda,'GradientTolerance',1e-8);
numNZCoeff = sum(Mdl.Beta~=0);
figure;
[h,hL1,hL2] = plotyy(log10(Lambda),log10(ce),...
    log10(Lambda),log10(numNZCoeff)); 
hL1.Marker = 'o';
hL2.Marker = 'o';
ylabel(h(1),'log_{10} classification error')
ylabel(h(2),'log_{10} nonzero-coefficient frequency')
xlabel('log_{10} Lambda')
title('Test-Sample Statistics')
hold off



[Mdl,FitInfo,HyperparameterOptimizationResults]  = fitclinear(trainingPredictors,trainingResponse,'Solver','sparsa','Verbose',0, ...
    'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'),'Learner','logistic')
