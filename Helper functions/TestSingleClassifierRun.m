%tests a single classifier run
perf=[];
for t=1:141
    tic
for i=1
   % SVMmdl=fitcsvm(predictors{1}(:,:,t),response{1});
    %SVMmdl=fitclinear(predictors{1}(:,:,t),response{1});
    SVMmdl=fitclinear(predictors{1}(:,:,t)',response{1},'ObservationsIn','columns','Learner','logistic','Lambda',1/60,...
        'Regularization','ridge','Solver','bfgs','FitBias',true,'PostFitBias',true);

% optstune = struct('Optimizer','bayesopt','ShowPlots',false,'Verbose',0);%,'CVPartition',c);
%     
%     [SVMmdl]  = fitclinear(predictors{1}(:,:,t)',response{1},'ObservationsIn','columns','Solver','bfgs','Verbose',0, ...
%         'Learner','logistic','OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
%         struct('AcquisitionFunctionName','expected-improvement-plus'),...
%         'HyperparameterOptimizationOptions',optstune,'FitBias',true,'PostFitBias',true);
%            
    outsvm=predict(SVMmdl,predictors{2}(:,:,1));
    outsvm(outsvm==4)=2;
    outsvm(outsvm==3)=1;
    perf(i,t)=sum(response{2}==outsvm)/length(outsvm);
end
toc
end
plot(perf,'g')