function residuals=GetArimaResiduals(X,Nar) % fits arima model and reports the residuals
series=X(:);

model=arima('ARLags',1:Nar);

fit=estimate(model,series,'Display','off');

residuals=infer(fit,series);

end