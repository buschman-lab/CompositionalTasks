

[beta_gls_obs,p_values_obs]=arrayfun(@(x) GeneralizedLeastSquaresTrendStatTest(TrlRng,ObservedMvMean(:,x)),1:116,'UniformOutput',0);
beta_gls_obs=cell2mat(beta_gls_obs);
p_values_obs=cell2mat(p_values_obs);

figure;
subplot(121)
plot(beta_gls_obs')

subplot(122)
plot(p_values_obs')


[beta_fgls_obs]=arrayfun(@(x) fgls(TrlRng,Observed(:,x)),1:116,'UniformOutput',0);

for i=1:50%NrepShuff
    fprintf('\n Cal fgls Shuff %i');
    beta_fgls_sh(:,:,i)=cell2mat(arrayfun(@(x) fgls(TrlRng,squeeze(Shuffle(i,:,x))'),1:NTim,'UniformOutput',0));
end

model = arima('ARLags', 10, 'D', 0, 'Constant', 0); % AR(1) model

[beta_arima_obs,p_arima_obs]=arrayfun(@(x) ARIMAtrendcorr(log(Observed(:,x)),TrlRng,model),1:116,'UniformOutput',0);
figure;
subplot(121)
plot(Time,cell2mat(beta_arima_obs))

subplot(122)
plot(Time,cell2mat(p_arima_obs))


