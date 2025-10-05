% checking the validity of modified Mann-Kendall test with a few examples
N=125; % number of samples
ABuff=rand(1,N);
num_stages=1:N;

SmoothSet=[1 5 10 50]; % smoothing factor
figure
n=0;
for M=SmoothSet
    A=movmean(ABuff,M,2);
    B=num_stages; 

    [tau_MMK, z_score, p_MMK, H] =  Modified_MannKendall_test(B', A', 0.001,0.001);
    [a,p]=corr(num_stages',A','type','Kendall');
    subplot(length(SmoothSet),1,1+n);
    plot(B,A);axis tight
    title(['K=' num2str(M) 'Kendall Corr w/o correction:a=' num2str(a,2), ' p=' num2str(p,2),...
        ' Corr with correction:a=' num2str(tau_MMK,2), ' p=' num2str(p_MMK,2)])
    xlabel('Time point');ylabel('value')
    n=n+1;
end

% sweep conditions for Mann-Kendall
NTrls=[16 50 125 500]; % number of samples
SmoothSet=[1 3 5 15 50]; % smoothing width
PvalSet=[0.1 0.05 0.01 0.001 0.0001]; % Pvalue threshold to detect autocorrelation lag
for Rep=1:1000  % repetition
for ntrl=NTrls
    ABuff=rand(1,ntrl);
    for SM=SmoothSet 
        for Pv=PvalSet
            fprintf('\nRunning Rep%i ntrl%i SM%i pval%0.4f',Rep,ntrl,SM,Pv)
            num_stages=1:ntrl;
            l1=find(ntrl==NTrls);l2=find(SM==SmoothSet);l3=Rep;l4=find(Pv==PvalSet);
            if ntrl<=SM
                PvalueSweep(l1,l2,Rep,l4)=NaN;
                AvalueSweep(l1,l2,Rep,l4)=NaN;
                PvalueSweep_NoCorr(l1,l2,Rep,l4)=NaN;
                AvalueSweep_NoCorr(l1,l2,Rep,l4)=NaN;

            else
                A=movmean(ABuff,SM,2);
                % test with Modiefied Man Kendall
                [AvalueSweep(l1,l2,Rep,l4),~,PvalueSweep(l1,l2,Rep,l4)]=Modified_MannKendall_test(num_stages', A', 0.001,Pv);
                % test with Kendall
                [PvalueSweep_NoCorr(l1,l2,Rep,l4),PvalueSweep_NoCorr(l1,l2,Rep,l4)]=corr(num_stages',A','type','Kendall');
            end
        end
    end
end
end

% find the number of times the pvalue has been less than 0.05
PvalueSweepSig=PvalueSweep<=0.001;
PvalueSweepSigAllReps=squeeze(sum(PvalueSweepSig,3))/1000;  % <-- Look at this var for final results

PvalueSweepSig_NoCorr=PvalueSweep_NoCorr<=0.001;
PvalueSweepSigAllReps_NoCorr=squeeze(sum(PvalueSweepSig_NoCorr,3))/1000;



% checking the validity of modified ARima test with a few examples
N=16; % number of samples
ABuff=rand(1,N);
num_stages=1:N;
beta=.1;
ABuff=ABuff+beta*num_stages;
SmoothSet=[1 5 10 ]; % smoothing factor
figure
n=0;
model = arima('ARLags', 3, 'D', 0, 'Constant', 0); % AR(1) model

for M=SmoothSet
    
    A=movmean(ABuff,M,2);
    B=num_stages; 

    [a_arima,p_arima] =  ARIMAtrendcorr(A',B',model,'Kendall');
    [a,p]=corr(num_stages',A','type','Kendall');
    subplot(length(SmoothSet),1,1+n);
    plot(B,A);axis tight
    title(['K=' num2str(M) 'Kendall Corr w/o correction:a=' num2str(a,2), ' p=' num2str(p,2),...
        ' Corr with correction:a=' num2str(a_arima,2), ' p=' num2str(p_arima,2)])
    xlabel('Time point');ylabel('value')
    n=n+1;
end


%%

% sweep conditions for Mann-Kendall
NTrls=16; % number of samples
SmoothSet=[1 5 10]; % smoothing width
LagSet=[1 2 3 4 5 7 10]; % Pvalue threshold to detect autocorrelation lag
PvalueSweep=[];PvalueSweep_NoCorr=[];
BetaSet=[0 3];
for b=1:length(BetaSet)
    for Rep=1:1000  % repetition
        tic
        fprintf('\nRunning Rep%i Beta %i',Rep,BetaSet(b))
        for ntrl=NTrls
            ABuff=rand(1,ntrl)+BetaSet(b)*(1:ntrl);
            for SM=SmoothSet
                for Pv=LagSet
                    model = arima('ARLags', Pv, 'D', 0, 'Constant', 0); % AR(1) model

                    num_stages=1:ntrl;
                    l1=find(ntrl==NTrls);l2=find(SM==SmoothSet);l3=Rep;l4=find(Pv==LagSet);
                    if ntrl<=SM
                        PvalueSweep(l1,l2,Rep,l4,b)=NaN;
                        AvalueSweep(l1,l2,Rep,l4,b)=NaN;
                        PvalueSweep_NoCorr(l1,l2,Rep,l4,b)=NaN;
                        AvalueSweep_NoCorr(l1,l2,Rep,l4,b)=NaN;
                    else
                        A=movmean(ABuff,SM,2);
                        % test with Modiefied Man Kendall
                        [AvalueSweep(l1,l2,Rep,l4,b),PvalueSweep(l1,l2,Rep,l4,b)]=ARIMAtrendcorr(A',num_stages',model);
                        % test with Kendall
                        [PvalueSweep_NoCorr(l1,l2,Rep,l4,b),PvalueSweep_NoCorr(l1,l2,Rep,l4,b)]=corr(num_stages',A','type','Kendall');
                    end
                end
            end
            toc
        end
    end
end

% find the number of times the pvalue has been less than 0.05
PvalueSweepSig=PvalueSweep<=0.05;
PvalueSweepSigAllReps=squeeze(sum(PvalueSweepSig,3))/1000;  % <-- Look at this var for final results

PvalueSweepSig_NoCorr=PvalueSweep_NoCorr<=0.05;
PvalueSweepSigAllReps_NoCorr=squeeze(sum(PvalueSweepSig_NoCorr,3))/1000;



%% Spearman 
% sweep conditions for Mann-Kendall
NTrls=16; % number of samples
SmoothSet=[1 5 10]; % smoothing width
LagSet=[1 2 3 4 5 7 10]; % Pvalue threshold to detect autocorrelation lag
PvalueSweepSpear=[];PvalueSweep_NoCorr=[];
BetaSet=[0 3];
for b=1:length(BetaSet)
    for Rep=1:1000  % repetition
        tic
        fprintf('\nRunning Rep%i Beta %i',Rep,BetaSet(b))
        for ntrl=NTrls
            ABuff=rand(1,ntrl)+BetaSet(b)*(1:ntrl);
            for SM=SmoothSet
                for Pv=LagSet
                    model = arima('ARLags', Pv, 'D', 0, 'Constant', 0); % AR(1) model

                    num_stages=1:ntrl;
                    l1=find(ntrl==NTrls);l2=find(SM==SmoothSet);l3=Rep;l4=find(Pv==LagSet);
                    if ntrl<=SM
                        PvalueSweepSpear(l1,l2,Rep,l4,b)=NaN;
                        AvalueSweepSpear(l1,l2,Rep,l4,b)=NaN;
                        PvalueSweep_NoCorrSpear(l1,l2,Rep,l4,b)=NaN;
                        AvalueSweep_NoCorrSpear(l1,l2,Rep,l4,b)=NaN;
                    else
                        A=movmean(ABuff,SM,2);
                        % test with Modiefied Man Kendall
                        [AvalueSweepSpear(l1,l2,Rep,l4,b),PvalueSweepSpear(l1,l2,Rep,l4,b)]=ARIMAtrendcorr(A',num_stages',model,'Spearman');
                        % test with Kendall
                        [PvalueSweep_NoCorrSpear(l1,l2,Rep,l4,b),PvalueSweep_NoCorrSpear(l1,l2,Rep,l4,b)]=corr(num_stages',A','type','Spearman');
                    end
                end
            end
            toc
        end
    end
end

% find the number of times the pvalue has been less than 0.05
PvalueSweepSpearSig=PvalueSweepSpear<=0.05;
PvalueSweepSpearSigAllReps=squeeze(sum(PvalueSweepSpearSig,3))/1000;  % <-- Look at this var for final results

PvalueSweepSig_NoCorrSpear=PvalueSweep_NoCorrSpear<=0.05;
PvalueSweepSigAllReps_NoCorrSpear=squeeze(sum(PvalueSweepSig_NoCorrSpear,3))/1000;

%% test correlations with shuffle 
ABuff=rand(1,16)+0.1*(1:16);SM=5;
MovABuff=movmean(ABuff,SM,2)';
Obs=corr(num_stages',MovABuff,'type','Spearman');
Shuff=arrayfun(@(x) corr(num_stages',MovABuff(randperm(16)),'type','Spearman'),1:1000);
sum(Shuff>=Obs)/1001


Shuff=arrayfun(@(x) ARIMAtrendcorr(MovABuff(randperm(16)),num_stages',model,'Spearman'),1:1000);




% checking the validity of ManKendallToolbox test with a few examples
N=16; % number of samples
ABuff=rand(1,N);
num_stages=1:N;
beta=0;
ABuff=ABuff+beta*num_stages;
SmoothSet=[1 5 10 ]; % smoothing factor
figure
n=0;
test_data=timetable(datetime(datevec([datenum([2000 01 01]):datenum([2000 01 16])])));
for M=SmoothSet
    
    test_data.Param=movmean(ABuff,M,2)';
    B=num_stages; 

    test_result=MK_tempAggr(test_data,0.001,'PW_method','TFPW_WS','alpha_MK',99,'alpha_ak',95, 'alpha_CL',95, 'alpha_Xhomo',90);
    p_Kendall=test_result.P;a_Kendall=test_result.slope;
    [a,p]=corr(num_stages',test_data.Param,'type','Kendall');
    subplot(length(SmoothSet),1,1+n);
    plot(B,test_data.Param);axis tight
    title(['K=' num2str(M) 'Kendall Corr w/o correction:a=' num2str(a,2), ' p=' num2str(p,2),...
        ' Corr with correction:a=' num2str(a_Kendall,2), ' p=' num2str(p_Kendall,2)])
    xlabel('Time point');ylabel('value')
    n=n+1;
end
