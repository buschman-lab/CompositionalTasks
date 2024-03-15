function W=TrainSOM(inputdata,varargin) 
% builds an SOM and trians it
% more to come 
% inputdata: X by N matirix denotes X samples with N dimensionas 
% Sina Tafazoli Ver 0: 23 July 2020

% set up options 
opts.dim1   = 10;   % number of first dimensions
opts.dim2   = 10;   % number of second dimensions 
opts.N      = size(inputdata,2);  % num of input dimensions 
opts.X      = size(inputdata,1);   % number of samples
opts.sigma  = 1;  % sigma of neigborhood fucntion ~ how many neighbors are being changed 
opts.alpha  = 0.03; % scalare factor that defines size of the wight correction 
opts.Niter  = 400;   % how many iterations
opts.verbose= 1;
opts        = ParseOptionalInputs(opts,varargin);
Col=colormap(jet(opts.Niter));

% prepare sheet of neurons 
Neu=zeros(opts.dim1,opts.dim2); % this the matrix for activity of each neuron 
NNeu=opts.dim1*opts.dim2;  % number of neurons
NeuSiz=[opts.dim1 opts.dim2]; 
[IND2SUB(1,:),IND2SUB(2,:)]=arrayfun(@(x) ind2sub(NeuSiz,x),1:NNeu,'UniformOutput',1); % jsut calculte all of our ind2subs in advance
% prepre the weight matrix 
% initialize it with random values
W=rand(NNeu,opts.N);

% setup neigborhood function 
hc=@(x) gauss2D(Neu,opts.sigma,x);
EucDist=@(x,y) EuclideanDist(x,y);

% start the training loop
W_all=[]; % all of our Ws just to keep track of them 

for iter=1:opts.Niter
    if opts.verbose;fprintf('\nTraining epoc %i',iter);end
   
    alpha(iter)=exp(-opts.alpha*iter);
    W_all(:,:,iter)=W; % keep track of evolution of weights over epochs
    % draw a random sample with replacement
    sampind=randsample(1:opts.X,1);
    sample=inputdata(sampind,:);
    
    % calculte Euclidean distance of this sample to all of the weights
    dist=arrayfun(@(n) EucDist(W(n,:),sample),1:NNeu);
    [mindist,mindist_ind]=min(dist); % find minimum distance
    i=IND2SUB(1,mindist_ind);j= IND2SUB(2,mindist_ind) ;
    
    % now update all of the weights
    neighborhood=hc([i j]);  % calculte all of the neighbors of the max node
    
    % update the new weights
    W=[cell2mat(arrayfun(@(x) transpose(W(x,:)+alpha(iter)*neighborhood(IND2SUB(1,x),IND2SUB(1,x))*(sample- W(x,:))),1:NNeu,'UniformOutput',0))]';
   
    if opts.verbose
        figure(1) 
        if iter==1;subplot(231);imagesc(neighborhood);axis square;colorbar;end
        subplot(232)
        plot(inputdata(:,1),inputdata(:,2),'.','markersize',7)
        axis square 
        subplot(233);plot(1:iter,alpha);drawnow
        subplot(234)
        plot(W(:,1),W(:,2),'.','color',Col(iter,:),'markersize',7)
        drawnow
        axis square    
    end
   % pause
end

   
    
    
end
    
    
   






