function [Clust_ind_Input,clust_ind_In_pheno] = SOMCluster(inputdata,varargin)
%SOMCLUSTER Summary of this function goes here
% clusters input data through a two level mechnism of SOM and Unsupervised
% clustring algorithm 
% Data->Decrease dimensionality with SOM-> plot data with TSNE, Sammons plot, curvilinear component analysis
% non metric multidimensional scalling -> cluster with Hierachical clusting, Phenograph or Kmeans 
% the input data is organized  XbyY where X are samples and Y are
% dimensions
% by Sina Tafazoli July 2020

% set up options 
opts.dim1    = 10;   % number of first dimensions
opts.dim2    = 10;   % number of second dimensions 
opts.verbose = 1;
opts.method  ='matlab'; % methods is matlab or my own code "Matlab" or "code"
opts.clustringMethod='Phenograph';  % can be 'HierarClust' or 'Phenograph'
opts.MaxClust=5;   % max cluster for hierarchical
opts         = ParseOptionalInputs(opts,varargin);
EucDist=@(x,y) EuclideanDist(x,y);

%
fp=fig_params;

% first train an SOM ( choose type of code we want to use to generate this)
if strcmpi(opts.method,'Matlab')
    [net,som_output] = ClusterSelfOrganizingMapMotif(inputdata');
    netW=net.IW{1};
elseif strcmpi(opts.method,'code')
    netW=TrainSOM(inputdata,'dim1',opts.dim1,'dim2',opts.dim2,'verbose',opts.verbose); % train SOM and get the weights   
end

% now apply phenograph on the weight matrix we have 
if strcmpi(opts.clustringMethod,'HierarClust')
      
    [clust_ind_In_pheno ,Linkage_In]  = HierarchicalClustring(inputdata,opts.MaxClust);
    [clust_ind_SOM,Linkage_SOM] = HierarchicalClustring(netW,opts.MaxClust);
else
    [clust_ind_In_pheno ] = Phenograph(inputdata,'Level',1);
    [clust_ind_SOM] = Phenograph(netW,'Level',1);
end
%% now find the cluster corresponding to each input for SOM
%% for each input find the neuron which has the minimum weight to it 
for i=1:size(inputdata,1)
    Dist=arrayfun(@(x) EucDist(inputdata(i,:),netW(x,:)),1:size(netW,1));
    [~,DistMinInd(i)]=min(Dist);
end
Clust_ind_Input=clust_ind_SOM(DistMinInd);

if opts.verbose
    figure(gcf)
    % show input data and som data together 
    [~,In_pca]=pca(inputdata);
    [~,SOM_pca]=pca(netW);

    [In_tsne]=tsne(inputdata);
    [SOM_tsne]=tsne(netW);
    subplot(241) ;PlotClusters(In_pca,clust_ind_In_pheno,   'pca input',fp) % pca of input
    subplot(242) ;PlotClusters(SOM_pca,clust_ind_SOM, 'pca SOM Weights',fp)% pca od SOM
    subplot(245) ;PlotClusters(In_tsne,clust_ind_In_pheno,  'tsne input',fp) % tsne of input
    subplot(246) ;PlotClusters(SOM_tsne,clust_ind_SOM,'tsne SOM Weights',fp)% tsne of SOM
    subplot(243) ;PlotClusters(In_pca,Clust_ind_Input,'pca input->SOM',fp) % SOM clustring of inputs
    subplot(247) ;PlotClusters(In_tsne,Clust_ind_Input ,'tsne input->SOM',fp)  % SOM clustring of inputs
    
    if strcmpi(opts.clustringMethod,'HierarClust')
        subplot(244);dendrogram(Linkage_In);title('dendogram input')
        subplot(248);dendrogram(Linkage_SOM);title('dendogram SOM')
    end
    
end

end

function PlotClusters(Data,ClustInd,Tit,fp)

Nclust=length(unique(ClustInd));
Col=distinguishable_colors(Nclust);

hold on
for i=1:Nclust
    ind=find(ClustInd==i); 
    plt1=arrayfun(@(x) plot(Data(x,1),Data(x,2),'marker',...
    '.','color',Col(i,:),'MarkerSize',7),ind,'UniformOutput',0); 
end
title([Tit '-' num2str(Nclust) ' clus'])
h=gca;
fp.FormatAxes(h);

end

