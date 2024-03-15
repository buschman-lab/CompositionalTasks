%% builds and autoencoder on motif images to reduce dimensionality 
function features=BuildAutoEncoderMotifs(Motifs,hiddenSize,varargin)

%% define important vars and classes
global AnalysisOpts AnalysisData

ManData=ManipulateData;
%% define options of autoencoder
opts.L2WeightReg=0.004;
opts.SparsityReg=4;
opts.SparsityProportion=0.15;
opts.hiddenSize=hiddenSize;

%Process optional inputs
if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs for options.'); end
for i = 1:2:length(varargin),
    try
        opts.(varargin{i}) = varargin{i+1};
    catch
        error('Couldn''t set option ''%s''.', varargin{2*i-1});
    end
end
%% reshape data to match autoencoder 
XTrain=ManData.Reshape2DMat2Cell(Motifs,AnalysisData.SizeW);

%% train autoencoder now 
features{1}=XTrain; % initialize it with our data
for i=1:length(opts.hiddenSize)
    fprintf('\nTraining encoder %i on %i data points',i,length(features{i}));
    autoenc{i} = trainAutoencoder(features{i},opts.hiddenSize(i),'L2WeightRegularization',opts.L2WeightReg,...
    'SparsityRegularization',opts.SparsityReg,'SparsityProportion',opts.SparsityProportion,'ShowProgressWindow',true);
    
    features{i+1}=encode(autoenc{i},features{i}); % evaluate the fueatures of this layer
end
% 
% features1=encode(autoenc1,XTrain);
% %% second autoencoder
% autoenc2 = trainAutoencoder(features1,opts.hiddenSize(2),'L2WeightRegularization',opts.L2WeightReg,...
% 'SparsityRegularization',opts.SparsityReg,...
% 'SparsityProportion',opts.SparsityProportion);
% 
% features2=encode(autoenc2,features1);
% %% third autoencoder
% autoenc3 = trainAutoencoder(features2,opts.hiddenSize(3),'L2WeightRegularization',opts.L2WeightReg,...
% 'SparsityRegularization',opts.SparsityReg,...
% 'SparsityProportion',opts.SparsityProportion);
% 
% features3=encode(autoenc3,features2);
% 
% %% forth autoencoder
% autoenc4 = trainAutoencoder(features3,opts.hiddenSize(4),'L2WeightRegularization',opts.L2WeightReg,...
% 'SparsityRegularization',opts.SparsityReg,...
% 'SparsityProportion',opts.SparsityProportion,'EncoderTransferFunction','purelin');
% 
% features4=encode(autoenc4,features2);
% 
% %% stack them together
% satckednet=stack(autoenc1,autoenc2,autoenc3,autoenc4);
% 
% Out=encode(satckednet,XTrain)
% end