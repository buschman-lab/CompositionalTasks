function [start_ind_all, end_ind_all, cluster_p_all, obs_clust_size_all,varargout] = ClusterCorrectionTJB(obs_data, shuff_data,SP, varargin)

ManData=ManipulateData;

% Threshold parameters
opts.Threshold = NaN;
opts.ThresholdPercentage = 5;
opts.UseSingleThreshold = 0;
opts.ZscoreData=0; % do zscore the observed and shuffle

% For double-sided tests
opts.DoubleSided = 0;
opts.LowerThreshold = NaN;


% How to calculate cluster size. Can 'Integral' for summing the function (equivalent to summming above threshold) or 'CountIndex' for counting # of bins
opts.ClusterSizeType = 'Integral';

%Plotting options
opts.PlotFigure = 1;
opts.PlotSignificantClusts=0;
opts.PlotSignificantClustsTH=0.06;
opts.Title = '';
opts.Time = [];

%Update from passsed name/value pairs
for i = 1:2:length(varargin),
    opts.(varargin{i}) = varargin{i+1};
end
 ExtTitle=['TH' num2str(opts.ThresholdPercentage)];
if opts.DoubleSided;ExtTitle=[ExtTitle 'DblSide'];end

% Data size variables
data_len = size(obs_data, 1);
if size(shuff_data, 1) ~= data_len,
    error('Shuffle data must match size of observed data.');
end
num_shuffs = size(shuff_data, 2);

% In case we are z-scoring the data 
if opts.ZscoreData
    [obs_data, shuff_data]=ManData.CalZscoreShuffle(obs_data, shuff_data);
end

% Determine the threshold
if opts.UseSingleThreshold, %Take a single value as our threshold
    if isnan(opts.Threshold), % Set threshold as percentage of the original data
        opts.Threshold = prctile(shuff_data(:), 100-opts.ThresholdPercentage);
    end
    if isnan(opts.LowerThreshold), %if lower threshold isn't specified need to set it
        if ~isnan(opts.Threshold),
            % either as inverse of threshold
            opts.LowerThreshold = -opts.Threshold;
        else 
            % or as percentage of the original data
            opts.LowerThreshold = prctile(shuff_data(:), opts.ThresholdPercentage);
        end
    end
    opts.Threshold = opts.Threshold*ones(data_len, 1);
    opts.LowerThreshold = opts.LowerThreshold*ones(data_len, 1);
else % Use a threshold that varies over time
    if isnan(opts.Threshold), % Set threshold as percentage of the original data
        opts.Threshold = prctile(shuff_data, 100-opts.ThresholdPercentage, 2);
    elseif isscalar(opts.Threshold),
        opts.Threshold = opts.Threshold*ones(data_len, 1);
    end
    if isnan(opts.LowerThreshold), % Set lower threshold as percentage of the original data
        opts.LowerThreshold = prctile(shuff_data, opts.ThresholdPercentage, 2);
    elseif isscalar(opts.LowerThreshold),
        opts.LowerThreshold = -opts.Threshold*ones(data_len, 1);
    end
end

% Apply threshold do the observed data and shuffled data
thresh_obs_data = (obs_data >= opts.Threshold);
thresh_shuff_data = (shuff_data >= repmat(opts.Threshold, [1 num_shuffs]));
if opts.DoubleSided,
    lower_thresh_obs_data = (obs_data <= opts.LowerThreshold);
    lower_thresh_shuff_data = (shuff_data <= repmat(opts.LowerThreshold, [1 num_shuffs]));
end

% Find times when the observed data started/ended above threshold
start_ind = find(diff(cat(1, 0, thresh_obs_data)) == 1);
end_ind = find(diff(cat(1, thresh_obs_data, 0)) == -1);
% Calculate cluster size
obs_clust_size = zeros(length(start_ind), 1);
for i = 1:length(start_ind),
    obs_clust_size(i) = CalcClustSize(obs_data, start_ind(i), end_ind(i), opts, 'above');
end

% Do the same for lower threshold, if requested
if opts.DoubleSided,    
    lower_start_ind = find(diff(cat(1, 0, lower_thresh_obs_data)) == 1);
    lower_end_ind = find(diff(cat(1, lower_thresh_obs_data, 0)) == -1);
    % Calculate cluster size
    obs_lower_clust_size = zeros(length(lower_start_ind), 1);
    for i = 1:length(lower_start_ind),
        obs_lower_clust_size(i) = CalcClustSize(obs_data, lower_start_ind(i), lower_end_ind(i), opts, 'below');
    end
end


%Initialize shuffle cluster size
shuff_clust_size = zeros(num_shuffs, 1);
shuff_lower_clust_size = zeros(num_shuffs, 1);
for cur_shuff = 1:num_shuffs,
    % Find start and end of clusters
    shuff_start_ind = find(diff(cat(1, 0, thresh_shuff_data(:, cur_shuff))) == 1);
    shuff_end_ind = find(diff(cat(1, thresh_shuff_data(:, cur_shuff), 0)) == -1);
    % If no cluster, its size is zero
    if ~isempty(shuff_start_ind),
        % Calculate the size of all of the found clusters
        temp_shuff_clust_size = zeros(length(shuff_start_ind), 1);
        for i = 1:length(shuff_start_ind),
            temp_shuff_clust_size(i) = CalcClustSize(shuff_data(:, cur_shuff), shuff_start_ind(i), shuff_end_ind(i), opts, 'above');
        end
        % Take the maximum to correct for multiple comparisons over time
        shuff_clust_size(cur_shuff) = max(temp_shuff_clust_size);
    end

    % Do the same for lower threshold?
    if opts.DoubleSided,        
        % Find start and end of clusters
        shuff_start_ind = find(diff(cat(1, 0, lower_thresh_shuff_data(:, cur_shuff))) == 1);
        shuff_end_ind = find(diff(cat(1, lower_thresh_shuff_data(:, cur_shuff), 0)) == -1);
        % If no cluster, its size is zero
        if ~isempty(shuff_start_ind),
            % Calculate the size of all of the found clusters
            temp_shuff_clust_size = zeros(length(shuff_start_ind), 1);
            for i = 1:length(shuff_start_ind),
                temp_shuff_clust_size(i) = CalcClustSize(shuff_data(:, cur_shuff), shuff_start_ind(i), shuff_end_ind(i), opts, 'below');
            end
            % Take the maximum to correct for multiple comparisons over time
            shuff_lower_clust_size(cur_shuff) = max(temp_shuff_clust_size);
        end
    end
end

% Calculate the probability of each observed cluster size
cluster_p = NaN*ones(length(obs_clust_size), 1);
for i = 1:length(obs_clust_size),
   % cluster_p(i) = 1-mean(obs_clust_size(i)>shuff_clust_size);
    cluster_p(i) = 1-mean(obs_clust_size(i)>[shuff_clust_size ;obs_clust_size(i)]);
end
if opts.DoubleSided,
    lower_cluster_p = NaN*ones(length(obs_lower_clust_size), 1);
    for i = 1:length(obs_lower_clust_size),
     %   lower_cluster_p(i) = mean(obs_lower_clust_size(i)>shuff_lower_clust_size);
        lower_cluster_p(i) = mean(obs_lower_clust_size(i)>[shuff_lower_clust_size;obs_lower_clust_size(i)]);
    end
end


% Plot figure
if opts.PlotFigure 
    if isempty(SP);figure;SP={[1 2 1],[1 2 2]};end
    if ~isempty(SP{1})
        subplot(SP{1}(1),SP{1}(2),SP{1}(3));
        histogram(shuff_clust_size, 50);
        v = axis;
        hold on;
        for i = 1:length(start_ind),
            plot(obs_clust_size(i)*[1 1], v(3:4), '-');
        end
        if opts.DoubleSided,
            histogram(shuff_lower_clust_size, 'FaceColor', 'red');
            v = axis;
            hold on;
            for i = 1:length(lower_start_ind),
                plot(obs_lower_clust_size(i)*[1 1], v(3:4), ':');
            end
        end
        title(['Distribution TJB code ' ExtTitle])
    end
    if ~isempty(SP{2})
        if opts.PlotSignificantClusts
            SigCl=cluster_p<=opts.PlotSignificantClustsTH;
            cluster_p=cluster_p(SigCl);
            start_ind=start_ind(SigCl);
            end_ind=end_ind(SigCl);
            if opts.DoubleSided
                SigCl=lower_cluster_p<=opts.PlotSignificantClustsTH;
                lower_cluster_p=lower_cluster_p(SigCl);
                lower_start_ind=lower_start_ind(SigCl);
                lower_end_ind=lower_end_ind(SigCl);
            end
        end
        subplot(SP{2}(1),SP{2}(2),SP{2}(3));
        if isempty(opts.Time),
            opts.Time = [1:length(obs_data)];
        end
        plot(opts.Time, obs_data);
        hold on;
        v = axis;
        for i = 1:length(start_ind),
            h = plot(opts.Time([start_ind(i) end_ind(i)]), 0.8*v(4)*[1 1], '-', 'LineWidth', 2);
            th = text(opts.Time(start_ind(i)), 0.85*v(4), ...
                sprintf('p=%f',cluster_p(i)), ...
                "Color", get(h, "Color"));
        end
        if opts.DoubleSided,
            for i = 1:length(lower_start_ind),
                h = plot(opts.Time([lower_start_ind(i) lower_end_ind(i)]), 0.2*v(4)*[1 1], ':', 'LineWidth', 2);
                th = text(opts.Time(lower_start_ind(i)), 0.15*v(4), ...
                    sprintf('p=%f',lower_cluster_p(i)), ...
                    "Color", get(h, "Color"));
            end
        end
        xlabel('Time')
        title(['Cluster correction TJB code' ExtTitle])

        if ~isempty(opts.Title),
            title(opts.Title);
        end
    end
end
% Output all both upper and lower 
start_ind_all=[lower_start_ind; start_ind];
end_ind_all=[lower_end_ind ;end_ind];
cluster_p_all=[lower_cluster_p; cluster_p];
obs_clust_size_all=[obs_lower_clust_size ;obs_clust_size];
%Output the lower threshold variables if tested
if opts.DoubleSided,
    varargout{1} = lower_start_ind;
    varargout{2} = lower_end_ind;
    varargout{3} = lower_cluster_p; 
    varargout{4} = obs_lower_clust_size;
end

end %Cluster correction function


%% Function to calculate cluster size
function clust_size = CalcClustSize(data, start_ind, end_ind, opts, direction),

if strcmpi(opts.ClusterSizeType, 'CountIndex'),
    clust_size = end_ind - start_ind + 1;
elseif strcmpi(opts.ClusterSizeType, 'Integral'),
    if strcmpi(direction, 'above'),
        clust_size = sum(data(start_ind:end_ind));
    elseif strcmpi(direction, 'below'),
        clust_size = sum(data(start_ind:end_ind));
    else
        error('Bad direction specified.');
    end
elseif strcmpi(opts.ClusterSizeType, 'IntegralAboveThreshold'),
    if strcmpi(direction, 'above'),
        clust_size = sum(data(start_ind:end_ind) - opts.Threshold(start_ind:end_ind));
    elseif strcmpi(direction, 'below'),
        clust_size = sum(data(start_ind:end_ind) - opts.LowerThreshold(start_ind:end_ind));
    else
        error('Bad direction specified.');
    end
else
    error('Bad cluster size calculation specified.');
end

end % CalcClustSize function