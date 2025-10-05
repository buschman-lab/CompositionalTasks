clc; clear; close all;

% Create a figure
figure;

% Generate a dummy image to apply colormap
imagesc(peaks(100)); % Example data
colormap(copper(22));    % Set colormap to 'copper'
colorbar_handle = colorbar; % Create colorbar

% Set the colorbar ticks from 10 to 73 with a step of 3
tick_values = 10:3:73;
colorbar_handle.Ticks = linspace(min(caxis), max(caxis), numel(tick_values));
colorbar_handle.TickLabels = arrayfun(@num2str, tick_values, 'UniformOutput', false);

% Label the colorbar
ylabel(colorbar_handle, 'Intensity');

% Set title
title('Copper Colormap with Custom Ticks');
