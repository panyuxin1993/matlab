close all,
clear,
clc;

M = im2double(imread('colorMap1.png'));
M = squeeze(M(ceil(size(M, 1)/2), :, :));

% figure(1);
% hold on;
% plot(M(:, 1), 'r', 'LineWidth', 2);
% plot(M(:, 2), 'g', 'LineWidth', 2);
% plot(M(:, 3), 'b', 'LineWidth', 2);
% title('RGB Curve of the ColorMap with Screenshot');
% xlim([0, size(M, 1)+1]);
% ylim([0, 1]);
% hold off;
% 
% colormap(M);
% colorbar;

%% Smoothing with polyfit
x = linspace(0, 1, size(M, 1))';
n = [10, 10, 8]; % control the fitting results
pr = polyfit(x, M(:, 1), n(1));
pg = polyfit(x, M(:, 2), n(2));
pb = polyfit(x, M(:, 3), n(3));

N = 2000; % control the precision
x2 = linspace(0, 1, N)';
M2 = [polyval(pr, x2), polyval(pg, x2), polyval(pb, x2)];
M2(M2 < 0) = 0;
M2(M2 > 1) = 1;

% figure(2);
% hold on;
% plot(M2(:, 1), 'r', 'LineWidth', 2);
% plot(M2(:, 2), 'g', 'LineWidth', 2);
% plot(M2(:, 3), 'b', 'LineWidth', 2);
% title('RGB Curve of the ColorMap with Polyfit Smoothing');
% xlim([0, size(M2, 1)+1]);
% ylim([0, 1]);
% hold off;
% 
% colormap(M2);
% colorbar;

%% 'M2' could be saved as '*.mat' in case of need. 
