function [outputcurve] = fPlotMean_CI(ts,neuralActivity,color,pSig,varargin)
%FPLOTMEAN_CI to plot mean and CI
%Input- 
%   ts- time stamps
%   neuralActiivty- n-by-m matrix, n trials, m=length(ts)
%   color- color for the plot
%   varargin- deciding whether or not plot additional lines, like
%   individual cases.

[ neuralActivityMean, neuralActivityCI ] = fMean_CI( neuralActivity,pSig );
if ~isempty(varargin)
    method=varargin{1};
    if contains(method,'cases')
        %plot individual traces
        color_case=(1+color)/2;
        plot(ts(1:size(neuralActivity,2)),neuralActivity,'Color',color_case,'linewidth',1);
        hold on;
    end
end

%shadow as se
xpatch=[ts(1:size(neuralActivity,2)), fliplr(ts(1:size(neuralActivity,2)))];
ypatch=[neuralActivityCI(1,:),fliplr(neuralActivityCI(2,:))];
indnan=isnan(xpatch)+isnan(ypatch);
if ~isempty(xpatch(~indnan))% remove nan values, and if all are nans, not plot se as shadow
    p=patch(xpatch(~indnan),ypatch(~indnan),color);%plot confidence interval
    p.FaceAlpha=0.3;
    p.EdgeColor='none';%color_mean_trace{nStim};%'none';
    hold on;
end

outputcurve=plot(ts(1:size(neuralActivity,2)),neuralActivityMean,'Color',color,'linewidth',2);
hold on;

end

