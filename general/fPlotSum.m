function [ outputcurve ] = fPlotSum( ts,neuralActivity,color,varargin )
%FPLOTMEAN_SE to plot mean and sem
%Input- 
%   ts- time stamps
%   neuralActiivty- n-by-m matrix, n trials, m=length(ts)
%   color- color for the plot
%   varargin- deciding whether or not plot additional lines, like
%   individual cases.

[ neuralActivityMean, neuralActivitySE ] = fMean_SE( neuralActivity );
if length(ts)~=size(neuralActivity,2)
    warning(['ts and neuralActivity have inconsistent size along time axis:',num2str(length(ts)),' vs. ', num2str(size(neuralActivity,2))]);
    if length(ts)<size(neuralActivity,2)
        neuralActivity=neuralActivity(:,1:length(ts));
        neuralActivitySE=neuralActivitySE(:,1:length(ts));
        neuralActivityMean=neuralActivityMean(:,1:length(ts));
    end
end
if ~isempty(varargin)
    method=varargin{1};
    if contains(method,'cases')
        %plot individual traces
        color_case=((1+color)/2+1)/2;
        plot(ts(1:size(neuralActivity,2)),neuralActivity','Color',color_case,'linewidth',1);
        hold on;
        if contains(method,'only cases')
            return
        end
    end
end

% %shadow as se
% xpatch=[ts(1:size(neuralActivity,2)), fliplr(ts(1:size(neuralActivity,2)))];
% ypatch=[neuralActivitySE(1,:),fliplr(neuralActivitySE(2,:))];
% indnan=isnan(xpatch)+isnan(ypatch);
% if ~isempty(xpatch(~indnan))% remove nan values, and if all are nans, not plot se as shadow
%     p=patch(xpatch(~indnan),ypatch(~indnan),color);%plot confidence interval
%     p.FaceAlpha=0.3;
%     p.EdgeColor='none';%color_mean_trace{nStim};%'none';
%     hold on;
% end

outputcurve=plot(ts(1:size(neuralActivity,2)),neuralActivityMean*size(neuralActivity,1),'Color',color,'linewidth',2);
hold on;

end



