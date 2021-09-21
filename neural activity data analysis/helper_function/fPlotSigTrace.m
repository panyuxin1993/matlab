function [baseline] = fPlotSigTrace(y,ind_baseline,dur_thresh)
%FPLOTSIGTRACE plot single trace of activities, with label of significance
%   Detailed explanation goes here
baseline=y(ind_baseline);
ind_sig=(y-mean(y)>std(baseline)*3);
if isempty(dur_thresh)
    dur_thresh=0;
end
if sum(ind_sig)>dur_thresh
    plot(y,'k-','LineWidth',1.5);
else
    plot(y,'k-','LineWidth',0.5);
end
end

