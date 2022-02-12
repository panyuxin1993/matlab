function [data_binned,ts_binned] = fBinActivity(data,binsize,binstep,ts_raw)
%FBINACTIVITY using binsize and binstep to calculate a binned data, usually
%used in order to reduce noise. 
%Input-
%   data- data to be binned, usually nROI-by-nFrame matrix
%   binsize- bin size of each bin
%   binstep- bin step to each bin to go through data
%   ts_raw- if data have corresponding time stamps, they should have same
%       scale at temporal dimension, a new ts is also calculated base on
%       same binning parameters
%Output-
%   data_binned- binned data, nROI-by-nFrame_binned matrix
%   ts_binned- binned data

nframe=floor((size(data,2)-1)/binstep)+1;
if isempty(ts_raw)
    warning('time stamp empty')
    ts_raw=1:size(data,2);
end
if length(ts_raw)~=size(data,2)
    warning('time stamp length not equal to frame number of activities')
    ts_raw=1:size(data,2);
end
ts_binned=ts_raw(1:binstep:end);%ts also need to be binned
ts_binned=ts_binned(1:nframe);
data_binned=[];
for it=1:nframe
    ind1it=max(1,(it-1)*binstep-floor(binsize/2)+1);
    ind2it=min(size(data,2),(it-1)*binstep+floor(binsize/2)+1);
    binActivity=data(:,ind1it:ind2it);
    binActivityMean=nanmean(binActivity,2);
    data_binned=[data_binned,binActivityMean];
end
end

