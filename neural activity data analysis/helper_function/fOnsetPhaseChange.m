function [onset,n_change] = fOnsetPhaseChange(series,threshold,varargin)
%FONSETPHASECHANGE find onset when a dominant stage changes(I call it phase
%change). The idea is to find the longest continuous duration of one stage
%and define the onset of that
%Input-
%   series- 0-1-vector to know the onset, to find the start of '1' series
%   leastduration- 5(default), a threshold for least continuous duration
%   threshold- confine the onset, e.g. not earlier than stimuli
%Output-
%   onset-a number of when start the longest continuous duration
%   n_change-number indicating how many times phase change happen, given
%   the criteria that duration is longer than least duration
if isempty(varargin)
    leastduration=5;
else
    leastduration=varargin{1};
end
diffSeries=diff(series);
diffSeries=reshape(diffSeries,1,[]);%to a row vector
diffSeries=[0,diffSeries];
tchange=find(diffSeries~=0);
if isempty(tchange)
    onset=[];
    n_change=0;
else
    if diffSeries(tchange(1))<0
        tchange=[1,tchange];
    end
    if diffSeries(tchange(end))>0
        tchange=[tchange,length(series)];
    end
end
%now, n-pairs, each up-down
ind_up=1:2:length(tchange);
ind_down=2:2:length(tchange);
duration=tchange(ind_down)-tchange(ind_up);
n_change=sum(duration>=leastduration);
ind_max=find((max(duration)==duration).*(duration>=leastduration));%this can result more than 1 vector
ind_1st=find(duration>=leastduration);
if isempty(ind_max)
    onset=length(series);%if not significant all the time, choose the end of the time series as starting point
else
%     onsets=tchange(2*ind_max(1)-1);
    onsets=tchange(2*ind_1st-1);
    onsets=onsets(onsets>threshold);
    if isempty(onsets)
        onset=length(series);
    else
        onset=onsets(1);
    end
end

end

