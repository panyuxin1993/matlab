function [ color ] = fMeanTraceColor( behrule,nStim )
%FMEANTRACECOLOR based on rule to decide which color used for mean trace 
%   Input is the rule of that behavior session
%   Output is color, a cell storing the color of individual traces
%   usually, left trials are blue,while right trials are red

color_left_low={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};%left-low--right-high
color_left_high=fliplr(color_left_low);%left-high--right-low
color=cell(1,nStim);
if mod(nStim,2)==1
    warning('Number of stimuli is odd');
elseif nStim>6
    warning('Number of stimuli is too much');
end
if strcmp(behrule,'low click rate-left')
    color_chosen=color_left_low;
else
    color_chosen=color_left_high;
end
for i=1:nStim/2
    color{1,i}=color_chosen{1,i};
end
for i=nStim/2+1:nStim
    color{1,i}=color_chosen{1,6-nStim+i};
end

end

