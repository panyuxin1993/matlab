function [indBaseline] = fIndexBaselineAroundMode(signal,varargin)
%FINDEXBASELINEAROUNDMODE calculate index of baseline base on signals.
%Input-
%   signal- n-by-m matrix of signals. n signals, each length is m
%Output-
%   indBaseline- n-by-m matrix of baseline, same size with input signal
if isempty(varargin)
    sigThreshSD=2;
else
    if mod(length(varargin),2)==0
        for i=1:2:length(varargin)
            switch varargin{i}
                case 'sigThreshSTD'
                    sigThreshSD=varargin{i+1};
            end
        end
    else
        warning('odd varargin input argument');
    end
end
indBaseline=false(size(signal));
data_std=nanstd(signal,0,2);
baselineStd=zeros(size(data_std));
for i=1:size(signal,1)
    [N,edges] = histcounts(signal(i,:),100);
    f0 = edges(N == max(N));
    if length(f0)>1
        f0=mean(f0);
    end
    indBaseline_temp=logical((signal(i,:)<f0+sigThreshSD*data_std(i)).*(signal(i,:)>f0-sigThreshSD*data_std(i)));
    [ segflagcell,segrawcell] = fSegment2Cell( indBaseline_temp,signal(i,:));
    baselineDur=cellfun(@length, segrawcell);
    baselineMaxDur_cell=segrawcell(baselineDur==max(baselineDur));
    baselineMaxDur=baselineMaxDur_cell{1};%maybe there more than one cell
    baselineStd(i,:)=nanstd(baselineMaxDur,0,2);
    indBaseline(i,:)=logical((signal(i,:)<f0+sigThreshSD*baselineStd(i,:)).*(signal(i,:)>f0-sigThreshSD*baselineStd(i,:)));
end
% scatter(data_std,baselineStd);
end

