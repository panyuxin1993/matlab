%local normalized, not for comparison between different trial types
function [dataout]=fLocalNormalized(datain, method)
%datain n-by-m matrix, each row a case, each column a time point
if strcmp(method,'minmax')
    disp('local normalized by (x-min)/(max-min)');
    minDataVal=min(datain,[],2);
    maxDataVal=max(datain,[],2);
    dataout=(datain-minDataVal)./(maxDataVal-minDataVal);
elseif strcmp(method,'z-score')
    disp('local normalized by z-score');
    dataout=(datain-nanmean(datain,2))./nanstd(datain,0,2);
end
end