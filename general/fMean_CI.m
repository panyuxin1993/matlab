function [ dataMean, dataCI ] = fMean_CI( matdata,p )
%FMEAN_CI Summary of this function goes here
%   calculate mean and confidence interval of a group of curves
%   Detailed explanation goes here
%   input matdata is m-n matrix, m samples with length of n; p is
%   {0.05,0.01,0.001}
%   output mean is 1-n vector, which is the mean of m samples; CI is 2-n
%   matrix, which is the confidence interval at n sampling point

dataMean=nanmean(matdata,1);
dataStd=nanstd(matdata,0,1);
% n=size(matdata,1);%how many example trace
matdataused=(~isnan(matdata));
n=sum(matdataused,1);%how many example trace
p=p*ones(size(n));
dataCI=[dataMean+tinv(p/2,n-1)./sqrt(n).*dataStd;dataMean+tinv(1-p/2,n-1)./sqrt(n).*dataStd];
end

