function [ dataMean,dataSE ] = fMean_SE( matdata )
%FMEAN_SE calculate mean and standard error of a group of curves
%   Detailed explanation goes here
%   input matdata is m-n matrix, m samples with length of n
%   output mean is 1-n vector, which is the mean of m samples; dataSE is 2-n
%   matrix, which is the standard error of data
dataMean=nanmean(matdata,1);
dataStd=nanstd(matdata,1,1);
matdataused=(~isnan(matdata));
n=sum(matdataused,1);%how many example trace
dataSE=[dataMean+dataStd./sqrt(n);dataMean-dataStd./sqrt(n)];
%fill nan-se data as mean
dataSE(1,isnan(dataSE(1,:)))=dataMean(isnan(dataSE(1,:)));
dataSE(2,isnan(dataSE(1,:)))=dataMean(isnan(dataSE(1,:)));
end

