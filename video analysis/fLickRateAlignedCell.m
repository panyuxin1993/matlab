function [ lickMat,lickRate ] = fLickRateAlignedCell( lickTime, base,binSize,trialLength,timeMin, binStep)
%FLICKRATEALIGNEDCELL Extended form of fLickRateAligned that support input
%data in the form of cell. 
%Input-
%   lickTime- vector of predicted licks
%   base- vector of event to be alinged, e.g. ground true licks, not
%   necessarily same size as lickTime
%   binSize- bin size used to calculated lick rate(moving average)
%   binStep- bin step used to calculated lick rate(moving average)
%   trialLength- length of single 'trial'(whole time window to show the
%   licking rate change)
%   timeMin- the start of each 'trial', relative to base as time 0.
%See also fLickRate and fLickRateAligned
if ~iscell(lickTime)||~iscell(base)
    warning('input is not cell, please use fLickRateAligned function');
end
nBaseElement=cellfun(@nansum, base);
nBaseElement=nansum(nBaseElement);%if base is ground true lickings, then nBaseElement is how many true licks 
sampleNum=floor((trialLength-binSize)/binStep)+1;%��λms
lickMat=[];%zeros(nBaseElement,trialLength);%��һάtrials,�ڶ�ά����trial����
lickRate=zeros(sampleNum,1);%��һά����trial����/��������
for i=1:length(base)
    lickMati=fLickRateAligned( lickTime{i}, base{i},binSize,trialLength,timeMin, binStep);
    lickMat=[lickMat;lickMati];
end
for j=1:sampleNum-ceil(binSize/binStep)
    lickRate(j+ceil(binSize/binStep))=nansum(nansum(lickMat(:,(binStep*(j+ceil(binSize/binStep)-1)-binSize+1):(binStep*(j+ceil(binSize/binStep)-1)))))*1000/binSize/length(base);
end
end


