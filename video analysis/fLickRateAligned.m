function [ lickMat,lickRate ] = fLickRateAligned( lickTime, base,binSize,trialLength,timeMin, binStep)
%fLickRateAligned calculate lickMat,a lick matrix, where lick =1, non lick=0; 
% and lickRate the corresponding licking rate vector when aigning to
% specific event. The primary scheme is to compare predicted lick against
% ground true. So the predicted licks should be aligned to ground true. 
%Input-
%   lickTime- vector of predicted licks
%   base- vector of event to be alinged, e.g. ground true licks, not
%   necessarily same size as lickTime
%   binSize- bin size used to calculated lick rate(moving average)
%   binStep- bin step used to calculated lick rate(moving average)
%   trialLength- length of single 'trial'(whole time window to show the
%   licking rate change)
%   timeMin- the start of each 'trial', relative to base as time 0.
%See also fLickRate
if iscell(lickTime)||iscell(base)
    warning('input is cell, please use fLickRateAlignedCell function');
end
sampleNum=floor((trialLength-binSize)/binStep)+1;%单位ms
lickMat=zeros(length(base),trialLength);%第一维trials,第二维单个trial长度
lickRate=zeros(sampleNum,1);%第一维代表trial长度/采样数量
for i=1:length(base)
    lickT=(lickTime-base(i))*1000;%单位s
    for j=1:length(lickT)
        if ~isnan(lickT(j)) && lickT(j)>-timeMin+1 && lickT(j)<trialLength-timeMin-1
            lickMat(i,round(lickT(j))+timeMin)=1;
        end
    end
end
for j=1:sampleNum-ceil(binSize/binStep)
    lickRate(j+ceil(binSize/binStep))=nansum(nansum(lickMat(:,(binStep*(j+ceil(binSize/binStep)-1)-binSize+1):(binStep*(j+ceil(binSize/binStep)-1)))))*1000/binSize/length(base);
end

end