function [ lickMat,lickRate ] = fLickRate( lickTime, binSize,trialLength,timeMin, trialInd, binStep)
%FLICKRATE Summary of this function goes here
% calculate a lick matrix, where lick =1, non lick=0;
%Input-
%   lickTime- vector of predicted licks
%   binSize- bin size used to calculated lick rate(moving average)
%   binStep- bin step used to calculated lick rate(moving average)
%   trialLength- length of single 'trial'(whole time window to show the
%   licking rate change)
%   timeMin- the start of each 'trial', relative to base as time 0.
%   trialInd- index of trial choosing to calculate licking rate.
%See also fLickRateAligned
sampleNum=floor((trialLength-binSize)/binStep)+1;
lickMat=zeros(length(lickTime),trialLength,3);%第一维trials,第二维单个trial长度，第三维代表总舔水，左右舔水
lickRate=zeros(sampleNum,3);%第一维代表trial长度/采样数量，第二维代表总舔水，左右舔水
for i=1:length(lickTime)
    for j=1:length(lickTime(i).lickleft)
        if ~isnan(lickTime(i).lickleft(j)) && lickTime(i).lickleft(j)-timeMin<size(lickMat,2) && lickTime(i).lickleft(j)-timeMin>0
            lickMat(i,lickTime(i).lickleft(j)-timeMin,2)=1;
        end
    end
    for j=1:length(lickTime(i).lickright)
        if ~isnan(lickTime(i).lickright(j)) && lickTime(i).lickright(j)-timeMin<size(lickMat,2) &&lickTime(i).lickright(j)-timeMin>0
            lickMat(i,lickTime(i).lickright(j)-timeMin,3)=1;
        end
    end
end
lickMat(:,:,1)=lickMat(:,:,2)+lickMat(:,:,3);
lickMat(~trialInd,:,:)=nan;%if trials not chosen, set to nan
for i=1:3
    for j=1:sampleNum-ceil(binSize/binStep)
        lickRate(j+ceil(binSize/binStep),i)=nansum(nansum(lickMat(:,(binStep*(j+ceil(binSize/binStep)-1)-binSize+1):(binStep*(j+ceil(binSize/binStep)-1)),i)))*1000/binSize/sum(trialInd);
    end
end


