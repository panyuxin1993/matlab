function [ ModeOfROI,AccumF,PreToneCaS,MeanF0s] = Find_mode_fun_pyx(Data_extract,SavedCaTrials )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
OnsetT = Data_extract.Tone_onset_time;
OnsetF = OnsetT/SavedCaTrials.FrameTime;

f_raw = SavedCaTrials.f_raw;
% FrameN=cellfun(@(x) size(x,2),f_raw);
% maxFramN=max(FrameN);
maxOnsetF=max(OnsetF);
if length(OnsetT) ~= length(SavedCaTrials.f_raw(:,1,1))
   warning('Length of behavior data ~= length of calsium signal data','color','r'); 
end

for nROI = 1:SavedCaTrials.nROIs
    PreToneCaS{nROI} = [];
    AccumF{nROI} = [];
    MeanF0s{nROI} = [];
    for i = 1:length(OnsetT)  %for each trials
%         temp = squeeze(f_raw(i,nROI,1:OnsetF(i)));
        temp = zeros(1,maxOnsetF);
        temp(1:OnsetF(i)) = f_raw{i}(nROI,1:OnsetF(i));
        temp(OnsetF(i):end) = nan;
        PreToneCaS{nROI} = [PreToneCaS{nROI} temp];
%         PreToneCaS{nROI} = [PreToneCaS{nROI} temp'];
        MeanF0s{nROI}(i) = nanmean(temp);
        if i == 1
            AccumF{nROI}(i) = 1;
        else
            AccumF{nROI}(i) = OnsetF(i)+AccumF{nROI}(i-1);
        end
    end
    ModeOfROI(nROI) = mode(PreToneCaS{nROI});
    if ModeOfROI(nROI) ~= 0 & nanmedian(PreToneCaS{nROI})/ModeOfROI(nROI)> 100
        ModeOfROI(nROI) = nanmedian(PreToneCaS{nROI});
    end
    if ModeOfROI(nROI)==0
        ModeOfROI(nROI) = nanmedian(PreToneCaS{nROI});
    end
    
end
end

