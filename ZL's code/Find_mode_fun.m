function [ ModeOfROI,AccumF,PreToneCaS,MeanF0s] = Find_mode_fun(Data_extract,SavedCaTrials )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
OnsetT = Data_extract.Tone_onset_time;
OnsetF = OnsetT/SavedCaTrials.FrameTime;

f_raw = SavedCaTrials.f_raw;
if length(OnsetT) ~= length(SavedCaTrials.f_raw(:,1,1))
   warning('Length of behavior data ~= length of calsium signal data','color','r'); 
end

for nROI = 1:SavedCaTrials.nROIs
    PreToneCaS{nROI} = [];
    AccumF{nROI} = [];
    MeanF0s{nROI} = [];
    for i = 1:length(OnsetT)
        temp = squeeze(f_raw(i,nROI,1:OnsetF(i)));
        PreToneCaS{nROI} = [PreToneCaS{nROI} temp'];
        MeanF0s{nROI}(i) = mean(temp);
        if i == 1
            AccumF{nROI}(i) = 1;
        else
            AccumF{nROI}(i) = OnsetF(i)+AccumF{nROI}(i-1);
        end
    end
    ModeOfROI(nROI) = mode(PreToneCaS{nROI});
    if ModeOfROI(nROI) ~= 0 & median(PreToneCaS{nROI})/ModeOfROI(nROI)> 100
        ModeOfROI(nROI) = median(PreToneCaS{nROI});
    end
    if ModeOfROI(nROI)==0
        ModeOfROI(nROI) = median(PreToneCaS{nROI});
    end
    
end
end

