function [ CorrTrial_out,WroTrial_out,MisTrial_out] = TrialType_fun( Data_extract,TestTrialNum ,f)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%  f:   0 Low_train Low_pro High_pro High_train         1 all frequecies
    if isempty(TestTrialNum)
        TestTrialNum = [1:length(Data_extract.Miss_Ind)];
    end
    Mis = Data_extract.Miss_Ind(TestTrialNum);
    Act = Data_extract.Action_choice(TestTrialNum);
    Tone = Data_extract.Tone_frequency(TestTrialNum);
    Fre = unique(Tone);
    for i = 1:length(Fre)
        if i <= length(Fre)/2
            CorrTrial(i,:) = Tone == Fre(i) & Act == 0;
            WroTrial(i,:) = Tone == Fre(i) & Act == 1;
        else
            CorrTrial(i,:) = Tone == Fre(i) & Act == 1;
            WroTrial(i,:) = Tone == Fre(i) & Act == 0;  
        end
        MisTrial(i,:) = Tone == Fre(i) & Mis;
    end
    if f ==1
        CorrTrial_out = CorrTrial;
        WroTrial_out = WroTrial;
        MisTrial_out = MisTrial;
    elseif f == 0
        CorrTrial_out = CorrTrial([1:2 7:8],:);
        WroTrial_out = WroTrial([1:2 7:8],:);
        MisTrial_out = MisTrial([1:2 7:8],:);     
        CorrTrial_out(2:3,:) = [mean(CorrTrial(2:4,:))>0;mean(CorrTrial(5:7,:))>0];
        WroTrial_out(2:3,:) = [mean(WroTrial(2:4,:))>0;mean(WroTrial(5:7,:))>0];
        MisTrial_out(2:3,:) = [mean(MisTrial(2:4,:))>0;mean(MisTrial(5:7,:))>0];  
    end

end

