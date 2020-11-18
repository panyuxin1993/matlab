function [ trialType] = TrialType_fun_pyx( Data_extract,TestTrialNum ,f)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%  f: 1 all click rates; 2 easy, hard, 2 times column
%     trialType.corInd=CorrTrial_out;
%     trialType.errInd=WroTrial_out;
%     trialType.missInd=MisTrial_out;
%     trialType.vioInd=VioTrial_out;
    if nargin ==1
        TestTrialNum = [1:length(Data_extract.Miss_Ind)];%get all trials in a session
        f=1;
    end
    if nargin == 2
        f=1;
    end
    if isempty(TestTrialNum)
        TestTrialNum = [1:length(Data_extract.Miss_Ind)];%get all trials in a session
    end
    Mis = Data_extract.Miss_Ind(TestTrialNum);
    Act = Data_extract.Action_choice(TestTrialNum);
    Stim = Data_extract.Stimulus(TestTrialNum);
    Stimuli = unique(Stim);
    animal=Data_extract.Animal_Name_File;
    if contains(animal,'pyx083')
        rule='low click rate-left';
    else
        rule='low click rate-right';
    end
    if strcmp(rule,'low click rate-left')
        for i = 1:length(Stimuli)
            if i <= length(Stimuli)/2
                CorrTrial(i,:) = Stim == Stimuli(i) & Act == 0;
                WroTrial(i,:) = Stim == Stimuli(i) & Act == 1;
            else
                CorrTrial(i,:) = Stim == Stimuli(i) & Act == 1;
                WroTrial(i,:) = Stim == Stimuli(i) & Act == 0;
            end
            MisTrial(i,:) = Stim == Stimuli(i) & Mis;
            VioTrial(i,:) = Stim == Stimuli(i) & Vio;
        end
    else
        for i = 1:length(Stimuli)
            if i <= length(Stimuli)/2
                CorrTrial(i,:) = Stim == Stimuli(i) & Act == 1;
                WroTrial(i,:) = Stim == Stimuli(i) & Act == 0;
            else
                CorrTrial(i,:) = Stim == Stimuli(i) & Act == 0;
                WroTrial(i,:) = Stim == Stimuli(i) & Act == 1;
            end
            MisTrial(i,:) = Stim == Stimuli(i) & Mis;
            VioTrial(i,:) = Stim == Stimuli(i) & Vio;
        end
    end
    if f ==1
        CorrTrial_out = CorrTrial;
        WroTrial_out = WroTrial;
        MisTrial_out = MisTrial;
        VioTrial_out = VioTrial;
    elseif f == 2  %2*n_stim column
        CorrTrial_out = logical(CorrTrial(1,:)+ CorrTrial(end,:));
        WroTrial_out = logical(WroTrial(1,:)+ WroTrial(end,:));
        MisTrial_out = logical(MisTrial(1,:)+ MisTrial(end,:));   
        VioTrial_out = logical(VioTrial(1,:)+ VioTrial(end,:));
        CorrTrial_out = [CorrTrial_out, logical(sum(CorrTrial(2:end-1,:),2))];
        WroTrial_out = [WroTrial_out, logical(sum(WroTrial(2:end-1,:),2))];
        MisTrial_out = [MisTrial_out, logical(sum(MisTrial(2:end-1,:),2))];
        VioTrial_out = [VioTrial_out, logical(sum(VioTrial(2:end-1,:),2))];
    end
    trialType.corInd=CorrTrial_out;
    trialType.errInd=WroTrial_out;
    trialType.missInd=MisTrial_out;
    trialType.vioInd=VioTrial_out;
end

