function [ trialType ,rule,trialtypestr] = fGetTrialType( Data_extract,TestTrialNum ,f ,outputType,varargin)
%FGETTRIALTYPE Summary of this function goes here
%   Detailed explanation goes here
%  f: 1 all click rates; 
%   -2 according to difficulty and sensory ipsi-easy, hard, contra-easy,hard
%   -3 according to choice(ipsi/contra first lick after stimOnset),
%   -4 according to sensory(low/high click rate)
%   to combine data from different sessions together, transformed the sign
%   according to fiber site and rule
%   outputType={'struct','matrix¡®}£»for matrix,
%   1d-cor/err/miss/vio,2d-stimuli,3d-trial,4d-opto/ctrl
%     trialType.corInd=CorrTrial_out;
%     trialType.errInd=WroTrial_out;
%     trialType.missInd=MisTrial_out;
%     trialType.vioInd=VioTrial_out;
switch nargin
    case 1
        TestTrialNum = [1:length(Data_extract.Miss_Ind)];%get all trials in a session
        f=1;
        outputType='matrix';
    case 2
        f=1;
        outputType='matrix';
    case 3
        outputType='matrix';
    otherwise
        %
end
combineCorErr='divideCorErr';%default settings
side='left';%default hemisphere, since for imaging, all neurons are in left hemisphere
divideOpto=0;
if f<=2 && nargin>4
    side=varargin{1};%={'left','right'}
    if nargin>5 
        combineCorErr=varargin{2};%={'combineCorErr','divideCorErr'}
    end
end
if f>=3 && nargin<5
    warning('lack of information which hemisphere to record');
elseif f>=3
    side=varargin{1};%={'left','right'}
    if nargin>5 
        combineCorErr=varargin{2};%={'combineCorErr','divideCorErr'}
    end     
    if nargin>6 && strcmp(varargin{3},'divideOpto')
        divideOpto=1;
    else
        divideOpto=0;
    end
else
    if nargin>4
        side=varargin{1};%={'left','right'}
    end
    if nargin>5 
        combineCorErr=varargin{2};%={'combineCorErr','divideCorErr'}
    end    
    if nargin>=6 && strcmp(varargin{3},'divideOpto')
        divideOpto=1;
    else
        divideOpto=0;
    end
end


if isempty(TestTrialNum)
    TestTrialNum = [1:length(Data_extract.Miss_Ind)];%get all trials in a session
end
Mis = Data_extract.Miss_Ind(TestTrialNum);
Vio = Data_extract.Vio_Ind(TestTrialNum);
Act = Data_extract.Action_choice(TestTrialNum);
Stim = Data_extract.Stimulus(TestTrialNum);
Stimuli = unique(Stim);
animal=Data_extract.Animal_Name_File;
if isfield(Data_extract,'rule')
    rule=Data_extract.rule;
else%if not extract rule from Data_extract variable, then extract here
    if contains(animal,'pyx159')||contains(animal,'PYX172')||contains(animal,'pyx095')||contains(animal,'pyx196')||contains(animal,'pyx237')||contains(animal,'pyx248')||contains(animal,'pyx249')||contains(animal,'pyx260')||contains(animal,'pyx255')||contains(animal,'pyx240')
        rule='low click rate-left';
    else
        rule='low click rate-right';
    end
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
if f <=2
    for i = 1:length(Stimuli)
        if contains(rule,side) %equal and simpler than this expression: (strcmp(rule,'low click rate-left') && strcmp(side,'left'))||(strcmp(rule,'low click rate-right') && strcmp(side,'right'))
            Stimuli_cat=Stimuli;
        else
            Stimuli_cat=flip(Stimuli);
        end
        CorrTrial_temp(i,:) = logical(sum(CorrTrial,1).*(Stim == Stimuli_cat(i)));
        WroTrial_temp(i,:) = logical(sum(WroTrial,1).*(Stim == Stimuli_cat(i)));
        MisTrial_temp(i,:) = logical(sum(MisTrial,1).*(Stim == Stimuli_cat(i)));
        VioTrial_temp(i,:) = logical(sum(VioTrial,1).*(Stim == Stimuli_cat(i)));
    end
    if f == 1
        CorrTrial_out=CorrTrial_temp;
        WroTrial_out=WroTrial_temp;
        MisTrial_out=MisTrial_temp;
        VioTrial_out=VioTrial_temp;
    elseif f == 2  %2*n_stim column
        CorrTrial_out = fSqueeze2(CorrTrial_temp);
        WroTrial_out = fSqueeze2(WroTrial_temp);
        MisTrial_out = fSqueeze2(MisTrial_temp);
        VioTrial_out = fSqueeze2(VioTrial_temp);
    end
end
if f==3%based on lick to ipsi/contra side related to brain region recorded
    if strcmp(side,'left')
        ipsiLick=Data_extract.Left_lick_time;
        contraLick=Data_extract.Right_lick_time;
    elseif strcmp(side,'right')
        contraLick=Data_extract.Left_lick_time;
        ipsiLick=Data_extract.Right_lick_time;
    end
    stimOnset=Data_extract.Stim_onset_time;
    lickFirst=zeros(1,length(ipsiLick));%1-left,2-right,nan-no lick
    for i=1:length(ipsiLick)
        ipsiLick{i}=double(ipsiLick{i});
        contraLick{i}=double(contraLick{i});
        ipsiLick{i}(ipsiLick{i}<stimOnset(i))=[];
        contraLick{i}(contraLick{i}<stimOnset(i))=[];
        if isempty(ipsiLick{i}) && ~isempty(contraLick{i})
            lickFirst(i)=2;
        elseif ~isempty(ipsiLick{i}) && isempty(contraLick{i})
            lickFirst(i)=1;
        elseif isempty(ipsiLick{i}) && isempty(contraLick{i})
            lickFirst(i)=nan;
        else %~isempty(ipsiLick{i}) && ~isempty(contraLick{i})
            if ipsiLick{i}(1)<contraLick{i}(1)
                lickFirst(i)=1;
            else
                lickFirst(i)=2;
            end
        end
    end
    for i=1:2%corresponding to ipsi, contra
        CorrTrial_out(i,:)=logical(sum(CorrTrial,1).*(lickFirst==i));
        WroTrial_out(i,:)=logical(sum(WroTrial,1).*(lickFirst==i));
        MisTrial_out(i,:)=logical(sum(MisTrial,1));%miss trials have no lick
        VioTrial_out(i,:)=logical(sum(VioTrial,1).*(lickFirst==i));
    end
end
if f==4%group by sensory
    nStim=floor(length(Stimuli)/2);
    boundary=(Stimuli(nStim)+Stimuli(nStim+1))/2;
    sensory={Stim<boundary,Stim>boundary};
    if contains(rule,side)
        sensory_cat=sensory;
    else
        sensory_cat=flip(sensory);
    end
    for i=1:2%corresponding to low,high freqency
        CorrTrial_out(i,:)=logical(sum(CorrTrial,1).*sensory_cat{i});
        WroTrial_out(i,:)=logical(sum(WroTrial,1).*sensory_cat{i});
        MisTrial_out(i,:)=logical(sum(MisTrial,1).*sensory_cat{i});%miss trials have no lick
        VioTrial_out(i,:)=logical(sum(VioTrial,1).*sensory_cat{i});
    end    
end

%using a variable output for easier knowing what the grouping method is
if f==1
    trialtypestr=strcat('stimuli_',combineCorErr);
elseif f==2
    trialtypestr=strcat('difficulty_',combineCorErr);
elseif f==3
    trialtypestr=strcat('first lick_',combineCorErr);%equal to 'choice'
elseif f==4
    trialtypestr=strcat('sensory_',combineCorErr);
end

%output by designed output type
if divideOpto==0
    if strcmp(outputType,'struct')
        if strcmp(combineCorErr,'combineCorErr')
            trialType.do=logical(CorrTrial_out+WroTrial_out);
        else
            trialType.corInd=CorrTrial_out;
            trialType.errInd=WroTrial_out;
        end
        trialType.missInd=MisTrial_out;
        trialType.vioInd=VioTrial_out;
    elseif strcmp(outputType,'matrix')
        if strcmp(combineCorErr,'combineCorErr')
            trialType=zeros(3,size(CorrTrial_out,1),size(CorrTrial_out,2));
            trialType(1,:,:)=logical(CorrTrial_out+WroTrial_out);
            trialType(2,:,:)=MisTrial_out;
            trialType(3,:,:)=VioTrial_out;
        else
            trialType=zeros(4,size(CorrTrial_out,1),size(CorrTrial_out,2));
            trialType(1,:,:)=CorrTrial_out;
            trialType(2,:,:)=WroTrial_out;
            trialType(3,:,:)=MisTrial_out;
            trialType(4,:,:)=VioTrial_out;
        end
    end
elseif divideOpto==1 %for simplicity, output type only support matrix
    if strcmp(combineCorErr,'combineCorErr')
        trialType=zeros(3,size(CorrTrial_out,1),size(CorrTrial_out,2),2);%last dimension is opto/ctrl
        for i2=1:size(trialType,2)
            for i4=1:size(trialType,4)
                trialType(1,i2,:,i4)=logical((CorrTrial_out(i2,:)+WroTrial_out(i2,:)).*(Data_extract.Opto_trial_index==(i4-1)));
                trialType(2,i2,:,i4)=logical(MisTrial_out(i2,:).*(Data_extract.Opto_trial_index==(i4-1)));
                trialType(3,i2,:,i4)=logical(VioTrial_out(i2,:).*(Data_extract.Opto_trial_index==(i4-1)));
            end
        end
    else
        trialType=zeros(4,size(CorrTrial_out,1),size(CorrTrial_out,2),2);%last dimension is opto/ctrl
        for i2=1:size(trialType,2)
            for i4=1:size(trialType,4)
                trialType(1,i2,:,i4)=logical(CorrTrial_out(i2,:).*(Data_extract.Opto_trial_index==(i4-1)));
                trialType(2,i2,:,i4)=logical(WroTrial_out(i2,:).*(Data_extract.Opto_trial_index==(i4-1)));
                trialType(3,i2,:,i4)=logical(MisTrial_out(i2,:).*(Data_extract.Opto_trial_index==(i4-1)));
                trialType(4,i2,:,i4)=logical(VioTrial_out(i2,:).*(Data_extract.Opto_trial_index==(i4-1)));
            end
        end
    end
end

end

function[out]=fSqueeze2(temp)
        ncol=size(temp,1);
        CorrTrial_out = zeros(4,size(temp,2));
        CorrTrial_out(1,:)=temp(1,:);
        CorrTrial_out(2,:)=sum(temp(2:floor(ncol/2),:));
        CorrTrial_out(3,:)=sum(temp(ceil(ncol/2)+1:ncol-1,:));
        CorrTrial_out(4,:)=temp(ncol,:);
        out=logical(CorrTrial_out);
end
