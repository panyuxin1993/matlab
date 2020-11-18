function [lickTime,maxdelay] =fGetLickTimes(sessionresult,align)%��õ���session��lick time��%align={'stim','go'};
    lickTime=struct('trialType',[],'go',[],'lickleft',[],'lickright',[],'choice',[]);
    maxdelay=0;
    for i=1:length(sessionresult)
        if strcmp(align,'go')
%             base=sessionresult{i}.Time_delayOffset;%offset���յ�vioӰ��
            base=sessionresult{i}.Time_delayOffset;
        else
            base=sessionresult{i}.Time_stimOnset;
        end
%         lickTime(i).delay=sessionresult{i}.Time_delayOffset-sessionresult{i}.Time_stimOnset;%delayOffset�����⣬violation trial��һ��violate��delay offset��
%       Ҳ��Ӧ����delayOnset+delay_duration�������violation������stimulus�׶���������취
        lickTime(i).go=sessionresult{i}.Time_delayOffset-base;
        lickTime(i).trialType=sessionresult{i}.Trial_Type;%��0��1
        left=strsplit(sessionresult{i}.Action_lickTimeLeft,'|');
        left(:,1)=[];
        lickTime(i).lickleft=str2double(left)-double(base);
        right=strsplit(sessionresult{i}.Action_lickTimeRight,'|');
        right(:,1)=[];
        lickTime(i).lickright=str2double(right)-double(base);
        lickTime(i).choice=double(sessionresult{i}.Action_choice);
    end     
end

