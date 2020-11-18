function [lickTime,maxdelay] =fGetLickTimes(sessionresult,align)%获得单个session的lick time；%align={'stim','go','answer'};
    lickTime=struct('trialType',[],'stimOn',[],'stimOff',[],'go',[],'lickleft',[],'lickright',[],'delay',[],'choice',[],'answer',[],'delayOff',[]);
    maxdelay=0;
    for i=1:length(sessionresult)
        if strcmp(align,'go')
%             base=sessionresult{i}.Time_delayOffset;%offset会收到vio影响;
            base=sessionresult{i}.Time_delayOnset+sessionresult{i}.Delay_duration;
        elseif strcmp(align,'answer')
            base=sessionresult{i}.Time_answer;
        elseif strcmp(align,'delayOff')%this may not equal to time of go cue in 1)violation trials, 2)pre-go-no-lick trials
            base=sessionresult{i}.Time_delayOffset;
        else
            base=sessionresult{i}.Time_stimOnset;
        end
%         lickTime(i).delay=sessionresult{i}.Time_delayOffset-sessionresult{i}.Time_stimOnset;%delayOffset有问题，violation trial中一旦violate就delay offset。
%       也许应该用delayOnset+delay_duration；但如果violation发生在stimulus阶段需另外想办法
        lickTime(i).delay=sessionresult{i}.Time_delayOnset-sessionresult{i}.Time_stimOnset+sessionresult{i}.Delay_duration;
        lickTime(i).stimOn=sessionresult{i}.Time_stimOnset-base;
        lickTime(i).stimOff=sessionresult{i}.Time_stimOffset-base;
        lickTime(i).go=sessionresult{i}.Time_delayOnset+sessionresult{i}.Delay_duration-base;%where go cue should happen
        lickTime(i).trialType=sessionresult{i}.Trial_Type;%左0右1
        lickTime(i).answer=sessionresult{i}.Time_answer-base;
        lickTime(i).delayOff= sessionresult{i}.Time_delayOffset;%this is the real time of delay
        left=strsplit(sessionresult{i}.Action_lickTimeLeft,'|');
        left(:,1)=[];
        lickTime(i).lickleft=str2double(left)-double(base);
        right=strsplit(sessionresult{i}.Action_lickTimeRight,'|');
        right(:,1)=[];
        lickTime(i).lickright=str2double(right)-double(base);
        lickTime(i).choice=double(sessionresult{i}.Action_choice);
        maxdelay=max(maxdelay,lickTime(i).delay);
    end     
end
