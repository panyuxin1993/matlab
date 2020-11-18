function [LCI,dataChoice,correct_ind,error_ind,currentSession, animal, date] = fGetLCI(npast,fileFolder,method,varargin)
%FGETLCI get LCI from raw mat file
%   input:
%   npast-set how many sessions data to show, 
%   fileFolder-set where data from, 
%   method-'proportion'(default)|'switch',decide the method to calculate
%   LCI, for 'proportion', LCI=max(left/right lick number)/sum(left and
%   right lick number); for 'switch', LCI= 1-number of switch/total lick
%   number
%   varargin-set which session to be show in detail(one subplot show changes 
%       within that session),if it is string, it decide filename of chosen
%       session, if it is num, it set index of chosen session
%   Output:LCI-LCI as a funciton of training day,(1d-LCI of all, correct, error
%   trials), dataChoice-LCI within a session, correct_ind,error_ind-ind of
%   representative session
dirmat=strcat(fileFolder,'*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
sort(filenames);
if isnumeric(varargin{1})
    i_session=varargin{1};
elseif ischar(varargin{1})
    indSession=cellfun(@(x) strcmp(x,varargin{1}),filenames);
    i_session=find(indSession);
    if isempty(i_session)% the specified session not exist
        i_session=length(filenames);%default is to plot the last session
        warning('Specified session not exist, using last session');
    end
end
if nargin==1
    i_session=length(filenames);%default is to plot the last session
end
LCI=zeros(npast, 3);%每行一个session，每一列LCI，左LCI,右LCI
if i_session-npast+1<=0%如果此值大于总数，则将所有的数据画出来；小于总数则画出指定数量
    istart=1;
else
    istart=i_session-npast+1;
end
for isession=istart:i_session
    load(strcat(fileFolder,filenames{isession}));  
    if isfield(SessionSettings{1},'clickRate_')
            clickright=SessionSettings{1}.clickRate_;%标记规则
    elseif isfield(SessionSettings{1},'clickRate_R')
        clickright=SessionSettings{1}.clickRate_R;%标记规则
    end
    stim=double(cellfun(@(x) x.Stim_clickRate, SessionResults));
    correct_ind=zeros(length(SessionResults),1);%cor not choice of left
    error_ind=zeros(length(SessionResults),1);
    actionchoice=cellfun(@(x) x.Action_choice, SessionResults);%0左，1右，2miss，3violation
    if clickright==20
        correct_ind=(stim<50).*(actionchoice==1)+(stim>50).*(actionchoice==0);
        error_ind=(stim<50).*(actionchoice==0)+(stim>50).*(actionchoice==1);
    else
        correct_ind=(stim<50).*(actionchoice==0)+(stim>50).*(actionchoice==1);
        error_ind=(stim<50).*(actionchoice==1)+(stim>50).*(actionchoice==0);
    end
    correct_ind=logical(correct_ind);
    error_ind=logical(error_ind);
    dataChoice=zeros(length(SessionResults),3);%第1-3列分别表示左边舔水次数ln，右边舔水次数rn，consistency index =|ln-rn|/(ln+rn)
    if strcmp(method,'proportion')
        for i=1:length(SessionResults)
            if SessionResults{1,i}.Action_choice~=2%0左，1右，2miss，3violation
                dataChoice(i,1)= SessionResults{1, i}.Action_numLickLeft;
                dataChoice(i,2)= SessionResults{1, i}.Action_numLickRight;
                dataChoice(i,3)=abs(dataChoice(i,1)-dataChoice(i,2))/(dataChoice(i,1)+dataChoice(i,2));
            else
                dataChoice(i,1)=0;
                dataChoice(i,2)=0;
                dataChoice(i,3)=nan;
            end
        end
    elseif strcmp(method,'switch')
        for i=1:length(SessionResults)
            if SessionResults{1,i}.Action_choice~=2%0左，1右，2miss，3violation
                left=strsplit(SessionResults{1, i}.Action_lickTimeLeft,'|');
                left(:,1)=[];
                right=strsplit(SessionResults{1, i}.Action_lickTimeRight,'|');
                right(:,1)=[];
                dataChoice(i,1)= fTimesSwitch(str2double(left),str2double(right));
                dataChoice(i,2)= SessionResults{1, i}.Action_numLickLeft+SessionResults{1, i}.Action_numLickRight;
                dataChoice(i,3)=1-dataChoice(i,1)/dataChoice(i,2);
            else
                dataChoice(i,1)=0;
                dataChoice(i,2)=0;
                dataChoice(i,3)=nan;
            end
        end
    end
    LCI(isession-istart+1, 1)=nanmean(dataChoice(:,3));
    LCI(isession-istart+1, 2)=nanmean(dataChoice(correct_ind,3));
    LCI(isession-istart+1, 3)=nanmean(dataChoice(error_ind,3));
end
currentFilename=filenames{i_session};
[currentSession, animal, date]=fFormatBehFilename(currentFilename,'animal-yyyymmdd');
if isempty(i_session)%if not specify a session, just calculate LCI of all input data
    warning('this date session not exist');
    dataChoice=nan;
    correct_ind=nan;
    error_ind=nan;
end
end

function [n_switch]=fTimesSwitch(a,b)
%calculated switch times between a and a (ascend)
if sum(~isnan(a))==0||sum(~isnan(b))==0
    n_switch=0;
else
    c=union(a,b);
    c=sort(c);
    d=ismember(c,a);
    e=diff(d);
    n_switch=length(find(e));
end
                    
end