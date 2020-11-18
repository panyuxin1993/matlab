function [animal_name,dataChoiceResult_final ] = fFindChoiceMat(fileFolder,varargin )
%FFINDCHOICEMAT Summary of this function goes here
%   Input: fileFolder that include all the files to be read
%   varargin={'session','animal'};which mean each cell of dataChoiceResult is a
%   sesion or an animal
%   dataChoiceResult{i}=zeros(length(SessionResults),6);
%   第1-6列分别表示
%   1-stimuli
%   2-Actionresults(1-correct,2-error,3-miss,4-violation,0-default,no data)
%   3-刺激对应的左右(用于旋转curve左右),uni也可以用
%   4-Delay_duration，
%   5-reaction time(from go cue)(using SessionResults{1, i}.Time_answer), 
%   6-experiment condition(1,3 规则左低右高，2,4左高右低；1,2打药左侧，3,4打药右侧
%   7-licking consistency index, each trials, for licking number, LCI=|left-right|/(left+right)
%   
dirmat=strcat(fileFolder,'\*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
sort(filenames);%将文件按照训练先后顺序排列起来
%dataChoiceResult{i}(:,3)对应规则
% click=[20,28,40,62,90,125];%下文直接从文件中总结出来
%stimlr=[1,2,3,4,5,6];
%四种uni-实验条件
curvetype=struct('rule',[],'site',[]);
curvetype.rule={'pyx012','pyx014','pyx019','pyx076','pyx071','pyx073';'pyx011','pyx015','pyx017','pyx070','pyx072','pyx079'};
curvetype.site={'left';'right'};

%读数据
%对于选择，对1，错2，miss为3，violation为4，没有数据就是0
%%对于actionchoice，0左，1右，2miss，3violation
dataChoiceResult=cell(length(filenames)+1,1);%最后一个cell表示总结果
clickrightrule=zeros(length(filenames),1);%总结果没有对应于某一边
animal_name=cell(length(filenames),1);% 最后通过end+1增加一个cell，将其最后一个命名为mean
n=0;%总trial数量
for i=1:length(filenames)
    load(strcat(fileFolder,'\',filenames{i}));
    dataChoiceResult{i}=zeros(length(SessionResults),7);%***%第1-5列分别表示当前刺激，Actionresults，刺激对应的左右,Delay_duration，RT from go cue,exp condition,LCI,
    if isfield(SessionSettings{1},'clickRate_')%老版python程序得出的结果
        clickright=SessionSettings{1}.clickRate_;%标记规则
    elseif isfield(SessionSettings{1},'clickRate_R')%基于v6最新版本Python得出结果，自19-12-15
        clickright=SessionSettings{1}.clickRate_R;
    end
    animal_name{i}=SessionSettings{1, 1}.animalName;
    for j=1:length(SessionResults)%分两次循环，首次只寻找stimulus
        %不能用exist('SessionResults{1,j}.Stim_clickRate','var');这种语句，field应该用isfield
        if  isfield(SessionResults{1,j},'Stim_clickRate') && ~isempty(SessionResults{j}.Stim_clickRate)
            if SessionResults{j}.Stim_clickRate<126%存在125170500这种类型，和后面volume合并的
                dataChoiceResult{i}(j,1)= SessionResults{j}.Stim_clickRate;
            else 
                dataChoiceResult{i}(j,1)= floor(SessionResults{j}.Stim_clickRate/10^6);
            end
        elseif isfield(SessionResults{1,j},'Stim_clickRate') && isempty(SessionResults{j}.Stim_clickRate)
            dataChoiceResult{i}(j,1)=nan;
        elseif ~isfield(SessionResults{1,j},'Stim_clickRate') %'clicks20'这种类型.或者20
            if isfield(SessionResults{1,j},'Stim_Type')
            dataChoiceResult{i}(j,1)=  str2double(SessionResults{j}.Stim_Type(7:end));
%             dataChoiceResult{i}(j,1)=  str2double(SessionResults{j}.Stim_Type);
            else
                dataChoiceResult{i}(j,1)=nan;
            end
        else
            dataChoiceResult{i}(j,1)=nan;
        end
        if isnan(dataChoiceResult{i}(j,1))%如果实在找不到stim，只能跳过这行了
            continue;
        end
    end
    clicktemp=unique(dataChoiceResult{i}(:,1));%适应多个probe
    click=[];
    for iclick=1:length(clicktemp)%删掉其中数量太少的点，提高程序容错能力
        if sum(clicktemp(iclick)==dataChoiceResult{i}(:,1))>10 %小于10个数据点认为是记录数据的乱码
            click=[click;clicktemp(iclick)];
        end
    end
        
    for j=1:length(SessionResults)
        freq_n=length(click);
        stimlr=sort(click);%asending sort clicks as left-right shifted stimuli
%         stimlr=1:freq_n;
        actionchoice=SessionResults{j}.Action_choice;%0左，1右，2miss，3violation
        ind=find(click==dataChoiceResult{i}(j,1));
        if isempty(ind)%如果找到stim,但不在设定的库里，跳过这个trial
            dataChoiceResult{i}(j,:)=nan;%该行设为nan
            continue;
        end
        %通过规则确定trial 为cor,err,miss,vio;以及是否需要改变刺激对应的side
        if clickright==20
            if (dataChoiceResult{i}(j,1)<50 && actionchoice==1) ||(dataChoiceResult{i}(j,1)>50 && actionchoice==0) 
                dataChoiceResult{i}(j,2)=1;
            elseif (dataChoiceResult{i}(j,1)<50 && actionchoice==0) ||(dataChoiceResult{i}(j,1)>50 && actionchoice==1) 
                dataChoiceResult{i}(j,2)=2;
            else
                dataChoiceResult{i}(j,2)=actionchoice+1;%actionchoice=2,3,即miss或者violation
            end
            dataChoiceResult{i}(j,3)=stimlr(freq_n+1-ind);
        else
            if (dataChoiceResult{i}(j,1)<50 && actionchoice==0) ||(dataChoiceResult{i}(j,1)>50 && actionchoice==1) 
                dataChoiceResult{i}(j,2)=1;
            elseif (dataChoiceResult{i}(j,1)<50 && actionchoice==1) ||(dataChoiceResult{i}(j,1)>50 && actionchoice==0) 
                dataChoiceResult{i}(j,2)=2;
            else
                dataChoiceResult{i}(j,2)=actionchoice+1;%actionchoice=2,3,即miss或者violation
            end
            dataChoiceResult{i}(j,3)=stimlr(ind);
        end
        dataChoiceResult{i}(j,4)=SessionResults{j}.Delay_duration;
        dataChoiceResult{i}(j,5)=SessionResults{j}.Time_answer-SessionResults{j}.Time_delayOffset;%对于violation trial， time_answer=0，该值为负值，不可能为0，因为lick=delayOffset时也是violation
        if isfield(SessionSettings{j},'ExpCondition')&&isfield(SessionSettings{j},'animalName')%早期版本的.mat没有ExpCondition等field。则默认为零
            prule=cellfun(@(x) strfind(SessionSettings{j}.animalName,x), curvetype.rule,'UniformOutput', false);
            psite=cellfun(@(x) strfind(SessionSettings{j}.ExpCondition,x), curvetype.site,'UniformOutput', false);
            prule=sum(~cellfun(@isempty,prule),2);%此行将3*2 cell变成了3*2 logical，并通过按照第二维加和，成为1*2数组
            psite=sum(~cellfun(@isempty,psite),2);
            if psite(1)>0  && prule(1)>0 
                dataChoiceResult{i}(j,6)=1;
            elseif psite(1)>0  && prule(2)>0 
                dataChoiceResult{i}(j,6)=2;
            elseif psite(2)>0  && prule(1)>0 
                dataChoiceResult{i}(j,6)=3;
            elseif psite(2)>0  && prule(2)>0 
                dataChoiceResult{i}(j,6)=4;
            else
                dataChoiceResult{i}(j,6)=0;
            end
        end
        %calculate licking consistency index(LCI)
        if SessionResults{j}.Action_choice~=2
            nll=double(SessionResults{j}.Action_numLickLeft);
            nrl=double(SessionResults{j}.Action_numLickRight);
            dataChoiceResult{i}(j,7)=abs(nll-nrl)/(nll+nrl);%LCI, need to change from int32 to double
        else
            dataChoiceResult{i}(j,7)=-1;%之前设置为nan，但是在处理时比较困难处理，计算容易出现nan；先设置为-1，倘有需要仍能够根据<0选出，并设为nan
        end
    end
    n=n+length(SessionResults);
end
dataChoiceResult{length(filenames)+1}=zeros(n,size(dataChoiceResult{1},2));
animal=unique(animal_name);
dataChoiceResult_animal=cell(length(animal)+1,1);%最后一个cell表示总结果
n=zeros(length(animal),1);%存每只动物trial数量
for i=1:length(filenames)
    temp=cellfun(@(x) strcmp(x,animal_name{i}),animal);
    animal_ind=find(temp);
    dataChoiceResult_animal{animal_ind}(n(animal_ind)+1:n(animal_ind)+size(dataChoiceResult{i},1),:)=dataChoiceResult{i}(:,:);
    n(animal_ind)=n(animal_ind)+size(dataChoiceResult{i},1);
end

if nargin==1%default, this way compatible with previous use
    dataChoiceResult_final=dataChoiceResult;
elseif strcmp(varargin{1},'session') 
    dataChoiceResult_final=dataChoiceResult;
elseif strcmp(varargin{1},'animal') 
    dataChoiceResult_final=dataChoiceResult_animal;
    animal_name=animal;
end    

n=0;%last cell save all the trial together into a large session
for i=1:length(dataChoiceResult_final)-1
%     for j=1:size(dataChoiceResult{i},1)
        dataChoiceResult_final{end}(n+1:n+size(dataChoiceResult_final{i},1),:)=dataChoiceResult_final{i}(:,:);
%     end
    n=n+size(dataChoiceResult_final{i},1);
end
animal_name{end+1}='mean';
    
end

