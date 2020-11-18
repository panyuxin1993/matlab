function [animal_name,dataChoiceResult_final ] = fFindChoiceMat(fileFolder,varargin )
%FFINDCHOICEMAT Summary of this function goes here
%   Input: fileFolder that include all the files to be read
%   varargin={'session','animal'};which mean each cell of dataChoiceResult is a
%   sesion or an animal
%   dataChoiceResult{i}=zeros(length(SessionResults),6);
%   ��1-6�зֱ��ʾ
%   1-stimuli
%   2-Actionresults(1-correct,2-error,3-miss,4-violation,0-default,no data)
%   3-�̼���Ӧ������(������תcurve����),uniҲ������
%   4-Delay_duration��
%   5-reaction time(from go cue)(using SessionResults{1, i}.Time_answer), 
%   6-experiment condition(1,3 ��������Ҹߣ�2,4����ҵͣ�1,2��ҩ��࣬3,4��ҩ�Ҳ�
%   7-licking consistency index, each trials, for licking number, LCI=|left-right|/(left+right)
%   
dirmat=strcat(fileFolder,'\*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
sort(filenames);%���ļ�����ѵ���Ⱥ�˳����������
%dataChoiceResult{i}(:,3)��Ӧ����
% click=[20,28,40,62,90,125];%����ֱ�Ӵ��ļ����ܽ����
%stimlr=[1,2,3,4,5,6];
%����uni-ʵ������
curvetype=struct('rule',[],'site',[]);
curvetype.rule={'pyx012','pyx014','pyx019','pyx076','pyx071','pyx073';'pyx011','pyx015','pyx017','pyx070','pyx072','pyx079'};
curvetype.site={'left';'right'};

%������
%����ѡ�񣬶�1����2��missΪ3��violationΪ4��û�����ݾ���0
%%����actionchoice��0��1�ң�2miss��3violation
dataChoiceResult=cell(length(filenames)+1,1);%���һ��cell��ʾ�ܽ��
clickrightrule=zeros(length(filenames),1);%�ܽ��û�ж�Ӧ��ĳһ��
animal_name=cell(length(filenames),1);% ���ͨ��end+1����һ��cell���������һ������Ϊmean
n=0;%��trial����
for i=1:length(filenames)
    load(strcat(fileFolder,'\',filenames{i}));
    dataChoiceResult{i}=zeros(length(SessionResults),7);%***%��1-5�зֱ��ʾ��ǰ�̼���Actionresults���̼���Ӧ������,Delay_duration��RT from go cue,exp condition,LCI,
    if isfield(SessionSettings{1},'clickRate_')%�ϰ�python����ó��Ľ��
        clickright=SessionSettings{1}.clickRate_;%��ǹ���
    elseif isfield(SessionSettings{1},'clickRate_R')%����v6���°汾Python�ó��������19-12-15
        clickright=SessionSettings{1}.clickRate_R;
    end
    animal_name{i}=SessionSettings{1, 1}.animalName;
    for j=1:length(SessionResults)%������ѭ�����״�ֻѰ��stimulus
        %������exist('SessionResults{1,j}.Stim_clickRate','var');������䣬fieldӦ����isfield
        if  isfield(SessionResults{1,j},'Stim_clickRate') && ~isempty(SessionResults{j}.Stim_clickRate)
            if SessionResults{j}.Stim_clickRate<126%����125170500�������ͣ��ͺ���volume�ϲ���
                dataChoiceResult{i}(j,1)= SessionResults{j}.Stim_clickRate;
            else 
                dataChoiceResult{i}(j,1)= floor(SessionResults{j}.Stim_clickRate/10^6);
            end
        elseif isfield(SessionResults{1,j},'Stim_clickRate') && isempty(SessionResults{j}.Stim_clickRate)
            dataChoiceResult{i}(j,1)=nan;
        elseif ~isfield(SessionResults{1,j},'Stim_clickRate') %'clicks20'��������.����20
            if isfield(SessionResults{1,j},'Stim_Type')
            dataChoiceResult{i}(j,1)=  str2double(SessionResults{j}.Stim_Type(7:end));
%             dataChoiceResult{i}(j,1)=  str2double(SessionResults{j}.Stim_Type);
            else
                dataChoiceResult{i}(j,1)=nan;
            end
        else
            dataChoiceResult{i}(j,1)=nan;
        end
        if isnan(dataChoiceResult{i}(j,1))%���ʵ���Ҳ���stim��ֻ������������
            continue;
        end
    end
    clicktemp=unique(dataChoiceResult{i}(:,1));%��Ӧ���probe
    click=[];
    for iclick=1:length(clicktemp)%ɾ����������̫�ٵĵ㣬��߳����ݴ�����
        if sum(clicktemp(iclick)==dataChoiceResult{i}(:,1))>10 %С��10�����ݵ���Ϊ�Ǽ�¼���ݵ�����
            click=[click;clicktemp(iclick)];
        end
    end
        
    for j=1:length(SessionResults)
        freq_n=length(click);
        stimlr=sort(click);%asending sort clicks as left-right shifted stimuli
%         stimlr=1:freq_n;
        actionchoice=SessionResults{j}.Action_choice;%0��1�ң�2miss��3violation
        ind=find(click==dataChoiceResult{i}(j,1));
        if isempty(ind)%����ҵ�stim,�������趨�Ŀ���������trial
            dataChoiceResult{i}(j,:)=nan;%������Ϊnan
            continue;
        end
        %ͨ������ȷ��trial Ϊcor,err,miss,vio;�Լ��Ƿ���Ҫ�ı�̼���Ӧ��side
        if clickright==20
            if (dataChoiceResult{i}(j,1)<50 && actionchoice==1) ||(dataChoiceResult{i}(j,1)>50 && actionchoice==0) 
                dataChoiceResult{i}(j,2)=1;
            elseif (dataChoiceResult{i}(j,1)<50 && actionchoice==0) ||(dataChoiceResult{i}(j,1)>50 && actionchoice==1) 
                dataChoiceResult{i}(j,2)=2;
            else
                dataChoiceResult{i}(j,2)=actionchoice+1;%actionchoice=2,3,��miss����violation
            end
            dataChoiceResult{i}(j,3)=stimlr(freq_n+1-ind);
        else
            if (dataChoiceResult{i}(j,1)<50 && actionchoice==0) ||(dataChoiceResult{i}(j,1)>50 && actionchoice==1) 
                dataChoiceResult{i}(j,2)=1;
            elseif (dataChoiceResult{i}(j,1)<50 && actionchoice==1) ||(dataChoiceResult{i}(j,1)>50 && actionchoice==0) 
                dataChoiceResult{i}(j,2)=2;
            else
                dataChoiceResult{i}(j,2)=actionchoice+1;%actionchoice=2,3,��miss����violation
            end
            dataChoiceResult{i}(j,3)=stimlr(ind);
        end
        dataChoiceResult{i}(j,4)=SessionResults{j}.Delay_duration;
        dataChoiceResult{i}(j,5)=SessionResults{j}.Time_answer-SessionResults{j}.Time_delayOffset;%����violation trial�� time_answer=0����ֵΪ��ֵ��������Ϊ0����Ϊlick=delayOffsetʱҲ��violation
        if isfield(SessionSettings{j},'ExpCondition')&&isfield(SessionSettings{j},'animalName')%���ڰ汾��.matû��ExpCondition��field����Ĭ��Ϊ��
            prule=cellfun(@(x) strfind(SessionSettings{j}.animalName,x), curvetype.rule,'UniformOutput', false);
            psite=cellfun(@(x) strfind(SessionSettings{j}.ExpCondition,x), curvetype.site,'UniformOutput', false);
            prule=sum(~cellfun(@isempty,prule),2);%���н�3*2 cell�����3*2 logical����ͨ�����յڶ�ά�Ӻͣ���Ϊ1*2����
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
            dataChoiceResult{i}(j,7)=-1;%֮ǰ����Ϊnan�������ڴ���ʱ�Ƚ����Ѵ����������׳���nan��������Ϊ-1��������Ҫ���ܹ�����<0ѡ��������Ϊnan
        end
    end
    n=n+length(SessionResults);
end
dataChoiceResult{length(filenames)+1}=zeros(n,size(dataChoiceResult{1},2));
animal=unique(animal_name);
dataChoiceResult_animal=cell(length(animal)+1,1);%���һ��cell��ʾ�ܽ��
n=zeros(length(animal),1);%��ÿֻ����trial����
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

