function [A] = fPlotDiffDelay(ndelaygroup,path,date,animal)
%FPLOTDIFFDELAY plot animal performance for different delay length in one
%session. ������ͬdelay���ȵ���ȷ�ʣ�miss,violation����,reaction time
%Input- 
%   ndelaygroup- 4��default), how many bins used to divide whole delay
%   duration
%   date- in the form 2020_04_11
%   animal- name of animal

%clear;

if contains(animal,'pyx17')
    date(date=='_')='';%���»���ȡ����
    filename=strcat(path,filesep,animal,'\',animal,'_',date,'_2AFC.mat');
elseif contains(animal,'pyx1') ||contains(animal,'pyx273') ||contains(animal,'pyx26')||contains(animal,'pyx23')||contains(animal,'pyx241')
    filename=strcat(path,filesep,animal,'\',date,'_',animal,'.mat');
elseif contains(animal,'pyx231')
    date(date=='_')='';%���»���ȡ����
    filename=strcat(path,filesep,animal,'\','PYX231',32,date,'.mat'); %here 32 means ' '
else %imaging
    date(date=='_')='';%���»���ȡ����
    filename=strcat(path,filesep,animal,'\',animal,'_',date,'.mat');
end
% filename='F:\FP\pyx215_20190613\2019_06_13_pyx215_FP.mat';
date(date=='_')='';%���»���ȡ����

dirmat=strcat(path,filesep,animal,'\*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
sort(filenames);
indSession=cellfun(@(x) strcmp(x,filename),filenames);
i_session=find(indSession);
if isempty(i_session)% the specified session not exist
    i_session=length(filenames);%default is to plot the last session
    warning(['Specified session not exist, using last session: ',filenames{i_session}]);
end
load([path,filesep,animal,filesep,filenames{i_session}]);
dataChoice=zeros(length(SessionResults),4);%��1-3�зֱ��ʾ��ǰ�̼���Action_choice��Delay_duration��reaction time
if isfield(SessionResults{1},'Trial_isRetractTrial') && isfield(SessionResults{end},'Trial_isRetractTrial') %may change settings
    ind_retract=cellfun(@(x) x.Trial_isRetractTrial==1,SessionResults);%logical, 0-1 vector
    ind_retract=ind_retract';
else
    ind_retract=zeros(length(SessionResults),1);
end
n=1; %n��trial
if isfield(SessionSettings{1},'clickRate_')
    clickright=SessionSettings{1}.clickRate_;
elseif isfield(SessionSettings{1},'clickRate_R')
    clickright=SessionSettings{1}.clickRate_R;
end
for i=1:length(SessionResults)
    if  isfield(SessionResults{i},'Stim_clickRate') && ~isempty(SessionResults{i}.Stim_clickRate)
        if SessionResults{i}.Stim_clickRate<126%����125170500�������ͣ��ͺ���volume�ϲ���
            dataChoice(i,1)= SessionResults{i}.Stim_clickRate;
        else
            dataChoice(i,1)= floor(SessionResults{i}.Stim_clickRate/10^6);%����125170��������
        end
    elseif isfield(SessionResults{i},'Stim_clickRate') && isempty(SessionResults{i}.Stim_clickRate)
        continue;
    elseif ~isfield(SessionResults{i}  ,'Stim_clickRate') %'clicks20'��������
        dataChoice(i,1)=  str2double(SessionResults{i}.Stim_Type(7:end));
    end
    dataChoice(i,2)=SessionResults{i}.Action_choice;%0��1�ң�2miss��3violation
    dataChoice(i,3)=SessionResults{i}.Delay_duration;
    dataChoice(i,4)=SessionResults{i}.Time_answer-SessionResults{i}.Time_delayOffset;
end
if dataChoice(1,1)>125 %��Щ�ļ���¼��Bug,��֮���volume������ֱ������click rate���棬���¸����ֶ�����λ��
    dataChoice(:,1)=floor(dataChoice(:,1)/1000);
end
maxdelay=max(dataChoice(:,3));
mindelay=min(dataChoice(:,3));
part=(maxdelay-mindelay)/ndelaygroup;
if part==0
    warning('invariable delay');
end
delay=dataChoice(:,3)<=(mindelay+part);
for i=1:ndelaygroup
    delaycol=(dataChoice(:,3)>(mindelay+i*part)).*(dataChoice(:,3)<=(mindelay+i*part+part));
    delay=[delay, delaycol];%���ܳɶ�ά���飬ÿһ�д���ͬ���ȵ�delay
end
boundary=50;
%boundary=200000;
if clickright==20
    left=(dataChoice(:,1)>boundary);
    right=(dataChoice(:,1)<boundary);
else
    left=(dataChoice(:,1)<boundary);
    right=(dataChoice(:,1)>boundary);
end
all=left+right;
side=[all,left,right];%���ܳɶ�ά���飬ÿһ�д���ͬ��
correct=left.*(dataChoice(:,2)==0)+right.*(dataChoice(:,2)==1);
error=right.*(dataChoice(:,2)==0)+left.*(dataChoice(:,2)==1);
miss=(dataChoice(:,2)==2);
violation=(dataChoice(:,2)==3);
do=correct+error;
nomiss=do+violation;
indRT=correct.*(dataChoice(:,4)>0);%ֻ����ȷ��trial

performDiffDelay=zeros(ndelaygroup,4,3);%��һά����delay���ȣ��ڶ�ά������ȷ��/miss/violation/RT������ά�����ܣ�����
for i=1:ndelaygroup
    for j=1:3
        performDiffDelay(i,1,j)=100*sum(double(correct.*side(:,j).*delay(:,i)))/sum(double(do.*side(:,j).*delay(:,i)));
        performDiffDelay(i,2,j)=100*sum(double(miss.*side(:,j).*delay(:,i)))/sum(double(side(:,j).*delay(:,i)));
        performDiffDelay(i,3,j)=100*sum(double(violation.*side(:,j).*delay(:,i).*(ind_retract==0)))/sum(nomiss.*double(side(:,j).*delay(:,i).*(ind_retract==0)));
        %         performDiffDelay(i,3,j)=100*sum(double(violation.*side(:,j).*delay(:,i)))/sum(nomiss.*double(side(:,j).*delay(:,i)));
        temp=double(dataChoice(:,4).*indRT.*side(:,j).*delay(:,i));
        temp2=temp;
        temp2(temp==0)=nan;
        performDiffDelay(i,4,j)=nanmedian(temp2);
    end
end

A=figure;
hold on;
%suptitle(animal_name);
%ha = tight_subplot(1,2,[1 1],[1 1],[1 1]);
set(gcf, 'position', [50 50 1600 400]);%����fig�ߴ�
%��ɫͳһ�������ɫ���ұߺ�ɫ������ɫ
%�����ע�����20-�ұ�125�����ĵ㣩�����125-�ұ�20��ʵ�ĵ㣩
% suptitle(strcat(animal,'-',date));
y_label={'% Correct';'% Miss';'%Violation';'reaction time(ms)'};%{}��cell���ͣ�������־��󳤶Ȳ�������⡣
% x_axis={strcat(num2str(mindelay),' - ',num2str(mindelay+part));...
%     strcat(num2str(mindelay+part),' - ',num2str(mindelay+2*part));...
%     strcat(num2str(mindelay+2*part),' - ',num2str(maxdelay));...
%     };
if mindelay==maxdelay
    x_plot=mindelay*ones(ndelaygroup,1);
else
    x_plot=mindelay+part/2:part:maxdelay-part/2;
end
for i=1:4
    subplot(1,4,i);
    H1=plot(x_plot,performDiffDelay(:,i,1),'k','LineWidth',2);	%����ȷ�ʣ�����miss
    %set(gca,'XTickLabel',x);     %������char�͵�ֵ��Ϊ�����꣬��Ӧ���Ȼ��������Ѻ�����ĳ��ַ���
    hold on;
    H2=plot(x_plot,performDiffDelay(:,i,2),'b','LineWidth',2);    %�����ȷ�ʣ�����miss
    hold on;
    H3=plot(x_plot,performDiffDelay(:,i,3),'r','LineWidth',2);    %�ұ���ȷ�ʣ�����miss
    hold on;
    H=[H1,H2,H3];
    if clickright==20 %��ʵ�ĵ�
        scatter(x_plot,performDiffDelay(:,i,1),50,'k','filled');%ʵ�ĵ�
        hold on;
        scatter(x_plot,performDiffDelay(:,i,2),50,'b','filled');
        hold on;
        scatter(x_plot,performDiffDelay(:,i,3),50,'r','filled');
        hold on;
    else
        scatter(x_plot,performDiffDelay(:,i,1),50,'k');%���ĵ�
        hold on;
        scatter(x_plot,performDiffDelay(:,i,2),50,'b');
        hold on;
        scatter(x_plot,performDiffDelay(:,i,3),50,'r');
        hold on;
    end
    if i<4
        plot([mindelay,maxdelay],[50,50],'--k','LineWidth',1);
        hold on;
        set(gca, 'YLim',[0 100]);
        ylimmax=100;
    else
        ylimrange=get(gca,'Ylim');
        ylimmax=ylimrange(end);
    end
    for j=1:ndelaygroup
        plot([mindelay+j*part,mindelay+j*part],[0,ylimmax],':k','LineWidth',1);
        hold on;
    end
    xlabel('Delay(ms)','FontName','Arial','FontSize',14);
    ylabel(y_label(i),'FontName','Arial','FontSize',14);
    box off; %ȡ��ͼ����߿�
    if mindelay==maxdelay
        set(gca, 'XLim',[mindelay-1 maxdelay+1]);
        set(gca,'XTick',mindelay,'XTickLabel',num2str(mindelay));%��������λ���ľ���
    else
        set(gca, 'XLim',[mindelay maxdelay]);
        set(gca,'XTick',mindelay:part:maxdelay);%x��ֻ��ע�����߽��
        set(gca,'XTickLabel',num2str(round(get(gca,'xTick')'/100)*100,'% .f'));%��������λ���ľ���
    end
    hl=legend(H,'Total','Left','Right','Location','Best');%legend��Ҫ��Plot���棬����֮���plot�����Զ�����legend
    set(hl,'Box','off');%�ص�legend�ı߿�
    set(gca,'FontName','Arial','FontSize',14);
end

end

