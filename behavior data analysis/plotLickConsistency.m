%plot lick consistency trial by trial in a session
%plot LCI across sessions until date chosen
%2 method of LCI- 1)proportion of one side licking divided by all licking;
%   2)times of change of licking side
clear;
date='2020_04_12';
animal='pyx239';
fileFolder=strcat('D:\xulab\behavior\',animal,'\');
npast=6;%plot 10 past LCI across day
if contains(animal,'pyx28')||strcmp(animal,'pyx196')||contains(animal,'pyx24')||contains(animal,'pyx25')||contains(animal,'pyx26')
    date(date=='_')='';%将下划线取消掉
    filename=strcat(animal,'_',date,'.mat');
elseif strcmp(animal,'PYX121') ||strcmp(animal,'PYX119')
    date(date=='_')='';%将下划线取消掉
    filename=strcat(animal,'_',date,'_2AFC.mat');
else
    filename=strcat(date,'_',animal,'.mat');
end
date(date=='_')='';%将下划线取消掉
strmethod={'proportion','switch'};
A=figure;    
set(gcf,'position', [100 100 600 400]);
for imethod=1:2
    %get LCI
    [LCI,dataChoice,correct_ind,error_ind,currentSession] = fGetLCI(npast,fileFolder,strmethod{imethod},filename);
    %plot LCI
    %ylabelstr=strcat('licking consistency index(',strmethod{imethod},')');
    ylabelstr=strcat('LCI(',strmethod{imethod},')');
    figure(A);
    subplot(2,2,imethod*2-1);%plot current LCI during session
    % plot(dataChoice(:,3),'k-');
    plot(find(correct_ind>0),dataChoice(correct_ind,3),'r-');
    hold on;
    plot(find(error_ind>0),dataChoice(error_ind,3),'b-');
    legend('correct','error');
    ylabel(ylabelstr);
    xlabel('trial number');
    title(currentSession);
    set(gca,'FontName','Arial','FontSize',14);
    subplot(2,2,imethod*2);%plot LCI across days,until date chosen
    plot(LCI(:,1),'-k');
    hold on;
    plot(LCI(:,2),'-r');
    hold on;
    plot(LCI(:,3),'-b');
    hold on;
    title(strcat('Recent LCI-',animal));
    legend('LCI','correct','error');
    ylabel(ylabelstr);
    xlabel('training day');
    set(gca,'FontName','Arial','FontSize',14);
end