function [A] = fPlotLickConsistency(path,date,animal,npast)
%FPLOTLICKCONSISTENCY plot licking consistency index(LCI) 
%Input-
%   npast- 6(default),how many LCI in the past would be plot in the figure
fileFolder=strcat(path,filesep,animal,'\');
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
%get LCI
[LCI,dataChoice,correct_ind,error_ind,currentSession] = fGetLCI(npast,fileFolder,'proportion',filename);
%plot LCI
A=figure;
set(gcf,'position', [100 100 800 300]);
subplot(1,2,1);%plot current LCI during session
% plot(dataChoice(:,3),'k-');
plot(find(correct_ind>0),dataChoice(correct_ind,3),'m-');
hold on;
plot(find(error_ind>0),dataChoice(error_ind,3),'g-');
legend('correct','error');
ylabel('licking consistency index');
xlabel('trial number');
title(currentSession);
set(gca,'FontName','Arial','FontSize',14);
subplot(1,2,2);%plot LCI across days,until date chosen
plot(LCI(:,1),'-k');
hold on;
plot(LCI(:,2),'-m');
hold on;
plot(LCI(:,3),'-g');
hold on;
title(strcat('Recent LCI-',animal));
legend('LCI','correct','error');
ylabel('licking consistency index');
xlabel('training day');
set(gca,'FontName','Arial','FontSize',14);
end

