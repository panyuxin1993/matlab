%plot lick density during delay when no violation(allow violation)
%get lick time
clear;
date='20180708';
animal='practice\LTY02';
fileFolder=strcat('D:\xulab\behavior\',animal);
% session=strcat('\',date,'_',animal,'.mat');
% dirmat=strcat(fileFolder,session);
dirmat=strcat(fileFolder,'\*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
sort(filenames);
lickTime=cell(length(filenames),1);%����cell��ʾһ��session������lick
figure;
% set(gcf, 'position', [0 0 1000 400]);%����fig�ߴ� 
ntotal=length(filenames);
nfig=1;
divideTrialType=2;%2��һ��session��������trialType������������,������һ��
resample=[];%'resample';%trial����̫��Ļ����Կɿ���resample
rasterSize=2;%ÿ��lick��ĺ��򳤶ȣ���Ϊ�ֱ���Ϊ1msС�󲻿�����ˮ��ô�ܼ�����ֵ�������Ե��ϰ�
for n=ntotal-nfig+1:ntotal
    load(strcat(fileFolder,'\',filenames{n}));
    [lickTime{n},maxdelay]=fGetLickTimes(SessionResults,'go');%align={'stim','go'};
    %����
    if divideTrialType==1%һ��һ��fig��ÿ��sessionһ��fig
        set(gcf,'position',[0,0,200*nfig,400]);
        subplot(2,nfig,n-ntotal+nfig);%�����趨�ܻ�ָ��������ͼ����violation��1000ms���״μ���300ms��violation
        fPlotLickRasterOneSession(lickTime{n},[0 1],'withVio','raw',filenames{n},animal,maxdelay,resample,rasterSize);
        subplot(2,nfig,nfig+1:nfig*2);
    else
        set(gcf,'position',[0,0,200*nfig,600]);
        subplot(3,nfig,n-ntotal+nfig);%��һ�л���ߵ�trial
        fPlotLickRasterOneSession(lickTime{n},0,'withVio','raw',filenames{n},animal,maxdelay,resample,rasterSize);
        subplot(3,nfig,n-ntotal+nfig*2);%�ڶ��л��ұߵ�trial
        fPlotLickRasterOneSession(lickTime{n},1,'withVio','raw',filenames{n},animal,maxdelay,resample,rasterSize);
        subplot(3,nfig,2*nfig+1:nfig*3);
    end
    [nchoiceL nchoiceR]=fGetAnsNum(SessionResults);
    [nlickL nlickR]=fGetLickNum(SessionResults);
    nchoice(i)=nchoiceL/nchoiceR;
    nlick(i)=nlickL/nlickR;
end

% fSavePPT(strcat(fileFolder,'\test.ppt'),filenames{n},'-f1');