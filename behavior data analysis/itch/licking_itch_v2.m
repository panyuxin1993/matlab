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
lickTime=cell(length(filenames),1);%单个cell表示一个session的所有lick
figure;
% set(gcf, 'position', [0 0 1000 400]);%控制fig尺寸 
ntotal=length(filenames);
nfig=1;
divideTrialType=2;%2则一个session区分左右trialType而画成两部分,否则画在一起
resample=[];%'resample';%trial数量太多的话可以可考虑resample
rasterSize=2;%每个lick点的横向长度，因为分辨率为1ms小鼠不可能舔水这么密集，此值甚至可以到上百
for n=ntotal-nfig+1:ntotal
    load(strcat(fileFolder,'\',filenames{n}));
    [lickTime{n},maxdelay]=fGetLickTimes(SessionResults,'go');%align={'stim','go'};
    %排序
    if divideTrialType==1%一共一行fig，每个session一张fig
        set(gcf,'position',[0,0,200*nfig,400]);
        subplot(2,nfig,n-ntotal+nfig);%这种设定能画指定数量张图，无violation的1000ms和首次加上300ms的violation
        fPlotLickRasterOneSession(lickTime{n},[0 1],'withVio','raw',filenames{n},animal,maxdelay,resample,rasterSize);
        subplot(2,nfig,nfig+1:nfig*2);
    else
        set(gcf,'position',[0,0,200*nfig,600]);
        subplot(3,nfig,n-ntotal+nfig);%第一行画左边的trial
        fPlotLickRasterOneSession(lickTime{n},0,'withVio','raw',filenames{n},animal,maxdelay,resample,rasterSize);
        subplot(3,nfig,n-ntotal+nfig*2);%第二行画右边的trial
        fPlotLickRasterOneSession(lickTime{n},1,'withVio','raw',filenames{n},animal,maxdelay,resample,rasterSize);
        subplot(3,nfig,2*nfig+1:nfig*3);
    end
    [nchoiceL nchoiceR]=fGetAnsNum(SessionResults);
    [nlickL nlickR]=fGetLickNum(SessionResults);
    nchoice(i)=nchoiceL/nchoiceR;
    nlick(i)=nlickL/nlickR;
end

% fSavePPT(strcat(fileFolder,'\test.ppt'),filenames{n},'-f1');