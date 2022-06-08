close all;

rootpath='E:\FP\pyx288_20200922';
temp=strsplit(rootpath,'\');
session_name=temp(end);
cd(rootpath);
SavingFolder=rootpath;
behaviorFile=dir([rootpath,'\*Virables.mat']);
FrameInfo = dir([rootpath,'\*.log']);
fileID = fopen([rootpath,'\',FrameInfo.name]);
C=textscan(fileID,'%d %d','HeaderLines',16);
%     ImagingSetup
ImagingSetup=FrameInfo.name(1:end-4);
disp(ImagingSetup);
ImagingSetup(ImagingSetup=='_')='-';
fclose(fileID);
TrialCount=C{1,1};
TrialStart_FrameCount=C{1,2};
fiberstr={'Left','Right'};

% compute dff or just load
if exist('dff_temp.mat','file')
    load('dff_temp');
else
    dff=fGetFPdff(SavingFolder);
end

% load Beh mat data and extract behavior events
load(behaviorFile.name);%load behavior data
[trialType,behrule] = fGetTrialType( Data_extract);%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
FrameRate=40;
FrameTime=1000/FrameRate;
frT=FrameTime*2;%2 channel, so framerate and frametime should be half
frameNumTime=[2,5];%from 2s before align point to 5s after align point
frameNum=double(round(frameNumTime*1000/frT));
frameRate=FrameRate/2;%2 channels
%settings for plotting
behEventAlign='delay onset';%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
masklick='no';
% yrange=[ -0.01 , 0.02 ];
strData={'1b2m','1m2b','470','410'};
iData=1;
ind_tr_1=1;
ind_1stFrame=round(double(TrialStart_FrameCount')/2);
f=3;
for nfiber=1:2
    IDstr=strcat('fiber-',fiberstr{nfiber});
    [figRaster,figMeanTrace] = fPlotRasterPSTHoptoVSctrl(Data_extract,behEventAlign, masklick, behEventSort,dff{1,iData}(nfiber,:),IDstr,frameNumTime,ind_1stFrame,frameRate,ind_tr_1,f);
    saveas(figRaster,[rootpath,filesep,strData{iData},'raster-', fiberstr{nfiber},'.pdf'],'pdf');
    saveas(figMeanTrace,[rootpath,filesep,strData{iData},'meantrace-', fiberstr{nfiber},'.pdf'],'pdf');
    
end