
%decide some global variable
behEventAlign='delay onset';%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
masklick='no';
frameNumTime=[1,3.5];%from 5s before align point to 5s after align point
% yrange=[ -0.5 , 1 ];
i_selectivity=4;%*********variable**************
selectivitystr={'stimuli','difficulty','sensory','choice'};%sensory means grouping difficulties;

cd 'D:\pyx340_20210205\im_data_reg\result_save';
CurrFolder = pwd;
rootname='pyx340_20210205';
savefolder=rootname;%'1-200trials';%'segNP';
% load([CurrFolder,'\',rootname, '_imaging_Virables.mat']);%load behavior data
dirmat=strcat(CurrFolder,'\*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
file_imaging=cellfun(@(x) contains(x,'Ca'), filenames);
i_file_imaging=find(file_imaging);
load(filenames{1,i_file_imaging});%load imaging data
file_beh=cellfun(@(x) contains(x,'imaging_Virables'), filenames);
i_file_beh=find(file_beh);
load(filenames{1,i_file_beh});%load behavior data

if ~exist(savefolder,'dir')
    mkdir(savefolder);
end
if exist('dff.mat','file')
    load([CurrFolder,'\','dff.mat']);%load dff
else
    dff=fnx_getDff(CurrFolder,rootname,'save figure');
end

ind_tr_1=1;%using data from trial 1
ntr = length(SavedCaTrials.f_raw); %use all data, or  just set the number of trials to use
frT = SavedCaTrials.FrameTime;
frameRate=1000/frT;
% align to behavior event
nFrameEachTrial=cellfun(@(x) size(x,2),SavedCaTrials.f_raw);
ind_1stFrame=zeros(1,length(nFrameEachTrial));
ind_1stFrame(1)=1;
ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
ind_1stFrame=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials

for roiNo = 1:size(SavedCaTrials.f_raw{1},1)
    IDstr=strcat('ROI-',num2str(roiNo));
    [figRaster,figMeanTrace] = fPlotRasterPSTHoptoVSctrl(Data_extract,behEventAlign, masklick, behEventSort,dff(roiNo,:),IDstr,frameNumTime,ind_1stFrame,frameRate,ind_tr_1,i_selectivity);
    saveas(figRaster,[CurrFolder,filesep,rootname,filesep,rootname,'ROI-',num2str(roiNo),'-algin to ',behEventAlign,'-sort ',behEventSort,'-cmpOptoCtrl-colorplot.png'],'png');
    saveas(figMeanTrace,[CurrFolder,filesep,rootname,filesep,rootname,'ROI-',num2str(roiNo),'-algin to ',behEventAlign,'-sort ',behEventSort,'-cmpOptoCtrl-mean trace.png'],'png');
%      saveas(figMeanTrace,[CurrFolder,filesep,rootname,filesep,rootname,'ROI-',num2str(roiNo),'-algin to ',behEventAlign,'-sort ',behEventSort,'-cmpOptoCtrl-mean trace.pdf'],'pdf');
    
    close all;
end