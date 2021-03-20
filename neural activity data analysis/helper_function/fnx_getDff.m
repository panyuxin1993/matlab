function [ dff ,ind_1stFr] = fnx_getDff( path,savefolder,varargin )
%FNX_GETDFF get dff data from Ca...f_raw data 
%   Input- path and saving folder of potential plots; th
%   Output- dff
if ~exist('path','var')
    path=rootname;%'1-200trials';%'segNP';
end
if ~exist('savefolder','var')
    savefolder=rootname;
end
savestr='dff';
if ~isempty(varargin)
    if strcmp(varargin{1},'save figure')
        figHow='save';
    else
        figHow='see';
    end
    if length(varargin)>1 && strcmp(varargin{2},'field_NPseg')
        savestr='dff_NPseg';
    end
end
        
CurrFolder = cd(path);
dirmat=strcat(CurrFolder,'\*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
file_imaging=cellfun(@(x) contains(x,'Ca'), filenames);
i_file_imaging=find(file_imaging);
load([path,filesep,filenames{1,i_file_imaging}]);%load imaging data

if ~exist(savefolder)
    mkdir(savefolder);
end

f_cat  = [];
ind_1stFr(1) = 1;
ind_tr_1=1;
ntr = length(SavedCaTrials.f_raw); %cancate all trials together, or  just set the number of trials to use
if strcmp(savestr,'dff')    
    for i = ind_tr_1:ntr
        f_cat = [f_cat SavedCaTrials.f_raw{i}];
        ind_1stFr(i+1) = size(f_cat,2) + 1;
    end
    ind_1stFr(i+1) = [];
else strcmp(savestr,'dff_NPseg')
    if ~isfield(SavedCaTrials,'SegNPdataAll')
        warning([savefolder,'no field named SegNPdataAll, so skip it and set dff empty']);
        dff=[];
        return;
    end
    for i = ind_tr_1:ntr
        f_cat = [f_cat nanmean(SavedCaTrials.SegNPdataAll{i}([6,7,10,11],:),1)];
        ind_1stFr(i+1) = size(f_cat,2) + 1;
    end
    ind_1stFr(i+1) = [];
end
dff=zeros(size(f_cat));
for roiNo = 1:size(f_cat,1) %SavedCaTrials.nROIs may be not true
%%
%substract changing baseline
    totalFr = size(f_cat,2);
    frT = SavedCaTrials.FrameTime;
    span = round(20*1000/frT/2); %every 40s, get a mean and substrate to compensate for baseline change
    ind_x = 0;
    x = [];
    for i = 1:totalFr
        %     ind_x = ind_x + 1;
        ind1 = i- span;
        ind2 = i+span;
        if ind1 <= 0
            ind1 = 1;
        end
        if ind2 > totalFr
            ind2 = totalFr;
        end
        x(i) = prctile(f_cat(roiNo,ind1:ind2),5);
    end
    f = f_cat(roiNo,:) - x + mean(x);
    % figure; histogram(f_cat_sub(roiNo,:),100)   
    %%
    %get f mode as f0
    [N,edges,bin] = histcounts(f,100);
    f0 = edges(N == max(N));
    if length(f0)==2
        f0=mean(f0);
    end
    %get dff
    dff(roiNo,:) = (f - f0)/f0;
    if exist('figHow','var')
        B=figure;
        set(gcf, 'position', [0 0 1500 500]);
        plot(f_cat(roiNo,:),'k-');
        hold on;
        plot(x,'c-');
        plot(f,'g-');
        plot(dff(roiNo,:),'b-');
        legend('f cat','moving f mean','f baseline correction','dff');
        set(gca,'FontName','Arial','FontSize',14);
        if strcmp(figHow,'save')
            saveas(B,[CurrFolder,'\',savefolder,'\',savefolder,'ROI-',num2str(roiNo),'-cat_f.jpg'],'jpg');
        end
        close all;
    end
end
save(savestr,'dff');
end

