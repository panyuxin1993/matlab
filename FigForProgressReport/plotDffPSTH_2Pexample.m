%continuing sampling, so link each trials one by one to get a 'long trial'
% load('F:\pyx083-20180522\im_data_reg_cpu\result_save\CaTrialsSIM_pyx083_20180522_920nm_power50_zoom4x_dftReg_.mat');
% load('D:\xulab\behavior\pyx083\Data_Virables\2018_05_22_pyx083-imaging_Virables.mat');
clear;
close all;
%plot data from summary xlsx file, batch analysis
%{
[num,txt,raw] =xlsread('D:\xulab\project\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
T=cell2table(raw(2:end,1:14));
T.Properties.VariableNames=strrep(raw(1,1:14),' ','_');%table variable name can't have ' ',so replace them
ind_session=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control');
ind_session=find(ind_session);
n_session=length(ind_session);
animal_unique=unique(T.animal(ind_session));
n_animal=length(animal_unique);
[sessionNamePool,pathPool,trial2includecell]=deal(cell(n_session,1));
for i=1:n_session
    sessionNamePool{i}=T.session{ind_session(i)};
    pathPool{i}=T.file_path{ind_session(i)};
    trial2includecell{i}=T.used_trial{ind_session(i)};
end
% trial2include=[1,64;69,82;93,215];
trial2exclude=[];
%}

%plot single session
%
% sessionNamePool={'pyx339_20210115'};setROI=22;%vlgut2 example
% sessionNamePool={'pyx314_20200907'};setROI=50;%vlgut2 example
% sessionNamePool={'pyx339_20210117'};setROI=34;%vlgut2 example
% sessionNamePool={'pyx338_20201219'};setROI=25;%vgat example, ipsi delay, contra lick
% sessionNamePool={'pyx311_20200806'};setROI=35;%vgat example, lick
sessionNamePool={'pyx338_20201219'};setROI=34;%vgat example, ipsi delay, contra lick

pathPool=cellfun(@(x) ['H:\2P\',x,'\im_data_reg\result_save'], sessionNamePool, 'UniformOutput',false);
trial2includecell={nan};
trial2exclude=[];
%}

behEventAlignPool={'delay onset','go cue'};
masklickPool={'yes','no'};
xlabelstr={'Time(s) from delay onset','Time(s) from go cue'};
% behEventAlignPool={'go cue'};
% masklickPool={'no'};

figCase=figure;
set(gcf,'Position',[100,100,200*length(behEventAlignPool),200]);
for i_align=1:length(behEventAlignPool)
    subplot(1,length(behEventAlignPool),i_align);
    hold on;
    %decide some global variable
    behEventAlign=behEventAlignPool{i_align};%align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
    behEventSort='go cue';% string can be in{'first lick','reward','go cue'};
    masklick=masklickPool{i_align};
    
    i_selectivity=4;%*********variable**************
    selectivitystr={'stimuli','sensory difficulty','sensory','choice'};%sensory means grouping difficulties;
    trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
    trialTypeStr=selectivitystr{i_selectivity};
    %plotting settings
    if strcmp(behEventAlign,'stim onset')
        frameNumTime=[0.5,2];%from 5s before align point to 5s after align point
    elseif strcmp(behEventAlign,'delay onset')
        frameNumTime=[1,1.4];%from 5s before align point to 5s after align point
    else
        frameNumTime=[1.8,3.5];%from 5s before align point to 5s after align point
    end
    ind_tr_1=1;%using data from trial 1
    
    %% plot individual cells' PSTH
    %
    for i_session=1:length(sessionNamePool)
        sessionName=sessionNamePool{i_session};
%         datapath=['D:\',sessionName];
        
        % yrange=[ ];
        clear CurrFolder;
%         cd([datapath,'\im_data_reg\result_save']);
        cd(pathPool{i_session});
        CurrFolder=pwd;
        savefolder=sessionName;%'1-200trials';%'segNP';
        % savefolder='F:\2P\example\';
        
        % load([CurrFolder,'\',sessionName, '_imaging_Virables.mat']);%load behavior data
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
        
        if ~exist(savefolder)
            mkdir(savefolder);
        end
        
        if exist('dff.mat','file')
            load([CurrFolder,'\','dff.mat']);%load dff
        else
            dff=fnx_getDff(CurrFolder,sessionName,'save figure');
        end
        for roiNo = setROI%1:size(SavedCaTrials.f_raw{1},1) %SavedCaTrials.nROIs may be not true
            disp(strcat('processing...',sessionName,'-ROI',num2str(roiNo)));
            savename_figdff=[CurrFolder,'\',savefolder,'\',savefolder,'ROI-',num2str(roiNo),'-cat_f_label_used_range.jpg'];
            savename_fig=[CurrFolder,'\',savefolder,'\',sessionName,'ROI-',num2str(roiNo),'-algin to ',behEventAlign,'-sort ',behEventSort,'colorplot.jpg'];
            curve_meanTrace=fPlotDffPSTH_ASession(dff(roiNo,:),ind_tr_1,Data_extract,...
                SavedCaTrials,frameNumTime,behEventAlign,masklick,...
                i_selectivity,trial2includecell{i_session},trial2exclude,...
                savename_figdff,['ROI-',num2str(roiNo)]);
        end
    end
    h=legend(curve_meanTrace,'ipsi','contra');
    set(h,'box','off');
    xlabel(xlabelstr{i_align});
    ylabel('\it\DeltaF/F');
    set(gca,'FontSize',10);
    %}
        %% plot field dff
    %{
    figure;
    ncol=4;
    nrow=ceil(length(sessionNamePool)/4);
    
%     savename_fig=['H:\2P\summary\field_dff\',sessionName,'field_dff-algin to ',behEventAlign,'-sort ',behEventSort,'colorplot.jpg'];

    for i_session=1:length(sessionNamePool)
        sessionName=sessionNamePool{i_session};
%         datapath=['D:\',sessionName];
        % yrange=[ ];
        clear CurrFolder;
%         cd([datapath,'\im_data_reg\result_save']);
        cd(pathPool{i_session});
        CurrFolder=pwd;
        savefolder='H:\2P\summary\';
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
        if exist('dff_NPseg.mat','file')
            load([CurrFolder,filesep,'dff_NPseg.mat']);%load dff
        else
            dff=fnx_getDff(CurrFolder,sessionName,'save figure','field_NPseg');
        end
        titlestr=strrep(sessionName,'_','-');
        savename_figdff=[savefolder,'field_dff-cat_f_label_used_range.jpg'];
        subplot(nrow,ncol,i_session);
        curve_meanTrace=fPlotDffPSTH_ASession(dff,ind_tr_1,Data_extract,...
            SavedCaTrials,frameNumTime,behEventAlign,masklick,...
            i_selectivity,trial2includecell{i_session},trial2exclude,...
            savename_figdff,titlestr);
        if ceil(i_session/ncol)==nrow
            xlabel(['Time (s) from ',behEventAlign],'FontSize',12);
        end
        if mod(i_session,ncol)==1
            ylabel('\it\DeltaF/F');
        end
    end

    subplot(nrow,ncol,i_session);
    if contains(trialTypeStr,'stimuli')
        %         h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','','contra hard','','contra easy'},'Location','best');
        h=legend(curve_meanTrace(:),num2str(Data_extract.Stimuli(:)),'Location','best');
    elseif contains(trialTypeStr,'difficulty')
        h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','contra hard','contra easy'},'Location','best');
    elseif contains(trialTypeStr,'choice')
        h=legend(curve_meanTrace(:),{'ipsi choice','contra choice'},'Location','best');
    elseif contains(trialTypeStr,'sensory')
        h=legend(curve_meanTrace(:),{'low click rate','high click rate'},'Location','best');
    end
    set(h,'box','off');
    clear curve_meanTrace curve_lick_meanTrace;
    %}
end


