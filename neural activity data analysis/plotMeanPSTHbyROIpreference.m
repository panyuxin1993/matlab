clear meanActivityByTrialType_combine cellActivityByTrialType_combine meanActivityByTrialType_combine_cor cellActivityByTrialType_combine_cor;
close all;
[num,txt,raw] =xlsread('C:\Users\PYX\Documents\DataSummary\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
T=cell2table(raw(2:end,1:15));
T.Properties.VariableNames=strrep(raw(1,1:15),' ','_');%table variable name can't have ' ',so replace them
%for soma 
ind_session=logical(strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*(~contains(T.behavior_performance,'probe')).*contains(T.ROI_type,'soma').*(~contains(T.field,'soma')));%not probe session
celltypePool={'M2','syn','vglut2','vgat'};%{'syn','vglut2','vgat'};
info_str='soma';

% %for dendrite
% ind_session=logical(strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*(~contains(T.behavior_performance,'probe')).*contains(T.ROI_type,'dendrite'));%not probe session
% celltypePool={'vglut2','vgat'};
% info_str='dendrite';

% %for soma corresponding to the dendrite
% ind_session=logical(strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*(~contains(T.behavior_performance,'probe')).*contains(T.ROI_type,'soma').*contains(T.field,'soma'));%not probe session
% celltypePool={'vglut2','vgat'};
% info_str='soma_near_dendrite';

% % %for spine
% ind_session=logical(strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*(~contains(T.behavior_performance,'probe')).*contains(T.ROI_type,'spine'));%not probe session
% celltypePool={'vglut2','vgat'};
% info_str='spine';


behEventAlignPool={'delay onset','delay onset','delay onset','go cue','first lick'};
epoch4preferPool={'delay','mid_delay','late_delay','response','lick'};%{'sound','delay','mid_delay','late_delay','response','lick'}
masklickPool={'yes','yes','yes','no','no'};
i_selectivity=3;%3-choice
behEventSortPool={'go cue','go cue','go cue','first lick','go cue'};
trial2include='all';
trial2exclude=[];
activity_typePool={'dff','spkr','zscored_dff','zscored_spkr'};%{'dff','spkr','zscored_dff','zscored_spkr'};

%
for i_celltype=1:length(celltypePool)
    ind_session2=logical(ind_session.*logical(strcmp(T.cell_type,celltypePool{i_celltype})+strcmp(T.cell_type,[celltypePool{i_celltype},'-flpo'])));
    Tchoose=T(ind_session2,:);
    % file_path='H:\2P\pyx349_20210418\im_data_reg\result_save';
    % session='pyx349_20210418';
    ind_ROI=[];
    %activity_typePool={'dff'};%{'dff','spkr'};
    frT=0;
    [cellActivityByTrialType_combine2show,meanActivityByTrialType_combine2show]=deal([]);
    for i_align=1:length(behEventAlignPool)
        behEventAlign=behEventAlignPool{i_align};
        epoch4prefer=epoch4preferPool{i_align};
        masklick=masklickPool{i_align};
        for i_activity_type=1:length(activity_typePool)
            activity_type=activity_typePool{i_activity_type};
            clear meanActivityByTrialType_combine;
            savenamestr=['E:\2P\summary\meanActivity_by_ROIpreference',filesep,celltypePool{i_celltype},'-',info_str,'-meanActivityByROI-align',behEventAlign,'-masklick',masklick,'-',activity_type,'.mat'];
            if false%exist(savenamestr,'file')
                load(savenamestr);
            else
                clear meanActivityByTrialType_combine_cor cellActivityByTrialType_combine_cor TActivityByTrialType_combine;
                for i=1:size(Tchoose,1)
                    file_path=Tchoose.file_path{i};
                    session=[Tchoose.session{i},'_',Tchoose.field{i}];
                    savepath=[file_path,filesep,session];
                    objsession=Session2P(session,file_path,trial2include,trial2exclude);
                    trialTypeStr='cor and err';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
                    AUCtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
                    AUCCorrectedMethod='SensoryChoiceOrthogonalSubtraction';%'None';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
                    [TAUC, Tmean,objsession] = objsession.mGetEpochAUC(Tchoose.cell_type,Tchoose.ROI_type{i},trialTypeStr,AUCtype,AUCCorrectedMethod);
                    [meanActivityByTrialType,cellActivityByTrialType,Tout]=objsession.mMeanPSTHbyROI(activity_type,[],behEventAlign,masklick,i_selectivity);
                    TAUC.date=cellfun(@(x) datestr(datetime(x,'InputFormat','yyyy/MM/dd'),'yyyymmdd'), TAUC.date,'UniformOutput',false);
                    TActivityByTrialType=join(Tout,TAUC);
                    %combine data together and plot population
                    if exist('meanActivityByTrialType_combine','var')
                        frT=max(objsession.metadata.frT,frT);
                        frN=max(size(meanActivityByTrialType_combine{1,1},2),size(meanActivityByTrialType{1,1},2));
                        meanActivityByTrialType_combine= fInter(meanActivityByTrialType_combine,frN);
                        meanActivityByTrialType=fInter(meanActivityByTrialType,frN);
                        meanActivityByTrialType_combine_cor=fInter(meanActivityByTrialType_combine_cor,frN);
                        disp(['frN:',num2str(frN),',meanActivityByTrialType_combine:',num2str(size(meanActivityByTrialType_combine{1,1},2)),',meanActivityByTrialType:',num2str(size(meanActivityByTrialType{1,1},2))]);
                        meanActivityByTrialType_combine=cellfun(@(x,y) vertcat(x,y),meanActivityByTrialType_combine,meanActivityByTrialType,'UniformOutput',false);
                        cellActivityByTrialType_combine=cellfun(@(x,y) vertcat(x,y),cellActivityByTrialType_combine,cellActivityByTrialType,'UniformOutput',false);
                        meanActivityByTrialType_combine_cor=cellfun(@(x,y) vertcat(x,y),meanActivityByTrialType_combine_cor,meanActivityByTrialType(:,1),'UniformOutput',false);
                        cellActivityByTrialType_combine_cor=cellfun(@(x,y) vertcat(x,y),cellActivityByTrialType_combine_cor,cellActivityByTrialType(:,1),'UniformOutput',false);
                    else
                        meanActivityByTrialType_combine=meanActivityByTrialType;
                        cellActivityByTrialType_combine=cellActivityByTrialType;
                        meanActivityByTrialType_combine_cor=meanActivityByTrialType(:,1);%choose only the correct trials
                        cellActivityByTrialType_combine_cor=cellActivityByTrialType(:,1);
                    end
                    %combine table together
                    if exist('TActivityByTrialType_combine','var')
                        TActivityByTrialType_combine=vertcat(TActivityByTrialType_combine,TActivityByTrialType);
                    else
                        TActivityByTrialType_combine=TActivityByTrialType;
                    end
                end
            end
            %plot raster and mean using ROIs from different sessions and
            %group them together
            save(savenamestr,'meanActivityByTrialType_combine','frT','cellActivityByTrialType_combine','meanActivityByTrialType_combine_cor','cellActivityByTrialType_combine_cor','TActivityByTrialType_combine');
%             %reorganized the activities based on ttest
%             if strcmp(behEventAlign,'stim onset')
%                 frameNumTime=[0.5,2];%from 5s before align point to 5s after align point
%                 time4prefer=[-0.8,1.5]; %time for comparing ipsi vs. contra
%             elseif strcmp(behEventAlign,'delay onset')
%                 frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
%                 time4prefer=[-0.3,1];
%             else
%                 frameNumTime=[1,2];%from 5s before align point to 5s after align point
%                 time4prefer=[0,1];
%             end
%             frameNum=floor(frameNumTime*1000/frT);
%             frame4prefer=floor(time4prefer*1000/frT);
%             ind4prefer=false(sum(frameNum)+1,1);
%             ind4prefer(frameNum(1)-frame4prefer(1):frameNum(1)+frame4prefer(2)+1)=true;
%             meanActivityReorganized = fReorganizeByPreference(cellActivityByTrialType_combine,meanActivityByTrialType_combine,ind4prefer);
            %reorganized the activities based on AUC
            if strcmp(behEventAlign,'stim onset')
                frameNumTime=[0.5,2];%from 5s before align point to 5s after align point
                %epoch4prefer='sound'; %time for comparing ipsi vs. contra
            elseif strcmp(behEventAlign,'delay onset')
                frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
                %epoch4prefer='delay';%{'mid_delay','late_delay'}
            elseif strcmp(behEventAlign,'go cue')
                frameNumTime=[1,2];%from 5s before align point to 5s after align point
                %epoch4prefer='response';
            elseif strcmp(behEventAlign,'first lick')
                frameNumTime=[1,2];%from 5s before align point to 5s after align point
                %epoch4prefer='lick';
            end
            meanActivityReorganized = fReorganizeByAUC(meanActivityByTrialType_combine,TActivityByTrialType_combine,epoch4prefer);
            %plot
            [figRaster,figMean]=fPlotRasterMean(meanActivityReorganized,meanActivityByTrialType_combine,behEventAlign,frameNumTime,frT);
            figure(figRaster);
            suptitle(strrep([celltypePool{i_celltype},'-align',behEventAlign,'-masklick',masklick,'-',activity_type,'-epoch4prefer_',epoch4prefer],'_','\_'));
            saveas(figRaster,['E:\2P\summary\meanActivity_by_ROIpreference',filesep,celltypePool{i_celltype},'-align',behEventAlign,'-masklick',masklick,'-',activity_type,'-epoch4prefer_',epoch4prefer,'_raster.pdf'],'pdf');
            saveas(figRaster,['E:\2P\summary\meanActivity_by_ROIpreference',filesep,celltypePool{i_celltype},'-align',behEventAlign,'-masklick',masklick,'-',activity_type,'-epoch4prefer_',epoch4prefer,'_raster.png'],'png');
            figure(figMean);
            suptitle(strrep([celltypePool{i_celltype},'-align',behEventAlign,'-masklick',masklick,'-',activity_type,'-epoch4prefer_',epoch4prefer],'_','\_'));
            saveas(figMean,['E:\2P\summary\meanActivity_by_ROIpreference',filesep,celltypePool{i_celltype},'-align',behEventAlign,'-masklick',masklick,'-',activity_type,'-epoch4prefer_',epoch4prefer,'_meanTrace.pdf'],'pdf');
            saveas(figMean,['E:\2P\summary\meanActivity_by_ROIpreference',filesep,celltypePool{i_celltype},'-align',behEventAlign,'-masklick',masklick,'-',activity_type,'-epoch4prefer_',epoch4prefer,'_meanTrace.png'],'png');
        end
    end
end
%}
%plot different alignment together, aiming to show consistent and
%shifted selectivities
%
behEventAlignPool2show={'delay onset','first lick'};
epoch4preferPool2show={'late_delay','lick'};%{'sound','delay','mid_delay','late_delay','response','lick'}
masklickPool={'yes','no'};
for i_celltype=1:length(celltypePool)
    ind_session2=logical(ind_session.*logical(strcmp(T.cell_type,celltypePool{i_celltype})+strcmp(T.cell_type,[celltypePool{i_celltype},'-flpo'])));
    Tchoose=T(ind_session2,:);
    % file_path='H:\2P\pyx349_20210418\im_data_reg\result_save';
    % session='pyx349_20210418';
    ind_ROI=[];
    activity_typePool={'dff','spkr','zscored_dff','zscored_spkr'};%{'dff','spkr'};
    frT=0;
    for i_activity_type=1:length(activity_typePool)
        activity_type=activity_typePool{i_activity_type};
        [cellActivityByTrialType_combine2show,meanActivityByTrialType_combine2show]=deal([]);
        for i_align=1:length(behEventAlignPool2show)
            behEventAlign=behEventAlignPool2show{i_align};
            masklick=masklickPool{i_align};

            clear meanActivityByTrialType_combine;
            savenamestr=['E:\2P\summary\meanActivity_by_ROIpreference',filesep,celltypePool{i_celltype},'-',info_str,'-meanActivityByROI-align',behEventAlign,'-masklick',masklick,'-',activity_type,'.mat'];

            if exist(savenamestr,'file')
                load(savenamestr);
            else
                clear meanActivityByTrialType_combine_cor cellActivityByTrialType_combine_cor TActivityByTrialType_combine;
                for i=1:size(Tchoose,1)
                    file_path=Tchoose.file_path{i};
                    session=[Tchoose.session{i},'_',Tchoose.field{i}];
                    savepath=[file_path,filesep,session];
                    objsession=Session2P(session,file_path,trial2include,trial2exclude);
                    [meanActivityByTrialType,cellActivityByTrialType]=objsession.mMeanPSTHbyROI(activity_type,[],behEventAlign,masklick,i_selectivity);
                    %combine data together and plot population
                    if exist('meanActivityByTrialType_combine_cor','var')
                        frT=max(objsession.metadata.frT,frT);
                        frN=max(size(meanActivityByTrialType_combine_cor{1,1},2),size(meanActivityByTrialType{1,1},2));
                        meanActivityByTrialType_combine_cor= fInter(meanActivityByTrialType_combine_cor,frN);
                        meanActivityByTrialType=fInter(meanActivityByTrialType,frN);
                        disp(['frN:',num2str(frN),',meanActivityByTrialType_combine_cor:',num2str(size(meanActivityByTrialType_combine_cor{1,1},2)),',meanActivityByTrialType:',num2str(size(meanActivityByTrialType{1,1},2))]);
                        meanActivityByTrialType_combine_cor=cellfun(@(x,y) vertcat(x,y),meanActivityByTrialType_combine_cor,meanActivityByTrialType(:,1),'UniformOutput',false);
                        cellActivityByTrialType_combine_cor=cellfun(@(x,y) vertcat(x,y),cellActivityByTrialType_combine_cor,cellActivityByTrialType(:,1),'UniformOutput',false);
                    else
                        meanActivityByTrialType_combine_cor=meanActivityByTrialType(:,1);%choose only the correct trials
                        cellActivityByTrialType_combine_cor=cellActivityByTrialType(:,1);
                    end
                end
            end
            save(savenamestr,'meanActivityByTrialType_combine','frT','cellActivityByTrialType_combine','meanActivityByTrialType_combine_cor','cellActivityByTrialType_combine_cor','TActivityByTrialType_combine');
             cellActivityByTrialType_combine2show=cat(2,cellActivityByTrialType_combine2show,cellActivityByTrialType_combine_cor);
            meanActivityByTrialType_combine2show=cat(2,meanActivityByTrialType_combine2show,meanActivityByTrialType_combine_cor);
        end
       
        %plot different alignment together, aiming to show consistent and
        %shifted selectivities
        %using ttest to determine the significance
%         [meanActivityReorganized,indIpsi,indContra,indNS] = fReorganizeByMultiPreference(cellActivityByTrialType_combine2show,meanActivityByTrialType_combine2show,behEventAlignPool2show,frT);
        %using AUC to determine the significance
        [meanActivityReorganized,indIpsi,indContra,indNS] = fReorganizeByMultiAUC(meanActivityByTrialType_combine2show,TActivityByTrialType_combine,behEventAlignPool2show,epoch4preferPool2show);
        [figRaster,figMean]=fRasterMeanMultiAlign(meanActivityReorganized,behEventAlignPool2show,frT,activity_type);
        figure(figRaster);
        suptitle(celltypePool{i_celltype});
        saveas(figRaster,['E:\2P\summary\meanActivity_by_ROIpreference',filesep,celltypePool{i_celltype},info_str,'-selectvitiy-',epoch4preferPool2show{1},'2',epoch4preferPool2show{2},'-masklick',masklick,'-',activity_type,'_raster2showConsistency.pdf'],'pdf');
        saveas(figRaster,['E:\2P\summary\meanActivity_by_ROIpreference',filesep,celltypePool{i_celltype},info_str,'-selectvitiy-',epoch4preferPool2show{1},'2',epoch4preferPool2show{2},'-masklick',masklick,'-',activity_type,'_raster2showConsistency.png'],'png');
        figure(figMean);
        suptitle(celltypePool{i_celltype});
        saveas(figMean,['E:\2P\summary\meanActivity_by_ROIpreference',filesep,celltypePool{i_celltype},info_str,'-selectvitiy-',epoch4preferPool2show{1},'2',epoch4preferPool2show{2},'-masklick',masklick,'-',activity_type,'_mean2showConsistency.pdf'],'pdf');
        saveas(figMean,['E:\2P\summary\meanActivity_by_ROIpreference',filesep,celltypePool{i_celltype},info_str,'-selectvitiy-',epoch4preferPool2show{1},'2',epoch4preferPool2show{2},'-masklick',masklick,'-',activity_type,'_mean2showConsistency.png'],'png');
    end
end
%loop again to see the individual cases from different sessions
%}
%saving figures from cases into a PPT using exportToPPT
pptFileName='E:\2P\summary\meanActivity_by_ROIpreference\cases.pptx';
isOpen  = exportToPPTX();
if ~isempty(isOpen),
    % If PowerPoint already started, then close first and then open a new one
    exportToPPTX('close');
end
slidesize=[16 9];
if exist(pptFileName,'file')
    exportToPPTX('open',pptFileName);
else
    exportToPPTX('new','Dimensions',slidesize, ...
        'Title','Comparing z-scored activities for different cell type neurons', ...
        'Author','PYX', ...
        'Comments','This file has been automatically generated by exportToPPTX');
end
behEventAlignPool={'delay onset','delay onset','delay onset','go cue','first lick'};
epoch4preferPool={'delay','mid_delay','late_delay','response','lick'};%{'sound','delay','mid_delay','late_delay','response','lick'}
masklickPool={'yes','yes','yes','no','no'};
for i_celltype=1:length(celltypePool)
    ind_session2=logical(ind_session.*logical(strcmp(T.cell_type,celltypePool{i_celltype})+strcmp(T.cell_type,[celltypePool{i_celltype},'-flpo'])));
    Tchoose=T(ind_session2,:);
    % file_path='H:\2P\pyx349_20210418\im_data_reg\result_save';
    % session='pyx349_20210418';
    ind_ROI=[];
    frT=0;
    for i=1:size(Tchoose,1)
        file_path=Tchoose.file_path{i};
        session=[Tchoose.session{i},'_',Tchoose.field{i}];
        savepath=[file_path,filesep,session];
        for i_activity_type=1:length(activity_typePool)
            activity_type=activity_typePool{i_activity_type};
            for i_align=1:length(behEventAlignPool)
                behEventAlign=behEventAlignPool{i_align};
                masklick=masklickPool{i_align};
                epoch4prefer=epoch4preferPool{i_align};
                objsession=Session2P(session,file_path,trial2include,trial2exclude);

                trialTypeStr='cor and err';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
                AUCtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
                AUCCorrectedMethod='SensoryChoiceOrthogonalSubtraction';%'None';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
                [TAUC, Tmean,objsession] = objsession.mGetEpochAUC(Tchoose.cell_type,Tchoose.ROI_type{i},trialTypeStr,AUCtype,AUCCorrectedMethod);
                [meanActivityByTrialType,cellActivityByTrialType,Tout]=objsession.mMeanPSTHbyROI(activity_type,[],behEventAlign,masklick,i_selectivity);
                TAUC.date=cellfun(@(x) datestr(datetime(x,'InputFormat','yyyy/MM/dd'),'yyyymmdd'), TAUC.date,'UniformOutput',false);
                TActivityByTrialType=join(Tout,TAUC);
                 %plot cases for single session
                %reorganized the activities based on AUC
                if strcmp(behEventAlign,'stim onset')
                    frameNumTime=[0.5,2];%from 5s before align point to 5s after align point
                    %epoch4prefer='sound'; %time for comparing ipsi vs. contra
                elseif strcmp(behEventAlign,'delay onset')
                    frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
                    %epoch4prefer='delay';%{'mid_delay','late_delay'}
                elseif strcmp(behEventAlign,'go cue')
                    frameNumTime=[1,2];%from 5s before align point to 5s after align point
                    %epoch4prefer='response';
                elseif strcmp(behEventAlign,'first lick')
                    frameNumTime=[1,2];%from 5s before align point to 5s after align point
                    %epoch4prefer='lick';
                end
                meanActivityReorganized = fReorganizeByAUC(meanActivityByTrialType,TActivityByTrialType,epoch4prefer);
                %plot
                [figRaster_case(i_align),figMean_case(i_align)]=fPlotRasterMean(meanActivityReorganized,meanActivityByTrialType,behEventAlign,frameNumTime,objsession.metadata.frT);
            end
            %save single case with different aligment in a PPT
            slideId = exportToPPTX('addslide');
            n_fig=length(behEventAlignPool);
            temp1=get(figRaster_case(1),'Position');
            temp2=get(figMean_case(1),'Position');
            widthRatio=temp1(3)/(temp1(3)+temp2(3));
            WR=[widthRatio,1-widthRatio];
            RH=[temp1(4)/temp1(3),temp2(4)/temp2(3)];
            for i_align=1:length(behEventAlignPool)
                exportToPPTX('addpicture',figRaster_case(i_align),'Position',[1+(slidesize(1)-1)*(i_align-1)/n_fig 1 slidesize(1)/n_fig/1.2*WR(1) slidesize(1)/n_fig/1.2*WR(1)*RH(1)]);
                exportToPPTX('addpicture',figMean_case(i_align),'Position',[1+(slidesize(1)-1)*(i_align-1)/n_fig+slidesize(1)/n_fig/1.2*WR(1) 1 slidesize(1)/n_fig/1.2*WR(2) slidesize(1)/n_fig/1.2*WR(2)*RH(2)]);
                exportToPPTX('addtext',session,'Position',[1 0 3 0.5],'Vert','top');
                exportToPPTX('addtext',celltypePool{i_celltype},'Position',[slidesize(1)-1 0 3 0.5],'Vert','top');
            end
            close(figRaster_case(:));
            close(figMean_case(:));
        end
    end
end
% Save and close (in one command)
newFile = exportToPPTX('saveandclose',pptFileName);
fprintf('New file has been saved: <a href="matlab:open(''%s'')">%s</a>\n',newFile,newFile);
%}

%assist function
function [A,figMean]=fPlotRasterMean(meanActivityReorganized,meanActivityByTrialType_combine,behEventAlign,frameNumTime,frT)
%% extract data and organized in a new way
%get preference of each ROI
frameNum=floor(frameNumTime*1000/frT);
%plot color plot,4 column(correct/error/miss/violation), 7 row(for each
%stimlus)+ mean trace
A=figure();
set(gcf, 'position', [50 50 150*size(meanActivityReorganized,2) 150*size(meanActivityReorganized,1)]);%3 rows, ipsi prefering, contra perfering, ns.; 2 columns, perfer, non-prefer
meanActivityRange=cell2mat(meanActivityByTrialType_combine);
meanActivityRange=reshape(meanActivityRange,[],1);
ColLimit = prctile(meanActivityRange,98);
titlestr={'Ipsi side','Contra side'};
ylabelstr={'ipsi preferring','contra preferring','n.s.'};

titlestr=strcat(titlestr);
if size(meanActivityReorganized,2)==6
    color_mean_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};
elseif size(meanActivityReorganized,2)==2
    color_mean_trace={[0 0 1],[1 0 0]};
end
%for plotting setting
generalMargin=[0.1,0.1];
fixMargin=0.05;
lenRatio=cellfun(@(x) size(x,1),meanActivityReorganized(:,1));
[prctPos,prctPos2] = fGetPrctPos(flip(lenRatio),fixMargin, generalMargin);
figMean=figure();
set(gcf, 'position', [50+150*size(meanActivityReorganized,2) 50 150 150*size(meanActivityReorganized,1)]);%3 rows, ipsi prefering, contra perfering, ns.

for nCelltype=1:size(meanActivityReorganized,1) %for each cell type, ie. ipsi prefering/contra prefering
    %sort peak of prefered side
    if nCelltype==2
        tempneuralActivity=meanActivityReorganized{nCelltype,2};
    else
        tempneuralActivity=meanActivityReorganized{nCelltype,1};
    end
    peaktemp=max(tempneuralActivity,[],2);
    peakind=zeros(size(tempneuralActivity,1),1);
    for i=1:size(tempneuralActivity,1)
        peakind(i)=min(find(tempneuralActivity(i,:)==peaktemp(i)));
    end
    [B,I]=sort(peakind);
    
    for nTrialtype=1:size(meanActivityReorganized,2) %for each stimulus/choice, ie. ipsi/contra
        neuralActivity=meanActivityReorganized{nCelltype,nTrialtype};
        figure(A);
        %ax=subplot(size(meanActivityReorganized,1),size(meanActivityReorganized,2),nTrialtype+size(meanActivityReorganized,2)*(nCelltype-1));%grid position
        [nrow,ncol]=size(meanActivityReorganized);
        pos_left=0.1+0.8/ncol*(nTrialtype-1)+0.8/ncol*0.15;
        pos_bottom=prctPos(nrow-nCelltype+1);
        pos_width=0.8/ncol*0.75;
        pos_height=prctPos2(nrow-nCelltype+1)-prctPos(nrow-nCelltype+1);
        pos_current=[pos_left,pos_bottom,pos_width,pos_height];
        subplot('Position',pos_current);
        hold on;
        imagesc(neuralActivity(I,:));
        plot([round(frameNum(1)),round(frameNum(1))],[0,size(neuralActivity,1)],'w--');%aligned event
        x_lim=[0,size(neuralActivity,2)];%get(gca,'Xlim');
        
        set(gca,'clim',[0 ColLimit]);
        set(gca,'ytick',size(neuralActivity,1),'yticklabel',size(neuralActivity,1),'ydir','normal');
        if nCelltype==size(meanActivityReorganized,1)
            set(gca,'xtick',[round(1000/frT*(frameNumTime(1)-floor(frameNumTime(1))))+1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
            %                 xlabel(['time(s) from ',behEventAlign],'FontName','Arial','FontSize',14);
        else
            set(gca,'xtick',[]);
        end
        if nCelltype==1 && nTrialtype==1
            ColBar = colorbar;
            set(ColBar,'position',[0.9 0.08 0.02 0.1445]);
        end
        if size(neuralActivity,2)>0 && ~isempty(neuralActivity)
            ylim([0.5 size(neuralActivity,1)+0.5]);
        end
        if nCelltype == 1
            title(titlestr{nTrialtype});
        end
        if nTrialtype==1
            ylabel(ylabelstr{nCelltype});
        end
        xlim(x_lim);
        if nCelltype==size(meanActivityReorganized,1) && nTrialtype==1
            xlabel(['Time (s) from ',behEventAlign]);
        end
        %plot mean trace
        figure(figMean);
        subplot(size(meanActivityReorganized,1),1,nCelltype);
        %[neuralActivityMean, neuralActivityCI]=fMean_CI(neuralActivity,0.05);
        ts=double((-frameNum(1):1:frameNum(2))*frT/1000);
        curve_meanTrace(nTrialtype)=fPlotMean_SE(ts,neuralActivity,color_mean_trace{nTrialtype});
        hold on;
        ylabel('\it\DeltaF/F');
        y_lim=get(gca,'Ylim');
        plot([0,0],y_lim,'k-');
        if strcmp(behEventAlign,'stim onset')
            plot([0.5,0.5],y_lim,'k-');
        elseif strcmp(behEventAlign,'delay onset')
            plot([-0.5,-0.5],y_lim,'k-');
        end
        if nCelltype==size(meanActivityReorganized,1) && nTrialtype==1
            xlabel(['Time (s) from ',behEventAlign]);
        end
        
    end
    if nCelltype == 1
        hl=legend(curve_meanTrace(:),{'ipsi','contra'});
        set(hl,'Box','off');
    end
end
end

function [figRaster,figMean]=fRasterMeanMultiAlign(meanActivityReorganized,behEventAlignPool,frT,activity_type)
%plot raster plot
figRaster=figure();
set(gcf, 'position', [50 50 150*size(meanActivityReorganized,3) 150*size(meanActivityReorganized,1)]);%3 rows, ipsi prefering, contra perfering, ns.; 2 columns, perfer, non-prefer
figMean=figure();
set(gcf, 'position', [50+150*size(meanActivityReorganized,3) 50 150*size(meanActivityReorganized,3)*2 150*size(meanActivityReorganized,1)]);%3 rows, ipsi prefering, contra perfering, ns.

tempcell=reshape(meanActivityReorganized,[],2);
meanActivityRange=cell2mat(tempcell);
meanActivityRange=reshape(meanActivityRange,[],1);
ColLimit = prctile(meanActivityRange,98);
titlestr={'Ipsi side','Contra side','Ipsi side','Contra side'};
%ylabelstr={'ipsi preferring','contra preferring','n.s.'};
ylabelstr={'ipsi','contra','n.s.'};
switch activity_type
    case 'dff'
        ylabel_str= '\it\DeltaF/F';
    case 'spkr'
        ylabel_str= 'spikes/s';
    case 'zscored_dff'
        ylabel_str= 'z-scored \it\DeltaF/F';
    case 'zscored_spkr'
        ylabel_str= 'z-scored spikes/s';
end
titlestr=strcat(titlestr);
color_trialtype={[0,0,1],[1,0,0],[0.5,0.5,0.5]};%ipsi, contra, n.s.
color_mean_trace={[0 0 1],[1 0 0]};

%for plotting setting
generalMargin=[0.1,0.1];
fixMargin=0.05;
lenRatio=cellfun(@(x) size(x,1),meanActivityReorganized(:,:,1));
lenRatio=sum(lenRatio,2);
[prctPos,prctPos2] = fGetPrctPos(flip(lenRatio),fixMargin, generalMargin);

for nCelltype=1:size(meanActivityReorganized,1) %for each cell type, ie. ipsi prefering/contra prefering
    ratio2ndEpoch=cellfun(@(x) size(x,1), meanActivityReorganized(nCelltype,:,1), 'UniformOutput', true);
    for nTrialtype=1:size(meanActivityReorganized,3) %for each stimulus/choice, ie. ipsi/contra
        neuralActivity_cell=reshape(meanActivityReorganized(nCelltype,:,nTrialtype),[],1);
        neuralActivity=cell2mat(neuralActivity_cell);
        figure(figRaster);
        %ax=subplot(size(meanActivityReorganized,1),size(meanActivityReorganized,2),nTrialtype+size(meanActivityReorganized,2)*(nCelltype-1));%grid position
        [nrow,~,ncol]=size(meanActivityReorganized);
        pos_left=0.1+0.8/ncol*(nTrialtype-1)+0.8/ncol*0.15;
        pos_bottom=prctPos(nrow-nCelltype+1);
        pos_width=0.8/ncol*0.75;
        pos_height=prctPos2(nrow-nCelltype+1)-prctPos(nrow-nCelltype+1);
        pos_current=[pos_left,pos_bottom,pos_width,pos_height];
        subplot('Position',pos_current);
        hold on;
        imagesc(neuralActivity);
        iAlign=ceil(nTrialtype/2);
        behEventAlign=behEventAlignPool{iAlign};
        if strcmp(behEventAlign,'stim onset')
            frameNumTime=[0.5,2];%from 5s before align point to 5s after align point
            time4prefer=[-0.8,1.5]; %time for comparing ipsi vs. contra
        elseif strcmp(behEventAlign,'delay onset')
            frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
            time4prefer=[-0.3,1];
        else
            frameNumTime=[1,2];%from 5s before align point to 5s after align point
            time4prefer=[0,1];
        end
        frameNum=floor(frameNumTime*1000/frT);
        plot([round(frameNum(1)),round(frameNum(1))],[0,size(neuralActivity,1)],'w--');%aligned event
        x_lim=[0,size(neuralActivity,2)];%get(gca,'Xlim');
        
        set(gca,'clim',[0 ColLimit]);
        set(gca,'ytick',size(neuralActivity,1),'yticklabel',size(neuralActivity,1),'ydir','normal');
        
        if nCelltype==size(meanActivityReorganized,1)
            set(gca,'xtick',[round(1000/frT*(frameNumTime(1)-floor(frameNumTime(1))))+1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
            %                 xlabel(['time(s) from ',behEventAlign],'FontName','Arial','FontSize',14);
        else
            set(gca,'xtick',[]);
        end
        if nCelltype==1 && nTrialtype==1
            ColBar = colorbar;
            set(ColBar,'position',[0.9 0.08 0.02 0.1445]);
        end
        if size(neuralActivity,2)>0 && ~isempty(neuralActivity)
            ylim([0.5 size(neuralActivity,1)+0.5]);
        end
        if nCelltype == 1
            title(titlestr{nTrialtype});
        end
        if nTrialtype==1
            ylabel(ylabelstr{nCelltype});
        end
        xlim(x_lim);
        if nCelltype==size(meanActivityReorganized,1) && mod(nTrialtype,2)==1
            xlabel(['Time (s) from ',behEventAlign]);
        end
        set(gca,'FontSize',10);
    end
    %plot bar to indicates second epoch preference
    posBarLeft=pos_left+pos_width;
    posBarWidth=0.8/ncol*0.05;
    posBarBottom=pos_bottom;
    posBarHeight=pos_height;
    pos_bar=[posBarLeft,posBarBottom,posBarWidth,posBarHeight];
    subplot('Position',pos_bar);
    hold on;
    ratio2ndEpoch=[0,ratio2ndEpoch];
    for i=1:size(meanActivityReorganized,2)
        patch([0,1,1,0],[sum(ratio2ndEpoch(1:i)),sum(ratio2ndEpoch(1:i)),sum(ratio2ndEpoch(1:i+1)),sum(ratio2ndEpoch(1:i+1))],color_trialtype{i});
    end
    set(gca,'Ylim',[0,sum(ratio2ndEpoch)],'ytick',[],'yticklabel',[]);
    %plot mean trace
    figure(figMean);
    n_col=size(meanActivityReorganized,2)*size(meanActivityReorganized,3)/2;
    n_row=size(meanActivityReorganized,1); 
    for i_col=1:n_col/2
        iAlign=ceil(i_col/3);
        behEventAlign=behEventAlignPool{iAlign};
        subplot(n_row,n_col,i_col+n_col*(nCelltype-1));
        ts=double((-frameNum(1):1:frameNum(2))*frT/1000);
        curve_meanTrace(1)=fPlotMean_SE(ts,meanActivityReorganized{nCelltype,i_col,1},color_mean_trace{1});%ipsi,1st align
        hold on;
        curve_meanTrace(2)=fPlotMean_SE(ts,meanActivityReorganized{nCelltype,i_col,2},color_mean_trace{2});%contra, 1st align
        if i_col==1
            ylabel({ylabelstr{nCelltype},ylabel_str});
        else
            ylabel(ylabel_str);
        end
        y_lim=get(gca,'Ylim');
        plot([0,0],y_lim,'k-');
        if strcmp(behEventAlign,'stim onset')
            plot([0.5,0.5],y_lim,'k-');
        elseif strcmp(behEventAlign,'delay onset')
            plot([-0.5,-0.5],y_lim,'k-');
        end
        if nCelltype==size(meanActivityReorganized,1) && i_col==2
            xlabel(['Time (s) from ',behEventAlign]);
        end
        if nCelltype ==1
            title([ylabelstr{i_col},'\_after\_',behEventAlignPool{2}]);
        end
    end
    for i_col=n_col/2+1:n_col
        iAlign=ceil(i_col/3);
        behEventAlign=behEventAlignPool{iAlign};
        subplot(n_row,n_col,i_col+n_col*(nCelltype-1));
        ts=double((-frameNum(1):1:frameNum(2))*frT/1000);
        curve_meanTrace(3)=fPlotMean_SE(ts,meanActivityReorganized{nCelltype,i_col-3,3},color_mean_trace{1});%ipsi,1st align
        hold on;
        curve_meanTrace(4)=fPlotMean_SE(ts,meanActivityReorganized{nCelltype,i_col-3,4},color_mean_trace{2});%contra, 1st align
        
        ylabel(ylabel_str);
        y_lim=get(gca,'Ylim');
        plot([0,0],y_lim,'k-');
        if strcmp(behEventAlign,'stim onset')
            plot([0.5,0.5],y_lim,'k-');
        elseif strcmp(behEventAlign,'delay onset')
            plot([-0.5,-0.5],y_lim,'k-');
        end
        if nCelltype==size(meanActivityReorganized,1) && i_col==5
            xlabel(['Time (s) from ',behEventAlign]);
        end
        if nCelltype ==1
            title([ylabelstr{i_col-3},'\_after\_',behEventAlignPool{2}]);
        end
    end
    if i_col == 1 && nCelltype==size(meanActivityReorganized,1)
        hl=legend(curve_meanTrace(1:2),{'ipsi','contra'});
        set(hl,'Box','off');
    end

end
end

function [meanActivityReorganized,indSigIpsi,indSigContra,indNS] = fReorganizeByAUC(meanActivityByTrialType_combine,TActivityByTrialType,epoch4prefer)
frT_v=cellfun(@(x) size(x{1,1}{1,1},2), TActivityByTrialType.cellActivityByTrialType,'UniformOutput',true);
frN=max(frT_v);
meanActivityByTrialType_combine= fInter(meanActivityByTrialType_combine,frN);
switch epoch4prefer
    case 'sound'
        epochAUC=TActivityByTrialType.sound;
        epochP=TActivityByTrialType.psound;
    case 'delay'
        epochAUC=TActivityByTrialType.delay;
        epochP=TActivityByTrialType.pdelay;        
    case 'mid_delay'
        epochAUC=TActivityByTrialType.mid_delay;
        epochP=TActivityByTrialType.pmid_delay;  
    case 'late_delay'
        epochAUC=TActivityByTrialType.late_delay;
        epochP=TActivityByTrialType.plate_delay;  
    case 'response'
        epochAUC=TActivityByTrialType.response;
        epochP=TActivityByTrialType.presponse;  
    case 'lick'
        epochAUC=TActivityByTrialType.lick;
        epochP=TActivityByTrialType.plick;  
end
epochP=1-abs(epochP-0.5)*2;
indSig=epochP<0.05;
indIpsi=epochAUC<0.5;
indContra=epochAUC>0.5;
indSigIpsi=logical(indSig.*indIpsi);
indSigContra=logical(indSig.*indContra);
indNS=logical((~indSigContra).*(~indSigIpsi));
%reorganized the data
meanActivityReorganized=cell(3,2);
meanActivityReorganized{1,1}=meanActivityByTrialType_combine{1}(indSigIpsi,:);
meanActivityReorganized{1,2}=meanActivityByTrialType_combine{2}(indSigIpsi,:);
meanActivityReorganized{2,1}=meanActivityByTrialType_combine{1}(indSigContra,:);
meanActivityReorganized{2,2}=meanActivityByTrialType_combine{2}(indSigContra,:);
meanActivityReorganized{3,1}=meanActivityByTrialType_combine{1}(indNS,:);
meanActivityReorganized{3,2}=meanActivityByTrialType_combine{2}(indNS,:);

end

function [meanActivityReorganized,indIpsi,indContra,indNS] = fReorganizeByMultiAUC(meanActivityByTrialType_combine,TActivityByTrialType,behEventAlignPool,epoch4preferPool)
%Input-
%   cellActivityByTrialType_combine- nTrialtyp-by-nAlignType cell array,
%   meanActivityByTrialType_combine- nTrialtyp-by-nAlignType cell array,
%   behEventAlignPool- 1-by-nAlignType cell array
nAlignType=length(behEventAlignPool);%only used for check, currently only 2 alignments
if (length(epoch4preferPool) == nAlignType) && (size(meanActivityByTrialType_combine,2) == nAlignType)
    for iAlign=1:nAlignType
        epoch4prefer=epoch4preferPool{iAlign};
        [~,indIpsi{iAlign},indContra{iAlign},indNS{iAlign}] = fReorganizeByAUC(meanActivityByTrialType_combine(:,iAlign),TActivityByTrialType,epoch4prefer);
     end
    %each pannel include one trial type using first alignment, and labeled
    %using second alignment
    %reorganized the data
    meanActivityReorganized=cell(3,3,4);%3d, assume ipsi/contra 2 trial types, and delay/lick 2 epochs
    indPrefer1=[indIpsi{1},indContra{1},indNS{1}];
    indPrefer2=[indIpsi{2},indContra{2},indNS{2}];
    for i=1:size(meanActivityReorganized,1)%preference during 1st epoch
        for j=1:size(meanActivityReorganized,2)%preference during 2nd epoch
            indPrefer=logical(indPrefer1(:,i).*indPrefer2(:,j));
            meanActivityReorganized{i,j,1}=meanActivityByTrialType_combine{1,1}(indPrefer,:);
            meanActivityReorganized{i,j,2}=meanActivityByTrialType_combine{2,1}(indPrefer,:);
            meanActivityReorganized{i,j,3}=meanActivityByTrialType_combine{1,2}(indPrefer,:);
            meanActivityReorganized{i,j,4}=meanActivityByTrialType_combine{2,2}(indPrefer,:);
            %sort based on the peak of the prefered activities during delay
            if i==2
                tempneuralActivity=meanActivityReorganized{i,j,2};
            else
                tempneuralActivity=meanActivityReorganized{i,j,1};
            end
            peaktemp=max(tempneuralActivity,[],2);
            peakind=zeros(size(tempneuralActivity,1),1);
            for k=1:size(tempneuralActivity,1)
                peakind(k)=find(tempneuralActivity(k,:)==peaktemp(k), 1 );
            end
            [B,I]=sort(peakind);
            meanActivityReorganized(i,j,:)=cellfun(@(x) x(I,:), meanActivityReorganized(i,j,:),'UniformOutput',false);
        end
    end
    
else
    error('Input size not consistent');
end
    
end

function [meanActivityReorganized,indIpsi,indContra,indNS] = fReorganizeByPreference(cellActivityByTrialType_combine,meanActivityByTrialType_combine,ind4prefer,varargin)
% extract data and organized in a new way
parameter4preference='max';
if ~isempty(varargin)
    if mod(length(varargin),2)~=0
        warning('Input Name-value pair for fReorganizeByPreference function not match');
    else
        for i=1:2:length(varargin)
            switch varargin{i}
                case 'parameter4preference'
                    parameter4preference=varargin{i+1};
                otherwise
                    warning(['Invalid name-value pair named ',varargin{i}]);
            end
        end
    end
end

%get preference of each ROI
cellActivityByTrialType_cor=cellActivityByTrialType_combine(:,1);%only plot correct trials

if length(ind4prefer)>=size(meanActivityByTrialType_combine{1,1},2)
    ind4prefer=ind4prefer(1:size(meanActivityByTrialType_combine{1,1},2));
else
    ind4prefer(length(ind4prefer)+1:size(meanActivityByTrialType_combine{1,1},2))=false;
end
[meanEpochActivity,maxEpochActivity,meanTrialEpochActivity,meanMaxTrialEpochActivity]=deal(cell(size(cellActivityByTrialType_cor)));
for i_cell=1:length(cellActivityByTrialType_cor)
    meanEpochActivity{i_cell}=cellfun(@(x) nanmean(x(:,ind4prefer),2), cellActivityByTrialType_cor{i_cell},'UniformOutput',false);
    maxEpochActivity{i_cell}=cellfun(@(x) max(x(:,ind4prefer),[],2), cellActivityByTrialType_cor{i_cell},'UniformOutput',false);
    meanTrialEpochActivity{i_cell}=cellfun(@nanmean, meanEpochActivity{i_cell},'UniformOutput',true);
    meanMaxTrialEpochActivity{i_cell}=cellfun(@nanmean, maxEpochActivity{i_cell},'UniformOutput',true);
end
indIpsi_byMean=meanTrialEpochActivity{1}>meanTrialEpochActivity{2};
indIpsi_byMax=meanMaxTrialEpochActivity{1}>meanMaxTrialEpochActivity{2};
indSig_byTtest_mean=cellfun(@(x,y) ttest2(x,y),meanEpochActivity{1},meanEpochActivity{2});%not perform multiple correction here
indSig_byTtest_max=cellfun(@(x,y) ttest2(x,y),maxEpochActivity{1},maxEpochActivity{2});
temp=indIpsi_byMean+indIpsi_byMax;%showing here, many neurons have inconsistent pereference based on mean acitivities or max activiities
figHist=figure();
histogram(temp);
indIpsi_byMean=logical(indIpsi_byMean.*indSig_byTtest_mean);
indContra_byMean=logical((~indIpsi_byMean).*indSig_byTtest_mean);
indNS_byMean=logical((~indIpsi_byMean).*(~indContra_byMean));
indIpsi_byMax=logical(indIpsi_byMax.*indSig_byTtest_max);
indContra_byMax=logical((~indIpsi_byMax).*indSig_byTtest_max);
indNS_byMax=logical((~indIpsi_byMax).*(~indContra_byMax));
if strcmp(parameter4preference,'mean')
    [indIpsi,indContra,indNS]=deal(indIpsi_byMean,indContra_byMean,indNS_byMean);
elseif strcmp(parameter4preference,'max')
    [indIpsi,indContra,indNS]=deal(indIpsi_byMax,indContra_byMax,indNS_byMax);
end

%reorganized the data
meanActivityReorganized=cell(3,2);
meanActivityReorganized{1,1}=meanActivityByTrialType_combine{1}(indIpsi,:);
meanActivityReorganized{1,2}=meanActivityByTrialType_combine{2}(indIpsi,:);
meanActivityReorganized{2,1}=meanActivityByTrialType_combine{1}(indContra,:);
meanActivityReorganized{2,2}=meanActivityByTrialType_combine{2}(indContra,:);
meanActivityReorganized{3,1}=meanActivityByTrialType_combine{1}(indNS,:);
meanActivityReorganized{3,2}=meanActivityByTrialType_combine{2}(indNS,:);
end

function [meanActivityReorganized,indIpsi,indContra,indNS] = fReorganizeByMultiPreference(cellActivityByTrialType_combine,meanActivityByTrialType_combine,behEventAlignPool,frT)
%Input-
%   cellActivityByTrialType_combine- nTrialtyp-by-nAlignType cell array,
%   meanActivityByTrialType_combine- nTrialtyp-by-nAlignType cell array,
%   behEventAlignPool- 1-by-nAlignType cell array
nAlignType=length(behEventAlignPool);
if (size(cellActivityByTrialType_combine,2) == nAlignType) && (size(meanActivityByTrialType_combine,2) == nAlignType)
    for iAlign=1:nAlignType
        behEventAlign=behEventAlignPool{iAlign};
        if strcmp(behEventAlign,'stim onset')
            frameNumTime=[0.5,2];%from 5s before align point to 5s after align point
            time4prefer=[-0.8,1.5]; %time for comparing ipsi vs. contra
        elseif strcmp(behEventAlign,'delay onset')
            frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
            time4prefer=[-0.3,1];
        else
            frameNumTime=[1,2];%from 5s before align point to 5s after align point
            time4prefer=[0,1];
        end
        frameNum=floor(frameNumTime*1000/frT);
        frame4prefer=floor(time4prefer*1000/frT);
        ind4prefer=false(sum(frameNum)+1,1);
        ind4prefer(frameNum(1)-frame4prefer(1):frameNum(1)+frame4prefer(2)+1)=true;
        [~,indIpsi{iAlign},indContra{iAlign},indNS{iAlign}] = fReorganizeByPreference(cellActivityByTrialType_combine(:,iAlign),meanActivityByTrialType_combine(:,iAlign),ind4prefer,'parameter4preference','mean');
    end
    %each pannel include one trial type using first alignment, and labeled
    %using second alignment
    %reorganized the data
    meanActivityReorganized=cell(3,3,4);%3d, assume ipsi/contra 2 trial types, and delay/lick 2 epochs
    indPrefer1=[indIpsi{1},indContra{1},indNS{1}];
    indPrefer2=[indIpsi{2},indContra{2},indNS{2}];
    for i=1:size(meanActivityReorganized,1)%preference during 1st epoch
        for j=1:size(meanActivityReorganized,2)%preference during 2nd epoch
            indPrefer=logical(indPrefer1(:,i).*indPrefer2(:,j));
            meanActivityReorganized{i,j,1}=meanActivityByTrialType_combine{1,1}(indPrefer,:);
            meanActivityReorganized{i,j,2}=meanActivityByTrialType_combine{2,1}(indPrefer,:);
            meanActivityReorganized{i,j,3}=meanActivityByTrialType_combine{1,2}(indPrefer,:);
            meanActivityReorganized{i,j,4}=meanActivityByTrialType_combine{2,2}(indPrefer,:);
            %sort based on the peak of the prefered activities during delay
            if i==2
                tempneuralActivity=meanActivityReorganized{i,j,2};
            else
                tempneuralActivity=meanActivityReorganized{i,j,1};
            end
            peaktemp=max(tempneuralActivity,[],2);
            peakind=zeros(size(tempneuralActivity,1),1);
            for k=1:size(tempneuralActivity,1)
                peakind(k)=find(tempneuralActivity(k,:)==peaktemp(k), 1 );
            end
            [B,I]=sort(peakind);
            meanActivityReorganized(i,j,:)=cellfun(@(x) x(I,:), meanActivityReorganized(i,j,:),'UniformOutput',false);
        end
    end
    
else
    error('Input size not consistent');
end
    
end

function [prctPos,prctPos2] = fGetPrctPos(lenRatio,fixMargin, generalMargin)
%from bottom to top, the height is
%generalMargin(1),lenRatio(1),fixMargin(1),lenRatio(2),fixMargin(2),...generalMargin(2)
prctPos=generalMargin(1);
unit=(1-sum(generalMargin)-length(lenRatio)*fixMargin)/sum(lenRatio);
prctPos2=prctPos+unit*lenRatio(1);
for i=2:length(lenRatio)
    prctPos=[prctPos, prctPos2(end)+fixMargin];
    prctPos2=[prctPos2, prctPos(end)+unit*lenRatio(i)];
end

end