clear meanActivityByTrialType_combine;
close all;
[num,txt,raw] =xlsread('C:\Users\PYX\Documents\DataSummary\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
T=cell2table(raw(2:end,1:15));
T.Properties.VariableNames=strrep(raw(1,1:15),' ','_');%table variable name can't have ' ',so replace them
ind_session=logical(strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*(~contains(T.behavior_performance,'probe')).*contains(T.ROI_type,'soma'));%not probe session

celltypePool={'M2','syn','vglut2','vgat'};%{'syn','vglut2','vgat'};
for i_celltype=1:length(celltypePool)
    ind_session2=logical(ind_session.*strcmp(T.cell_type,celltypePool{i_celltype}));
    Tchoose=T(ind_session2,:);
    % file_path='H:\2P\pyx349_20210418\im_data_reg\result_save';
    % session='pyx349_20210418';
    ind_ROI=[];
    behEventAlignPool={'delay onset','go cue','first lick'};
    masklickPool={'yes','no','no'};
    i_selectivity=3;%3-choice
    behEventSortPool={'go cue','first lick','go cue'};
    trial2include='all';
    trial2exclude=nan;
    activity_typePool={'spkr'};%{'dff','spkr'};
    frT=0;
    for i_align=1:length(behEventAlignPool)
        behEventAlign=behEventAlignPool{i_align};
        masklick=masklickPool{i_align};
        for i_activity_type=1:length(activity_typePool)
            activity_type=activity_typePool{i_activity_type};
            clear meanActivityByTrialType_combine;
            savenamestr=['E:\2P\summary\meanActivity_by_ROI',filesep,celltypePool{i_celltype},'-meanActivityByROI-align',behEventAlign,'-masklick',masklick,'-',activity_type,'.mat'];
            if exist(savenamestr,'file')
                load(savenamestr);
            else
                for i=1:size(Tchoose,1)

                    file_path=Tchoose.file_path{i};
                    session=[Tchoose.session{i},'_',Tchoose.field{i}];
                    savepath=[file_path,filesep,session];
                    objsession=Session2P(session,file_path,trial2include,trial2exclude);
                    meanActivityByTrialType=objsession.mMeanPSTHbyROI(activity_type,[],behEventAlign,masklick,i_selectivity);
                    if exist('meanActivityByTrialType_combine','var')
                        frT=max(objsession.metadata.frT,frT);
                        frN=max(size(meanActivityByTrialType_combine{1,1},2),size(meanActivityByTrialType{1,1},2));
                        meanActivityByTrialType_combine= fInter(meanActivityByTrialType_combine,frN);
                        meanActivityByTrialType=fInter(meanActivityByTrialType,frN);
                        disp(['frN:',num2str(frN),',meanActivityByTrialType_combine:',num2str(size(meanActivityByTrialType_combine{1,1},2)),',meanActivityByTrialType:',num2str(size(meanActivityByTrialType{1,1},2))]);
                        meanActivityByTrialType_combine=cellfun(@(x,y) vertcat(x,y),meanActivityByTrialType_combine,meanActivityByTrialType,'UniformOutput',false);
                    else
                        meanActivityByTrialType_combine=meanActivityByTrialType;
                    end
                end
            end
            save(savenamestr,'meanActivityByTrialType_combine','frT');
            %%
            %plot color plot,4 column(correct/error/miss/violation), 7 row(for each
            %stimlus)+ mean trace
            %plotting settings
            if strcmp(behEventAlign,'stim onset')
                frameNumTime=[0.5,2];%from 5s before align point to 5s after align point
            elseif strcmp(behEventAlign,'delay onset')
                frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
            else
                frameNumTime=[2,3.5];%from 5s before align point to 5s after align point
                meanActivityByTrialType_combine=meanActivityByTrialType_combine(:,1:2);%only plot correct and error trials
            end
            frameNum=frameNumTime*1000/frT;
            A=figure;
            set(gcf, 'position', [50 50 100*size(meanActivityByTrialType_combine,2) 100*(size(meanActivityByTrialType_combine,1)+1)]);%last row, mean trace
            meanActivityRange=cell2mat(meanActivityByTrialType_combine);
            meanActivityRange=reshape(meanActivityRange,[],1);
            ColLimit = prctile(meanActivityRange,98);
            titlestr={'Correct','Error','Miss','Violation'};
            suptitle([celltypePool{i_celltype},'-align',behEventAlign,'-masklick',masklick,'-',activity_type]);
            titlestr=strcat(titlestr);
            if size(meanActivityByTrialType_combine,1)==6
                color_mean_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};
            elseif size(meanActivityByTrialType_combine,1)==2
                color_mean_trace={[0 0 1],[1 0 0]};
            end
            %     if strcmp(rule,'low click rate-right')
            %         color_mean_trace=fliplr(color_mean_trace);
            %     end
            
            %sort peak
            tempneuralActivity=meanActivityByTrialType_combine{2,1};
            peaktemp=max(tempneuralActivity,[],2);
            peakind=zeros(size(tempneuralActivity,1),1);
            for i=1:size(tempneuralActivity,1)
                peakind(i)=min(find(tempneuralActivity(i,:)==peaktemp(i)));
            end
            [B,I]=sort(peakind);
            
            %[B,I]=sort(meanActivityByTrialType_combine{1,1}(:,round(frameNum(1))));
            for nStim=1:size(meanActivityByTrialType_combine,1) %for each stimulus
                
                for nResult=1:size(meanActivityByTrialType_combine,2) %4 column(correct/error/miss/violation),companied with 4 lick raster
                    neuralActivity=meanActivityByTrialType_combine{nStim,nResult};
                    figure(A);
                    subplot(size(meanActivityByTrialType_combine,1)+1,size(meanActivityByTrialType_combine,2),nResult+size(meanActivityByTrialType_combine,2)*(nStim-1));
                    hold on;
                    imagesc(neuralActivity(I,:));
                    plot([round(frameNum(1)),round(frameNum(1))],[0,size(neuralActivity,1)],'w--');
                    x_lim=[0,size(neuralActivity,2)];%get(gca,'Xlim');
                    
                    set(gca,'clim',[0 ColLimit]);
                    set(gca,'ytick',size(neuralActivity,1),'yticklabel',size(neuralActivity,1),'ydir','normal');
                    if nStim==size(meanActivityByTrialType_combine,1)
                        set(gca,'xtick',[round(1000/frT*(frameNumTime(1)-floor(frameNumTime(1))))+1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
                        %                 xlabel(['time(s) from ',behEventAlign],'FontName','Arial','FontSize',14);
                    else
                        set(gca,'xtick',[]);
                    end
                    if nStim==1 & nResult==1
                        ColBar = colorbar;
                        set(ColBar,'position',[0.93 0.1100 0.02 0.1445]);
                    end
                    if size(neuralActivity,2)>0
                        ylim([0.5 size(neuralActivity,1)+0.5]);
                    end
                    if nStim == 1
                        title(titlestr{nResult});
                    end
                    if nResult==1
                        if nStim==1
                            ylabel('ipsi');
                        else
                            ylabel('contra'); %' clicks/s'
                        end
                    end
                    xlim(x_lim);
                    %plot mean trace
                    subplot(size(meanActivityByTrialType_combine,1)+1,size(meanActivityByTrialType_combine,2),nResult+size(meanActivityByTrialType_combine,2)*size(meanActivityByTrialType_combine,1));
                    %[neuralActivityMean, neuralActivityCI]=fMean_CI(neuralActivity,0.05);
                    ts=double((-frameNum(1):1:frameNum(2))*frT/1000);
                    curve_meanTrace(nStim)=fPlotMean_SE(ts,neuralActivity,color_mean_trace{nStim});
                    y_lim=get(gca,'Ylim');
                    plot([0,0],y_lim,'k-');
                    if strcmp(behEventAlign,'stim onset')
                        plot([0.5,0.5],y_lim,'k-');
                    elseif strcmp(behEventAlign,'delay onset')
                        plot([-0.5,-0.5],y_lim,'k-');
                    end
                    if nStim==size(meanActivityByTrialType_combine,1) && nResult==1
                        xlabel(['Time (s) from ',behEventAlign]);
                    end
                    
                end
            end
            suptitle([celltypePool{i_celltype},'-align',behEventAlign,'-masklick',masklick,'-',activity_type]);
            saveas(A,['E:\2P\summary\meanActivity_by_ROI',filesep,celltypePool{i_celltype},'-align',behEventAlign,'-masklick',masklick,'-',activity_type,'.pdf'],'pdf');
        end
    end
end