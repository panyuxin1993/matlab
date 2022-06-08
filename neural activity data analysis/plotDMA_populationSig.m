%plot dynamic memory assembly, how neurons show significant activities
%across trials (show profile of neurons and of trials)
close all;
% dbstop if error;

behEventAlign='delay onset';
behEventSort='go cue';
prct_delay2calculate_pResp=1;

celltype_str={'M2','vglut2','vgat'};%{'M2','syn','vglut2','vgat'};
celltype_abbr_str={'M2','SC-E','SC-I'};
% activity_form_pool={'dff','spkr'};%%%%%%%%%%%%%%%%%%%%%%'dff' 2 std, spkr may >3std
% activity_smooth_binsize_pool={200,400};%ms, 200for dff, 400 for spkr
% threshold_dur_pool=[500,200];
activity_form_pool={'dff'};%%%%%%%%%%%%%%%%%%%%%%'dff' 2 std, spkr may >3std
activity_smooth_binsize_pool={200};%ms, 200for dff, 400 for spkr
threshold_dur_pool=[500];

trialTypeStr='cor';

saving_path='E:\2P\summary\population_significance';
[num,txt,raw] =xlsread('C:\Users\PYX\Documents\DataSummary\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
T=cell2table(raw(2:end,1:15));
T.Properties.VariableNames=strrep(raw(1,1:15),' ','_');%table variable name can't have ' ',so replace them

%-----hyper parameters----
selectivitySigEvent=zeros(length(activity_form_pool),4,length(celltype_str));%1dff/spkr,2d-none, ipsi,balence,contra,3d-cell type; value represent proportion of neurons
selectivitySigEventCounts=zeros(length(activity_form_pool),2,length(celltype_str));%2d ipsi, contra
labelstr_selectivity={'none','ipsi','balence','contra'};
labelstr_selectivity_sig_event={'ipsi','contra'};
pResponse_all=cell(length(celltype_str),length(activity_form_pool));
%-----loop each conditions----
for i_celltype=1:length(celltype_str)
    for i_activity_form=1:length(activity_form_pool)
        clear significantTrialNum4cell_all tbSigTrial4Cell;
        threshold_dur=threshold_dur_pool(i_activity_form);%continously 500ms significant activities to be judged as a significant event
        activity_form=activity_form_pool{i_activity_form};
        activity_smooth_binsize=activity_smooth_binsize_pool{i_activity_form};
        ind_session=logical(strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control')...
            .*logical(strcmp(T.cell_type,celltype_str{i_celltype})+strcmp(T.cell_type,[celltype_str{i_celltype},'-flpo']))...
            .*contains(T.ROI_type,'soma').*(~contains(T.field,'soma')));%SC ROIs, not probe session
        Tchoose=T(ind_session,:);
        save_name_mat=[saving_path,filesep,'TsigTrialNum4cell_all_',celltype_str{i_celltype},'_',activity_form,'_binSize-',activity_smooth_binsize,'_prctDelay2CalPResp-',num2str(prct_delay2calculate_pResp),'.mat'];

        %---------get table for each session---------
        %
        fig_n_event=figure;
        set(fig_n_event,'Position',[200,200,500,200]);
        fig_event_time=figure;
        set(fig_event_time,'Position',[100,0,1000,200*size(Tchoose,1)]);        
        for i_session=1:size(Tchoose,1)
            indrow=i_session;
            figure(fig_event_time);
            ax_temp1=subplot(size(Tchoose,1),4,i_session*4-3);
            ax_temp2=subplot(size(Tchoose,1),4,i_session*4-2);
            ax_temp3=subplot(size(Tchoose,1),4,i_session*4-1);
            ax_temp4=subplot(size(Tchoose,1),4,i_session*4);
            ax_input=[ax_temp1,ax_temp2,ax_temp3,ax_temp4];
            [tbSigTrial4Cell,figDMAsig,figSigCellNum,savingNameStr] = fGetDMAsigTableAsession(Tchoose.file_path{indrow},...
                Tchoose.session{indrow},Tchoose.animal{indrow},Tchoose.date{indrow},...
                Tchoose.field{indrow},celltype_str{i_celltype},...
                trialTypeStr,activity_form,activity_smooth_binsize,prct_delay2calculate_pResp,threshold_dur,...
                fig_n_event,ax_input);
            if ~isempty(figDMAsig)%if no figure, this is nan
                saveas(figDMAsig,[saving_path,filesep,Tchoose.session{indrow},'threshold_dur-',num2str(threshold_dur),'.pdf'],'pdf');
            end
            if exist('significantTrialNum4cell_all','var')
                significantTrialNum4cell_all=vertcat(significantTrialNum4cell_all,tbSigTrial4Cell);
            else
                significantTrialNum4cell_all=tbSigTrial4Cell;
            end
        end
        save(save_name_mat,'significantTrialNum4cell_all');
        figure(fig_n_event);
        subplot(1,2,1);
        title('All significant events');
        xlabel('Continues significant response times (ms)');
        ylabel('Trial counts per cell');
        subplot(1,2,2);
        title('Significant events during stim. and delay');
        xlabel('Continues significant response times (ms)');
        saveas(fig_n_event,[saving_path,filesep,'EventNum_',celltype_str{i_celltype},'_',activity_form_pool{i_activity_form},savingNameStr,'.png'],'png');
        saveas(fig_event_time,[saving_path,filesep,'EventTime_',celltype_str{i_celltype},'_',activity_form_pool{i_activity_form},savingNameStr,'.fig'],'fig');
        saveas(fig_event_time,[saving_path,filesep,'EventTime_',celltype_str{i_celltype},'_',activity_form_pool{i_activity_form},savingNameStr,'.png'],'png');
        set(fig_event_time,'PaperPosition',[0,0,2.5,size(Tchoose,1)]);      
        saveas(fig_event_time,[saving_path,filesep,'EventTime_',celltype_str{i_celltype},'_',activity_form_pool{i_activity_form},savingNameStr,'.pdf'],'pdf');
        %}
        
        %-------------further analyses-----------------
        load(save_name_mat);
        %--------- group neurons together---------------
        %
        figGroup=figure;
        set(gcf,'Position',[100,200,1000,200]);
        significantTrialNum4cell_all=sortrows(significantTrialNum4cell_all,{'contraMipsi','contra_sum'},{'descend','descend'});
        ind_row_show=logical((significantTrialNum4cell_all.ipsi_sum>0)+(significantTrialNum4cell_all.contra_sum>0));
        selectivitySigEvent(i_activity_form,1,i_celltype)=sum(~ind_row_show)/length(ind_row_show);
        selectivitySigEvent(i_activity_form,2,i_celltype)=sum(significantTrialNum4cell_all.ipsi_sum>significantTrialNum4cell_all.contra_sum)/length(ind_row_show);
        selectivitySigEvent(i_activity_form,4,i_celltype)=sum(significantTrialNum4cell_all.ipsi_sum<significantTrialNum4cell_all.contra_sum)/length(ind_row_show);
        selectivitySigEvent(i_activity_form,3,i_celltype)=1-sum(selectivitySigEvent(i_activity_form,[1,2,4],i_celltype));
        selectivitySigEventCounts(i_activity_form,1,i_celltype)=sum(significantTrialNum4cell_all.ipsi_sum);
        selectivitySigEventCounts(i_activity_form,2,i_celltype)=sum(significantTrialNum4cell_all.contra_sum);
        x=1:size(significantTrialNum4cell_all(ind_row_show,:),1);
        h=stem(x,[significantTrialNum4cell_all.ipsi_sum(ind_row_show),-significantTrialNum4cell_all.contra_sum(ind_row_show)],'filled','LineWidth',1);
        set(h(1),'Color','blue');
        set(h(2),'Color','red');
        xlabel('# cells');
        ylabel('# trials with significant responses');
        saveas(figGroup,[saving_path,filesep,'grouped_neuron_',celltype_str{i_celltype},'_',activity_form,'_threshold_dur-',num2str(threshold_dur),'.pdf'],'pdf');
        saveas(figGroup,[saving_path,filesep,'grouped_neuron_',celltype_str{i_celltype},'_',activity_form,'_threshold_dur-',num2str(threshold_dur),'.png'],'png');
        %}
        pResponse_all{i_celltype,i_activity_form}=significantTrialNum4cell_all.pSigResponse;

    end
end
%-----plot selectivity of significant event by cell----
%
fig_selectivitySigEvent=figure;
set(fig_selectivitySigEvent,'Position',[200,200,500,200]);
title_str_bar={'\it\DeltaF/F','Spiking rate'};

for i_subplot=1:2
    subplot(1,2,i_subplot);%dff/spkr data
    temp_selectivity=squeeze(selectivitySigEvent(i_subplot,:,:))';
    bar_color_map=[1,1,1;0,0,1;0.5,0.5,0.5;1,0,0];
    colormap(bar_color_map);%previous than R2017b
    set(groot,'defaultAxesColorOrder',bar_color_map);%after R2017b,and remember to remove the setting
    hb=bar(temp_selectivity,0.6,'stacked');
    set(gca,'Xlim',[0.5,size(temp_selectivity,1)+0.5],'XTick',1:size(temp_selectivity,1),'XTickLabel',celltype_abbr_str);
    if i_subplot==1
        ylabel('Proportion of cells');
    end
    if i_subplot==2
        hl=legend(labelstr_selectivity);
        set(hl,'Box','off');
    end
    title(title_str_bar{i_subplot});
    set(gca,'FontSize',12);
    box off;

end
set(fig_selectivitySigEvent,'PaperPosition',[1,1,4,1.5]);
saveas(fig_selectivitySigEvent,[saving_path,filesep,'selectivity_sig_event_threshold_dur-',num2str(threshold_dur),'.pdf'],'pdf');
%}
%-----plot selectivity of significant event counts----
%
fig_selectivitySigEventCounts=figure;
set(fig_selectivitySigEventCounts,'Position',[200,200,500,200]);
title_str_bar={'\it\DeltaF/F','Spiking rate'};

for i_subplot=1:2
    subplot(1,2,i_subplot);%dff/spkr data
    temp_selectivity=squeeze(selectivitySigEventCounts(i_subplot,:,:))';
    temp_selectivity=temp_selectivity./repmat(sum(temp_selectivity,2),1,2);
    bar_color_map=[0,0,1;1,0,0];
    colormap(bar_color_map);%previous than R2017b
    set(groot,'defaultAxesColorOrder',bar_color_map);%after R2017b,and remember to remove the setting
    hb=bar(temp_selectivity,0.6,'stacked');
    set(gca,'Xlim',[0.5,size(temp_selectivity,1)+0.5],'XTick',1:size(temp_selectivity,1),'XTickLabel',celltype_abbr_str);
    if i_subplot==1
        ylabel('Proportion of trials');
    end
    if i_subplot==2
        hl=legend(labelstr_selectivity_sig_event);
        set(hl,'Box','off');
    end
    title(title_str_bar{i_subplot});
    set(gca,'FontSize',12);
    box off;

end
set(fig_selectivitySigEventCounts,'PaperPosition',[1,1,4,1.5]);
saveas(fig_selectivitySigEventCounts,[saving_path,filesep,'selectivity_sig_event_threshold_dur-',num2str(threshold_dur),'.pdf'],'pdf');
%}
%%
%-----plot significant reponse probability of neurons------
fig_pRes=figure;
set(fig_pRes,'Position',[200,200,500,200]);
title_str_hist={'\it\DeltaF/F','Spiking rate'};
for i_subplot=1:2
    
    subplot(1,2,i_subplot);%dff/spkr data
    hold on;
    bar_color_map=[0,0,0;0,0,1;1,0,0];
    colormap(bar_color_map);%previous than R2017b
    set(groot,'defaultAxesColorOrder',bar_color_map);%after R2017b,and remember to remove the setting

    xlabel('P(Significant response)');
    for i_celltype=1:length(celltype_str)
        histogram(pResponse_all{i_celltype,i_subplot},'BinWidth',0.01,'DisplayStyle','stairs','Normalization','probability');
    end

    if i_subplot==1
%         ylabel('Proportion of trials');
        ylabel('Proportion of cells');
    end                               
%     if i_subplot==2
        hl=legend(celltype_abbr_str);
        set(hl,'Box','off');
%     end
    title(title_str_hist{i_subplot});
    set(gca,'FontSize',12);
    box off;

end
set(fig_pRes,'PaperPosition',[1,1,4,1.5]);
saveas(fig_pRes,[saving_path,filesep,'probability_sig_event_threshold_dur-',num2str(threshold_dur_pool(1)),'.pdf'],'pdf');