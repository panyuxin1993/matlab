close all;
clear;
[num,txt,raw] =xlsread('C:\Users\PYX\Documents\DataSummary\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
savepath='E:\2P\summary';
clear TAUC_combine Tmean_combine;
% trialTypeStr='cor and err';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
% AUCtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
% AUCCorrectedMethod='SensoryChoiceOrthogonalSubtraction';%'None';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
trialTypeStr='cor';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
AUCtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
AUCCorrectedMethod='None';%'None';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
manipulation='control';
%% calculate AUC
%{
T=cell2table(raw(2:end,1:15));
T.Properties.VariableNames=strrep(raw(1,1:15),' ','_');%table variable name can't have ' ',so replace them
ind_session1=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*(~contains(T.behavior_performance,'probe'));%not probe session
% ind_session2=logical(contains(T.cell_type,'syn').*contains(T.ROI_type,'soma'));
ind_session2=logical(strcmp(T.cell_type,'vglut2')+strcmp(T.cell_type,'vgat')+strcmp(T.cell_type,'M2')+strcmp(T.cell_type,'vglut2-flpo')+strcmp(T.cell_type,'vgat-flpo'));
% ind_session2=logical(contains(T.field,'soma')+contains(T.field,'dendrite'));
% ind_session2=logical(strcmp(T.cell_type,'M2'));
% ind_session2=logical(~contains(T.field,'soma').*contains(T.ROI_type,'soma'));%include all soma data, but exclude zoomin fields
ind_session= ind_session1 & ind_session2;
ind_trial_chosen4ScatterHist=logical(contains(T.ROI_type,'soma').*(~contains(T.field,'soma')));
ind_session= ind_session & ind_trial_chosen4ScatterHist;
ind_session=find(ind_session);
n_session=length(ind_session);
animal_unique=unique(T.animal(ind_session));
n_animal=length(animal_unique);

for i_session=1:n_session
    indrow=ind_session(i_session);
    [TAUC_currentSession, Tmean_currentSession] = fGetEpochAUCtableASession(T.file_path{indrow},T.animal{indrow},T.date{indrow},T.field{indrow},T.cell_type{indrow},T.ROI_type{indrow},trialTypeStr,AUCtype,AUCCorrectedMethod);
%     [TAUC_currentSession, Tmean_currentSession] = fGetEpochAUCtableASession_NPseg(T.file_path{indrow},T.animal{indrow},T.date{indrow},T.field{indrow},T.cell_type{indrow},trialTypeStr,AUCtype,AUCCorrectedMethod);    
    if strcmp(TAUC_currentSession.animal{1},'pyx298')
        a=1;
    end
    if isempty(TAUC_currentSession)
        continue;
    end
    if exist('TAUC_combine','var') && strcmp(AUCtype,'stimuli') %calculate choice probability
        if size(TAUC_currentSession.ITI,2)==2
            TAUC_combine=vertcat(TAUC_combine,TAUC_currentSession);
        elseif size(TAUC_currentSession.ITI,2)>2 %for sessions with probe sound, choose only the end sound
            TAUC_currentSession_formatted=TAUC_currentSession;
            for ivar=1:size(TAUC_currentSession_formatted,2)
                temp=table2array(TAUC_currentSession_formatted(:,ivar));
                if size(temp,2)>2
                    tempvarname=TAUC_currentSession_formatted.Properties.VariableNames{ivar};
                    TAUC_currentSession_formatted=removevars(TAUC_currentSession_formatted,TAUC_currentSession_formatted.Properties.VariableNames{ivar});
                    TAUC_currentSession_formatted = addvars(TAUC_currentSession_formatted,temp(:,[1,end]),'NewVariableNames',tempvarname,'Before',TAUC_currentSession_formatted.Properties.VariableNames{ivar});
                end
            end
            TAUC_combine=vertcat(TAUC_combine,TAUC_currentSession_formatted);
        end
    elseif exist('TAUC_combine','var') && ~strcmp(AUCtype,'stimuli')
        TAUC_combine=vertcat(TAUC_combine,TAUC_currentSession);
    else
        TAUC_combine=TAUC_currentSession;
    end
    if exist('Tmean_combine','var')
        Tmean_combine=vertcat(Tmean_combine,Tmean_currentSession);
    else
        Tmean_combine=Tmean_currentSession;
    end
end
save([savepath,filesep,'trialType',trialTypeStr,'-',AUCtype,'TepochAUC.mat'],'TAUC_combine');
save([savepath,filesep,'trialType',trialTypeStr,'-TepochMeanActivity.mat'],'Tmean_combine');
% save([savepath,filesep,'trialType',trialTypeStr,'-',AUCtype,'TepochAUC_NPseg.mat'],'TAUC_combine');
% save([savepath,filesep,'trialType',trialTypeStr,'-TepochMeanActivity_NPseg.mat'],'Tmean_combine');
%}
%% load data for further anaylyses
%
trialTypeStr='cor';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
AUCtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
load([savepath,filesep,'trialType',trialTypeStr,'-',AUCtype,'TepochAUC.mat']);
% load([savepath,filesep,'trialType',trialTypeStr,'-',AUCtype,'TepochAUC_NPseg.mat']);
% load('F:\2P\summary\trialTypecorTepochAUC.mat');
celltype='vglut2';
manipulation='control';
%T_AUC=TAUC_combine(strcmp(TAUC_combine.celltype,'syn')|contains(TAUC_combine.celltype,'+other'),:);%syn, ..+other
% T_AUC=TAUC_combine(strcmp(TAUC_combine.celltype,celltype),:);
% ind_trial_chosen=strcmp(TAUC_combine.animal,'pyx358') & strcmp(TAUC_combine.date,'2021/6/10');
ind_trial_chosen=logical(contains(TAUC_combine.ROItype,'soma').*contains(TAUC_combine.celltype,celltype));
% ind_trial_chosen=logical(contains(TAUC_combine.ROItype,'dendrite').*contains(TAUC_combine.celltype,celltype));
T_AUC=TAUC_combine(ind_trial_chosen,:);


%% plot histogram for choice AUC 
%{
% psig=[0.025,0.975];
psig=[0,0.025;0.975,1];
preference_threshold=0.5;
%plot histogram from ITI to lick of AUC
figAUC=fHistEpochAUC(T_AUC,psig,preference_threshold,'AUC');
a=1;

%indcate example as vertical line
%no such neuron
% Tsorted=sortrows(T_AUC,{'delay','lick'},{'descend','descend'});%vlgut2
% indCaseRow=5;%shift from ipsi to contra
% for i=1:5
%     subplot(1,5,i);
%     plot( table2array(Tsorted(indCaseRow,5+i))*ones(2,1),ylim,'k-');hold on;
%     box off;
% end
% text(0,1,[Tsorted.animal{indCaseRow},Tsorted.date{indCaseRow},'-',num2str(Tsorted.nROI(indCaseRow))],'Unit','Normalized');
% 
% Tsorted=sortrows(T_AUC,{'delay','lick'},{'ascend','descend'});%vglut2
% indCaseRow=14;%shift from ipsi to contra
% for i=1:5
%     subplot(1,5,i);
%     plot( table2array(Tsorted(indCaseRow,5+i))*ones(2,1),ylim,'b-');hold on;
%     box off;
% end
% text(0,0.9,[Tsorted.animal{indCaseRow},Tsorted.date{indCaseRow},'-',num2str(Tsorted.nROI(indCaseRow))],'Unit','Normalized','color','b');

% Tsorted=sortrows(T_AUC,{'delay','lick'},{'ascend','ascend'});%vgat
% indCaseRow=10;%shift from ipsi to contra
% for i=1:5
%     subplot(1,5,i);
%     plot( table2array(Tsorted(indCaseRow,5+i))*ones(2,1),ylim,'b-');hold on;
%     box off;
% end
% text(0,0.9,[Tsorted.animal{indCaseRow},Tsorted.date{indCaseRow},'-',num2str(Tsorted.nROI(indCaseRow))],'Unit','Normalized','color','b');
% 
% Tsorted=sortrows(T_AUC,{'lick','delay'},{'descend','descend'});%vgat
% indCaseRow=10;%shift from ipsi to contra
% for i=1:5
%     subplot(1,5,i);
%     plot( table2array(Tsorted(indCaseRow,5+i))*ones(2,1),ylim,'b-');hold on;
%     box off;
% end
% text(0,0.9,[Tsorted.animal{indCaseRow},Tsorted.date{indCaseRow},'-',num2str(Tsorted.nROI(indCaseRow))],'Unit','Normalized','color','b');

% Tsorted=sortrows(T_AUC,{'delay','response'},{'descend','descend'});
% indCaseRow=1;%shift from contra to contra
% for i=1:5
%     subplot(1,5,i);
%     plot( table2array(Tsorted(indCaseRow,5+i))*ones(2,1),ylim,'r-');
%     box off;
% end
% text(0,0.8,[Tsorted.animal{indCaseRow},Tsorted.date{indCaseRow},'-',num2str(Tsorted.nROI(indCaseRow))],'Unit','Normalized','color','r');
%save figure
saveas(figAUC,[savepath,filesep,celltype,'-',trialTypeStr,'-',AUCtype,'AUC at significance level of ',num2str(2*psig(1)),'.png'],'png');
set(figAUC,'PaperPosition',[0,0,8,2]);
saveas(figAUC,[savepath,filesep,celltype,'-',trialTypeStr,'-',AUCtype,'AUC at significance level of ',num2str(2*psig(1)),'.pdf'],'pdf');
%}

%% compare SC AUC and M2 AUC
%{
figCmpAUC=figure;
set(gcf,'Position',[100,100,400,200]);
load('E:\2P\summary\summary_mat\M2_delay_table.mat');
subplot(1,2,1);
if strcmp(AUCtype,'choice')
    histogram(M2_delay_table.Choice_auc,'FaceColor','w','BinWidth',0.05);hold on;
    histogram(M2_delay_table.Choice_auc(logical((M2_delay_table.Choice_p<0.05).*(M2_delay_table.Choice_auc<0.5))),'FaceColor','b','BinWidth',0.05);
    histogram(M2_delay_table.Choice_auc(logical((M2_delay_table.Choice_p<0.05).*(M2_delay_table.Choice_auc>0.5))),'FaceColor','r','BinWidth',0.05);
    xlabel('AUC');
    set(gca,'Xlim',[0,1]);
    set(gca,'FontSize',14);
    box off;
    subplot(1,2,2);
    histogram(M2_delay_table.Choice_auc,'EdgeColor','k','BinWidth',0.05,'Normalization','probability','DisplayStyle','stairs');
    hold on;
    histogram(T_AUC.delay,'EdgeColor',[1,0.5,0],'BinWidth',0.05,'Normalization','probability','DisplayStyle','stairs');
elseif strcmp(AUCtype,'sensory')
    histogram(M2_delay_table.Sensory_auc,'FaceColor','w','BinWidth',0.05);hold on;
    histogram(M2_delay_table.Sensory_auc(logical((M2_delay_table.Sensory_p<0.05).*(M2_delay_table.Sensory_auc<0.5))),'FaceColor','b','BinWidth',0.05);
    histogram(M2_delay_table.Sensory_auc(logical((M2_delay_table.Sensory_p<0.05).*(M2_delay_table.Sensory_auc>0.5))),'FaceColor','r','BinWidth',0.05);
    xlabel('AUC');
    set(gca,'Xlim',[0,1]);
    set(gca,'FontSize',14);
    box off;
    subplot(1,2,2);
    histogram(M2_delay_table.Sensory_auc,'EdgeColor','k','BinWidth',0.05,'Normalization','probability','DisplayStyle','stairs');
    hold on;
    histogram(T_AUC.delay,'EdgeColor',[1,0.5,0],'BinWidth',0.05,'Normalization','probability','DisplayStyle','stairs');
end
hl=legend('M2','SC');
set(hl,'box','off');
xlabel('AUC');
set(gca,'Xlim',[0,1]);
set(gca,'FontSize',14);
box off;
%}
%% compare M2 AUC and SC syn/vglut2/vgat AUC
%{
ind_trial_chosen_M2=strcmp(TAUC_combine.celltype,'M2');
T_AUC_M2=TAUC_combine(ind_trial_chosen_M2,:);
ind_trial_chosen_syn=strcmp(TAUC_combine.celltype,'syn');
T_AUC_syn=TAUC_combine(ind_trial_chosen_syn,:);
ind_trial_chosen_vglut2=strcmp(TAUC_combine.celltype,'vglut2');
T_AUC_vglut2=TAUC_combine(ind_trial_chosen_vglut2,:);
ind_trial_chosen_vgat=strcmp(TAUC_combine.celltype,'vgat');
T_AUC_vgat=TAUC_combine(ind_trial_chosen_vgat,:);
T_AUC_EI=TAUC_combine(logical(ind_trial_chosen_vgat+ind_trial_chosen_vglut2),:);
psig=[0,0.025;0.975,1];
preference_threshold=0.5;
bin_width=0.05;
figCmpAUC=figure;
set(figCmpAUC,'Position',[100,100,600,400]);
subplot(2,3,1);
histogram(T_AUC_M2.delay,'FaceColor','w','BinWidth',bin_width);
hold on;
histogram(T_AUC_M2.delay(findSig(T_AUC_M2.pdelay,psig) & T_AUC_M2.delay>preference_threshold),'FaceColor','r','BinWidth',bin_width);
histogram(T_AUC_M2.delay(findSig(T_AUC_M2.pdelay,psig) & T_AUC_M2.delay<preference_threshold),'FaceColor','b','BinWidth',bin_width);
ylabel('cell counts');
title('M2');
text(0.1,1,['n=',num2str(size(T_AUC_M2,1))],'Unit','Normalized');
set(gca,'Xlim',[0,1],'FontSize',14);
box off;

% subplot(2,3,2);
% histogram(T_AUC_syn.delay,'FaceColor','w','BinWidth',bin_width);
% hold on;
% histogram(T_AUC_syn.delay(findSig(T_AUC_syn.pdelay,psig) & T_AUC_syn.delay>preference_threshold),'FaceColor','r','BinWidth',bin_width);
% histogram(T_AUC_syn.delay(findSig(T_AUC_syn.pdelay,psig) & T_AUC_syn.delay<preference_threshold),'FaceColor','b','BinWidth',bin_width);
% title('SC syn');
% text(0.1,1,['n=',num2str(size(T_AUC_syn,1))],'Unit','Normalized');
% set(gca,'Xlim',[0,1]);
% box off;
% 
% subplot(2,3,3);
% histogram(T_AUC_M2.delay,'EdgeColor','k','BinWidth',0.05,'Normalization','probability','DisplayStyle','stairs','LineWidth',1);
% hold on;
% histogram(T_AUC_syn.delay,'EdgeColor',[1,0.5,0],'BinWidth',0.05,'Normalization','probability','DisplayStyle','stairs','LineWidth',1);
% hl=legend('M2','SC');
% set(hl,'box','off');
% xlabel('AUC');
% ylabel('Proportion');
% set(gca,'Xlim',[0,1]);
% box off;

subplot(2,3,2);
histogram(T_AUC_EI.delay,'FaceColor','w','BinWidth',bin_width);
hold on;
histogram(T_AUC_EI.delay(findSig(T_AUC_EI.pdelay,psig) & T_AUC_EI.delay>preference_threshold),'FaceColor','r','BinWidth',bin_width);
histogram(T_AUC_EI.delay(findSig(T_AUC_EI.pdelay,psig) & T_AUC_EI.delay<preference_threshold),'FaceColor','b','BinWidth',bin_width);
title('SC');
text(0.1,1,['n=',num2str(size(T_AUC_EI,1))],'Unit','Normalized');
set(gca,'Xlim',[0,1],'FontSize',14);
box off;
subplot(2,3,3);
histogram(T_AUC_M2.delay,'EdgeColor','k','BinWidth',0.05,'Normalization','probability','DisplayStyle','stairs','LineWidth',1);
hold on;
histogram(T_AUC_EI.delay,'EdgeColor',[1,0.5,0],'BinWidth',0.05,'Normalization','probability','DisplayStyle','stairs','LineWidth',1);
hl=legend('M2','SC');
set(hl,'box','off');
xlabel('AUC');
ylabel('Proportion');
set(gca,'Xlim',[0,1],'FontSize',14);
box off;

subplot(2,3,4);
histogram(T_AUC_vglut2.delay,'FaceColor','w','BinWidth',bin_width);
hold on;
histogram(T_AUC_vglut2.delay(findSig(T_AUC_vglut2.pdelay,psig) & T_AUC_vglut2.delay>preference_threshold),'FaceColor','r','BinWidth',bin_width);
histogram(T_AUC_vglut2.delay(findSig(T_AUC_vglut2.pdelay,psig) & T_AUC_vglut2.delay<preference_threshold),'FaceColor','b','BinWidth',bin_width);
title('SC vglut2');
text(0.1,1,['n=',num2str(size(T_AUC_vglut2,1))],'Unit','Normalized');
set(gca,'Xlim',[0,1],'FontSize',14);
box off;

subplot(2,3,5);
histogram(T_AUC_vgat.delay,'FaceColor','w','BinWidth',bin_width);
hold on;
histogram(T_AUC_vgat.delay(findSig(T_AUC_vgat.pdelay,psig) & T_AUC_vgat.delay>preference_threshold),'FaceColor','r','BinWidth',bin_width);
histogram(T_AUC_vgat.delay(findSig(T_AUC_vgat.pdelay,psig) & T_AUC_vgat.delay<preference_threshold),'FaceColor','b','BinWidth',bin_width);
title('SC vgat');
text(0.1,1,['n=',num2str(size(T_AUC_vgat,1))],'Unit','Normalized');
set(gca,'Xlim',[0,1],'FontSize',14);
box off;

subplot(2,3,6);
histogram(T_AUC_vglut2.delay,'EdgeColor',[0,0.5,1],'BinWidth',0.05,'Normalization','probability','DisplayStyle','stairs','LineWidth',1);
hold on;
histogram(T_AUC_vgat.delay,'EdgeColor',[1,0,0],'BinWidth',0.05,'Normalization','probability','DisplayStyle','stairs','LineWidth',1);
hl=legend('vglut2','vgat');
set(hl,'box','off');
xlabel('AUC');
ylabel('Proportion');
set(gca,'Xlim',[0,1],'FontSize',14);
box off;
%}
%% plot mean activities scatter with rotated histogram
%{
load([savepath,filesep,'trialType',trialTypeStr,'-TepochMeanActivity.mat']);
% load([savepath,filesep,'trialType',trialTypeStr,'-TepochMeanActivity_NPseg.mat']);
Tmean=Tmean_combine(ind_trial_chosen,:);
psig=[0,0.01];
preference_threshold=0;
Tipsi = fExtractTableField(Tmean,'ipsi');
Tcontra = fExtractTableField(Tmean,'contra');
epochstr={'ITI','sound','delay','response','lick'};
n_epoch=length(epochstr);
figScatterRotatedHist=figure;
set(gcf,'Position',[100,100,1000,180]);
figTemp=figure;
xy_lim=[0,0];
xy_lim_setted={[-1,1],[-0.5,0.5],[-0.5,0.5],[-1,2],[-2,4]};
data_diff=cell(1,n_epoch);
for i=1:n_epoch
    figure(figScatterRotatedHist);
    subplot(1,n_epoch,i);
    ind_tablecol=find(cellfun(@(x) strcmp(x,epochstr{i}), Tipsi.Properties.VariableNames));
    scatter(Tipsi{:,ind_tablecol},Tcontra{:,ind_tablecol},10,'k','filled');
    data_diff{i}=Tipsi{:,ind_tablecol}-Tcontra{:,ind_tablecol};
    xlabel('Mean activities ipsi choice');
    ylabel('Mean activities contra choice');
    x_lim_temp=get(gca,'Xlim');
    y_lim_temp=get(gca,'Ylim');
    if ~exist('xy_lim_setted','var')
        xy_lim(1)=min(min(x_lim_temp),min(y_lim_temp));
        xy_lim(2)=max(max(x_lim_temp),max(y_lim_temp));
    else
        xy_lim=xy_lim_setted{i};
    end
%     xy_lim(1)=min(min(min(xy_lim),min(x_lim_temp)),min(y_lim_temp));
%     xy_lim(2)=max(max(max(xy_lim),max(x_lim_temp)),max(y_lim_temp));
    title(epochstr{i});
    hold on;    
% end
% for i=1:n_epoch
    subplot(1,n_epoch,i);
    plot(xy_lim,[0,0],'k--');
    plot([0,0],xy_lim,'k--');
    plot(xy_lim,xy_lim,'k--');
    figure(figTemp);
    h=histogram(data_diff{i},'Normalization','probability','NumBins',20);
    [patch_array_x,patch_array_y]=hist2patch(h.BinEdges,h.Values*(xy_lim(2)-xy_lim(1))*1);
    [rotate_x,rotate_y] = cellfun(@(x,y) fRotateCurve(x,y,-45,0,0),patch_array_x,patch_array_y,'UniformOutput',false);
    final_x=cellfun(@(x) x+(xy_lim(2)-xy_lim(1))*0.8+xy_lim(1),rotate_x,'UniformOutput',false);
    final_y=cellfun(@(x) x+(xy_lim(2)-xy_lim(1))*0.8+xy_lim(1),rotate_y,'UniformOutput',false);
    line_x=mean(data_diff{i})*ones(1,2);
    line_y=[0,max(h.Values)*(xy_lim(2)-xy_lim(1))];
    [line_x_rotate,line_y_rotate]=fRotateCurve(line_x,line_y,-45,0,0);
    line_x_final=line_x_rotate+(xy_lim(2)-xy_lim(1))*0.8+xy_lim(1);
    line_y_final=line_y_rotate+(xy_lim(2)-xy_lim(1))*0.8+xy_lim(1);
    figure(figScatterRotatedHist);
    subplot(1,n_epoch,i);
    for i_patch=1:length(final_x)
        patch(final_x{i_patch},final_y{i_patch},[0.5,0.5,0.5],'EdgeColor','none');
    end
    plot(line_x_final,line_y_final,'k-');
    set(gca,'Xlim',xy_lim,'Ylim',xy_lim);
    axis equal;%so the x and y have same unit length
end
%}

% %plot histogram from ITI to lick of AUC
% figAUC=fHistEpochAUC(Tmean,psig,preference_threshold,'mean activities (\it\DeltaF/F\rm)',[-2,4],0.3);

%% plot bar
%{
n_cell= size(T_AUC,1);
psig=[0.025,0.975];
pAUCSigSC=zeros(3,4);%each row- contra, ipsi, n.s., each column- sound, delay, response, lick
pAUCSigSC(1,1)=sum(T_AUC.psound<psig(1))/n_cell;
pAUCSigSC(2,1)=sum(T_AUC.psound>psig(2))/n_cell;
pAUCSigSC(1,2)=sum(T_AUC.pdelay<psig(1))/n_cell;
pAUCSigSC(2,2)=sum(T_AUC.pdelay>psig(2))/n_cell;
pAUCSigSC(1,3)=sum(T_AUC.presponse<psig(1))/n_cell;
pAUCSigSC(2,3)=sum(T_AUC.presponse>psig(2))/n_cell;
pAUCSigSC(1,4)=sum(T_AUC.plick<psig(1))/n_cell;
pAUCSigSC(2,4)=sum(T_AUC.plick>psig(2))/n_cell;
pAUCSigSC(3,:)=1-sum(pAUCSigSC);
figBar=figure;
set(gcf,'Position',[100,100,400,160]);
subplot(1,2,2);
map=[1,0,0;0,0,1;0.5,0.5,0.5];%red-contra, ipsi-blue, n.s.-grey
colormap(map);%previous than R2017b
set(groot,'defaultAxesColorOrder',map);%after R2017b,and remember to remove the setting
bar(pAUCSigSC','stacked');

box off;
%set(gca,'XTick',1:4,'XTickLabel',{'sound','delay','response','lick'});
set(gca,'XTick',1:4,'XTickLabel',{'S','D','R','L'});
title(['SC ',celltype,' neurons']);
xlim([0,5]);
set(gca,'FontSize',12);
SC_choice_AUC_005=table(n_cell*pAUCSigSC(1,:)',n_cell*pAUCSigSC(2,:)',n_cell*ones(4,1),'VariableNames',{'contra_sig','ipsi_sig','total_n'},'RowNames',{'sound';'delay';'response';'lick'});
save('E:\2P\summary\SC_choice_AUC_005.mat','SC_choice_AUC_005');

subplot(1,2,1);%for SC projecting M2 neurons
if strcmp(AUCtype,'choice')
    load('E:\2P\summary\M2_choice_AUC_005.mat');
    pAUCSigM2=zeros(3,4);
    pAUCSigM2(1,:)=M2_choice_AUC_005.contra_sig./M2_choice_AUC_005.total_n;
    pAUCSigM2(2,:)=M2_choice_AUC_005.ipsi_sig./M2_choice_AUC_005.total_n;
    pAUCSigM2(3,:)=1-sum(pAUCSigM2);
elseif strcmp(AUCtype,'sensory')
    pAUCSigM2=zeros(3,4);
end
bar(pAUCSigM2','stacked');
box off;
%set(gca,'XTick',1:4,'XTickLabel',{'sound','delay','response','lick'});
set(gca,'XTick',1:4,'XTickLabel',{'S','D','R','L'});
title('SC-projecting M2 neurons');
set(gca,'FontSize',12);
xlim([0,5]);
h=legend('contra','ipsi','n.s.');
set(h,'box','off');
set(groot,'defaultAxesColorOrder','remove');
set(gcf,'PaperPosition',[0,0,3,1.2]);
saveas(figBar,[savepath,filesep,celltype,'-',trialTypeStr,'-',AUCtype,'-bar plot of AUC at significance level of ',num2str(2*psig(1)),'.pdf'],'pdf');
%}
%% compare early vs. late AUC
mkr_ScatterHist={'h','>','s'};
color_ScatterHist={'F16820','646464','46782D'};%vglut2,vgat,ALM
color_ScatterHist=fHex2RGB(color_ScatterHist);
celltype={'vglut2','vgat','M2'};

[figScatterHist]=fScatterHistCmpEpochAUCCellType(TAUC_combine,celltype,'sound','late_delay',color_ScatterHist,mkr_ScatterHist);
saveas(figScatterHist,'E:\2P\summary\AUC_cmp_epoch\scatterHist_EI_EarlyLate_AUC.pdf','pdf');
[figBarShift,R_celltype1,P_celltype1]=fBarShiftCellType(TAUC_combine,celltype,'sound','late_delay',color_ScatterHist);
saveas(figBarShift,'E:\2P\summary\AUC_cmp_epoch\bar_EI_EarlyLate_AUC.pdf','pdf');
%compare mid delay vs. lick AUC
[figScatterHist]=fScatterHistCmpEpochAUCCellType(TAUC_combine,celltype,'mid_delay','lick',color_ScatterHist,mkr_ScatterHist);
saveas(figScatterHist,'E:\2P\summary\AUC_cmp_epoch\scatterHist_EI_DelayLick_AUC.pdf','pdf');
[figBarShift,R_celltype2,P_celltype2]=fBarShiftCellType(TAUC_combine,celltype,'delay','lick',color_ScatterHist);
saveas(figBarShift,'E:\2P\summary\AUC_cmp_epoch\bar_EI_DelayLick_AUC.pdf','pdf');

%show individual session cases
%{
chosen_fields1={'sound','delay'};
chosen_fields2={'late_delay','lick'};

fig_corr_celltype=figure();
for i_celltype=1:length(celltype)
%     fScatterHistCmpEarlyLateBySession(TAUC_combine,celltype{i_celltype},color_ScatterHist{i_celltype},mkr_ScatterHist{i_celltype});
    [R_sessions,P_sessions]=fScatterHistCmpBySessions(TAUC_combine,celltype{i_celltype},chosen_fields1,chosen_fields2,color_ScatterHist{i_celltype},mkr_ScatterHist{i_celltype});
    n_session=size(R_sessions,1);
    indsig1=P_sessions(:,1)<0.05;
    indsig2=P_sessions(:,2)<0.05;
    figure(fig_corr_celltype);
    set(gcf,'Position',[100,100,300,300]);
    scatter((i_celltype-0.3)*ones(n_session-sum(indsig1),1),R_sessions(~indsig1,1),20,color_ScatterHist{i_celltype},mkr_ScatterHist{i_celltype});
    hold on;
    scatter((i_celltype-0.3)*ones(sum(indsig1),1),R_sessions(indsig1,1),20,color_ScatterHist{i_celltype},mkr_ScatterHist{i_celltype},'filled');
    scatter((i_celltype+0.3)*ones(n_session-sum(indsig2),1),R_sessions(~indsig2,2),20,color_ScatterHist{i_celltype},mkr_ScatterHist{i_celltype});
    scatter((i_celltype+0.3)*ones(sum(indsig2),1),R_sessions(indsig2,2),20,color_ScatterHist{i_celltype},mkr_ScatterHist{i_celltype},'filled');
    plot([i_celltype-0.3;i_celltype+0.3]*ones(1,n_session), R_sessions','Color',color_ScatterHist{i_celltype});
    scatter((i_celltype-0.35),R_celltype1(i_celltype),30,color_ScatterHist{i_celltype},mkr_ScatterHist{i_celltype});%group cells together and calculate correlation
    scatter((i_celltype+0.35),R_celltype2(i_celltype),30,color_ScatterHist{i_celltype},mkr_ScatterHist{i_celltype});%group cells together and calculate correlation
end
figure(fig_corr_celltype);
set(gca,'XTick',1:length(celltype),'XTickLabel',celltype);
ylabel('correlation coefficient');
saveas(fig_corr_celltype,'E:\2P\summary\AUC_cmp_epoch\corrcoef_SCEIM2_epochAUC.pdf','pdf');
%}
%% plot moving AUC by session
T=cell2table(raw(2:end,1:14));
T.Properties.VariableNames=strrep(raw(1,1:14),' ','_');%table variable name can't have ' ',so replace them
celltype={'M2','vglut2','vgat'};
for i_celltype=1:length(celltype)
    ind_session=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*strcmp(T.cell_type,celltype{i_celltype});
    ind_session=find(ind_session);
    n_session=length(ind_session);
    figure;
    ncol=4;
    nrow=ceil(n_session/ncol);
    set(gcf,'Position',[200,200,200*ncol,200*nrow]);
    for i_session=1:n_session
        indrow=ind_session(i_session);
        subplot(nrow,ncol,i_session);
        fPlotMovingAUC(TAUC_combine,T.animal{indrow},T.date{indrow},trialTypeStr);
    end
    suptitle([trialTypeStr,'-',AUCtype]);
end

%% choice probability
%{
%choose cells with significant delay/sound choice AUC
trialTypeStr_datapool='cor and err';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
AUCtype_datapool='sensory';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
load([savepath,filesep,'trialType',trialTypeStr_datapool,'-',AUCtype_datapool,'TepochAUC.mat']);
celltype='syn';
T_AUCs=TAUC_combine(strcmp(TAUC_combine.celltype,celltype),:);
psig=[0,0.025;0.975,1];
TsigDataInfosensory = fGetSigDataInfoTable(T_AUCs, psig,AUCtype_datapool);

AUCtype_datapool='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
load([savepath,filesep,'trialType',trialTypeStr_datapool,'-',AUCtype_datapool,'TepochAUC.mat']);
T_AUCc=TAUC_combine(strcmp(TAUC_combine.celltype,celltype),:);
psig=[0,0.025;0.975,1];
TsigDataInfochoice = fGetSigDataInfoTable(T_AUCc, psig,AUCtype_datapool);

% T_sig=join(TsigDataInfosensory,TsigDataInfochoice,'Keys',{'animal','date','nROI'});%choose neurons with either choice or sensory selectivity
T_sig=TsigDataInfochoice;%choose neurons with only choice selectivity

trialTypeStr='cor';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
AUCtype='stimuli';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
load([savepath,filesep,'trialType',trialTypeStr,'-',AUCtype,'TepochAUC.mat']);
T_AUC=TAUC_combine(strcmp(TAUC_combine.celltype,celltype),:);
% T_AUC.ID=strcat(T_AUC.animal,'-',T_AUC.date,'-',int2str(T_AUC.nROI));

T_join=join(T_AUC,T_sig,'Keys',{'animal','date','nROI'});

[flag1,flag2]=fSigFlag(T_AUCs.sound,T_AUCs.psound,T_AUCc.sound,T_AUCc.psound,psig);
T_join=addvars(T_join,flag1,flag2,'NewVariableNames',{'flags1','flags2'});
[flag1,flag2]=fSigFlag(T_AUCs.delay,T_AUCs.pdelay,T_AUCc.delay,T_AUCc.pdelay,psig);
T_join=addvars(T_join,flag1,flag2,'NewVariableNames',{'flagd1','flagd2'});
[flag1,flag2]=fSigFlag(T_AUCs.response,T_AUCs.presponse,T_AUCc.response,T_AUCc.presponse,psig);
T_join=addvars(T_join,flag1,flag2,'NewVariableNames',{'flagr1','flagr2'});
[flag1,flag2]=fSigFlag(T_AUCs.lick,T_AUCs.plick,T_AUCc.lick,T_AUCc.plick,psig);
T_join=addvars(T_join,flag1,flag2,'NewVariableNames',{'flagl1','flagl2'});

figAUC=figure;
set(figAUC,'Position',[100,100,600,200]);
subplot(1,3,1);
% histogram(T_join.sound,'BinWidth',0.05);
hold on;
histogram([T_join.sound(T_join.flagCAUCsound>0,2);T_join.sound(T_join.flagCAUCsound<0,1)],'FaceColor','r','BinWidth',0.05);
histogram(T_join.sound(T_join.flagCAUCsound<0,1),'FaceColor','b','BinWidth',0.05);
% xlabel('AUC error v.s.correct for prefered stimulus')
ylabel('Cell number');
set(gca,'Xlim',[0,1]);
set(gca,'FontSize',12);
title('sound');

subplot(1,3,2);
% histogram(T_join.delay(:,2),'BinWidth',0.05);
hold on;
histogram([T_join.delay(T_join.flagCAUCdelay>0,2);T_join.delay(T_join.flagCAUCdelay<0,1)],'FaceColor','r','BinWidth',0.05);
histogram(T_join.delay(T_join.flagCAUCdelay<0,1),'FaceColor','b','BinWidth',0.05);
xlabel('AUC error v.s.correct for prefered stimulus')
ylabel('Cell number');
set(gca,'Xlim',[0,1]);
title('delay');
set(gca,'FontSize',12);

subplot(1,3,3);
% histogram(T_join.delay(:,2),'BinWidth',0.05);
hold on;
histogram([T_join.lick(T_join.flagCAUClick>0,2);T_join.lick(T_join.flagCAUClick<0,1)],'FaceColor','r','BinWidth',0.05);
histogram(T_join.lick(T_join.flagCAUClick<0,1),'FaceColor','b','BinWidth',0.05);
% xlabel('AUC error v.s.correct for prefered stimulus')
ylabel('Cell number');
set(gca,'Xlim',[0,1]);
set(gca,'FontSize',12);
% legend('n.s.','contra','ipsi');
h=legend('contra','ipsi');
set(h,'box','off');
title('lick');
set(figAUC,'PaperPosition',[0,0,6,2]);
saveas(figAUC,[savepath,filesep,celltype,'-',trialTypeStr,'-',AUCtype,'AUC error vs correct prefered stimulus.pdf'],'pdf');

figAUC2=figure;
set(figAUC2,'Position',[100,100,600,200]);
subplot(1,3,1);
% histogram(T_join.sound,'BinWidth',0.05);
hold on;
histogram([T_join.sound(T_join.flags2,2);T_join.sound(T_join.flags1,1)],'FaceColor','w','BinWidth',0.05);
histogram([T_join.sound(logical(T_join.flags2.*findSig(T_join.psound(:,2),psig).*(T_join.sound(:,2)>0.5)),2);T_join.sound(logical(T_join.flags1.*findSig(T_join.psound(:,1),psig).*(T_join.sound(:,1)>0.5)),1)],'FaceColor','r','BinWidth',0.05);
histogram([T_join.sound(logical(T_join.flags2.*findSig(T_join.psound(:,2),psig).*(T_join.sound(:,2)<0.5)),2);T_join.sound(logical(T_join.flags1.*findSig(T_join.psound(:,1),psig).*(T_join.sound(:,1)<0.5)),1)],'FaceColor','b','BinWidth',0.05);
ylabel('Cell number');
set(gca,'Xlim',[0,1]);
title('sound');
set(gca,'FontSize',12);

subplot(1,3,2);
% histogram(T_join.sound,'BinWidth',0.05);
hold on;
histogram([T_join.delay(T_join.flagd2,2);T_join.delay(T_join.flagd1,1)],'FaceColor','w','BinWidth',0.05);
histogram([T_join.delay(logical(T_join.flagd2.*findSig(T_join.pdelay(:,2),psig).*(T_join.delay(:,2)>0.5)),2);T_join.delay(logical(T_join.flagd1.*findSig(T_join.pdelay(:,1),psig).*(T_join.delay(:,1)>0.5)),1)],'FaceColor','r','BinWidth',0.05);
histogram([T_join.delay(logical(T_join.flagd2.*findSig(T_join.pdelay(:,2),psig).*(T_join.delay(:,2)<0.5)),2);T_join.delay(logical(T_join.flagd1.*findSig(T_join.pdelay(:,1),psig).*(T_join.delay(:,1)<0.5)),1)],'FaceColor','b','BinWidth',0.05);
ylabel('Cell number');
set(gca,'Xlim',[0,1]);
title('delay');
xlabel('AUC correct v.s.error for prefered stimulus');
set(gca,'FontSize',12);

subplot(1,3,3);
% histogram(T_join.sound,'BinWidth',0.05);
hold on;
histogram([T_join.lick(T_join.flagl2,2);T_join.lick(T_join.flagl1,1)],'FaceColor','w','BinWidth',0.05);
histogram([T_join.lick(logical(T_join.flagl2.*findSig(T_join.plick(:,2),psig).*(T_join.lick(:,2)>0.5)),2);T_join.lick(logical(T_join.flagl1.*findSig(T_join.plick(:,1),psig).*(T_join.lick(:,1)>0.5)),1)],'FaceColor','r','BinWidth',0.05);
histogram([T_join.lick(logical(T_join.flagl2.*findSig(T_join.plick(:,2),psig).*(T_join.lick(:,2)<0.5)),2);T_join.lick(logical(T_join.flagl1.*findSig(T_join.plick(:,1),psig).*(T_join.lick(:,1)<0.5)),1)],'FaceColor','b','BinWidth',0.05);
set(gca,'FontSize',12);

ylabel('Cell number');
set(gca,'Xlim',[0,1]);
title('lick');
h=legend('n.s.','correct','error');
set(h,'box','off');

set(figAUC2,'PaperPosition',[0,0,6,2]);
saveas(figAUC2,[savepath,filesep,celltype,'-',trialTypeStr,'-',AUCtype,'AUC cor vs err prefered stimulus significance p0.05.pdf'],'pdf');

figAUC3=figure;%M2 delay AUC
set(figAUC3,'Position',[100,100,200,200]);
load('D:\download\M2_delay_table.mat');
psig=[0,0.05];
[flag1,flag2]=fSigFlag(M2_delay_table.Sensory_auc,M2_delay_table.Sensory_p,M2_delay_table.Choice_auc,M2_delay_table.Choice_p,psig);
M2_delay_table=addvars(M2_delay_table,flag1,flag2,'NewVariableNames',{'flagd1','flagd2'});
histogram([M2_delay_table.CP_auc_contra(M2_delay_table.flagd2);M2_delay_table.CP_auc_ipsi(M2_delay_table.flagd1)],'FaceColor','w','BinWidth',0.05);
hold on;
histogram([M2_delay_table.CP_auc_contra(logical(M2_delay_table.flagd2.*findSig(M2_delay_table.CP_p_contra,psig).*(M2_delay_table.CP_auc_contra>0.5)));M2_delay_table.CP_auc_ipsi(logical(M2_delay_table.flagd1.*findSig(M2_delay_table.CP_p_ipsi,psig).*(M2_delay_table.CP_auc_ipsi>0.5)))],'FaceColor','r','BinWidth',0.05);
histogram([M2_delay_table.CP_auc_contra(logical(M2_delay_table.flagd2.*findSig(M2_delay_table.CP_p_contra,psig).*(M2_delay_table.CP_auc_contra<0.5)));M2_delay_table.CP_auc_ipsi(logical(M2_delay_table.flagd1.*findSig(M2_delay_table.CP_p_ipsi,psig).*(M2_delay_table.CP_auc_ipsi<0.5)))],'FaceColor','b','BinWidth',0.05);
ylabel('Cell number');
set(gca,'Xlim',[0,1]);
title('delay');
xlabel('AUC correct v.s.error for prefered stimulus');
set(gca,'FontSize',12);
box off;

%plot distribution of delay AUC M2 and SC together
figCmpM2SCdelay=figure;
set(figCmpM2SCdelay,'Position',[100,100,400,200]);
subplot(1,2,1);
histogram([T_join.delay(T_join.flagd2,2);T_join.delay(T_join.flagd1,1)],'FaceColor','r','BinWidth',0.05);
hold on;
histogram([M2_delay_table.CP_auc_contra(M2_delay_table.flagd2);M2_delay_table.CP_auc_ipsi(M2_delay_table.flagd1)],'FaceColor','b','BinWidth',0.05);
h=legend('SC','M2');
set(h,'box','off');
ylabel('Cell number');
set(gca,'FontSize',12);
set(gca,'Xlim',[0,1]);
box off;
subplot(1,2,2);
histogram([T_join.delay(T_join.flagd2,2);T_join.delay(T_join.flagd1,1)],'EdgeColor','r','BinWidth',0.05,'Normalization','probability','DisplayStyle','stairs');
hold on;
histogram([M2_delay_table.CP_auc_contra(M2_delay_table.flagd2);M2_delay_table.CP_auc_ipsi(M2_delay_table.flagd1)],'EdgeColor','b','BinWidth',0.05,'Normalization','probability','DisplayStyle','stairs');
plot([0.5,0.5],[0,0.3],'k--');
h=legend('SC','M2');
set(h,'box','off');
ylabel('Probabilty');
set(gca,'FontSize',12);
set(gca,'Xlim',[0,1]);
box off;
%}
%% functions
function [figAUC] =fHistEpochAUC(T_AUC,psig,threshold,xlabel_str,varargin)
%plot epoch AUC/activities amplitude with psig indicating significance
%level
if strcmp(xlabel_str,'AUC')
    bin_width=0.05;
else
    if ~isempty(varargin) && length(varargin)>=2
        bin_width=varargin{2};
    else
        bin_width=round((max(T_AUC.lick)-min(T_AUC.lick))/20,2);
    end
end

figAUC=figure;
n_animal=length(unique(T_AUC.animal));
set(figAUC,'Position',[100,100,800,200]);
subplot(1,5,1);
histogram(T_AUC.ITI,'FaceColor','w','BinWidth',bin_width);
hold on;
histogram(T_AUC.ITI(findSig(T_AUC.pITI,psig) & T_AUC.ITI>threshold),'FaceColor','r','BinWidth',bin_width);
histogram(T_AUC.ITI(findSig(T_AUC.pITI,psig) & T_AUC.ITI<threshold),'FaceColor','b','BinWidth',bin_width);
ylabel('cell counts');
title('ITI');
text(0.1,1,['n=',num2str(size(T_AUC,1)),'cells from ',num2str(n_animal),'animals'],'Unit','Normalized');

box off;
subplot(1,5,2);
histogram(T_AUC.sound,'FaceColor','w','BinWidth',bin_width);
hold on;
histogram(T_AUC.sound(findSig(T_AUC.psound,psig) & T_AUC.sound>threshold),'FaceColor','r','BinWidth',bin_width);
histogram(T_AUC.sound(findSig(T_AUC.psound,psig) & T_AUC.sound<threshold),'FaceColor','b','BinWidth',bin_width);
title('sound');

subplot(1,5,3);
histogram(T_AUC.delay,'FaceColor','w','BinWidth',bin_width);
hold on;
histogram(T_AUC.delay(findSig(T_AUC.pdelay,psig) & T_AUC.delay>threshold),'FaceColor','r','BinWidth',bin_width);
histogram(T_AUC.delay(findSig(T_AUC.pdelay,psig) & T_AUC.delay<threshold),'FaceColor','b','BinWidth',bin_width);
title('delay');

subplot(1,5,4);
histogram(T_AUC.response,'FaceColor','w','BinWidth',bin_width);
hold on;
histogram(T_AUC.response(findSig(T_AUC.presponse,psig) & T_AUC.response>threshold),'FaceColor','r','BinWidth',bin_width);
histogram(T_AUC.response(findSig(T_AUC.presponse,psig) & T_AUC.response<threshold),'FaceColor','b','BinWidth',bin_width);
title('response');

subplot(1,5,5);
histogram(T_AUC.lick,'FaceColor','w','BinWidth',bin_width);
hold on;
histogram(T_AUC.lick(findSig(T_AUC.plick,psig) & T_AUC.lick>threshold),'FaceColor','r','BinWidth',bin_width);
histogram(T_AUC.lick(findSig(T_AUC.plick,psig) & T_AUC.lick<threshold),'FaceColor','b','BinWidth',bin_width);
title('lick');

if strcmp(xlabel_str,'AUC')
    for i=1:5
        subplot(1,5,i);
        set(gca,'Xlim',[0,1]);
        xlabel(xlabel_str);
        if i==1
            h=legend('n.s.','contra','ipsi','AutoUpdate','off');
            set(h,'box','off');
        end
        box off;
        set(gca,'FontSize',12);
    end  
else
    if ~isempty(varargin)
        x_lim=varargin{1};
    else
        x_lim=[0,1];
        for i=1:5
            subplot(1,5,i);
            tempx_lim=get(gca,'Xlim');
            x_lim(1)=min(x_lim(1),tempx_lim(1));
            x_lim(2)=max(x_lim(2),tempx_lim(2));
        end
    end
    for i=1:5
        subplot(1,5,i);
        xlabel(xlabel_str);
        if i==1
            h=legend('n.s.','activated','inhibited','AutoUpdate','off');
            set(h,'box','off');
        end
        set(gca,'Xlim',x_lim);
        box off;
        set(gca,'FontSize',12);
    end
end

end

function [Tout] = fExtractTableField(Tin,fieldName)
size_in=size(table2cell(Tin));
Cout=cell(size_in);

%extract the field name of some variables in a table
n_var=length(Tin.Properties.VariableNames);
for i_var=1:n_var
    if isstruct(Tin{1,i_var})
        temp=Tin{:,i_var};
        switch fieldName
            case 'ipsi'
                Cout(:,i_var)=arrayfun(@(x) x.ipsi, temp,'UniformOutput',false);
            case 'contra'
                Cout(:,i_var)=arrayfun(@(x) x.contra, temp,'UniformOutput',false);
            case 'mean'
                Cout(:,i_var)=arrayfun(@(x) x.mean, temp,'UniformOutput',false);
        end
    else
        Cout(:,i_var)=table2cell(Tin(:,i_var));
    end
end
Tout=cell2table(Cout,'VariableNames',Tin.Properties.VariableNames);
end

function [patch_array_x,patch_array_y]=hist2patch(binEdges,values)
patch_array_x=cell(1,length(values));
patch_array_y=cell(1,length(values));
for i_cell=1:length(values)
    patch_array_x{i_cell}=[binEdges(i_cell),binEdges(i_cell+1),binEdges(i_cell+1),binEdges(i_cell)];
    patch_array_y{i_cell}=[0,0,values(i_cell),values(i_cell)];
end
end

function [flag1,flag2]=fSigFlag(SAUC,pSAUC,CAUC,pCAUC,p_range)
%p_range- the range of significance p value
if length(SAUC)~=length(CAUC) 
    warning('unequal rows of input table');
    return;
end
[indsigS,indsigC,flag1,flag2]=deal(zeros(length(SAUC),1));
for i=1:size(p_range,1) %每一行一个区间段，各行之间取并集
    indsigS=logical((pSAUC<=p_range(i,2)).*(pSAUC>=p_range(i,1))+indsigS);
    indsigC=logical((pCAUC<=p_range(i,2)).*(pCAUC>=p_range(i,1))+indsigC);
end
for irow=1:length(SAUC)
    if (indsigS(irow) && ~indsigC(irow)) || (indsigS(irow) && indsigC(irow) && abs(SAUC(irow)-0.5)>=abs(CAUC(irow)-0.5)) 
        if SAUC(irow)>0.5
            flag2(irow)=1;
        else
            flag1(irow)=1;
        end
    elseif (~indsigS(irow) && indsigC(irow))|| (indsigS(irow) && indsigC(irow) && abs(SAUC(irow)-0.5)<abs(CAUC(irow)-0.5)) 
        if CAUC(irow)>0.5
            flag2(irow)=1;
        else
            flag1(irow)=1;
        end
    end
end
flag1=logical(flag1);
flag2=logical(flag2);
end

function [TsigDataInfo] = fGetSigDataInfoTable(Tin, p_range,strAUCtype)
% summary the significance table using Tin and setted range of significant
% p value
[indsigsound,indsigdelay,indsigresponse,indsiglick]=deal(zeros(size(Tin,1),1));
for i=1:size(p_range,1) %每一行一个区间段，各行之间取并集
    indsigsound=logical((Tin.psound<=p_range(i,2)).*(Tin.psound>=p_range(i,1))+indsigsound);
    indsigdelay=logical((Tin.pdelay<=p_range(i,2)).*(Tin.pdelay>=p_range(i,1))+indsigdelay);
    indsigresponse=logical((Tin.presponse<=p_range(i,2)).*(Tin.presponse>=p_range(i,1))+indsigresponse);
    indsiglick=logical((Tin.plick<=p_range(i,2)).*(Tin.plick>=p_range(i,1))+indsiglick);
end
indsigContraSound=logical(indsigsound.*(Tin.sound>0.5));
indsigIpsiSound=logical(indsigsound.*(Tin.sound<0.5));
indNSSound= ~indsigsound;
indsigContraDelay=logical(indsigdelay.*(Tin.delay>0.5));
indsigIpsiDelay=logical(indsigdelay.*(Tin.delay<0.5));
indNSDelay= ~indsigdelay;
indsigContraResponse=logical(indsigresponse.*(Tin.response>0.5));
indsigIpsiResponse=logical(indsigresponse.*(Tin.response<0.5));
indNSResponse= ~indsigresponse;
indsigContraLick=logical(indsiglick.*(Tin.lick>0.5));
indsigIpsiLick=logical(indsiglick.*(Tin.lick<0.5));
indNSLick= ~indsiglick;
TsigDataInfo=Tin(:,{'animal','date','nROI'});
if strcmp(strAUCtype,'sensory')
    TsigDataInfo.flagSAUCsound(indsigContraSound)=1;
    TsigDataInfo.flagSAUCsound(indsigIpsiSound)=-1;
    TsigDataInfo.flagSAUCsound(indNSSound)=0;
    TsigDataInfo.flagSAUCdelay(indsigContraDelay)=1;
    TsigDataInfo.flagSAUCdelay(indsigIpsiDelay)=-1;
    TsigDataInfo.flagSAUCdelay(indNSDelay)=0;
    TsigDataInfo.flagSAUCresponse(indsigContraResponse)=1;
    TsigDataInfo.flagSAUCresponse(indsigIpsiResponse)=-1;
    TsigDataInfo.flagSAUCresponse(indNSResponse)=0;
    TsigDataInfo.flagSAUClick(indsigContraLick)=1;
    TsigDataInfo.flagSAUClick(indsigIpsiLick)=-1;
    TsigDataInfo.flagSAUClick(indNSLick)=0;
elseif strcmp(strAUCtype,'choice')
    TsigDataInfo.flagCAUCsound(indsigContraSound)=1;
    TsigDataInfo.flagCAUCsound(indsigIpsiSound)=-1;
    TsigDataInfo.flagCAUCsound(indNSSound)=0;
    TsigDataInfo.flagCAUCdelay(indsigContraDelay)=1;
    TsigDataInfo.flagCAUCdelay(indsigIpsiDelay)=-1;
    TsigDataInfo.flagCAUCdelay(indNSDelay)=0;
    TsigDataInfo.flagCAUCresponse(indsigContraResponse)=1;
    TsigDataInfo.flagCAUCresponse(indsigIpsiResponse)=-1;
    TsigDataInfo.flagCAUCresponse(indNSResponse)=0;
    TsigDataInfo.flagCAUClick(indsigContraLick)=1;
    TsigDataInfo.flagCAUClick(indsigIpsiLick)=-1;
    TsigDataInfo.flagCAUClick(indNSLick)=0;
end
end

function [figScatterHist,R_out,P_out]=fScatterHistCmpASession(Tin,chosen_field1,chosen_field2,color,mkr)
data1=Tin{:,chosen_field1};
data2=Tin{:,chosen_field2};
figScatterHist=figure;
set(gcf,'Position',[100,100,400,400]);
h1=subplot(3,3,[2,3,5,6]);
hold on;
plot([0,1],[0.5,0.5],'k--');
plot([0.5,0.5],[0,1],'k--');
patch_color=fHex2RGB({'646464'});
patch([0,0.5,0.5,0],[0.5,0.5,1,1],patch_color{1},'FaceAlpha',0.3,'EdgeColor','none');
patch([0.5,0.5,1,1],[0,0.5,0.5,0],patch_color{1},'FaceAlpha',0.3,'EdgeColor','none');
%h = scatterhist(TAUC.AUCearly,TAUC.AUClate,'Group',TAUC.celltype,'Color',color,'Marker',mkr);
curve_scatter=scatter(h1,data1,data2,30,color,mkr);hold on;
%plot linear fitting showing relationship of x and y
% f=fit(data1,data2,'poly1');
% plot(f);
[R,P]=corrcoef(data1,data2);
R_out=R(2,1);
P_out=P(2,1);
%calculate how much ROI shift preference
p_shift=sum((data1-0.5).*(data2-0.5)<0)/length(data1);

set(gca,'Xlim',[0,1],'Ylim',[0,1],'FontSize',14,'FontName','Arial');
set(h1,'position',[0.35,0.35,0.6,0.6]);
text(0.1,0.9,{['ipsi ',num2str(Tin.ipsi_performance{1})];['contra ',num2str(Tin.contra_performance{1})]},'Unit','Normalized');
text(0.7,0.9,{['rho=',num2str(R(2,1))];['p=',num2str(P(2,1))]},'Unit','Normalized');
text(0.1,0.1,['P(shift preference)=',num2str(p_shift)],'Unit','Normalized');
% scatter(h1,TAUC.sound(indExample),TAUC.late_delay(indExample),30,color,mkr,'filled');

h2=subplot(3,3,[8,9]);
boxplot(data1,Tin.celltype,'orientation','horizontal','color',color,'Notch','on');
xlabel([strrep(chosen_field1,'_','\_'),' AUC']);
[h,p1]=ttest(data1-0.5);
text(h2,0.9,0.5,plabelsymbol(p1),'Unit','Normalized');
set(h2,'Xlim',[0,1],'FontSize',14,'FontName','Arial','box','off','XTickLabel',{' '},'YTickLabel',{' '});
set(h2,'position',[0.35,0.15,0.6,0.1]);

h3=subplot(3,3,[1,4]);
boxplot(data2,Tin.celltype,'orientation','vertical','color',color,'Notch','on');
[h,p2]=ttest(data2-0.5);
text(h3,0.5,0.9,plabelsymbol(p2),'Unit','Normalized');
set(h3,'Ylim',[0,1],'FontSize',14,'FontName','Arial','box','off','XTickLabel',{' '},'YTickLabel',{' '});
ylabel([strrep(chosen_field2,'_','\_'),' AUC']);
set(h3,'position',[0.15,0.35,0.1,0.6]);

suptitle(strrep(char(Tin.session(1)),'_','\_'));
end

function [R_sessions,P_sessions]=fScatterHistCmpBySessions(Tin,chosen_celltype,chosen_fields1,chosen_fields2,color,mkr)
%plot scatter of early vs. late delay activities session by session
%chosen_celltype- string indicating the cell type
%chosen_fields1/2- cellarrays indicating the fields to be compared
celltype=categorical(Tin.celltype);
Tin.celltype=[];
Tin=addvars(Tin,celltype,'After','field');
ind_shown=logical((Tin.celltype==chosen_celltype));
TAUC_celltype=Tin(ind_shown,:);
TAUC_celltype.session=categorical(strcat(TAUC_celltype.field,'_',strrep(TAUC_celltype.date,'/','-')));
session_pool=unique(TAUC_celltype.session);

if length(chosen_fields1)~=length(chosen_fields2)
    warning('different lengths of chosen fields to be compared');
    return;
end
[R_sessions,P_sessions]=deal(zeros(length(session_pool),length(chosen_fields1)));
for i=1:length(session_pool)
    %choose data from one session
    TAUC=TAUC_celltype(TAUC_celltype.session==session_pool(i),:);
    for i_field=1:length(chosen_fields1)
        field1=chosen_fields1{i_field};
        field2=chosen_fields2{i_field};
        [figScatterHist(i_field),R_sessions(i,i_field),P_sessions(i,i_field)]=fScatterHistCmpASession(TAUC,field1,field2,color,mkr);
        saveas(figScatterHist(i_field),['E:\2P\summary\AUC_cmp_epoch\cases\',char(session_pool(i)),'_',field1,'_vs_',field2,'.png'],'png');
    end
    %save figure to PPT using exportToPPT
    pptFileName='E:\2P\summary\AUC_cmp_epoch\cases\cases.pptx';
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
            'Title','Comparing AUC from different epochs', ...
            'Author','PYX', ...
            'Comments','This file has been automatically generated by exportToPPTX');
    end
    slideId = exportToPPTX('addslide');
    n_fig=length(chosen_fields1);
    for i_field=1:n_fig
        exportToPPTX('addpicture',figScatterHist(i_field),'Position',[1+(slidesize(1)-1)*(i_field-1)/n_fig 1 slidesize(1)/n_fig/2 slidesize(1)/n_fig/2]);
        exportToPPTX('addtext',char(session_pool(i)),'Position',[1+(slidesize(1)-1)*(i_field-1)/n_fig 0 3 0.5],'Vert','top');
    end
    exportToPPTX('addtext',char(chosen_celltype),'Position',[slidesize(1)-1 0 3 0.5],'Vert','top');

    % Save and close (in one command)
    newFile = exportToPPTX('saveandclose',pptFileName);
    fprintf('New file has been saved: <a href="matlab:open(''%s'')">%s</a>\n',newFile,newFile);
    close(figScatterHist(:));
    
end
end

function [R_sessions,P_sessions]=fScatterHistCmpEarlyLateBySession(Tin,chosen_celltype,color,mkr)
%plot scatter of early vs. late delay activities session by session
%chosen_celltype- string indicating the cell type

celltype=categorical(Tin.celltype);
Tin.celltype=[];
Tin=addvars(Tin,celltype,'After','field');
ind_shown=logical((Tin.celltype==chosen_celltype));
TAUC_celltype=Tin(ind_shown,:);

TAUC_celltype.session=categorical(strcat(TAUC_celltype.field,'_',strrep(TAUC_celltype.date,'/','-')));
session_pool=unique(TAUC_celltype.session);
[R_sessions,P_sessions]=deal(zeros(size(session_pool)));
for i=1:length(session_pool)
    
    %choose data from one session
    TAUC=TAUC_celltype(TAUC_celltype.session==session_pool(i),:);
    [figScatterHist,R_sessions(i),P_sessions(i)]=fScatterHistCmpASession(TAUC,'sound','late_delay',color,mkr);
    saveas(figScatterHist,['E:\2P\summary\AUC_early_vs_late\cases\',char(session_pool(i)),'.png'],'png');
    %save figure to PPT
    saveppt('E:\2P\summary\AUC_early_vs_late\cases\cases',char(session_pool(i)),figScatterHist);
    close(figScatterHist);
end
end

function [figScatterHist]=fScatterHistCmpEpochAUCCellType(Tin,chosen_celltype,chosen_fields1,chosen_fields2,color,mkr)
figScatterHist=figure;
set(gcf,'Position',[100,100,600,600]);
h1=subplot(3,3,[2,3,5,6]);
hold on;
plot([0,1],[0.5,0.5],'k--');
plot([0.5,0.5],[0,1],'k--');
patch([0,0.5,0.5,0],[0.5,0.5,1,1],color{2},'FaceAlpha',0.3,'EdgeColor','none');
patch([0.5,0.5,1,1],[0,0.5,0.5,0],color{2},'FaceAlpha',0.3,'EdgeColor','none');
ind_all_shown=zeros(size(Tin,1),1);
for icelltype=length(chosen_celltype):-1:1
    celltype=categorical(Tin.celltype);
    Tin.celltype=[];
    Tin=addvars(Tin,celltype,'After','field');
    ind_shown=logical((Tin.celltype==chosen_celltype{icelltype}));
    TAUC=Tin(ind_shown,:);
    ind_all_shown=ind_all_shown+ind_shown;
    data1=TAUC{:,chosen_fields1};
    data2=TAUC{:,chosen_fields2};
        
%     h = scatterhist(TAUC.AUCearly,TAUC.AUClate,'Group',TAUC.celltype,'Color',color{icelltype},'Marker',mkr{icelltype});
    curve_scatter(icelltype)=scatter(h1,data1,data2,10,color{icelltype},mkr{icelltype});hold on;
    set(gca,'Xlim',[0,1],'Ylim',[0,1],'FontSize',14,'FontName','Arial');
    set(h1,'position',[0.35,0.35,0.6,0.6]);

%     scatter(h1,TAUC.sound(indExample),TAUC.late_delay(indExample),30,color{icelltype},mkr{icelltype},'filled');
    h2=subplot(3,3,[8,9]);
    boxplot(data1,char(TAUC.celltype),'orientation','horizontal','color',color{icelltype},'Notch','on');
%     hold on;
    set(h2,'Xlim',[0,1],'FontSize',14,'FontName','Arial','box','off','XTickLabel',{' '},'YTickLabel',{' '});
    set(h2,'position',[0.35,0.05,0.6,0.2]);
    h3=subplot(3,3,[1,4]);
    boxplot(data2,char(TAUC.celltype),'orientation','vertical','color',color{icelltype},'Notch','on');
%     hold on;
    set(h3,'Ylim',[0,1],'FontSize',14,'FontName','Arial','box','off','XTickLabel',{' '},'YTickLabel',{' '});
    set(h3,'position',[0.05,0.35,0.2,0.6]);
end
%replot boxplot together
TAUC=Tin(logical(ind_all_shown),:);
data1=TAUC{:,chosen_fields1};
data2=TAUC{:,chosen_fields2};
if length(chosen_celltype)==2
    boxplot(h2,data1,TAUC.celltype,'orientation','horizontal','Notch','on','GroupOrder',chosen_celltype,'colors',[color{1};color{2}]);
    early1=data1(TAUC.celltype==chosen_celltype{1});
    early2=data1(TAUC.celltype==chosen_celltype{2});
    [h,p1]=ttest2(early1,early2);
    text(h2,0.9,0.5,plabelsymbol(p1),'Unit','Normalized');
    set(h2,'Xlim',[0,1],'FontSize',14,'FontName','Arial','box','off','XTickLabel',{' '},'YTickLabel',{' '});
    set(h2,'position',[0.35,0.05,0.6,0.2]);
    
    boxplot(h3,data2,TAUC.celltype,'orientation','vertical','Notch','on','GroupOrder',chosen_celltype,'colors',[color{1};color{2}]);
    late1=data2(TAUC.celltype==chosen_celltype{1});
    late2=data2(TAUC.celltype==chosen_celltype{2});
    [h,p2]=ttest2(late1,late2);
    text(h3,0.5,0.9,plabelsymbol(p2),'Unit','Normalized');
    set(h3,'Ylim',[0,1],'FontSize',14,'FontName','Arial','box','off','XTickLabel',{' '},'YTickLabel',{' '});
    set(h3,'position',[0.05,0.35,0.2,0.6]);
    
    % hl=legend(curve_scatter(:),'vglut2' , 'vgat' , 'ALM terminal','AutoUpdate','off');
    hl=legend(curve_scatter(:), 'vglut2','vgat' ,'AutoUpdate','off','Position',[0.1,0.1,0.1,0.1]);
    set(hl,'Box','Off');
    % text(0,1,['pearson correlation correct=',num2str(corr(TAUC.AUCearly,TAUC.AUClate))]);
elseif length(chosen_celltype)==3
    boxplot(h2,data1,TAUC.celltype,'orientation','horizontal','Notch','on','GroupOrder',chosen_celltype,'colors',[color{1};color{2};color{3}]);
    early1=data1(TAUC.celltype==chosen_celltype{1});
    early2=data1(TAUC.celltype==chosen_celltype{2});
    [h,p1]=ttest2(early1,early2);
    text(h2,0.9,0.5,plabelsymbol(p1),'Unit','Normalized');
%     set(h2,'Xlim',[0,1],'FontSize',14,'FontName','Arial','box','off','XTickLabel',{' '},'YTickLabel',{' '});
    set(h2,'Xlim',[0,1],'FontSize',14,'FontName','Arial','box','off');
    set(h2,'position',[0.35,0.05,0.6,0.2]);
    
    boxplot(h3,data2,TAUC.celltype,'orientation','vertical','Notch','on','GroupOrder',chosen_celltype,'colors',[color{1};color{2};color{3}]);
    late1=data2(TAUC.celltype==chosen_celltype{1});
    late2=data2(TAUC.celltype==chosen_celltype{2});
    [h,p2]=ttest2(late1,late2);
    text(h3,0.5,0.9,plabelsymbol(p2),'Unit','Normalized');
%     set(h3,'Ylim',[0,1],'FontSize',14,'FontName','Arial','box','off','XTickLabel',{' '},'YTickLabel',{' '});
    set(h3,'Ylim',[0,1],'FontSize',14,'FontName','Arial','box','off');
    set(h3,'position',[0.05,0.35,0.2,0.6]);
    
    % hl=legend(curve_scatter(:),'vglut2' , 'vgat' , 'ALM terminal','AutoUpdate','off');
    hl=legend(curve_scatter(:), 'vglut2','vgat' ,'AutoUpdate','off','Position',[0.1,0.1,0.1,0.1]);
    set(hl,'Box','Off');
    % text(0,1,['pearson correlation correct=',num2str(corr(TAUC.AUCearly,TAUC.AUClate))]);
end

end

function [figBarShift,R_celltype,P_celltype]=fBarShiftCellType(Tin,chosen_celltype,chosen_fields1,chosen_fields2,color)
figBarShift=figure;
set(gcf,'Position',[100,100,300,300]);
plot([0,1],[0.5,0.5],'k--');
plot([0.5,0.5],[0,1],'k--');
patch([0,0.5,0.5,0],[0.5,0.5,1,1],color{2},'FaceAlpha',0.3,'EdgeColor','none');
patch([0.5,0.5,1,1],[0,0.5,0.5,0],color{2},'FaceAlpha',0.3,'EdgeColor','none');
[R_celltype,P_celltype]=deal(zeros(size(chosen_celltype)));
xstr=cell(size(chosen_celltype));
for icelltype=1:length(chosen_celltype)
    celltype=categorical(Tin.celltype);
    Tin.celltype=[];
    Tin=addvars(Tin,celltype,'After','field');
    ind_shown=logical((Tin.celltype==chosen_celltype{icelltype}));
    TAUC=Tin(ind_shown,:);
    data1=TAUC{:,chosen_fields1};
    data2=TAUC{:,chosen_fields2};
    %plot linear fitting showing relationship of x and y
    [R,P]=corrcoef(data1,data2);
    R_celltype(icelltype)=R(2,1);
    P_celltype(icelltype)=P(2,1);
    %calculate how much ROI shift preference
    n_shift=sum((data1-0.5).*(data2-0.5)<0);
    p_shift=n_shift/length(data1);
    bar(icelltype,1-p_shift,'FaceColor',color{icelltype},'EdgeColor',color{icelltype});
    text(icelltype,0.1,[num2str(n_shift),'/',num2str(length(data1))]);
    hold on;
end

plot([0,length(chosen_celltype)+1],[0.5,0.5],'k--');
set(gca,'Ylim',[0,1],'XTick',1:length(chosen_celltype),'XTickLabel',chosen_celltype);
ylabel('% consistent sessions');
title_str=['AUC change ',chosen_fields1,' vs. ',chosen_fields2];
title(strrep(title_str,'_','\_'));
for icelltype=1:length(chosen_celltype)
    ind_i=logical(Tin.celltype==chosen_celltype{icelltype});
    TAUC=Tin(ind_i,:);
    data1=TAUC{:,chosen_fields1};
    data2=TAUC{:,chosen_fields2};
    n_shift=sum((data1-0.5).*(data2-0.5)<0);
    n_i=sum(ind_i);
    for jcelltype=icelltype+1:length(chosen_celltype)
        ind_j=logical(Tin.celltype==chosen_celltype{jcelltype});
        n_j=sum(ind_j);
        TAUCj=Tin(ind_j,:);
        data1j=TAUCj{:,chosen_fields1};
        data2j=TAUCj{:,chosen_fields2};
        n_shiftj=sum((data1j-0.5).*(data2j-0.5)<0);
        [h,p, chi2stat,df] = prop_test([n_shift,n_shiftj] , [n_i,n_j], false);
        plot([icelltype+0.1,jcelltype-0.1],(0.8+(jcelltype-icelltype)*0.05)*ones(1,2),'k-');
        text((icelltype+jcelltype)/2,0.8+(jcelltype-icelltype+1)*0.05,plabelsymbol(p));
    end
end
box off;
end

function [figMovingAUC] = fPlotMovingAUC(TAUC_combine,animal,date,trialTypeStr,figMovingAUC)
indrow=logical(strcmp(TAUC_combine.animal,animal).*strcmp(TAUC_combine.date,date));
TAUC=TAUC_combine(indrow,:);
nROI=size(TAUC,1);
frT =  34.2163;%SavedCaTrials.FrameTime;usually 3X
frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
frameNum=double(round(frameNumTime*1000/frT));

hold on;
for i=1:nROI
    if strcmp(trialTypeStr,'cor')
        AUC=TAUC.delayMovingAUC{i,1}.cor;
    elseif strcmp(trialTypeStr,'cor and err')
        AUC=TAUC.delayMovingAUC{i,1}.do;
    end
    ts=double((-frameNum(1):frameNum(2))*frT/1000);
    if length(ts)==length(AUC)
        plot(ts,AUC);
    else
        lim=min(length(ts),length(AUC));
        plot(ts(1:lim),AUC(1:lim));
        disp(strcat(animal,date,'-ROI-',num2str(i),'-ts-',num2str(length(ts)),'-AUC-',num2str(length(AUC))));
    end
end
title([animal,date,'-',num2str(i),'ROIs']);
ylim([0,1]);
plot([-frameNumTime(1),frameNumTime(2)],[0.5,0.5],'k-');
xlabel('Time (s) from delay onset');
ylabel('AUC');
box off;
set(gca,'FontSize',12);

end


function [ind]=findSig(var,pSig)
ind=zeros(length(var),1);
for i=1:size(pSig,1) %每一行一个区间段，各行之间取并集
    ind=logical((var<=pSig(i,2)).*(var>=pSig(i,1))+ind);
end
end