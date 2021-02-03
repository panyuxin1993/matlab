close all;
[num,txt,raw] =xlsread('D:\xulab\project\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
savepath='H:\2P\summary';
clear TAUC_combine;
trialTypeStr='cor and err';%{'cor and err','do','cor'}
AUCtype='sensory';%{'choice','sensory'};
%{
T=cell2table(raw(2:end,1:11));
T.Properties.VariableNames=strrep(raw(1,1:11),' ','_');%table variable name can't have ' ',so replace them
ind_session=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control');
ind_session=find(ind_session);
n_session=length(ind_session);
animal_unique=unique(T.animal(ind_session));
n_animal=length(animal_unique);

for i_session=n_session:-1:1
    indrow=ind_session(i_session);
    TAUC_currentSession = fGetEpochAUCtableASession(T.file_path{indrow},T.animal{indrow},T.date{indrow},T.field{indrow},T.cell_type{indrow},trialTypeStr,AUCtype);
    if exist('TAUC_combine','var')
        TAUC_combine=vertcat(TAUC_combine,TAUC_currentSession);
    else
        TAUC_combine=TAUC_currentSession;
    end
end
save([savepath,filesep,'trialType',trialTypeStr,'-',AUCtype,'TepochAUC.mat'],'TAUC_combine');
%}
%% plot histogram
load([savepath,filesep,'trialType',trialTypeStr,'-',AUCtype,'TepochAUC.mat']);
% load('F:\2P\summary\trialTypecorTepochAUC.mat');
celltype='syn';
manipulation='control';
%T_AUC=TAUC_combine(strcmp(TAUC_combine.celltype,'syn')|contains(TAUC_combine.celltype,'+other'),:);%syn, ..+other
T_AUC=TAUC_combine(strcmp(TAUC_combine.celltype,celltype),:);
n_animal=length(unique(T_AUC.animal));
psig=[0.025,0.975];
figAUC=figure;
set(figAUC,'Position',[100,100,800,200]);
subplot(1,5,1);
histogram(T_AUC.ITI,'BinWidth',0.05);
hold on;
histogram(T_AUC.ITI(T_AUC.pITI<psig(1)),'FaceColor','r','BinWidth',0.05);
histogram(T_AUC.ITI(T_AUC.pITI>psig(2)),'FaceColor','m','BinWidth',0.05);
ylabel('cell counts');
xlabel('AUC');
title('ITI');
set(gca,'Xlim',[0,1]);

box off;
subplot(1,5,2);
histogram(T_AUC.sound,'BinWidth',0.05);
hold on;
histogram(T_AUC.sound(T_AUC.psound<psig(1)),'FaceColor','r','BinWidth',0.05);
histogram(T_AUC.sound(T_AUC.psound>psig(2)),'FaceColor','m','BinWidth',0.05);
title('sound');
xlabel('AUC');
set(gca,'Xlim',[0,1]);
box off;
subplot(1,5,3);
histogram(T_AUC.delay,'BinWidth',0.05);
hold on;
histogram(T_AUC.delay(T_AUC.pdelay<psig(1)),'FaceColor','r','BinWidth',0.05);
histogram(T_AUC.delay(T_AUC.pdelay>psig(2)),'FaceColor','m','BinWidth',0.05);
title('delay');
xlabel('AUC');
set(gca,'Xlim',[0,1]);
box off;
subplot(1,5,4);
histogram(T_AUC.response,'BinWidth',0.05);
hold on;
histogram(T_AUC.response(T_AUC.presponse<psig(1)),'FaceColor','r','BinWidth',0.05);
histogram(T_AUC.response(T_AUC.presponse>psig(2)),'FaceColor','m','BinWidth',0.05);
title('response');
xlabel('AUC');
set(gca,'Xlim',[0,1]);
box off;
subplot(1,5,5);
histogram(T_AUC.lick,'BinWidth',0.05);
hold on;
histogram(T_AUC.lick(T_AUC.plick<psig(1)),'FaceColor','r','BinWidth',0.05);
histogram(T_AUC.lick(T_AUC.plick>psig(2)),'FaceColor','m','BinWidth',0.05);
title('lick');
xlabel('AUC');
text(0.1,0.9,['n=',num2str(size(T_AUC,1)),'cells from ',num2str(n_animal),'animals'],'Unit','Normalized');
set(gca,'Xlim',[0,1]);
box off;
%indcate example as vertical line
Tsorted=sortrows(TAUC_combine,{'delay','response'},{'ascend','descend'});
% indCaseRow=49;%shift from ipsi to contra
% for i=1:5
%     subplot(1,5,i);
%     plot( table2array(Tsorted(indCaseRow,5+i))*ones(2,1),ylim,'k-');hold on;
%     box off;
% end
% Tsorted=sortrows(TAUC_combine,{'delay','response'},{'descend','descend'});
% indCaseRow=1;%shift from contra to contra
% for i=1:5
%     subplot(1,5,i);
%     plot( table2array(Tsorted(indCaseRow,5+i))*ones(2,1),ylim,'r-');
%     box off;
% end
%save figure
saveas(figAUC,[savepath,filesep,celltype,'-',trialTypeStr,'-',AUCtype,'AUC at significance level of ',num2str(2*psig(1)),'.png'],'png');
set(figAUC,'PaperPosition',[0,0,8,2]);
saveas(figAUC,[savepath,filesep,celltype,'-',trialTypeStr,'-',AUCtype,'AUC at significance level of ',num2str(2*psig(1)),'.pdf'],'pdf');

%% plot bar
n_cell= size(T_AUC,1);
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
set(gcf,'PaperPosition',[0,0,3,1.2]);
subplot(1,2,2);
bar(pAUCSigSC','stacked');
map=[1,0,0;0,0,1;0.5,0.5,0.5];%red-contra, ipsi-blue, n.s.-grey
colormap(map);
box off;
%set(gca,'XTick',1:4,'XTickLabel',{'sound','delay','response','lick'});
set(gca,'XTick',1:4,'XTickLabel',{'S','D','R','L'});
title('SC neurons');
xlim([0,5]);
%set(gca,'FontSize',12);
SC_choice_AUC_005=table(pAUCSigSC(1,:)',pAUCSigSC(2,:)',n_cell*ones(size(pAUCSigSC,2),1),'VariableNames',{'contra_sig','ipsi_sig','total_n'},'RowNames',{'sound';'delay';'response';'lick'});
save('H:\2P\summary\SC_choice_AUC_005.mat','SC_choice_AUC_005');

subplot(1,2,1);%for SC projecting M2 neurons
load('H:\2P\summary\M2_choice_AUC_005.mat');
pAUCSigM2=zeros(3,4);
pAUCSigM2(1,:)=M2_choice_AUC_005.contra_sig./M2_choice_AUC_005.total_n;
pAUCSigM2(2,:)=M2_choice_AUC_005.ipsi_sig./M2_choice_AUC_005.total_n;
pAUCSigM2(3,:)=1-sum(pAUCSigM2);
bar(pAUCSigM2','stacked');
box off;
%set(gca,'XTick',1:4,'XTickLabel',{'sound','delay','response','lick'});
set(gca,'XTick',1:4,'XTickLabel',{'S','D','R','L'});
title('SC-projecting M2 neurons');
%set(gca,'FontSize',12);
xlim([0,5]);
h=legend('contra','ipsi','n.s.');
set(h,'box','off');
saveas(figBar,[savepath,filesep,celltype,'-',trialTypeStr,'-',AUCtype,'-bar plot of AUC at significance level of ',num2str(2*psig(1)),'.pdf'],'pdf');