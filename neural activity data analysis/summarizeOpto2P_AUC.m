close all;
clear;
[num,txt,raw] =xlsread('C:\Users\PYX\Documents\DataSummary\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
savepath='E:\2P\summary';
clear TAUC_combine Tmean_combine;
trialTypeStr='cor and err';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
AUCtype='sensory';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
AUCCorrectedMethod='SensoryChoiceOrthogonalSubtraction';%'None';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
% trialTypeStr='cor';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
% AUCtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
% AUCCorrectedMethod='None';%'None';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
manipulation='control';
%% calculate AUC
%
T=cell2table(raw(2:end,1:15));
T.Properties.VariableNames=strrep(raw(1,1:15),' ','_');%table variable name can't have ' ',so replace them
% ind_session1=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*(~contains(T.behavior_performance,'probe'));%not probe session
ind_session1=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'opto');
%  ind_session2=logical(contains(T.cell_type,'syn').*contains(T.ROI_type,'soma'));
ind_session2=logical(strcmp(T.cell_type,'vglut2')+strcmp(T.cell_type,'vgat')+strcmp(T.cell_type,'M2')+strcmp(T.cell_type,'vglut2-flpo')+strcmp(T.cell_type,'vgat-flpo'));
% ind_session2=logical(strcmp(T.cell_type,'vglut2')+strcmp(T.cell_type,'vgat')+strcmp(T.cell_type,'vglut2-flpo')+strcmp(T.cell_type,'vgat-flpo'));
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

filename_AUC=[savepath,filesep,'trialType',trialTypeStr,'-',AUCtype,'TepochAUC.mat'];
filename_MeanAcitivity=[savepath,filesep,'trialType',trialTypeStr,'-TepochMeanActivity.mat'];
%{
for i_session=1:n_session
    indrow=ind_session(i_session);
    trial2include=T.used_trial{indrow};
    trial2exclude=[];
    session=[T.session{indrow},'_',T.field{indrow}];
    file_path=T.file_path{indrow};
    savepath=[file_path,filesep,session];
    objsession=SessionOpto2P(session,file_path,trial2include,trial2exclude);
    [TAUC_ctrl, Tmean_ctrl,objsession] = objsession.mGetEpochAUC(T.cell_type{indrow},T.ROI_type{indrow},trialTypeStr,AUCtype,AUCCorrectedMethod,'ctrl');
    [TAUC_opto, Tmean_opto,objsession] = objsession.mGetEpochAUC(T.cell_type{indrow},T.ROI_type{indrow},trialTypeStr,AUCtype,AUCCorrectedMethod,'opto');
    if exist('TAUC_combine_opto','var') && strcmp(AUCtype,'stimuli') %calculate choice probability
        if size(TAUC_currentSession.ITI,2)==2
            TAUC_combine_opto=vertcat(TAUC_combine_opto,TAUC_opto);
            TAUC_combine_ctrl=vertcat(TAUC_combine_ctrl,TAUC_ctrl);
        end
    elseif exist('TAUC_combine_opto','var') && ~strcmp(AUCtype,'stimuli')
        TAUC_combine_opto=vertcat(TAUC_combine_opto,TAUC_opto);
        TAUC_combine_ctrl=vertcat(TAUC_combine_ctrl,TAUC_ctrl);
    else
        TAUC_combine_opto=TAUC_opto;
        TAUC_combine_ctrl=TAUC_ctrl;
    end
    if exist('Tmean_combine_opto','var')
        Tmean_combine_opto=vertcat(Tmean_combine_opto,Tmean_opto);
        Tmean_combine_ctrl=vertcat(Tmean_combine_ctrl,Tmean_ctrl);
    else
        Tmean_combine_opto=Tmean_opto;
        Tmean_combine_ctrl=Tmean_ctrl;
    end
end
save(filename_AUC,'TAUC_combine_opto','TAUC_combine_ctrl');
save(filename_MeanAcitivity,'Tmean_combine_opto','Tmean_combine_ctrl');
%}
%plot mean activities and selectivity
load(filename_AUC);
load(filename_MeanAcitivity);
figSummary=figure;
set(gcf,'Position',[100,100,600,600]);
color_scatter={[0.5,0.5,0.5],[0,0,1],[1,0,0]};%n.s.,ipsi, contra
legend_str={'n.s.','ipsi','contra'};
indSig1= (1-abs(TAUC_combine_ctrl.pdelay-0.5)*2)<0.05;
subplot(2,2,1);%plot mean activities change, label different selectivity (AUC) neurons
ctrl_mean=arrayfun(@(x) x.mean, Tmean_combine_ctrl.delay);
ind_ipsi=logical(indSig1.*(TAUC_combine_ctrl.delay<0.5));
ind_contra=logical(indSig1.*(TAUC_combine_ctrl.delay>=0.5));
opto_mean=arrayfun(@(x) x.mean, Tmean_combine_opto.delay);
curve=fScatterGroupPlot(ctrl_mean,opto_mean,[~indSig1,ind_ipsi,ind_contra],color_scatter);
hl=legend(curve,legend_str);
set(hl,'Box','off','Location','southeast');
title('Delay activities');
xlabel('Control');
ylabel('Opto');

subplot(2,2,2);%plot activity differences change
ctrl_contra=arrayfun(@(x) x.contra, Tmean_combine_ctrl.delay);
ctrl_ipsi=arrayfun(@(x) x.ipsi, Tmean_combine_ctrl.delay);
opto_contra=arrayfun(@(x) x.contra, Tmean_combine_opto.delay);
opto_ipsi=arrayfun(@(x) x.ipsi, Tmean_combine_opto.delay);
ctrl_diff=ctrl_contra-ctrl_ipsi;
opto_diff=opto_contra-opto_ipsi;
fScatterGroupPlot(ctrl_diff,opto_diff,[~indSig1,ind_ipsi,ind_contra],color_scatter);
title('Delay activities (contra-ipsi)');
xlabel('Control');
ylabel('Opto');

subplot(2,2,3);%plot AUC change
indSig2= (1-abs(TAUC_combine_opto.pdelay-0.5)*2)<0.05;
indSig=logical(indSig1+indSig2);
fScatterGroupPlot(TAUC_combine_ctrl.delay,TAUC_combine_opto.delay,[~indSig1,ind_ipsi,ind_contra],color_scatter,indSig);
title([AUCtype,' AUC by ',trialTypeStr,' trials']);
xlabel('Control');
ylabel('Opto');

subplot(2,2,4);%plot AUC change
fScatterGroupPlot(abs(TAUC_combine_ctrl.delay-0.5)*2,abs(TAUC_combine_opto.delay-0.5)*2,[~indSig1,ind_ipsi,ind_contra],color_scatter,indSig);
title('Absolute AUC');
xlabel('Control');
ylabel('Opto');

%% helper function
function [curve]=fScatterPlot(x,y,color,varargin)
if isempty(varargin)
    curve=scatter(x,y,10,color);hold on;
else
    indHighlight=varargin{1};
    curve=scatter(x(~indHighlight),y(~indHighlight),10,color);hold on;
    curve=scatter(x(indHighlight),y(indHighlight),10,color,'filled');
end
temp1=min(min(x),min(y));
temp2=max(max(x),max(y));
plot([temp1,temp2],[temp1,temp2]);

set(gca,'Xlim',[temp1,temp2],'Ylim',[temp1,temp2]);
if jbtest(x) && jbtest(y)%adtest(x) && adtest(y)
    [h,p]=ttest(x,y);
    text(0,0.9,['t test p=',num2str(p)],'Unit','Normalized');
else
    p=signrank(x,y);
    text(0,0.9,['signed rank test p=',num2str(p)],'Unit','Normalized');
end
set(gca,'FontSize',12,'FontName','Arial');
end

function [curve]=fScatterGroupPlot(x_in,y_in,ind,color,varargin)
if size(ind,2)~=length(color)
    error('Input size not consistent');
end
for i=1:length(color)
    if isempty(varargin)
        x=x_in(ind(:,i));
        y=y_in(ind(:,i));
        curve(i)=scatter(x,y,10,color{i});hold on;
    else
        indHighlight=varargin{1};
        x=x_in(ind(:,i));
        y=y_in(ind(:,i));
        x_norm=x_in(logical(ind(:,i).*(~indHighlight)));
        y_norm=y_in(logical(ind(:,i).*(~indHighlight)));
        x_high=x_in(logical(ind(:,i).*indHighlight));
        y_high=y_in(logical(ind(:,i).*indHighlight));
        curve(i)=scatter(x_norm,y_norm,10,color{i});hold on;
        curve_opto=scatter(x_high,y_high,10,color{i},'filled');
    end
    temp1=min(min(x_in),min(y_in));
    temp2=max(max(x_in),max(y_in));
    plot([temp1,temp2],[temp1,temp2]);
    
    set(gca,'Xlim',[temp1,temp2],'Ylim',[temp1,temp2]);
    if jbtest(x) && jbtest(y)%adtest(x) && adtest(y)
        [h,p]=ttest(x,y);
        text(0,1-0.1*i,['t test p=',num2str(p)],'Unit','Normalized');
    else
        p=signrank(x,y);
        text(0,1-0.1*i,['signed rank test p=',num2str(p)],'Unit','Normalized');
    end
end
set(gca,'FontSize',12,'FontName','Arial');
end
