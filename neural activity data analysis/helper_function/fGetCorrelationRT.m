function [Tcorrelation] = fGetCorrelationRT(filepath,animal,date,field,celltype,trialTypeStr,selectivityType,correlationMethod,varargin)
%FGETCORRELATIONRT using data from one session to calculate the correlation
%between RT and some parameter by neurons(activities) or by sessions(SVM
%predictions)
%Input-
%	filepath, animal, date of the chosen session
%   trialTypeStr={'cor and err','cor'} decide {'combineCorErr','divideCorErr'}
%   selectivityType-choose from selectivitystr to decide the label of
%       trial type
%   correlationMethod={'cells','SVM'}, etc. decide which to be correlated
%       correlated with RT, eithr by cell or by population(like SVM)

cd(filepath);
dirmat=strcat(filepath,filesep,'*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
file_imaging=cellfun(@(x) contains(x,'Ca'), filenames);
load(filenames{file_imaging});%load imaging data
load('dff.mat');%load dff
file_beh=cellfun(@(x) contains(x,'Virables'), filenames);
load(filenames{file_beh});

datestr=strrep(date,'/','-');
if strcmp(trialTypeStr,'cor')
    combineCorErr='divideCorErr';%and only use correct trial
elseif strcmp(trialTypeStr,'cor and err')
    combineCorErr='combineCorErr';%{'combineCorErr','divideCorErr'}
end
selectivitystr={'stimuli','sensory difficulty','sensory','choice'};%sensory means grouping difficulties;
trialTypeVar=[1,2,4,3];%corresponding to variable 'selectivitystr',decide what trial type means
i_selectivity=cellfun(@(x) strcmp(x,selectivityType),selectivitystr);%*********variable**************
i_selectivity=find(i_selectivity);
str_nFrames='1s';%'500ms';%'1s'
if strcmp(correlationMethod,'SVM')
    nRepeat=100;%default
    pTraining=0.9;%default
    if ~isempty(varargin)
        nRepeat=varargin{1}(1);
        pTraining=varargin{1}(2);
    end
elseif strcmp(correlationMethod,'cells')%test correlation between RT and individual cell activities
    nROI=size(SavedCaTrials.f_raw{1},1);
    [delayMovingAUC,pdelayMovingAUC]=deal(cell(nROI,1));
    [varanimal{1:nROI}]=deal(animal);
    varanimal=reshape(varanimal,[],1);
    [vardate{1:nROI}]=deal(date);
    vardate=reshape(vardate,[],1);
    ind_tr_1=1;
    ntr=length(SavedCaTrials.f_raw);
    frT = SavedCaTrials.FrameTime;
    % align to behavior event
    nFrameEachTrial=cellfun(@(x) size(x,2),SavedCaTrials.f_raw);
    ind_1stFrame=zeros(1,length(nFrameEachTrial));
    ind_1stFrame(1)=1;
    ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
    ind_1stFrame=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials
    frameNumTime=[1,1.5];%from 5s before align point to 5s after align point
    frameNum=double(round(frameNumTime*1000/frT));
    [behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( Data_extract, ind_1stFrame, SavedCaTrials.FrameTime );%get behavior event time
    [trialType,~,~] = fGetTrialType( Data_extract,[],trialTypeVar(i_selectivity),'matrix','left',combineCorErr);
    RT=double(Data_extract.Answer_time-Data_extract.Go_time);
    RT(RT<0)=nan;
    %calculate for both side, since not know which is the preferred side
    if strcmp(selectivityType,'choice')|| strcmp(selectivityType,'sensory')
        label_choice = fTrialType2Label(trialType,2);
        [ITI, sound,delay, response, lick, pITI, psound, pdelay, presponse, plick]= deal(nan(nROI,size(trialType,2)));
        if contains(trialTypeStr,'cor') % 'cor'| 'cor and err'
            ind_trial=logical(reshape(sum(trialType(1,:,:),2),[],1));%only correct trials
        end
        n_col=5;
        n_row=7;
        ind_fig=0;
        for roiNo = 1:size(SavedCaTrials.f_raw{1},1) %SavedCaTrials.nROIs may be not true
            if mod(roiNo,n_row)==1
                ind_fig=ind_fig+1;
                figCorrCase(ind_fig)=figure;
                set(gcf,'PaperPosition',[0,0,8,11]);
            end
            disp([animal,date,'rioNo',num2str(roiNo)]);
            [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,dff(roiNo,:),frT,str_nFrames);
            indiIpsi=logical((label_choice==1).*ind_trial);
            indiContra=logical((label_choice==2).*ind_trial);
            ax=subplot(n_row,n_col,1+mod(roiNo-1,n_row)*n_col);
            [ITI(roiNo,:),pITI(roiNo,:)]=fCorrcoefIpsiContra(RT',T_SigbyEpoch.ITI,indiIpsi,indiContra,ax);
            ylabel(ax,'RT(ms)');
            ax=subplot(n_row,n_col,2+mod(roiNo-1,n_row)*n_col);
            [sound(roiNo,:),psound(roiNo,:)]=fCorrcoefIpsiContra(RT',T_SigbyEpoch.sound,indiIpsi,indiContra,ax);
            ax=subplot(n_row,n_col,3+mod(roiNo-1,n_row)*n_col);
            [delay(roiNo,:),pdelay(roiNo,:)]=fCorrcoefIpsiContra(RT',T_SigbyEpoch.delay,indiIpsi,indiContra,ax);
            ax=subplot(n_row,n_col,4+mod(roiNo-1,n_row)*n_col);
            [response(roiNo,:),presponse(roiNo,:)]=fCorrcoefIpsiContra(RT',T_SigbyEpoch.response,indiIpsi,indiContra,ax);
            ax=subplot(n_row,n_col,5+mod(roiNo-1,n_row)*n_col);
            [lick(roiNo,:),plick(roiNo,:)]=fCorrcoefIpsiContra(RT',T_SigbyEpoch.lick,indiIpsi,indiContra,ax);
            if mod(roiNo,n_row)==0
                xlabel(ax,'dff');
            end
        end
%         mkdir('correlation_RT_activity');
%         for i=1:ind_fig
%             saveas(figCorrCase(i),[filepath,filesep,'correlation_RT_activity',filesep,'case-',num2str(i),'.pdf'],'pdf');
%             saveas(figCorrCase(i),[filepath,filesep,'correlation_RT_activity',filesep,'case-',num2str(i),'.png'],'png');
%         end
        close all;
    end
    [varfield{1:nROI}]=deal(field);
    varfield=reshape(varfield,[],1);
    [varcelltype{1:nROI}]=deal(celltype);
    varcelltype=reshape(varcelltype,[],1);
    varROI=(1:size(SavedCaTrials.f_raw{1},1));
    varROI=reshape(varROI,[],1);
    Tcorrelation=table(varanimal,vardate,varfield,varcelltype,varROI,ITI, sound,delay, response, ...
        lick, pITI, psound, pdelay, presponse, plick,'VariableNames',...
        {'animal','date','field','celltype','nROI','ITI','sound','delay','response',...
        'lick','pITI','psound','pdelay','presponse','plick'});
end

end

function [rho,p]=fCorrcoefIpsiContra(RT,activity,indiIpsi,indiContra,ax)
[struct_rho.ipsi,struct_p.ipsi]=corrcoef(activity(indiIpsi),RT(indiIpsi),'Rows','pairwise');
[struct_rho.ipsi,struct_p.ipsi]=deal(struct_rho.ipsi(2,1),struct_p.ipsi(2,1));
[struct_rho.contra,struct_p.contra]=corrcoef(activity(indiContra),RT(indiContra),'Rows','pairwise');
[struct_rho.contra,struct_p.contra]=deal(struct_rho.contra(2,1),struct_p.contra(2,1));
rho=[struct_rho.ipsi,struct_rho.contra];
p=[struct_p.ipsi,struct_p.contra];
% scatter(ax,activity(indiIpsi),RT(indiIpsi),20,'b');hold on;
% scatter(ax,activity(indiContra),RT(indiContra),20,'r');
% set(gca,'Ylim',[0,1500]);
% if struct_p.ipsi<0.05
%     text(ax,0.1,0.7,['rho=',num2str(round(struct_rho.ipsi,2)),',',plabelsymbol(struct_p.ipsi)],'Color','b','Units','normalized');
% else
%     text(ax,0.1,0.7,['rho=',num2str(round(struct_rho.ipsi,2)),',',plabelsymbol(struct_p.ipsi)],'Color','k','Units','normalized');
% end
% if struct_p.contra<0.05
%     text(ax,0.1,0.9,['rho=',num2str(round(struct_rho.contra,2)),',',plabelsymbol(struct_p.contra)],'Color','r','Units','normalized');
% else
%     text(ax,0.1,0.9,['rho=',num2str(round(struct_rho.contra,2)),',',plabelsymbol(struct_p.contra)],'Color','k','Units','normalized');
% end

end