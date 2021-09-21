function [Tcorrelation] = fGetCorrelationRsc(filepath,animal,date,field,celltype,trialTypeStr,selectivityType,varargin)
%fGetCorrelationSRC using data from one session to calculate the spike
%count correlation (rsc or noise correlation) for contralateral correct
%trials among neurons.
%Input-
%	filepath, animal, date of the chosen session
%   trialTypeStr={'cor and err','cor'} decide {'combineCorErr','divideCorErr'}
%   selectivityType-choose from selectivitystr to decide the label of
%       trial type
%   AUCtype- load corresponding AUC file to correlated with RSC


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
str_nFrames='500ms';%'500ms';%'1s'
%load AUC data
[TAUC] = fGetEpochAUCtableASession(filepath,animal,date,field,celltype,'cor and err',selectivityType,'SensoryChoiceOrthogonalSubtraction');

fileNameT=[filepath,filesep,animal,'-',datestr,'trialType',trialTypeStr,'-',selectivityType,'-timeBin',str_nFrames,'-Trsc.mat'];
if exist(fileNameT,'file')
    load(fileNameT);
    npair=size(TAUC,1);
    disp(['Table exist, use ',animal,date,';nPair=',num2str(npair)]);   
else
    nROI=size(SavedCaTrials.f_raw{1},1);
    npair=nchoosek(nROI,2);
    [varanimal{1:npair}]=deal(animal);
    varanimal=reshape(varanimal,[],1);
    [vardate{1:npair}]=deal(date);
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
    
    %calculate for only contralateral choice correct side
    if strcmp(selectivityType,'choice')|| strcmp(selectivityType,'sensory')
        label_choice = fTrialType2Label(trialType,2);
        [ITI, sound,delay, response, lick, pITI, psound, pdelay, presponse, plick]= deal(nan(npair,size(trialType,2)));
        [corrAUCITI,corrAUCsound, corrAUCdelay,corrAUCresponse,corrAUClick]= deal(nan(npair,1));
        varROI=[];
        varROI2=[];
        if contains(trialTypeStr,'cor') % 'cor'| 'cor and err'
            ind_trial=logical(reshape(sum(trialType(1,:,:),2),[],1));%only correct trials
        end
        indiIpsi=logical((label_choice==1).*ind_trial);
        indiContra=logical((label_choice==2).*ind_trial);
        n_col=5;
        n_row=7;
        ind_pair=1;
        %     ind_fig=0;
        for roiNo = 1:size(SavedCaTrials.f_raw{1},1) %SavedCaTrials.nROIs may be not true
            [T_SigbyEpoch,str_nFrames] = fGetSigBehEpoch(behEventFrameIndex,dff(roiNo,:),frT,str_nFrames);
            for roiNo2=roiNo+1:size(SavedCaTrials.f_raw{1},1)
                %             if mod(ind_pair,n_row)==1
                %                 ind_fig=ind_fig+1;
                %                 figCorrCase(ind_fig)=figure;
                %                 set(gcf,'PaperPosition',[0,0,8,11]);
                %             end
                disp([animal,date,'rioNo',num2str(roiNo),'-',num2str(roiNo2)]);
                [T_SigbyEpoch2,str_nFrames2] = fGetSigBehEpoch(behEventFrameIndex,dff(roiNo2,:),frT,str_nFrames);
                %             ax=subplot(n_row,n_col,1+mod(ind_pair-1,n_row)*n_col);
                ax=0;
                [ITI(ind_pair,:),pITI(ind_pair,:)]=fCorrcoefIpsiContra(T_SigbyEpoch.ITI,T_SigbyEpoch2.ITI,indiIpsi,indiContra,ax);
                %             ylabel(ax,'dff');
                %             ax=subplot(n_row,n_col,2+mod(ind_pair-1,n_row)*n_col);
                [sound(ind_pair,:),psound(ind_pair,:)]=fCorrcoefIpsiContra(T_SigbyEpoch.sound,T_SigbyEpoch2.sound,indiIpsi,indiContra,ax);
                %             ax=subplot(n_row,n_col,3+mod(ind_pair-1,n_row)*n_col);
                [delay(ind_pair,:),pdelay(ind_pair,:)]=fCorrcoefIpsiContra(T_SigbyEpoch.delay,T_SigbyEpoch2.delay,indiIpsi,indiContra,ax);
                %             ax=subplot(n_row,n_col,4+mod(ind_pair-1,n_row)*n_col);
                [response(ind_pair,:),presponse(ind_pair,:)]=fCorrcoefIpsiContra(T_SigbyEpoch.response,T_SigbyEpoch2.response,indiIpsi,indiContra,ax);
                %             ax=subplot(n_row,n_col,5+mod(ind_pair-1,n_row)*n_col);
                [lick(ind_pair,:),plick(ind_pair,:)]=fCorrcoefIpsiContra(T_SigbyEpoch.lick,T_SigbyEpoch2.lick,indiIpsi,indiContra,ax);
                %             if mod(ind_pair,n_row)==0
                %                 xlabel(ax,'dff');
                %             end
                varROI=[varROI;roiNo];
                varROI2=[varROI2;roiNo2];
                %calculate AUC correlation
                corrAUCITI(ind_pair)=fcorrAUC(TAUC.ITI(roiNo),TAUC.ITI(roiNo2));
                corrAUCsound(ind_pair)=fcorrAUC(TAUC.sound(roiNo),TAUC.sound(roiNo2));
                corrAUCdelay(ind_pair)=fcorrAUC(TAUC.delay(roiNo),TAUC.delay(roiNo2));
                corrAUCresponse(ind_pair)=fcorrAUC(TAUC.response(roiNo),TAUC.response(roiNo2));
                corrAUClick(ind_pair)=fcorrAUC(TAUC.lick(roiNo),TAUC.lick(roiNo2));
                ind_pair=ind_pair+1;
            end
        end
        %     mkdir('Rsc');
        %     for i=1:ind_fig
        %         saveas(figCorrCase(i),[filepath,filesep,'Rsc',filesep,'case-',num2str(i),'.pdf'],'pdf');
        %         saveas(figCorrCase(i),[filepath,filesep,'Rsc',filesep,'case-',num2str(i),'.png'],'png');
        %     end
        %     close all;
    end
    [varfield{1:npair}]=deal(field);
    varfield=reshape(varfield,[],1);
    [varcelltype{1:npair}]=deal(celltype);
    varcelltype=reshape(varcelltype,[],1);
    
    Tcorrelation=table(varanimal,vardate,varfield,varcelltype,varROI,varROI2,...
        ITI,sound,delay,response,lick, pITI, psound, pdelay, presponse,plick,...
        corrAUCITI,corrAUCsound, corrAUCdelay,corrAUCresponse,corrAUClick,...
        'VariableNames',  {'animal','date','field','celltype','nROI','nROI2',...
        'ITI','sound','delay','response', 'lick','pITI','psound','pdelay',...
        'presponse','plick','corrAUCITI','corrAUCsound','corrAUCdelay',...
        'corrAUCresponse','corrAUClick'});
end
%plot correlation of AUC and RSC
figAUCRSC=figure;
set(gcf,'position',[100,100,500,100]);
for i=1:5
    subplot(1,5,i);
    temp=table2array(Tcorrelation(:,6+i));
    scatter(temp(:,2),abs(table2array(Tcorrelation(:,16+i))),20,'k');
    hold on;
end
title([animal,'-',datestr,'-',Tcorrelation.celltype{1}]);
save(fileNameT,'Tcorrelation');
saveas(figAUCRSC,[filepath,filesep,'figAUCRSC.jpg'],'jpg');
set(figAUCRSC,'paperposition',[0,0,8,1.5]);
saveas(figAUCRSC,['H:\2P\summary',filesep,'corr_Rsc_AUC',filesep,animal,'-',datestr,'-',Tcorrelation.celltype{1},'.pdf'],'pdf');
end

function corrAUC= fcorrAUC(AUC1,AUC2)
corrAUC=abs(AUC1-AUC2);
end
function [rho,p]=fCorrcoefIpsiContra(data1,data2,indiIpsi,indiContra,varargin)
[struct_rho.ipsi,struct_p.ipsi]=corrcoef(data2(indiIpsi),data1(indiIpsi),'Rows','pairwise');
[struct_rho.ipsi,struct_p.ipsi]=deal(struct_rho.ipsi(2,1),struct_p.ipsi(2,1));
[struct_rho.contra,struct_p.contra]=corrcoef(data2(indiContra),data1(indiContra),'Rows','pairwise');
[struct_rho.contra,struct_p.contra]=deal(struct_rho.contra(2,1),struct_p.contra(2,1));
rho=[struct_rho.ipsi,struct_rho.contra];
p=[struct_p.ipsi,struct_p.contra];
% if ~isempty(varargin)
%     ax=varargin{1};
%     scatter(ax,data2(indiIpsi),data1(indiIpsi),20,'b');hold on;
%     scatter(ax,data2(indiContra),data1(indiContra),20,'r');
%     if struct_p.ipsi<0.05
%         text(ax,0.1,0.7,['rho=',num2str(round(struct_rho.ipsi,2)),',',plabelsymbol(struct_p.ipsi)],'Color','b','Units','normalized');
%     else
%         text(ax,0.1,0.7,['rho=',num2str(round(struct_rho.ipsi,2)),',',plabelsymbol(struct_p.ipsi)],'Color','k','Units','normalized');
%     end
%     if struct_p.contra<0.05
%         text(ax,0.1,0.9,['rho=',num2str(round(struct_rho.contra,2)),',',plabelsymbol(struct_p.contra)],'Color','r','Units','normalized');
%     else
%         text(ax,0.1,0.9,['rho=',num2str(round(struct_rho.contra,2)),',',plabelsymbol(struct_p.contra)],'Color','k','Units','normalized');
%     end
% end
end