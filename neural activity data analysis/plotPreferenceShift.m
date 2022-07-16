close all;
[num,txt,raw] =xlsread('C:\Users\PYX\Documents\DataSummary\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
T=cell2table(raw(2:end,1:15));
T.Properties.VariableNames=strrep(raw(1,1:15),' ','_');%table variable name can't have ' ',so replace them
%for soma 
ind_session=logical(strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*(~contains(T.behavior_performance,'probe')).*contains(T.ROI_type,'soma').*(~contains(T.field,'soma')));%not probe session
celltypePool={'M2','vglut2','vgat'};%{'syn','vglut2','vgat'};
info_str='soma';
epochPool={'sound','delay','response','lick'};
epochStrPool={'S','D','R','L'};
% epochPool={'delay','mid_delay','late_delay'};
% epochStrPool={'D','MD','LD'};

for i_celltype=1:length(celltypePool)
    ind_session2=logical(ind_session.*logical(strcmp(T.cell_type,celltypePool{i_celltype})+strcmp(T.cell_type,[celltypePool{i_celltype},'-flpo'])));
    Tchoose=T(ind_session2,:);
    % file_path='H:\2P\pyx349_20210418\im_data_reg\result_save';
    % session='pyx349_20210418';
    ind_ROI=[];
    %activity_typePool={'dff'};%{'dff','spkr'};
    frT=0;
    [cellActivityByTrialType_combine2show,meanActivityByTrialType_combine2show]=deal([]);
    savename_str_epoch=cell2mat(epochStrPool);
%     trialTypeStr='cor and err';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
%     AUCtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
%     AUCCorrectedMethod='SensoryChoiceOrthogonalSubtraction';%'None';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
    trialTypeStr='cor';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
    AUCtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
    AUCCorrectedMethod='None';%'None';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
    clear meanActivityByTrialType_combine;
    savenamestr=['E:\2P\summary\preferenceShift',filesep,celltypePool{i_celltype},'-',info_str,'-epochs_',savename_str_epoch,'-',trialTypeStr,'-',AUCtype,'.mat'];
    if false%exist(savenamestr,'file')
        load(savenamestr);
    else
        clear TAUC_combine;
        for i=1:size(Tchoose,1)
            file_path=Tchoose.file_path{i};
            session=[Tchoose.session{i},'_',Tchoose.field{i}];
            savepath=[file_path,filesep,session];
            trial2include='all';
            trial2exclude=[];
            objsession=Session2P(session,file_path,trial2include,trial2exclude);

            [TAUC, Tmean,objsession] = objsession.mGetEpochAUC(Tchoose.cell_type,Tchoose.ROI_type{i},trialTypeStr,AUCtype,AUCCorrectedMethod);
            %combine table together
            if exist('TAUC_combine','var')
                TAUC_combine=vertcat(TAUC_combine,TAUC);
            else
                TAUC_combine=TAUC;
            end
        end
        nEpochs=length(epochPool);
        [varTypes{1:nEpochs}] = deal('categorical');
        T_preference=table('Size',[size(TAUC,1),length(epochPool)],'VariableTypes',varTypes,'VariableNames',epochPool);
        for i=1:nEpochs
            switch epochPool{i}
                case 'sound'
                    indIpsi= (TAUC.sound<0.5);
                    indContra=(TAUC.sound>=0.5);
                    indSig=(1-abs(TAUC.psound-0.5)*2<0.05);
                case 'delay'
                    indIpsi= (TAUC.delay<0.5);
                    indContra=(TAUC.delay>=0.5);
                    indSig=(1-abs(TAUC.pdelay-0.5)*2<0.05);
                case 'response' 
                    indIpsi= (TAUC.response<0.5);
                    indContra=(TAUC.response>=0.5);
                    indSig=(1-abs(TAUC.presponse-0.5)*2<0.05);
                case 'lick'
                    indIpsi= (TAUC.lick<0.5);
                    indContra=(TAUC.lick>=0.5);
                    indSig=(1-abs(TAUC.plick-0.5)*2<0.05);
                case 'mid_delay'
                    indIpsi= (TAUC.mid_delay<0.5);
                    indContra=(TAUC.mid_delay>=0.5);
                    indSig=(1-abs(TAUC.pmid_delay-0.5)*2<0.05);
                case 'late_delay'
                    indIpsi= (TAUC.late_delay<0.5);
                    indContra=(TAUC.late_delay>=0.5);
                    indSig=(1-abs(TAUC.plate_delay-0.5)*2<0.05);
            end
            indSigIpsi= logical(indIpsi.*indSig);
            indSigContra= logical(indContra.*indSig);
            indNS= ~indSig;
            indNSIpsi=logical(indIpsi.*indNS);
            indNSContra = logical(indContra.*indNS);
            ind_col=zeros(1,size(T_preference,2));
            ind_col(i)=1;
            ind_col=logical(ind_col);
            T_preference=fFillTableVal(T_preference,indSigIpsi, ind_col,strcat(epochStrPool{i},'-SI'));
            T_preference=fFillTableVal(T_preference,indSigContra, ind_col,strcat(epochStrPool{i},'-SC'));
            T_preference=fFillTableVal(T_preference,indNSIpsi, ind_col,strcat(epochStrPool{i},'-NI'));
            T_preference=fFillTableVal(T_preference,indNSContra, ind_col,strcat(epochStrPool{i},'-NC'));
        end
        varNames={'left','n','right'};
        [varL,varN,varR]=deal([]);
        for i=1:nEpochs-1
            [C,ia,ic]=unique(T_preference(:,i:i+1));
            a_counts = accumarray(ic,1);
            varL=[varL;table2array(C(:,1))];
            varR=[varR;table2array(C(:,2))];
            varN=[varN;a_counts];
        end
        Tsankey=table(varL,varN,varR,'VariableNames',varNames);   
        Lst_sankey=table2cell(Tsankey);
        save(savenamestr,'TAUC_combine','T_preference','Tsankey','Lst_sankey');
        
    end
    cat_item=unique([Tsankey.left;Tsankey.right]);
    colorList_sankey=zeros(length(cat_item),3);%each category one color
    for i=1:length(cat_item)
        temp=strsplit(string(cat_item(i)),'-');
        switch temp{2}
            case 'SI'
                colorList_sankey(i,:)=[0,0,1];
            case 'SC'
                colorList_sankey(i,:)=[1,0,0];
            case 'NI'
                colorList_sankey(i,:)=[0.7,0.7,1];
            case 'NC'
                colorList_sankey(i,:)=[1,0.7,0.7];
        end
    end
    figSankey=figure();
    sankeyHdl=fsankey(figSankey,'XLim',[0,2],'YLim',[0,1],'PieceWidth',0.15,'Table',Tsankey,'Color',colorList_sankey,'layer_str',epochStrPool);   
    title(celltypePool{i_celltype});
    saveas(figSankey,['E:\2P\summary\preferenceShift',filesep,celltypePool{i_celltype},'-',info_str,'-epochs_',savename_str_epoch,'_sankey.pdf'],'pdf');
    saveas(figSankey,['E:\2P\summary\preferenceShift',filesep,celltypePool{i_celltype},'-',info_str,'-epochs_',savename_str_epoch,'_sankey.png'],'png');
    
end

function [Tout] =fFillTableVal(Tin,ind_row, ind_col,val)
if sum(ind_row) ==0 || sum(ind_col)==0
    Tout=Tin;
else
    Tout=Tin;
    tempCell=cell(sum(ind_row),sum(ind_col));
    [tempCell{:}]=deal(val);
    Tout(ind_row,ind_col)=tempCell;
end
end