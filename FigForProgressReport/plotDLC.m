%plot coordinates and likelihood in using deeplabcut
%coordinates of which body part
%for coordinates
%bodyparts={'Tongue','Tongue','LeftHandFingerTip','LeftHandFingerTip','RightHandFingerTip','RightHandFingerTip','Nose','Nose','LeftWhiskerTip','LeftWhiskerTip','RightWhiskerTip','RightWhiskerTip','LeftLickPort','LeftLickPort','RightLickPort','RightLickPort'};%{Tongue,LeftHandFingerTip,LeftHandFingerRoot,LeftLickPort,Nose,LeftWhiskerTip,LeftWhiskerRoot, etc.}
%coordinates={'x','y','x','y','x','y','x','y','x','y','x','y','x','y','x','y'};%{x,y,likelihood};
% %for likelihood
% bodyparts={'Tongue','LeftHandFingerTip','RightHandFingerTip','LeftHandFingerRoot','RightHandFingerRoot','Nose','LeftWhiskerTip','RightWhiskerTip','LeftLickPort','RightLickPort'};
% coordinates={'likelihood','likelihood','likelihood','likelihood','likelihood','likelihood','likelihood','likelihood','likelihood','likelihood'};

dbstop if error;

close all;
clear;
%% plot SC FP video
[num,txt,raw] =xlsread('D:\xulab\project\fiber_photometry_data_summary.xlsx');%criteria to choose sessions come from this file
savepath='H:\FP';
summaryFile=[savepath,filesep,'imaging_video_data_summary.xlsx'];
trialTypeStr='cor';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
AUCtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
AUCCorrectedMethod='balencedCorErrTrialNum';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
manipulation='control';
celltype_Str='SC vgat';
T=cell2table(raw(2:end,1:14));
T.Properties.VariableNames=strrep(raw(1,1:14),' ','_');%table variable name can't have ' ',so replace them
ind_session=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*strcmp(T.behavior_video,'DLC tracked').*strcmp(T.experiment,celltype_Str);
ind_session=find(ind_session);
n_session=length(ind_session);
animal_unique=unique(T.animal(ind_session));
n_animal=length(animal_unique);
[var_path,var_animal,var_session,var_video,var_DLC,var_OLED,var_iter]=deal(cell(n_session,1));
for i_session=1:n_session
    indrow=ind_session(i_session);
    var_path{i_session}=T.file_path{indrow};
    cd(T.file_path{indrow});
    CurrFolder=pwd;
    dirmat=strcat(T.file_path{indrow},filesep,'*.mat');
    dirs=dir(dirmat);
    dircell=struct2cell(dirs);
    filenames=dircell(1,:);
    file_beh_ind=cellfun(@(x) contains(x,'_Virables.mat'), filenames);
    file_beh_str=filenames{file_beh_ind};
    file_beh=strsplit(file_beh_str,'_Virables');
    var_animal{i_session}=T.animal{i_session};
    var_session{i_session}=file_beh{1};
    diravi=strcat(T.file_path{indrow},filesep,'video',filesep,'*.avi');
    dirs=dir(diravi);
    dircell=struct2cell(dirs);
    filenames=dircell(1,:);
    file_avi_ind=cellfun(@(x) ~contains(x,'crop'), filenames);
    var_video{i_session}=strrep(filenames{file_avi_ind},'.avi','');
    dircsv=strcat(T.file_path{indrow},filesep,'video',filesep,'*.csv');
    dirs=dir(dircsv);
    dircell=struct2cell(dirs);
    filenames=dircell(1,:);
    file_oled_ind=cellfun(@(x) contains(x,'OLED'), filenames);
    file_DLC_ind=cellfun(@(x) contains(x,'DLC'), filenames);
    var_DLC{i_session}=strrep(filenames{file_DLC_ind},'.csv','');
    var_OLED{i_session}=strrep(filenames{file_oled_ind},'.csv','');
    var_iter{i_session}='iteration-1';
end
Tout1=table(var_path,var_animal,var_session,var_video,var_DLC,var_OLED);
Tout1.Properties.VariableNames={'rootpath','animal','session','videoName','DLCFileName1','OLEDFileName'};
Tout2=table(var_session,var_iter,var_iter,var_iter,var_iter,var_iter,var_iter);
Tout2.Properties.VariableNames={'session','Tongue','LeftHand','RightHand','Nose','LeftLickPort','RightLickPort'};
writetable(Tout1,summaryFile,'Sheet',1);
writetable(Tout2,summaryFile,'Sheet',2);

filenameID=['check delay movement of SC FP results-',celltype_Str];
%% plot SC imaging video
%summarized the data in an excel file
%{
[num,txt,raw] =xlsread('D:\xulab\project\imaging_data_summary.xlsx');%criteria to choose sessions come from this file
savepath='H:\2P';
summaryFile=[savepath,filesep,'imaging_video_data_summary.xlsx'];
trialTypeStr='cor';%{'cor and err','do','cor'}; should be 'cor' if AUCtype is stimuli, otherwise may merge cor and err together
AUCtype='choice';%{'choice','sensory'};'stimuli' means comparing cor and err for each stimuli
AUCCorrectedMethod='balencedCorErrTrialNum';%'balencedCorErrTrialNum';%'SensoryChoiceOrthogonalSubtraction';
manipulation='control';
celltype_Str='vgat';
T=cell2table(raw(2:end,1:17));
T.Properties.VariableNames=strrep(raw(1,1:17),' ','_');%table variable name can't have ' ',so replace them
ind_session=strcmp(T.used_as_data,'yes').*strcmp(T.manipulation,'control').*strcmp(T.behavior_video,'DLC tracked').*strcmp(T.cell_type,celltype_Str).*(~strcmp(T.animal,'pyx290'));
ind_session=find(ind_session);
n_session=length(ind_session);
animal_unique=unique(T.animal(ind_session));
n_animal=length(animal_unique);
[var_path,var_animal,var_session,var_video,var_DLC,var_OLED,var_iter]=deal(cell(n_session,1));
for i_session=1:n_session
    indrow=ind_session(i_session);
    var_path{i_session}=T.file_path{indrow};
    cd(T.file_path{indrow});
    CurrFolder=pwd;
    dirmat=strcat(T.file_path{indrow},filesep,'*.mat');
    dirs=dir(dirmat);
    dircell=struct2cell(dirs);
    filenames=dircell(1,:);
    file_beh_ind=cellfun(@(x) contains(x,'_Virables.mat'), filenames);
    file_beh_str=filenames{file_beh_ind};
    file_beh=strsplit(file_beh_str,'_Virables');
    var_animal{i_session}=T.animal{i_session};
    var_session{i_session}=file_beh{1};
    diravi=strcat(T.file_path{indrow},filesep,'video',filesep,'*.avi');
    dirs=dir(diravi);
    dircell=struct2cell(dirs);
    filenames=dircell(1,:);
    file_avi_ind=cellfun(@(x) ~contains(x,'crop'), filenames);
    var_video{i_session}=strrep(filenames{file_avi_ind},'.avi','');
    dircsv=strcat(T.file_path{indrow},filesep,'video',filesep,'*.csv');
    dirs=dir(dircsv);
    dircell=struct2cell(dirs);
    filenames=dircell(1,:);
    file_oled_ind=cellfun(@(x) contains(x,'OLED'), filenames);
    file_DLC_ind=cellfun(@(x) contains(x,'DLC'), filenames);
    var_DLC{i_session}=strrep(filenames{file_DLC_ind},'.csv','');
    var_OLED{i_session}=strrep(filenames{file_oled_ind},'.csv','');
    var_iter{i_session}='iteration-1';
end
Tout1=table(var_path,var_animal,var_session,var_video,var_DLC,var_OLED);
Tout1.Properties.VariableNames={'rootpath','animal','session','videoName','DLCFileName1','OLEDFileName'};
Tout2=table(var_session,var_iter,var_iter,var_iter,var_iter,var_iter,var_iter);
Tout2.Properties.VariableNames={'session','Tongue','LeftHand','RightHand','Nose','LeftLickPort','RightLickPort'};
writetable(Tout1,summaryFile,'Sheet',1);
writetable(Tout2,summaryFile,'Sheet',2);

filenameID=['check delay movement of SC imaging results-',celltype_Str];
%}
%% for data strored distributedly
%
bodyparts={'Tongue','Tongue','LeftHand','LeftHand';
    'Tongue','Tongue','RightHand','RightHand';
    'LeftLickPort','LeftLickPort','Nose','Nose'};%{Tongue,LeftHandFingerTip,LeftHandFingerRoot,LeftLickPort,Nose,LeftWhiskerTip,LeftWhiskerRoot, etc.}
coordinates={'x','x','x','x';
    'y','y','x','x';    
    'x','x','x','x'};
behEventAlign={'delay onset','first lick','delay onset','first lick';
    'delay onset','first lick','delay onset','first lick';
    'delay onset','first lick','delay onset','first lick'};%align to which event(string can be in {'stim onset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
titlestr={'Tongue (x)','Tongue (x)','Left forepaw (x)','Left forepaw (x)';
    'Tongue (y)','Tongue (y)','Right forepaw (x)','Right forepaw (x)';
    'Left lick port (x)','Left lick port (x)','Nose (x)','Nose (x)'};
yrange={[-40,40],[-40,40],[-100,100],[-100,100];
    [-40,40],[-40,40],[-40,40],[-40,40];
    [-10,10],[-10,10],[-10,10],[-10,10]};
nonOverlappingBin=1;

pSigTtest=0.01;
treshold4likelihood=0.1;
fCorrelationBehaviorDistributedSession( savepath,bodyparts',coordinates',pSigTtest,treshold4likelihood,behEventAlign',filenameID ,yrange',titlestr',nonOverlappingBin);
%}
%% plot pyx290 as example of licking during delay
%{
filepath='H:\video tracking\SC_imaging_video';
bodyparts={'Tongue','Tongue','LeftHand','LeftHand';
    'Tongue','Tongue','RightHand','RightHand';
    'LeftLickPort','LeftLickPort','Nose','Nose'};%{Tongue,LeftHandFingerTip,LeftHandFingerRoot,LeftLickPort,Nose,LeftWhiskerTip,LeftWhiskerRoot, etc.}
coordinates={'x','x','x','x';
    'y','y','x','x';    
    'x','x','x','x'};
behEventAlign={'delay onset','first lick','delay onset','first lick';
    'delay onset','first lick','delay onset','first lick';
    'delay onset','first lick','delay onset','first lick'};%align to which event(string can be in {'stim onset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
filenameID='see SC imaging results pyx290 whether have movement';
titlestr={'Tongue (x)','Tongue (x)','Left forepaw (x)','Left forepaw (x)';
    'Tongue (y)','Tongue (y)','Right forepaw (x)','Right forepaw (x)';
    'Left lick port (x)','Left lick port (x)','Nose (x)','Nose (x)'};
yrange={[-40,40],[-40,40],[-100,100],[-100,100];
    [-40,40],[-40,40],[-40,40],[-40,40];
    [-10,10],[-10,10],[-10,10],[-10,10]};
nonOverlappingBin=1;
%}

%% plot M2 imaging data
%{
% filepath='H:\video tracking\M2 imaging video';
% bodyparts={'Tongue'};%{Tongue,LeftHandFingerTip,LeftHandFingerRoot,LeftLickPort,Nose,LeftWhiskerTip,LeftWhiskerRoot, etc.}
% coordinates={'x'};
% behEventAlign={'delay onset'};%align to which event(string can be in {'stim onset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
% filenameID='try align to delay-non overlapping bin-iteration-1+2';
% titlestr={'Tongue (x)'};
% yrange={[-15,30]};
% nonOverlappingBin=2;

% filepath='H:\video tracking\M2 imaging video';
% bodyparts={'Tongue','Tongue','LeftLickPort';
%     'Tongue','Tongue','LeftLickPort';
%     'Tongue','Tongue','LeftLickPort'};%{Tongue,LeftHandFingerTip,LeftHandFingerRoot,LeftLickPort,Nose,LeftWhiskerTip,LeftWhiskerRoot, etc.}
% coordinates={'likelihood','likelihood','likelihood';
%     'x','x','x';
%     'y','y','y'};
% behEventAlign={'stim onset','first lick','first lick';
%     'stim onset','first lick','first lick';
%     'stim onset','first lick','first lick'};%align to which event(string can be in {'stim onset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
% filenameID='method-confirmation';

%for iteration-1~3
% filepath='H:\video tracking\M2 imaging video';
% bodyparts={'Tongue','Tongue','LeftHandFingerTip','LeftHandFingerTip';
%     'Tongue','Tongue','RightHandFingerTip','RightHandFingerTip';
%     'LeftLickPort','LeftLickPort','Nose','Nose'};%{Tongue,LeftHandFingerTip,LeftHandFingerRoot,LeftLickPort,Nose,LeftWhiskerTip,LeftWhiskerRoot, etc.}
% % for iteration-4~5
% % bodyparts={'Tongue','Tongue','LeftHand','LeftHand';
% %     'Tongue','Tongue','RightHand','RightHand';
% %     'RightLickPort','RightLickPort','Nose','Nose'};%{Tongue,LeftHand,LeftHandFingerRoot,LeftLickPort,Nose,LeftWhiskerTip,LeftWhiskerRoot, etc.}
% coordinates={'x','x','x','x';
%     'y','y','x','x';    
%     'x','x','x','x'};
% behEventAlign={'delay onset','first lick','delay onset','first lick';
%     'delay onset','first lick','delay onset','first lick';
%     'delay onset','first lick','delay onset','first lick'};%align to which event(string can be in {'stim onset', 'go cue','first lick','first left lick','first right lick', 'answer','reward','start'},
% filenameID='iteration-1+2-non overlapping bin-2-low likelihood2nan-align delay';
% titlestr={'Tongue (x)','Tongue (x)','Left forepaw (x)','Left forepaw (x)';
%     'Tongue (y)','Tongue (y)','Right forepaw (x)','Right forepaw (x)';
%     'Left lick port (x)','Left lick port (x)','Nose (x)','Nose (x)'};
% yrange={[-40,40],[-40,40],[-200,200],[-200,200];
%     [-40,40],[-40,40],[-200,200],[-200,200];
%     [-20,20],[-20,20],[-20,20],[-20,20]};
% nonOverlappingBin=2;
%}
%% main code for data strored in one folder
%{
pSigTtest=0.01;
treshold4likelihood=0.1;
fCorrelationBehaviorMultiSession( filepath,bodyparts',coordinates',pSigTtest,treshold4likelihood,behEventAlign',filenameID ,yrange',titlestr',nonOverlappingBin);
%}
%% plot likelihood histogram
%{
filepath='F:\video tracking\M2 imaging video';
savepath=[filepath,filesep,'summary figures'];
summaryFile=[filepath,filesep,'imaging_video_data_summary.xlsx'];
[~,~,temp]=xlsread(summaryFile,1);
dataSummaryT=cell2table(temp(2:end,:));
dataSummaryT.Properties.VariableNames =temp(1,:);
[~,datasource,~]=xlsread(summaryFile,2);
fr=24;%frame rate of the video
bodyparts={'Tongue','LeftHandFingerTip','RightHandFingerTip','Nose','LeftLickPort'};
coordinates={'likelihood','likelihood','likelihood','likelihood','likelihood'};
colorhex={'76FD75','0081FA','7EFF79','FF9500','7D0002'};
colorhist=fHex2RGB(colorhex);
threshold4likelihood=0.1;
histLikelihood=figure;

for iBody=1:length(bodyparts)
    likelihoodall=[];
    for iSession=1:size(dataSummaryT,1)%can be a loop
        col_datasource=cellfun(@(x) strcmp(bodyparts{iBody},x),datasource(1,:));
        row_datasource=cellfun(@(x) strcmp(dataSummaryT.session{iSession},x),datasource(:,1));
        DLCiteration=datasource{row_datasource,col_datasource};
        switch DLCiteration %DLC model from which iteration
            case 'iteration-1'
                file_trace=[filepath,filesep,'iteration1',filesep,dataSummaryT.DLCFileName1{iSession},'.csv'];
            case 'iteration-2'
                file_trace=[filepath,filesep,'iteration2',filesep,dataSummaryT.DLCFileName2{iSession},'.csv'];
            case 'iteration-3'
                file_trace=[filepath,filesep,'iteration3',filesep,dataSummaryT.DLCFileName3{iSession},'.csv'];
        end
        name_file_trace=strsplit(file_trace,'.');
        if ~exist([name_file_trace{1},'.mat'],'file')
            [dcnum,dctxt,~] = xlsread(file_trace,1);%this process is time consuming, so not do very often
            save([name_file_trace{1},'.mat'],'dcnum','dctxt');
        else
            load([name_file_trace{1},'.mat'])
        end
        
        indcolbody=cellfun(@(x) strcmp(bodyparts{iBody},x),dctxt(2,:));
        indcol_likelihood=cellfun(@(x) strcmp('likelihood',x),dctxt(3,:));
        indcol_like=find(indcolbody.*indcol_likelihood);
        bodycoli=dcnum(:,indcol_like);
        likelihoodall=[likelihoodall;bodycoli];
    end
    figure(histLikelihood);
    histcurve(iBody)=histogram(likelihoodall,'BinWidth',0.05,'DisplayStyle','stairs','EdgeColor',colorhist{iBody},'LineWidth',1.5);hold on;
end
lh=legend(histcurve(:),bodyparts);
set(lh,'box','off');
xlabel('likelihood');
ylabel('Frame count');
set(gcf, 'PaperPosition', [0 0 2.5 2]);
set(gca,'FontSize',12,'FontName','Arial','box','off');
set(gca,'Xlim',[-0.1,1.1]);
saveas(histLikelihood,[savepath,filesep,'histogram-likelihood.pdf'],'pdf');
%}
function [colorhist]=fHex2RGB(colorhex)
colorhist=cell(size(colorhex));
for i=1:length(colorhex)
    rgb=zeros(1,3);
    for iRGB=1:3
        strRGB=colorhex{i}(2*iRGB-1:2*iRGB);
        rgb(iRGB)=(hex2dec(strRGB))/255;
    end
    colorhist{i}=rgb;
end
end
