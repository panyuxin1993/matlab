function [figBehavior] = fOptoPerformance(rootpath,date,animal)
%FOPTOPERFORMANCE Summary of this function goes here
%   Detailed explanation goes here
tempdate=strrep(date,'_','');
dirmat=strcat(rootpath,filesep,animal,filesep,animal,'_',tempdate,'*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
flag=cellfun(@(x) contains(x,'opto-2P.mat'), filenames);
figBehavior=figure;
if sum(flag)==0
    figure(figBehavior);
    text(0,1,'no opto behavior data');
    warning('no opto session, not plot opto performance');
    return;
end
file_beh=cellfun(@(x) contains(x,'Virables'), filenames);
if sum(file_beh)==0 %if still no Data_extract variable
    filename_behdata =fDataExtract(1,[rootpath,filesep,animal],[animal,'_',tempdate,'*opto-2P.mat']);
else
    filename_behdata = filenames{file_beh};
end
path_behdata=strcat(rootpath,filesep,animal,filesep,filename_behdata);
load(path_behdata);
f=1;%all click rates
[trialType,rule] = fGetTrialType( Data_extract,[],f,'matrix','left','divideCorErr','divideOpto');%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials
%delete after using the file
delete(path_behdata);

nTrial=zeros(size(trialType,1),size(trialType,2),size(trialType,4));% 3d-opto/non-opto,2d-stimuli(usually 2 or 4)
pRightChoice=zeros(size(trialType,2),size(trialType,4));% 2d-opto/non-opto,1d-stimuli(usually 2 or 4)
for nResult=1:size(trialType,1) %4 column(correct/error/miss/violation),companied with 4 lick raster
    for nStim=1:size(trialType,2) %for each stimulus
        for nOpto=1:size(trialType,4)
            selectedTrialInd=trialType(nResult,nStim,:,nOpto);
            selectedTrialInd=logical(squeeze(selectedTrialInd))';
            temp=find(selectedTrialInd);
            nTrial(nResult,nStim,nOpto)=length(temp);
        end
    end
end
for nOpto=1:size(trialType,4)
    for nStim=1:size(trialType,2)/2 %for each stimulus, first half
        switch rule
            case 'low click rate-left'
                pRightChoice(nStim,nOpto)=nTrial(2,nStim,nOpto)/(nTrial(1,nStim,nOpto)+nTrial(2,nStim,nOpto));
            case 'low click rate-right'
                pRightChoice(nStim,nOpto)=nTrial(1,nStim,nOpto)/(nTrial(1,nStim,nOpto)+nTrial(2,nStim,nOpto));
        end
    end
    for nStim=size(trialType,2)/2+1:size(trialType,2) %for each stimulus, first half
        switch rule
            case 'low click rate-left'
                pRightChoice(nStim,nOpto)=nTrial(1,nStim,nOpto)/(nTrial(1,nStim,nOpto)+nTrial(2,nStim,nOpto));
            case 'low click rate-right'
                pRightChoice(nStim,nOpto)=nTrial(2,nStim,nOpto)/(nTrial(1,nStim,nOpto)+nTrial(2,nStim,nOpto));
        end
    end
end

figure(figBehavior);
curve(1)=scatter(1:size(pRightChoice,2),pRightChoice(:,1),'k');hold on;%control
curve(2)=scatter(1:size(pRightChoice,2),pRightChoice(:,2),'r');%opto
legend(curve(:),'ctrl','opto');
set(gca,'Xlim',[0,size(pRightChoice,2)+1],'Ylim',[0,1]);
ylabel('P(Right Choice)');
xlabel('stimuli');
set(gca,'FontSize',14);

end

