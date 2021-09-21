function [score,score_shuffle] = fSVM(Tin,Y, nRepeat, pTraining,binsize,binstep)
%FSVM automaticly select balanced number of trials(rows) for 2 categories
%to fit SVM model and return the accuracy of the model
%Input-
%   Tin-table of data, m-n-t matrix, m trials,n neurons, t frames/epochs
%   Y- class labels made of 0/1
%   nRepeat- 100(default),times of trainig-testing iteration
%   pTraining- 0.8(default),proportion of trials used for training
%Output-
%   score- accuracy of training, nRepeat-by-t matrix
if size(Tin,1)~=length(Y)
    warning('input not same length');
    return;
end
vY=unique(Y);
if length(vY)~=2
    warning('input not 2 category, quit SVM');
    return;
end
ind1=find(vY(1)==Y);
ind2=find(vY(2)==Y);


nTrial=min(length(ind1),length(ind2));%using minimal number of trials for one category, thus balance the trial number
nframe=floor((size(Tin,3)-binsize)/binstep)+1;
[score, score_shuffle]=deal(zeros(nRepeat,nframe));
for it=1:nframe
%     %check input
%     clf;
%     histogram(Tin(ind1,:,it));hold on;
%     histogram(Tin(ind2,:,it));
    
    ind1it=max(1,(it-1)*binstep-floor(binsize/2)+1);
    ind2it=min(size(Tin,3),(it-1)*binstep+floor(binsize/2)+1);
    binActivity=Tin(:,:,ind1it:ind2it);
    binMeanActivity=nanmean(binActivity,3);%n-by-1 vector
    Tit=squeeze(binMeanActivity);
    parfor iRepeat=1:nRepeat
        indused1=randperm(length(ind1),nTrial);
        indused2=randperm(length(ind2),nTrial);
        indtrain1=indused1(1:floor(pTraining*length(indused1)));
        indtest1=indused1(floor(pTraining*length(indused1))+1:end);
        indtrain2=indused1(1:floor(pTraining*length(indused2)));
        indtest2=indused1(floor(pTraining*length(indused2))+1:end);
        if sum(~isnan(Tit(indtrain1,1)))<3||sum(~isnan(Tit(indtrain2,1)))<3
            score(iRepeat,it)=nan;
        else
            indtrain=[indtrain1,indtrain2];
            indtest=[indtest1,indtest2];
            Ttrain=Tit(indtrain,:);
            Ytrain=Y(indtrain);
            Ttrain=array2table(Ttrain);
            Mdl=fitcsvm(Ttrain,Ytrain,'Standardize',true);
            label=predict(Mdl,Tit(indtest,:));
            Ytrue=Y(indtest);
            score(iRepeat,it)=nansum(Ytrue==label)/length(Ytrue);
            indTrain_temp=randperm(length(indtrain));
            indTrain_shuffle=indtrain(indTrain_temp);
            Ytrain_shuffle=Y(indTrain_shuffle);
            indTest_temp=randperm(length(indtest));
            indTest_shuffle=indtest(indTest_temp);
            Mdl_shuffle=fitcsvm(Ttrain,Ytrain_shuffle,'Standardize',true);
            label_shuffle=predict(Mdl_shuffle,Tit([indtest1,indtest2],:));
            Yshuffle=Y(indTest_shuffle);
            score_shuffle(iRepeat,it)=nansum(Yshuffle==label_shuffle)/length(Yshuffle);
        end
    end
end
% %plot accuracy of each session
% colorScore={[1,0,0],[0,1,0],[0.5,0.5,0.5]};
% pSig=0.05;
% fPlotMean_CI(1:size(score,2),score,colorScore{1},pSig);
% fPlotMean_CI(1:size(score,2),score_shuffle,colorScore{3},pSig);
end

