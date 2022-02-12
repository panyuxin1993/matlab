function [score_tt, score_shuffle_tt,score,score_shuffle,figMeanSig] = fCrossTimeSVM(Tin,Y, nRepeat, pTraining,binsize,binstep,ts_raw)
%fCrossTimeSVM automaticly select balanced number of trials(rows) for 2 categories
%to fit SVM model and return the accuracy of the model, similar to fSVM,
%but can train model at one time point and test at another time point
%Input-
%   Tin-table of data, m-n-t matrix, m trials,n neurons, t frames/epochs
%   Y- class labels made of 0/1
%   nRepeat- 100(default),times of trainig-testing iteration
%   pTraining- 0.8(default),proportion of trials used for training
%Output-
%   score- accuracy of training, nRepeat X t_train X t_test matrix
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
ts=ts_raw(ceil(binsize/2):binstep:end);%ts also need to be binned
ts=ts(1:nframe);
[score, score_shuffle]=deal(zeros(nRepeat,nframe,nframe));
[score_tt, score_shuffle_tt]=deal(zeros(nRepeat,nframe));
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
    for it_test=1:nframe
        ind1it_test=max(1,(it_test-1)*binstep-floor(binsize/2)+1);
        ind2it_test=min(size(Tin,3),(it_test-1)*binstep+floor(binsize/2)+1);
        binActivity_test=Tin(:,:,ind1it_test:ind2it_test);
        binMeanActivity_test=nanmean(binActivity_test,3);%n-by-1 vector
        Tit_test=squeeze(binMeanActivity_test);
        parfor iRepeat=1:nRepeat
            indused1=randperm(length(ind1),nTrial);
            indused2=randperm(length(ind2),nTrial);
            indtrain1=indused1(1:floor(pTraining*length(indused1)));
            indtest1=indused1(floor(pTraining*length(indused1))+1:end);
            indtrain2=indused1(1:floor(pTraining*length(indused2)));
            indtest2=indused1(floor(pTraining*length(indused2))+1:end);
            if sum(~isnan(Tit(indtrain1,1)))<3||sum(~isnan(Tit(indtrain2,1)))<3
                score(iRepeat,it,it_test)=nan;
            else
                indtrain=[indtrain1,indtrain2];
                indtest=[indtest1,indtest2];
                Ttrain=Tit(indtrain,:);
                Ytrain=Y(indtrain);
                Ttrain=array2table(Ttrain);
                Mdl=fitcsvm(Ttrain,Ytrain,'Standardize',true);
                label=predict(Mdl,Tit_test(indtest,:));
                Ytrue=Y(indtest);
                score(iRepeat,it,it_test)=nansum(Ytrue==label)/length(Ytrue);
                indTrain_temp=randperm(length(indtrain));
                indTrain_shuffle=indtrain(indTrain_temp);
                Ytrain_shuffle=Y(indTrain_shuffle);
                indTest_temp=randperm(length(indtest));
                indTest_shuffle=indtest(indTest_temp);
                Mdl_shuffle=fitcsvm(Ttrain,Ytrain_shuffle,'Standardize',true);
                label_shuffle=predict(Mdl_shuffle,Tit_test([indtest1,indtest2],:));
                Yshuffle=Y(indTest_shuffle);
                score_shuffle(iRepeat,it,it_test)=nansum(Yshuffle==label_shuffle)/length(Yshuffle);
            end
        end
    end
end
%plot accuracy of each session
figMeanSig=figure;
set(gcf,'position',[200,200,700,200]);
colorScore={[1,0,0],[0,1,0],[0.5,0.5,0.5]};
pSig=0.05;
for iRepeat=1:nRepeat
    for iframe=1:nframe
        score_tt(iRepeat,iframe)=score(iRepeat,iframe,iframe);
        score_shuffle_tt(iRepeat,iframe)=score_shuffle(iRepeat,iframe,iframe);
    end
end
subplot(1,3,1);
fPlotMean_CI(1:size(score_tt,2),score_tt,colorScore{1},pSig);
fPlotMean_CI(1:size(score_tt,2),score_shuffle_tt,colorScore{3},pSig);
%color plot of the cross temporal decoding
ax2=subplot(1,3,2);
[ meandata, sigmat,pmat] = fPlotMeanSig2D_SVM(ax2, score,ts,pSig, 0.5,'no bar');
xlabel('Training time');
ylabel('Testing time');
title('Real data');
ax3=subplot(1,3,3);
[ meandata, sigmat,pmat] = fPlotMeanSig2D_SVM(ax3, score_shuffle,ts,pSig, 0.5,'colorbar');
xlabel('Training time');
ylabel('Testing time');
title('Shuffled data');
end

%helper function
function [ meandata, sigmat,pmat] = fPlotMeanSig2D_SVM(ax,data,ts,pSig, baseline,strbar)
%FPLOTMEANSIG2D plot mean in color and draw contour lines ot represent
%significant range; Currently mainly used for showing SVM decoding accuracy
%Input-
%   data- r-by-m-by-n matrix, m rows, n columns and r repeats, m=n
%   pSig- p value that decide whether data was significant from baseline
%   baseline- double, compared with data
%Output-
%   ax- axis to plot
%   meandata- m-by-n matrix store the mean of each point
%   ts- time stamps 
%   sigmat- m-by-n matrix made of 0-1, which indicates which position is
%   significant
%   pmat- m-by-n matrix with each point indicates the p value of that
%   posision

meandata = nanmean(data,1);
pmat=ones(size(data,2),size(data,3));
% %this method seems to slow
% for i=1:size(data,2)
%     for j=1:size(data,3)
%         test_mean=bootstrp(1000,@mean,data(:,i,j));
%         temp_p=nansum(test_mean>baseline)/length(test_mean);
%         pmat(i,j)=1-2*abs(temp_p-0.5);
%     end
% end
for i=1:size(data,2)
    for j=1:size(data,3)
        [h,temp_p]=ttest(data(:,i,j)-0.5);
        pmat(i,j)=temp_p;
    end
end
sigmat= (pmat<pSig);
%plot mean
axes(ax);
clims=[0,1];
meandata=squeeze(meandata);
imagesc(ts,ts,meandata,clims);
hold on;
if strcmp(strbar,'colorbar')
    colormap(jet);
    chandle=colorbar;
    set(chandle,'Position',[0.91,0.1,0.02,0.8]);
    chandle.Label.String= 'Decoding accuracy';
end
%plot significance contour
contour(ts,ts,sigmat,1,'--','LineColor','w');


end


