function [score,score_shuffle,p_ct,p_shuffle,ts_binned] = fCrossTimeSVM(Tin,Y, nRepeat, pTraining,binsize,binstep,ts_raw,varargin)
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
%   score_shuffle- shuffled version
%   p_ct- cross time p value of the svm decoding accuracy
%   p_shuffle- cross time p value of shuffled decoder accuracy
CTtype='CT';%default, means that cross time SVM, also may choose 'TT', that only calculate diag of the matrix
if ~isempty(varargin)
    CTtype=varargin{1};
end
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
ts_binned=ts(1:nframe);
[score, score_shuffle]=deal(zeros(nRepeat,nframe,nframe));
[p_ct, p_shuffle]=deal(ones(nframe, nframe));
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
        if ischar(CTtype) && (strcmp(CTtype,'TT')||strcmp(CTtype,'EE')) && it ~= it_test
            score(:,it,it_test)=nan;
            score_shuffle(:,it,it_test)=nan;
            p_ct(it,it_test)=nan;
            p_shuffle(it,it_test)=nan;
        elseif isnumeric(CTtype) && (it~=CTtype || it_test ~= CTtype)
            score(:,it,it_test)=nan;
            score_shuffle(:,it,it_test)=nan;
            p_ct(it,it_test)=nan;
            p_shuffle(it,it_test)=nan;
        elseif ischar(CTtype) && strcmp(CTtype,'nan')
            score(:,it,it_test)=nan;
            score_shuffle(:,it,it_test)=nan;
            p_ct(it,it_test)=nan;
            p_shuffle(it,it_test)=nan;
        else
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
            %calculate the p value for each point
            nboot=1000;
            [~,p_ct(it,it_test)]=fBootstrpMeanP(score(:,it,it_test),nboot,0.5);
            [~,p_shuffle(it,it_test)]=fBootstrpMeanP(score_shuffle(:,it,it_test),nboot,0.5);
        end
    end
end

end
%helper function
function [meanout,pout] = fBootstrpMeanP(activity,nboot,baseline)
%calculate the p value for each point
%Input- activity, vector or 2d matrix, 1st d is repeat number
bootmean=bootstrp(nboot,@nanmean,activity);
meanout=nanmean(bootmean);
pout=sum(bootmean<baseline)./sum(~isnan(bootmean));
pout_temp=pout;
pout(pout_temp<0.5)=pout_temp(pout_temp<0.5)*2;
pout(pout_temp>=0.5)=(1-pout_temp(pout_temp>=0.5))*2;
end