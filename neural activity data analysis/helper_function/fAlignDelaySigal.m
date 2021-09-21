function [ dff_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( dff, behEventFrameIndex, frameNum,varargin )
%FALIGNDELAYSIGAL Summary of this function goes here align Cal2+ signal to
%delay onset, and mask activity after delay offset(go cue ) as nan to look
%at delay activity
%   varargin- 'BaselineCorrected'(default)|'raw'


% base=behEventFrameIndex.stimOnset;
base=behEventFrameIndex.stimOffset;%delay onset
% base=behEventFrameIndex.go;
%delayEnd=behEventFrameIndex.go;
behEvent_aligned =structfun(@(x) x-base+frameNum(1)+1,behEventFrameIndex,'UniformOutput',false);
nlength=frameNum(2)+frameNum(1)+1;
if max(behEvent_aligned.go)>frameNum(2)+frameNum(1)+1
    warning('delay longer than shown frame number');
%     return;
end
dff_aligned=nan(length(base),sum(frameNum)+1); %each row means one trial
for i=1:length(base) %caution, Subscript indices must either be real positive integers or logicals.
    if base(i)<=frameNum(1) 
        dff_aligned(i,1:nlength-(base(i)+frameNum(2)))=nan;
        dff_aligned(i,nlength+1-(base(i)+frameNum(2)):end)=dff(1:base(i)+frameNum(2));
    elseif  base(i)+frameNum(2)>length(dff)
        dff_aligned(i,1:sum(frameNum)+1-(base(i)+frameNum(2)-size(dff,2)))=dff(base(i)-frameNum(1):end);
        dff_aligned(i,sum(frameNum)+1-(base(i)+frameNum(2)-size(dff,2))+1:end)=nan;
    else
        if ~isnan(base(i)) %when align to some event that is nan,this trial remain zero
            dff_aligned(i,:)=dff(base(i)-frameNum(1):base(i)+frameNum(2));
        end
    end
    if behEvent_aligned.go(i)<frameNum(2)+frameNum(1)+1 %mask dff after go cue as nan
        dff_aligned(i,behEvent_aligned.go(i)+1:end)=nan;
    end
    if behEvent_aligned.lickFirst(i)<frameNum(2)+frameNum(1)+1 %mask dff after first lick as nan; useful for violation trials
        dff_aligned(i,behEvent_aligned.lickFirst(i)+1:end)=nan;
    end
    licking_aligned.leftLick{i} = [];
    licking_aligned.rightLick{i} =[];
end
if ~isempty(varargin) && strcmp(varargin{1},'raw')
    dff_aligned=dff_aligned;
else
    baseline=zeros(size(dff_aligned,1),1);
    for i=1:size(dff_aligned,1)
        if isnan(behEvent_aligned.stimOnset(i))
            baseline(i)=nan;
        else
            baseline(i)=mean(dff_aligned(i,1:behEvent_aligned.stimOnset(i)),2);%for each trial(each row),pre-aligned point baseline
        end
    end
    dff_aligned=dff_aligned-repmat(baseline,1,size(dff_aligned,2));%substract baseline for each trial(row)
end
end

