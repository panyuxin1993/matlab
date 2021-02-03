function [ dff_aligned, behEvent_aligned, licking_aligned ] = fAlignSigBehEvent( dff, behEventFrameIndex, lickingFrameIndex, alignTo,frameNum,varargin )
%FALIGNSIGBEHEVENT Align signals to different behavior events
%   Detailed explanation 
%Input is dff( in a continues whole-session-long trial,is a vector for a ROI), 
%   behavior event frame index( struct, including stimOnset, answer time, 
%   first lick and so on), 
%   align to which event(string can be in {'stimOnset', 'go cue','first lick',...
%   'first left lick','first right lick', 'answer','reward'}, and
%   frame number( [frame number before 0, frame number after 0])
%   varargin- 'BaselineCorrected'(default)|'raw'
%Output is dff_aligned( matrix of dff, each row a trial),
%   behEvent_aligned(struct, each event stored in one field, which is a vector)
%  
switch alignTo
    case 'stim onset'
        base=behEventFrameIndex.stimOnset;
    case 'delay onset'
        base=behEventFrameIndex.stimOffset;%delay onset
    case 'go cue'
        base=behEventFrameIndex.go;
    case 'first lick'
        base=behEventFrameIndex.lickFirst;
    case 'first left lick'
        base=behEventFrameIndex.lickFirst_left;
    case 'first right lick'
        base=behEventFrameIndex.lickFirst_right;
    case 'last lick'
        base=behEventFrameIndex.lickLast;
    case 'answer'
        base=behEventFrameIndex.ansTime;
    case 'reward'
        base=behEventFrameIndex.rewTime;
    otherwise
        base=behEventFrameIndex.start;%using the 1st frame in a trial
end
behEvent_aligned =structfun(@(x) x-base+frameNum(1)+1,behEventFrameIndex,'UniformOutput',false);
dff_aligned=nan(length(base),sum(frameNum)+1); %each row means one trial
for i=1:length(base) %caution, Subscript indices must either be real positive integers or logicals.
    if base(i)<=frameNum(1) 
        dff_aligned(i,1:frameNum(1)-base(i)+1)=nan;
        dff_aligned(i,frameNum(1)-base(i)+2:end)=dff(1:base(i)+frameNum(2));
    elseif  base(i)+frameNum(2)>size(dff,2)
%         dff_aligned(i,1:base(i)+frameNum(2)-size(dff,2))=dff(base(i)-frameNum(1):end);
%         dff_aligned(i,base(i)+frameNum(2)-size(dff,2)+1:end)=nan;
        dff_aligned(i,1:sum(frameNum)+1-(base(i)+frameNum(2)-size(dff,2)))=dff(base(i)-frameNum(1):end);
        dff_aligned(i,sum(frameNum)+1-(base(i)+frameNum(2)-size(dff,2))+1:end)=nan;
    else
        if ~isnan(base(i)) %when align to some event that is nan,this trial remain zero
            dff_aligned(i,:)=dff(1,base(i)-frameNum(1):base(i)+frameNum(2));
        end
    end
    % align licking raster
    licking_aligned.leftLick{i} = lickingFrameIndex.leftLick{i} - base(i)+frameNum(1)+1;
    licking_aligned.rightLick{i} = lickingFrameIndex.rightLick{i}- base(i)+frameNum(1)+1;
end
if ~isempty(varargin) && strcmp(varargin{1},'raw')
    dff_aligned=dff_aligned;
elseif strcmp(alignTo,'stim onset')|| strcmp(alignTo,'delay onset') %and other condition
    for i=1:size(dff_aligned,1)
        if isnan(behEvent_aligned.stimOnset(i))
            baseline(i)=nan;
        else
            baseline(i)=nanmean(dff_aligned(i,1:behEvent_aligned.stimOnset(i)),2);%for each trial(each row),pre-aligned point baseline
        end
    end
    baseline=reshape(baseline,[],1);
    dff_aligned=dff_aligned-repmat(baseline,1,size(dff_aligned,2));%substract baseline for each trial(row)
end
end

