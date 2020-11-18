function [ timepoint,endDff,neuralActivityGrouped ] = fGetEndDff( neuralActivity, smooth,varargin )
%FGETENDDFF using neural activity data to get end dff with correspoinding
%timepoint
%   Input- neuralActivity, nTrial-by-nFrame
%   Output- endDff, vector of dff values; timepoint, vector of
%   corresponding timepoint of end point; neuralActivityGrouped, grouped
%   neural activities according to varargin
timepoint=zeros(size(neuralActivity,1),1);
endDff=zeros(size(neuralActivity,1),1);
for i=1:size(neuralActivity,1)
    temp=find(isnan(neuralActivity(i,:)));
    if ~isempty(temp)
        timepoint(i)=temp(1)-1;
    else
        timepoint(i)=size(neuralActivity,2);
    end
    if timepoint(i)==0
        endDff(i)=nan;
        timepoint(i)=nan;
    else
        endDff(i)=neuralActivity(i,timepoint(i));
    end
end
if ~strcmp(smooth,'smooth')%smooth ='raw'
    neuralActivityGrouped=neuralActivity;
else
    [ timepoint,endDff,neuralActivityGrouped ] =fSmooth(timepoint,endDff,neuralActivity,varargin);
end
end

function [ xout,yout ,neuralActivityGrouped] = fSmooth( xin,yin,neuralActivity,varargin )
%FSMOOTH smooth curve
%   input- x and y; varargin indicated which methods to smooth
%   output- x and y, note the length may change;

if isempty(varargin)
    method='bin';
    bin=range(xin)/10;
else
    method=varargin{1}{1};
    if strcmp(method,'bin')
        bin=varargin{1}{2};
    elseif strcmp(method,'span')
        span=varargin{1}{2};
    end
end
[x,I]=sort(xin);
y=yin(I);
%here I need to leave out some data points where like x=[12 27 27 28 28
%...] the first one must be abnormal
diffx=diff(x);
for i=1:length(diffx)/100%only scan for the 1st one third where error may happen; later when delay is longer, it is usually happen,so do not scan that
    if diffx(i)>1%abnormal
        x(1:i)=nan;
        y(1:i)=nan;
    end
end       
%%%
neuralActivity=neuralActivity(I,:);
if strcmp(method,'bin')
    nbin=ceil(range(x)/bin);
    xout=zeros(nbin,1);
    yout=zeros(nbin,1);
    neuralActivityGrouped=zeros(nbin,size(neuralActivity,2));
    limspace=nanmin(x):bin:nanmax(x)+bin-mod(range(x),bin);%here plus bin-mod.. is to ensure limspace cover range(x) 
    %limspace=x(2):bin:max(x)+bin-mod(range(x),bin);%using x(1) may have some unknown small value
%     disp(limspace);
    for i=1:nbin
        ind=logical((limspace(i)<x).*(x<limspace(i+1)));
        xout(i)=mean(x(ind));
        yout(i)=mean(y(ind));
        neuralActivityGrouped(i,:)=nanmean(neuralActivity(ind,:),1);%here, rather than use nanmean, I garantee that all the mean data are from not-nan data,reduce variability;however, some data has a nan row,so still use nanmean
    end 
elseif strcmp(method,'span')
    n=length(x)-span+1;
    xout=zeros(n,1);
    yout=zeros(n,1);
    neuralActivityGrouped=zeros(n,size(neuralActivity,2));
    for i=1:n
        xout(i)=mean(x(i:i+span-1));
        yout(i)=mean(y(i:i+span-1));
        neuralActivityGrouped(i,:)=nanmean(neuralActivity(i:i+span-1,:),1);
    end
end
end