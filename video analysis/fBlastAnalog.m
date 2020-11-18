function [ indexesRefInData ] = fBlastAnalog( refvector,datavector,varargin )
%FBLASTANALOG using similar idea to treat analog data format
%Input-
%   refvector- usually the ground true, a vector to be aligned
%   datavector- usually a predicted data to be test, a vector not
%   neccessarily same length as datavector
%   criteria- [0.9,1.1](default), a range decide whether element in
%   datavector match any in refvector.
%Output-
%   indexesRefInData- vector of same length as refvector, each point to the
%   index of element in datavector that match the element in refvector; if
%   multiple possible candidates, choose the most likely one(e.g. least
%   abs(delta(ref,data)) as the criteria); if non corresponding candidates,
%   set as nan.

if isempty(varargin)
    criteria=[0.99,1.01];
else
    criteria=varargin{1};
end
matMatch=false(length(refvector),length(datavector));
matRatio=zeros(length(refvector),length(datavector));
for i=1:length(refvector)
    matMatch(i,:)=(refvector(i)*criteria(2)>=datavector).*(refvector(i)*criteria(1)<=datavector);
    matRatio(i,:)=abs(datavector/refvector(i)-1);
end
figBlast=figure;
set(gcf,'Position',[50,50,700,200]);
subplot(1,3,1);
imagesc(matMatch);
title('matched element');
xlabel('input2');
ylabel('input1');
subplot(1,3,2);
clim=[0,criteria(2)-1];
imagesc(matRatio,clim);
title('value ratio');
xlabel('input2');
ylabel('input1');
colorbar;
%% blast, find longest continues match
matchStartCell=cell(1,length(refvector));
matchLengthCell=cell(1,length(refvector));
for i=1:length(refvector)
    templength=zeros(1,sum(matMatch(i,:)));%these can be start point
    if isempty(templength)
        matchLengthCell{i}=nan;
        matchStartCell{i}=nan;
    else
        tempStart=templength;
        iStart=0;
        for j=1:length(datavector)
            startpoint=matMatch(i,j);
            if startpoint
                iStart=iStart+1;
                tempStart(iStart)=j;
                templength(iStart)=1;
                for k=1:min(length(datavector)-j,length(refvector)-i)
                    if matMatch(i+k,j+k)
                        templength(iStart)=k+1;
                    else
                        break;
                    end
                end     
            end
        end
        matchLengthCell{i}=templength;
        matchStartCell{i}=tempStart;
    end
end
nUnmatched=length(refvector);  
indUnmatched=true(1,length(refvector));
matchLength=nan(1,length(refvector));
matchStart=nan(1,length(refvector));
while nUnmatched>0
    maxMatchLength=cellfun(@max, matchLengthCell);
    maxMatchLength(~indUnmatched)=nan;
    valueMatchLength=max(maxMatchLength(indUnmatched));
    temp=find(maxMatchLength==valueMatchLength);
    if isempty(temp)
        break;
    end
    indMatchLength=temp(1);%if has many blocks
    indtemp0=find(matchLengthCell{indMatchLength}==max(matchLengthCell{indMatchLength}));
    indtemp=indtemp0;
    if indMatchLength>1 && ~isnan(max(matchStart(1:indMatchLength-1)))
        indtemp(matchStartCell{indMatchLength}(indtemp0)<=max(matchStart(1:indMatchLength-1)))=nan;
    end
    if indMatchLength<length(refvector) && ~isnan(min(matchStart(indMatchLength+1:end)))
        indtemp(matchStartCell{indMatchLength}(indtemp0)>=min(matchStart(indMatchLength+1:end)))=nan;
    end
    indtemp(isnan(indtemp))=[];%rule out inappropriate dots
    if isempty(indtemp)
        matchLengthCell{indMatchLength}(indtemp0)=nan;%if the max-length dots are not appropriate, then set to nan and continue loops
        matchStartCell{indMatchLength}(indtemp0)=nan;
    else
        tempvar=1:valueMatchLength;
        matchStart(indMatchLength:indMatchLength+valueMatchLength-1)=matchStartCell{indMatchLength}(indtemp(1))+tempvar-1;
        matchLength(indMatchLength:indMatchLength+valueMatchLength-1)=valueMatchLength-tempvar+1;
        indUnmatched(indMatchLength:indMatchLength+valueMatchLength-1)=false;
        nUnmatched=nUnmatched-valueMatchLength;
    end
end
figure(figBlast);
subplot(1,3,3);
for i=1:length(refvector)
    curveraw=plot([i,i],[matchStart(i),matchStart(i)+matchLength(i)],'k-','LineWidth',2);hold on;
end
%% rule out some single points
indstayMatchStart=false(size(matchStart));
indstayMatchStart(1)=true;
indstayMatchStart(end)=true;
for i=2:length(refvector)-1
    if matchStart(i)>matchStart(i-1)&&matchStart(i)<matchStart(i+1)
        indstayMatchStart(i)=true;
    end
end
for i=2:length(refvector)-1
    if indstayMatchStart(i-1) && ~indstayMatchStart(i) && matchStart(i)==matchStart(i-1)+1
        indstayMatchStart(i)=true;
    elseif indstayMatchStart(i+1) && ~indstayMatchStart(i) && matchStart(i)==matchStart(i+1)-1
        indstayMatchStart(i)=true;
    end
end
matchStart(~indstayMatchStart)=nan;
figure(figBlast);
subplot(1,3,3);
for i=1:length(refvector)
    curverefined=plot([i,i],[matchStart(i),matchStart(i)+matchLength(i)],'r-','LineWidth',2);hold on;
end
legend([curveraw,curverefined],'raw','refined');
xlabel('start trial');
ylabel('covered trial range');
indexesRefInData=matchStart;
% close(figBlast);
end

