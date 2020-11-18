function [C] = fExtractConsecutiveVector(A,option)
%FEXTRACTCONSECUTIVEVECTOR extract consecutive vector from a given vector
%inspired by https://www.ilovematlab.cn/thread-187896-1-1.html
%   Input: A-vector to be analyzed; option-{'first','longest'} decides
%   which of the segmented result will be chose
B = [0 find(diff(A)~=1) length(A)];
segV=cell(1,length(B));
for i = 1:length(B)-1
    segV{i}=A((B(i)+1):B(i+1));
end
if strcmp(option,'first')
    C=segV{1,1};
elseif strcmp(option,'longest')
    temp=cellfun(@length,segV);
    C=segV{1,temp==max(temp)};
end

