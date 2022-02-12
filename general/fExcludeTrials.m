function [outvector] = fExcludeTrials(indTrial2exclude,data, varargin)
%FEXCLUDETRIALS exclude specified trials as nans
%Intput- 
%indTrial2exclude-2d matrix of start and end of excluded trials as blocks
%data- raw data to be processed
%varargin{1}- 'raw','logical'(output is a true-false vector)
if size(indTrial2exclude,2)==1 %here a column verctor, each row a trial to exclude
    indTrial2exclude=repmat(indTrial2exclude,1,2);
end
if length(data)==1 %here data means the trial number
    data=1:data;
end
if isempty(varargin) || (~isempty(varargin) && strcmp(varargin{1},'raw'))
    outvector=data;
    for i=1:size(indTrial2exclude,1)
        for j=indTrial2exclude(i,1):indTrial2exclude(i,2)
            if iscell(outvector)
                outvector{j}=[];
            else
                outvector(j)=nan;
            end
        end
    end   
elseif (~isempty(varargin) && strcmp(varargin{1},'logical'))
    n=length(data);
    outvector=ones(n,1);
    for i=1:size(indTrial2exclude,1)
        outvector(indTrial2exclude(i,1):indTrial2exclude(i,2))=0;
    end
    outvector=logical(outvector);
end
strout='Trial included is No.';
for i=1:size(indTrial2exclude,1)
    strout=strcat(strout,num2str(indTrial2exclude(i,1)),'-',num2str(indTrial2exclude(i,2)),',');
end       
disp(strout);
end

