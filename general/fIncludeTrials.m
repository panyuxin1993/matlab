function [outvector] = fIncludeTrials(indTrial2include,data,varargin)
%FINCLUDETRIALS Summary of this function goes here
%FEXCLUDETRIALS exclude specified trials as nans
%Intput- 
%indTrial2exclude-2d matrix of start and end of excluded trials as blocks
%data- raw data to be processed
%varargin{1}- 'raw','logical'(output is a true-false vector)
if length(data)==1 %here data means the trial number
    data=1:data;
end
if strcmp(indTrial2include,'all')
    indTrial2include=[1,length(data)];
elseif ischar(indTrial2include)
    indTrial2include=str2num(indTrial2include);
elseif size(indTrial2include,2)==1 %here a column verctor, each row a trial to include
    indTrial2include=repmat(indTrial2include,1,2);
end
if isempty(varargin) || (~isempty(varargin) && strcmp(varargin{1},'raw'))
    if iscell(data)
        outvector=cell(size(data));
    else
        outvector=nan(size(data));
    end       
    for i=1:size(indTrial2include,1)
        for j=indTrial2include(i,1):indTrial2include(i,2)
            if iscell(outvector)
                outvector{j}=data{j};
            else
                outvector(j)=data(j);
            end
        end
    end   
elseif (~isempty(varargin) && strcmp(varargin{1},'logical'))
    n=length(data);
    outvector=zeros(n,1);
    for i=1:size(indTrial2include,1)
        outvector(indTrial2include(i,1):indTrial2include(i,2))=1;
    end
    outvector=logical(outvector);
end
strout='Trial included is No.';
for i=1:size(indTrial2include,1)
    strout=strcat(strout,num2str(indTrial2include(i,1)),'-',num2str(indTrial2include(i,2)),',');
end       
disp(strout);
end

