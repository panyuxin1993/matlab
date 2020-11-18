function [animal_name,history] = fGetLearningHistory(folder,option)
%FGETLEARNINGHISTORY get learning history of one animal
%   input are folder-folder of data from one animal;
%   option-{'num','txt}(default is num), specified output is number/string
%   output are animal_name; history-a vector 
dirmat=strcat(folder,'\*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
sort(filenames);%将文件按照训练先后顺序排列起来
nfile=length(filenames);
history=cell(nfile,1);
if strcmp(option,'txt')
    dict={'2AFC','retract','delay','violation'};
elseif strcmp(option,'num')
    dict={1,2,3,4};
end
for i=1:nfile
    load(strcat(folder,'\',filenames{i}));
    animal_name=SessionSettings{1, 1}.animalName;%here for simplicity, in one folder just one animal
    if contains(SessionSettings{1, 1}.ExpCondition,'retract')%maybe for some cases, default is 'retract'
        history{i,1}=dict{2};
    elseif isfield(SessionResults,'Trial_isRetractTrial') && SessionResults{1, 1}.Trial_isRetractTrial>0
        history{i,1}=dict{2};
    elseif contains(SessionSettings{1, 1}.ExpCondition,'delay')
        history{i,1}=dict{3};
    elseif contains(SessionSettings{1, 1}.ExpCondition,'vio')
        history{i,1}=dict{4};
    else
        history{i,1}=dict{1};
    end
        
end

