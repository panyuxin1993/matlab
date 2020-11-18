%automatic reading folder information and fill a summary table in appended
%way, which can be used in both 2P and FP dataset

%parameters for 2P
%path='F:\2P';
filepath='F:\video tracking\M2 imaging video';
filename = [filepath, filesep,'imaging_data_summary.xlsx'];
% % parameters for FP
% filepath='F:\FP';
% filename = [path, filesep,'FP_data_summary.xlsx'];

files = dir(strcat(filepath,'\*2AFC.beh*'));%select files with certain motif in file name
varTypes = {'datetime','string','string'};
T=table('Size',[length(files),3],'VariableTypes',varTypes,'VariableNames',{'date','animal','session'});
formatOut = 'yyyy/mm/dd';
for n_session=1:length(files)
    name_session=strsplit(files(n_session).name,'_');
    T.animal(n_session)=name_session(1);
    [Y,M,D]=datevec(name_session(2),'yyyymmdd');
    T.date(n_session)=datetime([Y,M,D],'InputFormat','yyyymmdd');%datestr(date_vector);
    T.session(n_session)=files(n_session).name;
end
if ~exist(filename,'file')    
    writetable(T,filename,'Sheet',1);
else
    [num,txt,raw] =xlsread(filename,1);
    row_ind=size(raw,1)+1;
    session_list=raw(:,3);%'session' analog to ID of each row
    for indT=1:size(T,1)
        if sum(sum(strcmp(T.session(indT),session_list)))==0%not included in existed file, then appended a row
            xlrange=['A',num2str(row_ind),':C',num2str(row_ind)];
            writetable(T(indT,:),filename,'Sheet',1,'Range',xlrange,'WriteVariableNames',false);
            row_ind=row_ind+1;
        end
    end    
end