function [ animal_name,dataChoice ] = fFindChoice(  fileFolder,fileNum )
%fFINDCHOICE Summary of this function goes here
%   read .txt/.beh file, only two column, stimulus and
%   result(1-correct,2-error,3-miss,4-violation,0-default,no data)
%   fileNum means number of files to open, if larger than total number of
%   file in this folder, then open all files.
%   One fileFolder means one animal, so combine the data into one variable
%   Detailed explanation goes here
dirbeh=strcat(fileFolder,'\*.beh');
%dirbeh=strcat(fileFolder,'\*.txt');

dirs=dir(dirbeh);
for i=length(dirs):-1:1%ȥ�����ļ�
    if dirs(i).bytes==0
        dirs(i)=[];
    end
end

dircell=struct2cell(dirs);
    filenames=dircell(1,:);
sort(filenames);%���ļ�����ѵ���Ⱥ�˳����������
if fileNum<length(filenames)%���ָ����ȡn���ļ����ȡ����˳��������n���ļ�
    filenames=filenames(length(filenames)-fileNum+1:end);
end
namep=strfind(fileFolder,'\');
animal_name = fileFolder(namep(end)+1:end);
%���ڴ��ĳֻС��ĳ��ĳ��trial��click_rate��ѡ������
%dataChoice(x,1001,1)����������ݣ���4λ����ʾ,dataChoice(x,1001,2)��ʾanswer period
%dataChoice(x,1002,1)���'clickRate_L'��dataChoice(x,1002,2)���clickRate_R
dataChoice=zeros(length(filenames),1002,2);
%dataChoice=zeros(length(filenames),2002,2);%�����ļ�������
for i=1:length(filenames)
    file=cell2mat(strcat(fileFolder,'\',filenames(i)));
    fid=fopen(file,'r');
    n=1;%��¼�ǵڼ���trial
   while ~feof(fid)%�����ļ���β
       s=fgetl(fid);%ÿ�ζ���һ��
       if strfind(s,'answerPeriod')==1     %�ҵ�answerperiod
          ap=strfind(s,'=');
          answerPeriod=s(ap(1)+1:end);
          dataChoice(i,1001,2)=str2num(answerPeriod);
       end
       if strfind(s,'clickRate_L')==1   %�ҵ���session��Ӧ��click_rate����Ϊ�Ĺ���
           cp=strfind(s,'=');
           clickRate_L=str2num(s(cp(1)+1:end));
       end
       if strfind(s,'clickRate_R')==1   %�ҵ���session��Ӧ��click_rate����Ϊ�Ĺ���
           cp=strfind(s,'=');
           clickRate_R=str2num(s(cp(1)+1:end));
       end
       if contains(s,'Click Rate' )  %���ڸ��ַ�����λ����1
           p=[0 strfind(s,':') length(s)+1];%�ҵ��ո��λ�ã���λƵ��ֵλ��
           freq=s(p(2)+1:end);
           dataChoice(i,n,1)=str2num(freq);
       end 
       if contains(s,'o/')
           j=1;
           while j<2
               s=fgetl(fid);
               if ~isempty(s)
                   j=j+1;
               end
           end
%            for j=1:3 %�Զ������j=1:6���ֶ������j=1:3  %�������trial����ȷ��񣻴���Ǳ�ڵ�bug�ǵ���¼�ļ���������ʱ����������һ�У�����Ҫj=1:4��������
%                s=fgetl(fid);
%            end
%            if isempty(s)  %Ҳ���Ǹ��취����Щ�ĵ�����������������õİ취�ǰѿ���ȫ��ȥ��
%                for j=1:3
%                    s=fgetl(fid);
%                end
%            end
           if strfind(s,'ERROR')%����ѡ�񣬶�1����2��missΪ3��violationΪ4��û�����ݾ���0
               dataChoice(i,n,2)=2;
           elseif strfind(s,'CORRECT')
               dataChoice(i,n,2)=1;
           elseif strfind(s,'MISS')
               dataChoice(i,n,2)=3;
           elseif strfind(s,'VIOLATION')
               dataChoice(i,n,2)=4;
           end
           n=n+1;
       end
   end
   filename=cell2mat(filenames(i));
   date_p=cell2mat(strfind(filenames(i),'_'));  %ʼ��ע�����ݸ�ʽ��Ҫ��cell�ĳ�mat
   if length(date_p)==3 %�ϰ��ļ�������ʽ
        date_m=filename(date_p(1)+1:date_p(2)-1);    %�·�
        date_d=filename(date_p(2)+1:date_p(3)-1);    %����
        date=str2double(date_m)*100+str2double(date_d);
   elseif length(date_p)==1%�°��ļ�������ʽ
       date=str2double(filename(date_p(1)+5:date_p(1)+8));
   end
   
   %dataChoice(i,1001,1)=date;
   dataChoice(i,1002,1)=clickRate_L;
   dataChoice(i,1002,2)=clickRate_R;
   fclose(fid);
end


end

