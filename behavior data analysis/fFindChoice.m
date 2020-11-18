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
for i=length(dirs):-1:1%去掉空文件
    if dirs(i).bytes==0
        dirs(i)=[];
    end
end

dircell=struct2cell(dirs);
    filenames=dircell(1,:);
sort(filenames);%将文件按照训练先后顺序排列起来
if fileNum<length(filenames)%如果指定读取n个文件则读取按照顺序最后的那n个文件
    filenames=filenames(length(filenames)-fileNum+1:end);
end
namep=strfind(fileFolder,'\');
animal_name = fileFolder(namep(end)+1:end);
%用于存放某只小鼠某天某个trial的click_rate和选择结果，
%dataChoice(x,1001,1)存放日期数据，以4位数表示,dataChoice(x,1001,2)表示answer period
%dataChoice(x,1002,1)存放'clickRate_L'，dataChoice(x,1002,2)存放clickRate_R
dataChoice=zeros(length(filenames),1002,2);
%dataChoice=zeros(length(filenames),2002,2);%两个文件合起来
for i=1:length(filenames)
    file=cell2mat(strcat(fileFolder,'\',filenames(i)));
    fid=fopen(file,'r');
    n=1;%记录是第几个trial
   while ~feof(fid)%读到文件结尾
       s=fgetl(fid);%每次读入一行
       if strfind(s,'answerPeriod')==1     %找到answerperiod
          ap=strfind(s,'=');
          answerPeriod=s(ap(1)+1:end);
          dataChoice(i,1001,2)=str2num(answerPeriod);
       end
       if strfind(s,'clickRate_L')==1   %找到该session对应的click_rate和行为的规则
           cp=strfind(s,'=');
           clickRate_L=str2num(s(cp(1)+1:end));
       end
       if strfind(s,'clickRate_R')==1   %找到该session对应的click_rate和行为的规则
           cp=strfind(s,'=');
           clickRate_R=str2num(s(cp(1)+1:end));
       end
       if contains(s,'Click Rate' )  %存在该字符，且位置是1
           p=[0 strfind(s,':') length(s)+1];%找到空格的位置，定位频率值位置
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
%            for j=1:3 %自动存的是j=1:6，手动存的是j=1:3  %读出这个trial的正确与否；存在潜在的bug是当记录文件出现乱码时，可能跳过一行，即需要j=1:4方到达结果
%                s=fgetl(fid);
%            end
%            if isempty(s)  %也不是个办法，有些文档会出现其它情况，最好的办法是把空行全部去掉
%                for j=1:3
%                    s=fgetl(fid);
%                end
%            end
           if strfind(s,'ERROR')%对于选择，对1，错2，miss为3，violation为4，没有数据就是0
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
   date_p=cell2mat(strfind(filenames(i),'_'));  %始终注意数据格式，要把cell改成mat
   if length(date_p)==3 %老版文件命名方式
        date_m=filename(date_p(1)+1:date_p(2)-1);    %月份
        date_d=filename(date_p(2)+1:date_p(3)-1);    %日期
        date=str2double(date_m)*100+str2double(date_d);
   elseif length(date_p)==1%新版文件命名方式
       date=str2double(filename(date_p(1)+5:date_p(1)+8));
   end
   
   %dataChoice(i,1001,1)=date;
   dataChoice(i,1002,1)=clickRate_L;
   dataChoice(i,1002,2)=clickRate_R;
   fclose(fid);
end


end

