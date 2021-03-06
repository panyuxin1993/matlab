function [ defaultSide ] = fDefaultSide( fileFolder )
%FDEFAULTSIDE Summary of this function goes here
%   Detailed explanation goes here
%   将出声音前舔的多的认为是default side
%   一种方法是将一天总体来看，用总的在给光前舔的左边数量和右边数量进行分析；
%   另一种方法是每一个trial分析出是左边还是右边，再将每天的trial进行统计
dirbeh=strcat(fileFolder,'\*.beh');
dirs=dir(dirbeh);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
sort(filenames);%将文件按照训练先后顺序排列起来
defaultSide=zeros(7,length(filenames));    %三行分别存日期，左边总舔的次数和右边总舔的次数,error_stay,answerPeriod，左边default的次数和右边default的次数
for i=1:length(filenames)
    file=cell2mat(strcat(fileFolder,'\',filenames(i)));
    fid=fopen(file,'r');
    filename=cell2mat(filenames(i));
    date_p=cell2mat(strfind(filenames(i),'_'));  %始终注意数据格式，要把cell改成mat
    date_m=filename(date_p(1)+1:date_p(2)-1);    %月份
    date_d=filename(date_p(2)+1:date_p(3)-1);    %日期
    date=str2double(date_m)*100+str2double(date_d);
    defaultSide(1,i)=date;                       %第一行存日期
    while ~feof(fid)%读到文件结尾
       s=fgetl(fid);%每次读入一行
       if strfind(s,'answerPeriod')==1     %找到answerperiod
          ap=strfind(s,'=');
          answerPeriod=s(ap(1)+1:end);
          defaultSide(5,i)=str2num(answerPeriod);
       end
       if strfind(s,'Error_stay_number')==1
           esnp=strfind(s,'=');
           defaultSide(4,i)=str2num(s(esnp(1)+1:end));
       end
       t_beforeTimeLeft=[];
       t_beforeTimeRight=[];
       if strfind(s,'o/')==1   %存在该字符，且位置是1
           p=strfind(s,'/');                %每个指标都用/隔开
           tsop=strfind(s,'Time_stimOnset');%找到刺激开始时间
           pp=find(p==(tsop-1));
           s_stimOnset=s(p(pp)+1:p(pp+1)-1);  %记录刺激开始时间
           q=strfind(s_stimOnset,'=');
           t_stimOnset=str2num(s_stimOnset(q(1)+1:end));    %将原始数据转化为数值存储
           %左边舔水分析
           lp=strfind(s,'Action_numLickLeft');%找到刺激开始时间
           lpp=find(p==(lp-1));
           s_numLickLeft=s(p(lpp)+1:p(lpp+1)-1);    %记录左侧舔动作发生次数
           q=strfind(s_numLickLeft,'=');
           n_numLickLeft=str2num(s_numLickLeft(q(1)+1:end)); %将左侧舔的原始数据转化为数值存储
           if  n_numLickLeft~=0   %舔水次数不为零，下面才去记录舔水时刻，否则来声音之前必然没有舔水，可以直接进入下一个trial
               s_lickTimeLeft=s(p(lpp+1)+1:p(lpp+2)-1); %记录左边舔的时刻
               q=strfind(s_lickTimeLeft,'|');
               n_numLickLeft=length(q)-1;            %***有的数据中舔水数量统计和实际数量不一致，为了不是数组越界，此处选用实测的舔水时刻数组大小作为舔水数量**
               t_lickTimeLeft=zeros(1,n_numLickLeft);%定义行向量，存储每次舔的时刻
               for j=1:n_numLickLeft-1              %用向量记录每一次舔水时刻,防止越界，最后一个值单独处理
                   if ~isempty(str2num(s_lickTimeLeft(q(j+1)+1:q(j+2)-1)))   %竟然有这种情况|||连续的|，只好加一个判断是否为空的了
                        t_lickTimeLeft(j)=t_stimOnset-str2num(s_lickTimeLeft(q(j+1)+1:q(j+2)-1));
                   end
               end
               if ~isempty(str2num(s_lickTimeLeft(q(n_numLickLeft+1)+1:end))) %竟然有这种情况|||连续的|，只好加一个判断是否为空的了
                    t_lickTimeLeft(n_numLickLeft)=t_stimOnset-str2num(s_lickTimeLeft(q(n_numLickLeft+1)+1:end));%最后一个值单独处理                        
               end
               t_beforeTimeLeft=t_lickTimeLeft(t_lickTimeLeft>0); %大于零表示在出声音之前的舔；t_beforeTimeRight=t_lickTimeRight>0;只会得到一个逻辑矩阵，不是从中选出一些元素组成新矩阵                        
               if ~isempty(t_beforeTimeLeft)               %在出声音前有舔的数据，否则将值设为0
                    defaultSide(2,i)=defaultSide(2,i)+sum(t_beforeTimeLeft*1);
               end
           end
            %右边舔水分析
           rp=strfind(s,'Action_numLickRight');%找到刺激开始时间
           rpp=find(p==(rp-1));
           s_numLickRight=s(p(rpp)+1:p(rpp+1)-1);    %记录右侧舔动作发生次数
           q=strfind(s_numLickRight,'=');
           n_numLickRight=str2num(s_numLickRight(q(1)+1:end)); %将右侧舔的原始数据转化为数值存储
           if  n_numLickRight~=0
               s_lickTimeRight=s(p(rpp+1)+1:p(rpp+2)-1); %记录右边舔的时刻
               q=strfind(s_lickTimeRight,'|');
               n_numLickRight=length(q)-1;             %***有的数据中舔水数量统计和实际数量不一致，为了不是数组越界，此处选用实测的舔水时刻数组大小作为舔水数量**
               t_lickTimeRight=zeros(1,n_numLickRight);%定义行向量，存储每次舔的时刻
               for jr=1:n_numLickRight-1              %用向量记录每一次舔水时刻
                   if ~isempty(str2num(s_lickTimeRight(q(jr+1)+1:q(jr+2)-1))) %竟然有这种情况|||连续的|，只好加一个判断是否为空的了
                        t_lickTimeRight(jr)=t_stimOnset-str2num(s_lickTimeRight(q(jr+1)+1:q(jr+2)-1));%******数据文件有误，统计舔总数与实际有记录的各次舔的时刻并非一致。所以只能够用实际数到的有个数的舔时刻来定义数组大小
                   end
               end
               if  ~isempty(str2num(s_lickTimeRight(q(n_numLickRight+1)+1:end)))%竟然有这种情况|||连续的|，只好加一个判断是否为空的了
                    t_lickTimeRight(n_numLickRight)=t_stimOnset-str2num(s_lickTimeRight(q(n_numLickRight+1)+1:end));%最后一个单独对待           
               end
               t_beforeTimeRight=t_lickTimeRight(t_lickTimeRight>0); %大于零表示在出声音之前的舔；t_beforeTimeRight=t_lickTimeRight>0;只会得到一个逻辑矩阵，不是从中选出一些元素组成新矩阵
               if ~isempty(t_beforeTimeRight)               %在出声音前有舔的数据，否则将值设为0
                   defaultSide(3,i)=defaultSide(3,i)+sum(t_beforeTimeRight*1);
               end
           end
           %trial的左右判断
          if   sum(t_beforeTimeLeft)>sum(t_beforeTimeRight)
              defaultSide(6,i)=defaultSide(6,i)+1;
          elseif sum(t_beforeTimeLeft)<sum(t_beforeTimeRight)
              defaultSide(7,i)=defaultSide(7,i)+1;
          end
       end
    end
    fclose(fid);
end
end

