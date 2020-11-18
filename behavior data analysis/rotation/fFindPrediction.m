function [ dataPrediction ] = fFindPrediction( fileFolder )
%FFINDPREDICTION Summary of this function goes here
%   Detailed explanation goes here
%   寻找到出声音前最近一次的两侧舔的时间，用先于出声音的微秒数量表示（两侧都记录）。
%   若没有在出声音前舔，用-1表示，舔了则用出声音时刻减去最后一次舔的时刻。
%   此外同样需要记录某次trial 的应选方向以及正确与否
%   一般上，stimDuration = 300，responseDelay = 300，故Time_stimOnset之后600方开始记录选择
%   这600ms的选择比较难以分析，暂且舍去

%   读目录下的所有数据文件.beh
dirbeh=strcat(fileFolder,'\*.beh');
dirs=dir(dirbeh);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
sort(filenames);%将文件按照训练先后顺序排列起来
%***************
%     第三维4个值分别表示频率，选择正确与否，左侧出声音前最后一次舔距离出声音时间间隔，右侧时间间隔
%***************
%   另外第二维1001表示日期,error_stay
dataPrediction=zeros(length(filenames),1001,4);
for i=1:length(filenames)
    file=cell2mat(strcat(fileFolder,'\',filenames(i)));
    fid=fopen(file,'r');
    n=1;%记录是第几个trial
   while ~feof(fid)%读到文件结尾
       s=fgetl(fid);%每次读入一行
       if strfind(s,'Error_stay_number = 0')==1              %找到error_stay
           dataPrediction(i,1001,2)=0;
       else
           dataPrediction(i,1001,2)=1;
       end
       if strfind(s,'Tone Freq')==1   %存在该字符，且位置是1
           p=[0 strfind(s,' ') length(s)+1];%找到空格的位置，定位频率值位置
           freq=s(p(3)+1:p(4)-1);
           dataPrediction(i,n,1)=str2num(freq);
           for j=1:2                %读出这个trial的prediction
               s=fgetl(fid);
           end
           p=strfind(s,'/');                %每个指标都用/隔开
           if ~isempty(p)                   %防止最后一个trial只有start数据，没有具体舔的数据
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
               if  n_numLickLeft==0
                   dataPrediction(i,n,3)=0;
               else
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
                   if ~isempty(t_beforeTimeLeft)               %在出声音前有舔的数据，否则认为没有预测，也将值设为0
                       t_lastTimeLeft=min(t_beforeTimeLeft);%最后次在出声音前舔，也就是说最小的正值
                       dataPrediction(i,n,3)=t_lastTimeLeft;  
                   else
                       dataPrediction(i,n,3)=0;
                   end
               end
                %右边舔水分析
               rp=strfind(s,'Action_numLickRight');%找到刺激开始时间
               rpp=find(p==(rp-1));
               s_numLickRight=s(p(rpp)+1:p(rpp+1)-1);    %记录右侧舔动作发生次数
               q=strfind(s_numLickRight,'=');
               n_numLickRight=str2num(s_numLickRight(q(1)+1:end)); %将右侧舔的原始数据转化为数值存储
               if  n_numLickRight==0
                   dataPrediction(i,n,4)=0;
               else
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
                   if ~isempty(t_beforeTimeRight)               %在出声音前有舔的数据，否则认为没有预测，也将值设为0
                       t_lastTimeRight=min(t_beforeTimeRight);%最后次在出声音前舔，也就是说最小的正值
                       dataPrediction(i,n,4)=t_lastTimeRight;   
                   else
                       dataPrediction(i,n,4)=0;
                   end
               end
           end
           s=fgetl(fid);            %读出这个trial的正确与否
           if strfind(s,'ERROR')    %对于选择，对1，错2，miss为3，没有数据就是0
               dataPrediction(i,n,2)=2;
           elseif strfind(s,'CORRECT')
               dataPrediction(i,n,2)=1;
           else
               dataPrediction(i,n,2)=3;
           end
           n=n+1;
       end
       
   end
   filename=cell2mat(filenames(i));
   date_p=cell2mat(strfind(filenames(i),'_'));  %始终注意数据格式，要把cell改成mat
   date_m=filename(date_p(1)+1:date_p(2)-1);    %月份
   date_d=filename(date_p(2)+1:date_p(3)-1);    %日期
   date=str2double(date_m)*100+str2double(date_d);
   dataPrediction(i,1001,1)=date;
   fclose(fid);
end

end

