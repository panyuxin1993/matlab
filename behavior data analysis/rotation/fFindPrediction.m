function [ dataPrediction ] = fFindPrediction( fileFolder )
%FFINDPREDICTION Summary of this function goes here
%   Detailed explanation goes here
%   Ѱ�ҵ�������ǰ���һ�ε��������ʱ�䣬�����ڳ�������΢��������ʾ�����඼��¼����
%   ��û���ڳ�����ǰ����-1��ʾ���������ó�����ʱ�̼�ȥ���һ�����ʱ�̡�
%   ����ͬ����Ҫ��¼ĳ��trial ��Ӧѡ�����Լ���ȷ���
%   һ���ϣ�stimDuration = 300��responseDelay = 300����Time_stimOnset֮��600����ʼ��¼ѡ��
%   ��600ms��ѡ��Ƚ����Է�����������ȥ

%   ��Ŀ¼�µ����������ļ�.beh
dirbeh=strcat(fileFolder,'\*.beh');
dirs=dir(dirbeh);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
sort(filenames);%���ļ�����ѵ���Ⱥ�˳����������
%***************
%     ����ά4��ֵ�ֱ��ʾƵ�ʣ�ѡ����ȷ�����������ǰ���һ������������ʱ�������Ҳ�ʱ����
%***************
%   ����ڶ�ά1001��ʾ����,error_stay
dataPrediction=zeros(length(filenames),1001,4);
for i=1:length(filenames)
    file=cell2mat(strcat(fileFolder,'\',filenames(i)));
    fid=fopen(file,'r');
    n=1;%��¼�ǵڼ���trial
   while ~feof(fid)%�����ļ���β
       s=fgetl(fid);%ÿ�ζ���һ��
       if strfind(s,'Error_stay_number = 0')==1              %�ҵ�error_stay
           dataPrediction(i,1001,2)=0;
       else
           dataPrediction(i,1001,2)=1;
       end
       if strfind(s,'Tone Freq')==1   %���ڸ��ַ�����λ����1
           p=[0 strfind(s,' ') length(s)+1];%�ҵ��ո��λ�ã���λƵ��ֵλ��
           freq=s(p(3)+1:p(4)-1);
           dataPrediction(i,n,1)=str2num(freq);
           for j=1:2                %�������trial��prediction
               s=fgetl(fid);
           end
           p=strfind(s,'/');                %ÿ��ָ�궼��/����
           if ~isempty(p)                   %��ֹ���һ��trialֻ��start���ݣ�û�о����������
                tsop=strfind(s,'Time_stimOnset');%�ҵ��̼���ʼʱ��
                pp=find(p==(tsop-1));
               s_stimOnset=s(p(pp)+1:p(pp+1)-1);  %��¼�̼���ʼʱ��
               q=strfind(s_stimOnset,'=');
               t_stimOnset=str2num(s_stimOnset(q(1)+1:end));    %��ԭʼ����ת��Ϊ��ֵ�洢
               %�����ˮ����
               lp=strfind(s,'Action_numLickLeft');%�ҵ��̼���ʼʱ��
               lpp=find(p==(lp-1));
               s_numLickLeft=s(p(lpp)+1:p(lpp+1)-1);    %��¼���������������
               q=strfind(s_numLickLeft,'=');
               n_numLickLeft=str2num(s_numLickLeft(q(1)+1:end)); %��������ԭʼ����ת��Ϊ��ֵ�洢
               if  n_numLickLeft==0
                   dataPrediction(i,n,3)=0;
               else
                   s_lickTimeLeft=s(p(lpp+1)+1:p(lpp+2)-1); %��¼������ʱ��
                   q=strfind(s_lickTimeLeft,'|');
                   n_numLickLeft=length(q)-1;            %***�е���������ˮ����ͳ�ƺ�ʵ��������һ�£�Ϊ�˲�������Խ�磬�˴�ѡ��ʵ�����ˮʱ�������С��Ϊ��ˮ����**
                   t_lickTimeLeft=zeros(1,n_numLickLeft);%�������������洢ÿ�����ʱ��
                   for j=1:n_numLickLeft-1              %��������¼ÿһ����ˮʱ��,��ֹԽ�磬���һ��ֵ��������
                       if ~isempty(str2num(s_lickTimeLeft(q(j+1)+1:q(j+2)-1)))   %��Ȼ���������|||������|��ֻ�ü�һ���ж��Ƿ�Ϊ�յ���
                            t_lickTimeLeft(j)=t_stimOnset-str2num(s_lickTimeLeft(q(j+1)+1:q(j+2)-1));
                       end
                   end
                   if ~isempty(str2num(s_lickTimeLeft(q(n_numLickLeft+1)+1:end))) %��Ȼ���������|||������|��ֻ�ü�һ���ж��Ƿ�Ϊ�յ���
                        t_lickTimeLeft(n_numLickLeft)=t_stimOnset-str2num(s_lickTimeLeft(q(n_numLickLeft+1)+1:end));%���һ��ֵ��������                        
                   end
                   t_beforeTimeLeft=t_lickTimeLeft(t_lickTimeLeft>0); %�������ʾ�ڳ�����֮ǰ����t_beforeTimeRight=t_lickTimeRight>0;ֻ��õ�һ���߼����󣬲��Ǵ���ѡ��һЩԪ������¾���                        
                   if ~isempty(t_beforeTimeLeft)               %�ڳ�����ǰ��������ݣ�������Ϊû��Ԥ�⣬Ҳ��ֵ��Ϊ0
                       t_lastTimeLeft=min(t_beforeTimeLeft);%�����ڳ�����ǰ��Ҳ����˵��С����ֵ
                       dataPrediction(i,n,3)=t_lastTimeLeft;  
                   else
                       dataPrediction(i,n,3)=0;
                   end
               end
                %�ұ���ˮ����
               rp=strfind(s,'Action_numLickRight');%�ҵ��̼���ʼʱ��
               rpp=find(p==(rp-1));
               s_numLickRight=s(p(rpp)+1:p(rpp+1)-1);    %��¼�Ҳ�������������
               q=strfind(s_numLickRight,'=');
               n_numLickRight=str2num(s_numLickRight(q(1)+1:end)); %���Ҳ����ԭʼ����ת��Ϊ��ֵ�洢
               if  n_numLickRight==0
                   dataPrediction(i,n,4)=0;
               else
                   s_lickTimeRight=s(p(rpp+1)+1:p(rpp+2)-1); %��¼�ұ����ʱ��
                   q=strfind(s_lickTimeRight,'|');
                   n_numLickRight=length(q)-1;             %***�е���������ˮ����ͳ�ƺ�ʵ��������һ�£�Ϊ�˲�������Խ�磬�˴�ѡ��ʵ�����ˮʱ�������С��Ϊ��ˮ����**
                   t_lickTimeRight=zeros(1,n_numLickRight);%�������������洢ÿ�����ʱ��
                   for jr=1:n_numLickRight-1              %��������¼ÿһ����ˮʱ��
                       if ~isempty(str2num(s_lickTimeRight(q(jr+1)+1:q(jr+2)-1))) %��Ȼ���������|||������|��ֻ�ü�һ���ж��Ƿ�Ϊ�յ���
                            t_lickTimeRight(jr)=t_stimOnset-str2num(s_lickTimeRight(q(jr+1)+1:q(jr+2)-1));%******�����ļ�����ͳ����������ʵ���м�¼�ĸ������ʱ�̲���һ�¡�����ֻ�ܹ���ʵ���������и�������ʱ�������������С
                       end
                   end
                   if  ~isempty(str2num(s_lickTimeRight(q(n_numLickRight+1)+1:end)))%��Ȼ���������|||������|��ֻ�ü�һ���ж��Ƿ�Ϊ�յ���
                        t_lickTimeRight(n_numLickRight)=t_stimOnset-str2num(s_lickTimeRight(q(n_numLickRight+1)+1:end));%���һ�������Դ�           
                   end
                   t_beforeTimeRight=t_lickTimeRight(t_lickTimeRight>0); %�������ʾ�ڳ�����֮ǰ����t_beforeTimeRight=t_lickTimeRight>0;ֻ��õ�һ���߼����󣬲��Ǵ���ѡ��һЩԪ������¾���
                   if ~isempty(t_beforeTimeRight)               %�ڳ�����ǰ��������ݣ�������Ϊû��Ԥ�⣬Ҳ��ֵ��Ϊ0
                       t_lastTimeRight=min(t_beforeTimeRight);%�����ڳ�����ǰ��Ҳ����˵��С����ֵ
                       dataPrediction(i,n,4)=t_lastTimeRight;   
                   else
                       dataPrediction(i,n,4)=0;
                   end
               end
           end
           s=fgetl(fid);            %�������trial����ȷ���
           if strfind(s,'ERROR')    %����ѡ�񣬶�1����2��missΪ3��û�����ݾ���0
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
   date_p=cell2mat(strfind(filenames(i),'_'));  %ʼ��ע�����ݸ�ʽ��Ҫ��cell�ĳ�mat
   date_m=filename(date_p(1)+1:date_p(2)-1);    %�·�
   date_d=filename(date_p(2)+1:date_p(3)-1);    %����
   date=str2double(date_m)*100+str2double(date_d);
   dataPrediction(i,1001,1)=date;
   fclose(fid);
end

end

