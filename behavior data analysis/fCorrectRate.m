function [ performance ] = fCorrectRate( dataChoice )
%FCORRECTRATE Summary of this function goes here
%   Detailed explanation goes here
% 输入选择矩阵数据，计算某日的正确率（包括两侧），相对偏好bias.
n_day=size(dataChoice,1);     %获得训练天数
%每一列分别表示总（7000Hz侧/28000Hz侧）trial数量，正确率（不包含miss），miss_rate，violation
%rate(不包括miss）
%以及bias,日期
performance=zeros(n_day,14);    
for i=1:n_day
    n_trial=find(dataChoice(i,:,1)==0); 
    if isempty(n_trial)%做完1000trial时n_tial=null;
        performance(i,1)=1000;
    else
        performance(i,1)=n_trial(1)-1;  %总trial数
    end
%     m_7k=double(dataChoice(i,:,1)==7000);  %同时要把逻辑矩阵转化为数值矩阵
%     m_28k=double(dataChoice(i,:,1)==28000);
    m_20=double((dataChoice(i,:,1)<50).*(dataChoice(i,:,1)>0));%行向量用50做标准的话，就能使得probe也进行计算。注意防止出现0为没做的trial
    m_125=double(dataChoice(i,:,1)>50);
    m_correct=double(dataChoice(i,:,2)==1);
    m_error=double(dataChoice(i,:,2)==2);
    m_miss=double(dataChoice(i,:,2)==3);
    m_violation=double(dataChoice(i,:,2)==4);
    m_nomiss=m_correct+m_error+m_violation;
    m_total=m_nomiss+m_miss;
    performance(i,2)=sum(m_correct)/sum(m_correct+m_error);   %正确率
    performance(i,3)=sum(m_miss)/sum(m_total);     %miss_rate
    performance(i,4)=sum(m_violation)/sum(m_nomiss);   %不包括miss的violation
    performance(i,5)=sum(m_20);     %总20数    
    performance(i,6)=sum(m_correct*m_20')/((m_correct+m_error)*m_20');      %20正确率, .'才是转置
    performance(i,7)=sum(m_miss*m_20')/performance(i,5);        %20miss_rate
    performance(i,8)=sum(m_violation*m_20')/sum(m_nomiss*m_20');   %不包括miss的violation
    performance(i,9)=sum(m_125);    %总125数
    performance(i,10)=sum(m_correct*m_125')/((m_correct+m_error)*m_125');     %125正确率
    performance(i,11)=sum(m_miss*m_125')/performance(i,9);       %125miss_rate
    performance(i,12)=sum(m_violation*m_125')/sum(m_nomiss*m_125');   %不包括miss的violation
%     performance(i,4)=sum(m_7k);     %总7k数
%     performance(i,7)=sum(m_28k);    %总28k数
%     performance(i,5)=sum(m_correct*m_7k')/performance(i,4);      %7k正确率, .'才是转置
%     performance(i,6)=sum(m_error*m_7k')/performance(i,4);        %7k错误率
%     performance(i,8)=sum(m_correct*m_28k')/performance(i,7);     %28k正确率
%     performance(i,9)=sum(m_error*m_28k')/performance(i,7);       %28k错误率
    %归一化的bias 指数
    %performance(i,10)=(performance(i,5)-performance(i,8))/(performance(i,5)+performance(i,8));  %bias;
    %简单相减得到的bias
    performance(i,13)=(performance(i,5)-performance(i,8));
    performance(i,14)=dataChoice(i,1001,1);             %日期
end

end

