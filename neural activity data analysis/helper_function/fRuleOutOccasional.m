function [out]= fRuleOutOccasional(in,threshold)
%rule out data points with short duration(<threshold) 
%in- a vector, nan/num
out=in;
indsig=~isnan(in);
diffSeries=diff(indsig);
diffSeries=reshape(diffSeries,1,[]);%to a row vector
diffSeries=[0,diffSeries];
tchange=find(diffSeries~=0);
if isempty(tchange)
    onset=[];
    n_change=0;
else
    if diffSeries(tchange(1))<0
        tchange=[1,tchange];
    end
    if diffSeries(tchange(end))>0
        tchange=[tchange,length(in)];
    end
end
%now, n-pairs, each up-down
ind_up=1:2:length(tchange);
ind_down=2:2:length(tchange);
duration=tchange(ind_down)-tchange(ind_up);
indRuleOut=find(duration<threshold);
for i=1:length(indRuleOut)
    out(tchange(ind_up(indRuleOut(i))):tchange(ind_down(indRuleOut(i))))=nan;
end
end


