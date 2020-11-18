clear;
date='20200710';
animal='1';
date(date=='_')='';%将下划线取消掉
filename=strcat('D:\xulab\behavior\practice\DWQ\',animal,'\',date,'_#',animal,'.mat');
date(date=='_')='';%将下划线取消掉  
load(filename); 
pL=cellfun(@(x) x.leftProb,SessionSettings);
choiceL=cellfun(@(x) x.Action_choice==0, SessionResults);
choiceR=cellfun(@(x) x.Action_choice==1, SessionResults);
choiceMiss=cellfun(@(x) x.Action_choice==2, SessionResults);
missrate=sum(choiceMiss)/length(choiceMiss);
nleftlick=cellfun(@(x) x.Action_numLickLeft,SessionResults);
nrightlick=cellfun(@(x) x.Action_numLickRight,SessionResults);
for i=1:length(choiceL)-20
    pChoiceR(i)=sum(choiceR(i:i+20))/sum(choiceL(i:i+20)+choiceR(i:i+20));
    pLickR(i)=sum(nrightlick(i:i+20))/(sum(nrightlick(i:i+20))+sum(nleftlick(i:i+20)));
end
pChoiceR_value=sum(choiceR)/sum(choiceL+choiceR);
pLickR_value=sum(nrightlick)/(sum(nrightlick)+sum(nleftlick));
figure;
set(gcf,'position',[0,0,900,300]);
subplot(1,3,1);
plot(pChoiceR);hold on;
titlestr={'Probability of Right Choice when ';strcat('Probability of Right Reward is ',num2str(100-pL(1)))}';
title(titlestr);
text(0,0.9,['p(right choice)=',num2str(pChoiceR_value)],'Unit','Normalized','FontSize',14);
plot([1,length(pChoiceR)],[pChoiceR_value,pChoiceR_value],'k--');
xlabel('Trials');
ylim([0,1]);

subplot(1,3,2);
plot(pLickR);hold on;
title('Probability of Right Lick');
text(0,0.9,['p(right lick)=',num2str(pLickR_value)],'Unit','Normalized','FontSize',14);
plot([1,length(pLickR)],[pLickR_value,pLickR_value],'k--');
xlabel('Trials');
ylim([0,1]);

subplot(1,3,3);
plot(pChoiceR./pLickR);hold on;
title('Ratio of Right Choice and Right Lick');
xlabel('Trials');
text(0,0.9,['Miss rate=',num2str(missrate)],'Unit','Normalized','FontSize',14);
