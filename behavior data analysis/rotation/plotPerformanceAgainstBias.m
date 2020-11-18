dataChoice=fFindChoice('E:\�ҵ��ļ�\�Ϻ�����\��ת\������\ѵ�����\pyx03');
performance=fCorrectRate(dataChoice);
date=performance(:,11);
answerPeriod=find(performance(:,12)==300);
%c1=[0:255/(size(performance,1)-1):255;0:255/(size(performance,1)-1):255;0:255/(size(performance,1)-1):255]';
plot(performance(:,10),performance(:,2),'k','LineWidth',2); %����ȷ�ʺ���bias
hold on;
leftCorrectRateNoMiss=performance(:,5)./(performance(:,5)+performance(:,6));
rightCorrectRateNoMiss=performance(:,8)./(performance(:,8)+performance(:,9));

% normalized bias
% biasNoMiss=(leftCorrectRateNoMiss-rightCorrectRateNoMiss)./(leftCorrectRateNoMiss+rightCorrectRateNoMiss);
% simple bias
biasNoMiss=(leftCorrectRateNoMiss-rightCorrectRateNoMiss);

%c2=zeros(size(performance,1),3);
%c2(:,1)=(0:255/(size(performance,1)-1):255)';
plot(biasNoMiss,performance(:,2)./(performance(:,2)+performance(:,3)),'r','LineWidth',2);%����miss����ȷ�ʺ�bias
hold on;
scatter(performance(:,10),performance(:,2),40,'k');%���Ƶ�ͼ
hold on;
scatter(biasNoMiss,performance(:,2)./(performance(:,2)+performance(:,3)),40,'r');%���Ƶ�ͼ
hold on;
xlabel('Bias score','FontName','Arial','FontSize',14);
ylabel('Performance','FontName','Arial','FontSize',14);
title('Performance versus bias of mouse','FontName','Arial','FontWeight','Bold','FontSize',16)
set(gca,'FontName','Arial','FontSize',14)%����������̶��������ƣ���С
xmin=min(min(performance(:,2)),min(biasNoMiss));
xmax=max(max(performance(:,2)),max(biasNoMiss));
plot([0,0],[0,1],'-k','LineWidth',1);
hold on;
plot([xmin,xmax],[0.8,0.8],'--k','LineWidth',1);
%legend('Performance vs. bias','Performance vs. bias without miss','Location', 'NorthWestOutside');
