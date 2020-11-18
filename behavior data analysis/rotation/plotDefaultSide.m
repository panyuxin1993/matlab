defaultSide=fDefaultSide( 'E:\�ҵ��ļ�\�Ϻ�����\��ת\������\ѵ�����\pyx04' );

%�ⲿ����Ϊ��ͬһ��ͼ�ϱȽ���ȷ�ʺ�default side֮��Ĺ�ϵ
% choice=fFindChoice( 'E:\�ҵ��ļ�\�Ϻ�����\��ת\������\ѵ�����\pyx09' );
% performance=fCorrectRate( choice );
% plot(performance(:,2),'r','LineWidth',2);   %����ȷ�ʣ���miss
% hold on;

date=defaultSide(1,:);
leftLick=defaultSide(2,:);
rightLick=defaultSide(3,:);
leftDefault=defaultSide(6,:);
rightDefault=defaultSide(7,:);
errorStay=find(defaultSide(4,:)==1);
bias=(leftLick-rightLick)./(leftLick+rightLick);
biasByTrial=(leftDefault-rightDefault)./(leftDefault+rightDefault);
answerPeriod=find(defaultSide(5,:)==300);
plot(bias,'k','LineWidth',2); 
hold on;
plot(biasByTrial,'r','LineWidth',2);
hold on;
plot([errorStay(1),errorStay(1)],[0,1],'--b','LineWidth',1);
plot([errorStay(end),errorStay(end)],[0,1],'--b','LineWidth',1);
plot([0,size(defaultSide,2)],[0,0],'-k','LineWidth',1);
xlabel('Training day','FontName','Arial','FontSize',14);
ylabel('Default side score','FontName','Arial','FontSize',14);
title('Default side score of mouse','FontName','Arial','FontWeight','Bold','FontSize',16)
set(gca,'FontName','Arial','FontSize',14)%����������̶��������ƣ���С
plot([answerPeriod(1),answerPeriod(1)],[0,1],'--k','LineWidth',1);
plot([answerPeriod(end),answerPeriod(end)],[0,1],'--k','LineWidth',1);
%legend('Default side score by lick','Default side score by trial','Location', 'NorthWestOutside');
legend('Default side score by lick','Performance','Location', 'NorthWestOutside');
