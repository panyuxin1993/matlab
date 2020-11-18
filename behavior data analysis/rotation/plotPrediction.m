a=fFindPrediction( 'E:\我的文件\上海神经所\轮转\徐宁龙\训鼠分析\pyx09' );
b=fPredictionCorrectRate(a);
data=b(1,:);
correctRate=b(2,:);
errorStay=find(b(8,:)==1);
correctRateWithoutNull=b(9,:);
staySameSidePredictRateAfterWrong=b(10,:);
% staySameSidePredictRateAfterWrongWithoutNull=b(11,:);
% correctPredictRateAfterWrong=b(12,:);
% wrongPredictRateAfterWrong=b(13,:);
predictionAltSideRateAfterWrong=b(14,:);
plot(correctRateWithoutNull,'r','LineWidth',2);
hold on;
plot(correctRate,'k','LineWidth',2);
hold on;
plot(staySameSidePredictRateAfterWrong,'g','LineWidth',2);
hold on;
% plot(staySameSidePredictRateAfterWrongWithoutNull,'b','LineWidth',2);
% hold on;
% plot(correctPredictRateAfterWrong,'c','LineWidth',2);
% hold on;
% plot(wrongPredictRateAfterWrong,'y','LineWidth',2);
% hold on;
plot(predictionAltSideRateAfterWrong,'b','LineWidth',2);
xlabel('Training day','FontName','Arial','FontSize',14);
ylabel('Prediction correct rate','FontName','Arial','FontSize',14);
title('Prediction of mouse','FontName','Arial','FontWeight','Bold','FontSize',16)
set(gca,'FontName','Arial','FontSize',14)%设置坐标轴刻度字体名称，大小
plot([errorStay(1),errorStay(1)],[0,1],'--b','LineWidth',1);
hold on;
plot([errorStay(end),errorStay(end)],[0,1],'--b','LineWidth',1);
hold on;
plot([0,size(b,2)],[0.5,0.5],'--k','LineWidth',1);%指示随机猜测值
hold on;
% legend('Correct rate in prediction trial ','Correct rate in all trial',...
% 'Stay same side rate after Wrong Trial','Stay same side rate after Wrong Trial without miss','Correct rate after wrong',...
% 'Error rate after wrong','Alt side rate after wrong trial',...
% 'Location');
% legend('Correct rate in prediction trial ','Correct rate in all trial',...
% 'Stay same side rate after wrong trial','Alt side rate after wrong trial',...
% 'Location','BestOutside');