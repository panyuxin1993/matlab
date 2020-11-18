function [ predictionRate ] = fPredictionCorrectRate( prediction )
%FPREDICTIONCORRECTRATE Summary of this function goes here
%   Detailed explanation goes here
%   ����С��prediction��ȷ�����Ҳ����Է�Ӧ����ǳ̶�
%   
%14�зֱ��¼ʱ�䣬��Ԥ����ȷ�ʣ������ʣ������ȷ�ʣ������ʣ��ұ���ȷ�ʴ�����,...
%  error_stay,�޳�miss����ȷ��,��һ��trial����֮�����Ȼѡ�ò���(�������߲�����miss�������ڴ���error_stay=1,����һ�λ���ͬһ�ֻ࣬��ѡ�뱾�δ���ѡ��ĶԲ෽������
%   ��������ȷ�ʣ������ʣ������ѡ�ò���
predictionRate=zeros(14,size(prediction,1));     

for i=1:size(prediction,1)
    predictionRate(1,i)=prediction(i,1001,1);     %1�м�¼ʱ��
    predictionRate(8,i)=prediction(i,1001,2);     %��¼error_stay
    n_trial=find(prediction(i,:,1)==0);
    if isempty(n_trial)
        totalTrial=1000;
    else
        totalTrial=n_trial(1)-1;            %�����trial��
    end
    leftTrial=sum(find(prediction(i,:,1)==28000));%������trial��
    rightTrial=sum(find(prediction(i,:,1)==7000));%����ұ�trial��
    predictionLeft=(prediction(i,:,3)<prediction(i,:,4)).*(0<prediction(i,:,3))+...  %��0<a<bʱȡa
                    (prediction(i,:,4)==0).*(0<prediction(i,:,3));              %��a>0��b=0ʱ
    predictionRight=(prediction(i,:,3)>prediction(i,:,4)).*(0<prediction(i,:,4))+... 
                    (prediction(i,:,3)==0).*(0<prediction(i,:,4));
    predictionNone=1-predictionLeft-predictionRight;
    realLeft=(prediction(i,:,1)==28000);
    realRight=(prediction(i,:,1)==7000);
    predictionLeftCorrect=predictionLeft.*realLeft;
    predictionRightCorrect=predictionRight.*realRight;
    predictionLeftError=predictionRight.*realLeft;
    predictionRightError=predictionLeft.*realRight;
    WrongPredictionTrialNextLeft=[0,predictionLeftError];   %�ҵ��������һ��trial
    WrongPredictionTrialNextRight=[0,predictionRightError];    
    predictionSameSideAfterWrong=WrongPredictionTrialNextLeft.*[predictionLeft,0]+WrongPredictionTrialNextRight.*[predictionRight,0];%��֤.*��ά�����
    predictionAltSideAfterWrong=WrongPredictionTrialNextLeft.*[predictionRight,0]+WrongPredictionTrialNextRight.*[predictionLeft,0];
    %*************ע�⣬�����correct��ʾ����Ԥ��Զ��Ǹ�trialѡ��***************
    correctAfterWrong=[(predictionLeftCorrect+predictionRightCorrect),0].*(WrongPredictionTrialNextLeft+WrongPredictionTrialNextRight);
    wrongAfterWrong=[(predictionLeftError+predictionRightError),0].*(WrongPredictionTrialNextLeft+WrongPredictionTrialNextRight);
    correctRateAfterWrong=sum(correctAfterWrong)/sum(WrongPredictionTrialNextLeft+WrongPredictionTrialNextRight);
    predictionTrial=(prediction(i,:,3)>0)+(prediction(i,:,4)>0)-(prediction(i,:,3)>0).*(prediction(i,:,4)>0);%����Ԥ���trial
    n_predictionTrial=sum(predictionTrial);             %����Ԥ���trial����
    predictionCorrectRate=sum(predictionLeftCorrect+predictionRightCorrect)/totalTrial;
    predictionErrorRate=sum(predictionLeftError+predictionRightError)/totalTrial;
    predictionCorrectWithoutNull=sum(predictionLeftCorrect+predictionRightCorrect)/n_predictionTrial;
    predictionSameSideRateAfterWrong=sum(predictionSameSideAfterWrong)/sum(predictionLeftError+predictionRightError);     %����ȷ�ʣ�ֻҪ���ڴ���֮���trial
    predictionAltSideRateAfterWrong=sum(predictionAltSideAfterWrong)/sum(predictionLeftError+predictionRightError); 
    predictionSameSideRateAfterWrongWithoutNuLL=sum(predictionSameSideAfterWrong)/sum(predictionSameSideAfterWrong+predictionAltSideAfterWrong); %������miss����ȷ��
    predictionRate(2,i)=predictionCorrectRate;
    predictionRate(9,i)=predictionCorrectWithoutNull;
    predictionRate(10,i)=predictionSameSideRateAfterWrong;
    predictionRate(11,i)=predictionSameSideRateAfterWrongWithoutNuLL;
    predictionRate(12,i)=correctRateAfterWrong;
    predictionRate(13,i)=sum(wrongAfterWrong)/sum(predictionLeftError+predictionRightError);
    predictionRate(14,i)=predictionAltSideRateAfterWrong;
   % predictionRate(2:7,i)=[predictionCorrectRate; predictionErrorRate; predictionLeftCorrect/leftTrial;predictionLeftError/leftTrial;predictionRightCorrect/rightTrial;predictionRightError/rightTrial];
end

end

