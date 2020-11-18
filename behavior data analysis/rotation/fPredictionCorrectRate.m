function [ predictionRate ] = fPredictionCorrectRate( prediction )
%FPREDICTIONCORRECTRATE Summary of this function goes here
%   Detailed explanation goes here
%   计算小鼠prediction正确情况，也许可以反应其机智程度
%   
%14行分别记录时间，总预测正确率，错误率，左边正确率，错误率，右边正确率错误率,...
%  error_stay,剔除miss的正确率,上一个trial错误之后的仍然选该侧率(包含或者不包含miss）（由于错误，error_stay=1,故下一次还是同一侧，只能选与本次错误选择的对侧方能做对
%   错误后的正确率，错误率，错误后不选该侧率
predictionRate=zeros(14,size(prediction,1));     

for i=1:size(prediction,1)
    predictionRate(1,i)=prediction(i,1001,1);     %1行记录时间
    predictionRate(8,i)=prediction(i,1001,2);     %记录error_stay
    n_trial=find(prediction(i,:,1)==0);
    if isempty(n_trial)
        totalTrial=1000;
    else
        totalTrial=n_trial(1)-1;            %求出总trial数
    end
    leftTrial=sum(find(prediction(i,:,1)==28000));%求出左边trial数
    rightTrial=sum(find(prediction(i,:,1)==7000));%求出右边trial数
    predictionLeft=(prediction(i,:,3)<prediction(i,:,4)).*(0<prediction(i,:,3))+...  %当0<a<b时取a
                    (prediction(i,:,4)==0).*(0<prediction(i,:,3));              %当a>0且b=0时
    predictionRight=(prediction(i,:,3)>prediction(i,:,4)).*(0<prediction(i,:,4))+... 
                    (prediction(i,:,3)==0).*(0<prediction(i,:,4));
    predictionNone=1-predictionLeft-predictionRight;
    realLeft=(prediction(i,:,1)==28000);
    realRight=(prediction(i,:,1)==7000);
    predictionLeftCorrect=predictionLeft.*realLeft;
    predictionRightCorrect=predictionRight.*realRight;
    predictionLeftError=predictionRight.*realLeft;
    predictionRightError=predictionLeft.*realRight;
    WrongPredictionTrialNextLeft=[0,predictionLeftError];   %找到做错的下一个trial
    WrongPredictionTrialNextRight=[0,predictionRightError];    
    predictionSameSideAfterWrong=WrongPredictionTrialNextLeft.*[predictionLeft,0]+WrongPredictionTrialNextRight.*[predictionRight,0];%保证.*的维数相等
    predictionAltSideAfterWrong=WrongPredictionTrialNextLeft.*[predictionRight,0]+WrongPredictionTrialNextRight.*[predictionLeft,0];
    %*************注意，这里的correct表示的是预测对而非该trial选对***************
    correctAfterWrong=[(predictionLeftCorrect+predictionRightCorrect),0].*(WrongPredictionTrialNextLeft+WrongPredictionTrialNextRight);
    wrongAfterWrong=[(predictionLeftError+predictionRightError),0].*(WrongPredictionTrialNextLeft+WrongPredictionTrialNextRight);
    correctRateAfterWrong=sum(correctAfterWrong)/sum(WrongPredictionTrialNextLeft+WrongPredictionTrialNextRight);
    predictionTrial=(prediction(i,:,3)>0)+(prediction(i,:,4)>0)-(prediction(i,:,3)>0).*(prediction(i,:,4)>0);%进行预测的trial
    n_predictionTrial=sum(predictionTrial);             %进行预测的trial数量
    predictionCorrectRate=sum(predictionLeftCorrect+predictionRightCorrect)/totalTrial;
    predictionErrorRate=sum(predictionLeftError+predictionRightError)/totalTrial;
    predictionCorrectWithoutNull=sum(predictionLeftCorrect+predictionRightCorrect)/n_predictionTrial;
    predictionSameSideRateAfterWrong=sum(predictionSameSideAfterWrong)/sum(predictionLeftError+predictionRightError);     %总正确率，只要是在错误之后的trial
    predictionAltSideRateAfterWrong=sum(predictionAltSideAfterWrong)/sum(predictionLeftError+predictionRightError); 
    predictionSameSideRateAfterWrongWithoutNuLL=sum(predictionSameSideAfterWrong)/sum(predictionSameSideAfterWrong+predictionAltSideAfterWrong); %不包括miss的正确率
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

