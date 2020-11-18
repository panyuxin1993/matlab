function [answerNum] = fGetAnsNum(sessionresult)
choice=cellfun(@(x) x.Action_choice, sessionresult);
answerNum=[sum(choice==0) sum(choice==1) ];
end