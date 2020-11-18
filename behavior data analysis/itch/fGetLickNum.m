function [lickNum]= fGetLickNum(sessionresult)
leftLick=cellfun(@(x) x.Action_numLickLeft,sessionresult);
rightLick=cellfun(@(x) x.Action_numLickRight,sessionresult);
lickNum=[sum(leftLick) sum(rightLick)];
end

