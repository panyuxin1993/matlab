function [str]=plabelsymbol(pvalue)
%PLABELSYMBOL transfer p value to * label
%   input pvalue, a number
%   output str, a string
if pvalue<0.05 && pvalue>=0.01
    str=' *';
elseif pvalue<0.01 && pvalue>=0.001
    str=' **';
elseif pvalue<0.001
    str=' ***';
elseif pvalue>=0.05 && pvalue<0.01
    pvalue=round(pvalue,2);%保留两位即可
    str=strcat('p=',num2str(pvalue));
elseif isnan(pvalue)
    str='';
else
    pvalue=round(pvalue,1);%保留1位即可
    str=strcat('p=',num2str(pvalue));
end
end



