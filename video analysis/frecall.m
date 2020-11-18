function [ratio,hit]= frecall(data,template,window)
%input: data and template are both vectors,use the second as
%baseline,window means the time window that tolerate time shifts,[a,b]
%output: ratio of elements in template, for which data recalled
% upgrade 2020.1.12 to compatible with data that organized by trials(each
% trial a cell)
if (~iscell(data))&&(~iscell(template))
    hit=zeros(1,length(template));
    for i=1:length(template)
        hit(i)=sum((data<template(i)+window(2)).*(data>template(i)-window(1)));
    end
    ratio=1-sum(hit==0)/length(hit);
elseif iscell(data)&&iscell(template)&&(length(data)==length(template))%here data and template should both be cells with same length
    hit=cell(1,length(data)); 
    for n_trial=1:length(data)
        if isnan(template{1,n_trial})
            hit{1,n_trial}=[];
        else
            hit{1,n_trial}=zeros(1,length(template{1,n_trial}));
            for i=1:length(template{1,n_trial})
                hit{1,n_trial}(i)=sum((data{1,n_trial}<template{1,n_trial}(i)+window(2)).*(data{1,n_trial}>template{1,n_trial}(i)-window(1)));
            end
        end
    end
    n_template=cellfun(@length , hit);
    n_template=sum(n_template);
    n_hit=cellfun(@(x) sum(x~=0), hit);
    n_hit=sum(n_hit);
    ratio=n_hit/n_template;
else
    warning('Two input variables are not same type');
end
end