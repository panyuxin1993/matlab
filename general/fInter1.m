function [cellout] = fInter1(cellin,varargin)
%FINTER1 before using cell2mat, to ensure that each cell have same size,
%this function help to re-estimate values based on interp1
%   Detailed explanation goes here
flag=cellfun(@(x) ~isvector(x), cellin,'UniformOutput',true);
if sum(sum(flag))>0
    warning('input cell array contain not only vector data, so not process');
    cellout=cellin;
    return;
end
if isempty(varargin)
    lengthmax=cellfun(@(x) max(size(x)), cellin,'UniformOutput',false);
else
    lengthmax=varargin{1};
end
matlength=cell2mat(lengthmax);
xx=cellfun(@(x) (1:x)/x,lengthmax,'UniformOutput',false);
xq=(1:max(matlength,[],'all'))/max(matlength,[],'all');
cellout=cellfun(@(x,y) interp1(x,y,xq), xx,cellin,'UniformOutput',false);
end

