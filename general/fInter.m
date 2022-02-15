function [cellout] = fInter(cellin,varargin)
%FINTER2 before using cell2mat, to ensure that each cell have same size,
%this function help to re-estimate values based on inter1 and interp2
%   Detailed explanation goes here
flag_matrix=cellfun(@(x) ismatrix(x), cellin,'UniformOutput',true);
flag=cellfun(@(x) isvector(x), cellin,'UniformOutput',true);
if sum(sum(flag_matrix))==sum(sum(ones(size(cellin))))
    lengthmax=cellfun(@(x) max(size(x,2)), cellin,'UniformOutput',false);
    if isempty(varargin)
        matlength=cell2mat(lengthmax);
        lengthmaxmax=max(matlength,[],'all');
    else
        lengthmaxmax=varargin{1};
    end
    
    xx=cellfun(@(x) (1:x)/x,lengthmax,'UniformOutput',false);
    xq=(1:lengthmaxmax)/lengthmaxmax;
    yy=1:size(cellin{1,1},1);
    [Xq,Yq]=meshgrid(xq,yy);
    cellout=cellfun(@(x,y) interp2(x,yy,y,Xq,Yq), xx,cellin,'UniformOutput',false);
elseif sum(sum(flag))==sum(sum(ones(size(cellin))))
    lengthmax=cellfun(@(x) max(size(x)), cellin,'UniformOutput',false);
    if isempty(varargin)
        matlength=cell2mat(lengthmax);
        lengthmaxmax=max(matlength,[],'all');
    else
        lengthmaxmax=varargin{1};
    end
    
    xx=cellfun(@(x) (1:x)/x,lengthmax,'UniformOutput',false);
    xq=(1:lengthmaxmax)/lengthmaxmax;
    cellout=cellfun(@(x,y) interp1(x,y,xq), xx,cellin,'UniformOutput',false);
else
    warning('input cell array contain not only vector data or matrix data, so not process');
    cellout=cellin;
    return;
end

end