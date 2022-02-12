function [ segflagcell,segrawcell] = fSegment2Cell( inflag,inraw )
%FSEGMENT2CELL change a vector into cell array based on flag
%Input-
%   inflag- logical vector, indicating each element in/not in output data
%   inraw- vector, raw data that will be segmented, same length as inflag
%Output-
%   segrawcell- n-by-1 cell array, each store the data of a segment
%   segflagcell- n-by-1 cell array, each store the flag of a segment
% tic;
nelement=length(inflag);
if nelement~=length(inraw)
    warning('input not same length');
end
if sum(inflag)==0 %no data, no segment
    segflagcell=[];
    segrawcell=[];
else
    nseg=1;%at least one segment
    segflagcell=cell(nseg,1);
    segflagcell{1}=false(1,nelement);
    for i=1:nelement-1
        if inflag(i) && inflag(i+1)
            segflagcell{nseg}(i)=true;
        elseif inflag(i) && ~inflag(i+1)
            segflagcell{nseg}(i)=true;
            nseg=nseg+1;
        elseif ~inflag(i) && inflag(i+1)
            segflagcell{nseg}=false(1,nelement);
        end
    end
    segflagcell=reshape(segflagcell,[],1);
%     [inrawcell{1:length(segflagcell)}]=deal(inraw);
%     inrawcell=reshape(inrawcell,[],1);
    segrawcell=cellfun(@(x) inraw(x),segflagcell,'UniformOutput',false);
end
% toc;
end

