function ind_chosen = fInd_chosen(chosen_cellarray, defined_cellarray)
%FIND_CHOSEN change the chosen_cellarray from cell array to matrix of same
%size;
%Input-
%   defined_cellarray- vector of cells, each defined the index of content
%   if happened in chosen_cellarray
ind_chosen=zeros(size(chosen_cellarray));
for i=1:length(defined_cellarray)
    temp=cellfun(@(x) strcmp(x,defined_cellarray{i}), chosen_cellarray,'UniformOutput',1);
    ind_chosen=ind_chosen+i*temp;
end

end

