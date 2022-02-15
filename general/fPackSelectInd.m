function [packed_ind] = fPackSelectInd(ind_all,n_show)
%fPackSelectInd using n_show to pack index to cell array
%   Detailed explanation goes here
n_cell=ceil(sum(ind_all)/n_show);
packed_ind=cell(1,n_cell);
ind_start=1;
for i=1:n_cell
    packed_ind{i}=ind_all;
    if ind_start<=length(ind_all) 
        if ind_start>1
            packed_ind{i}(1:ind_start)=0;
        end
        for ind_ind=ind_start:length(ind_all)
            if sum(ind_all(ind_start:ind_ind))>=n_show
                break;
            end
        end
        packed_ind{i}(ind_ind+1:end)=0;
        ind_start=ind_ind+1;
    end

end
end

