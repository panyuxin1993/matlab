f_cat  = [];
ind_1stFr(1) = 1;
ntr = length(SavedCaTrials.f_raw);
for i = 1:ntr
    f_cat = [f_cat SavedCaTrials.f_raw{i}];
    ind_1stFr(i+1) = size(f_cat,2) + 1;
end
ind_1stFr(i+1) = [];

%%
roiNo = 3;
totalFr = size(f_cat,2);
frT = SavedCaTrials.FrameTime;
span = round(40*1000/frT/2);
ind_x = 0;
x = [];
for i = 1:totalFr
%     ind_x = ind_x + 1;
    ind1 = i- span;
    ind2 = i+span;
    if ind1 <= 0
        ind1 = 1;
    end
    if ind2 > totalFr
        ind2 = totalFr;
    end 
    x(i) = prctile(f_cat(roiNo,ind1:ind2),5);
end
%%
f = f_cat(roiNo,:) - x + mean(x);
% figure; histogram(f_cat_sub(roiNo,:),100)
[N,edges,bin] = histcounts(f,100);
f0 = edges(N == max(N));
dff = (f - f0)/f0*100;
figure; plot(dff)