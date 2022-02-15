function [ meandata, sigmat,pmat] = fPlotMeanSig2D(ax,data,ts,pSig, baseline)
%FPLOTMEANSIG2D plot mean in color and draw contour lines ot represent
%significant range; Currently mainly used for showing SVM decoding accuracy
%Input-
%   data- m-by-n-by-r matrix, m rows, n columns and r repeats, m=n
%   pSig- p value that decide whether data was significant from baseline
%   baseline- double, compared with data
%Output-
%   ax- axis to plot
%   meandata- m-by-n matrix store the mean of each point
%   ts- time stamps 
%   sigmat- m-by-n matrix made of 0-1, which indicates which position is
%   significant
%   pmat- m-by-n matrix with each point indicates the p value of that
%   posision

meandata = nanmean(data,3);
pmat=ones(size(data,1),size(data,2));
for i=1:size(data,1)
    for j=1:size(data,2)
        test_mean=bootstrp(2000,@mean,data(i,j,:));
        temp_p=nansum(test_mean>baseline)/length(test_mean);
        pmat(i,j)=1-2*abs(temp_p-0.5);
    end
end
sigmat= (pmat<pSig);
%plot mean

axes(ax);
clims=[0,100];
imagesc(ts,ts,meandata,clims);
hold on;
colormap(jet);
chandle=colorbar;
set(chandle,'Position',[0.9,0.1,0.05,0.8]);
chandle.Label.String= 'Decoding accuracy';
%plot significance contour
contour(ts,ts,sigmat,1,'--','LineColor','w');

end

