function [ meandata, sigmat,pmat] = fPlotMeanSig2D_SVM(ax,data,p_data,ts,pSig,strbar)
%FPLOTMEANSIG2D plot mean in color and draw contour lines ot represent
%significant range; Currently mainly used for showing SVM decoding accuracy
%Input-
%   data- r-by-m-by-n matrix, m rows, n columns and r repeats, m=n
%   pdata- m-by-n matrix, p value of each data point
%   pSig- p value that decide whether data was significant from baseline
%Output-
%   ax- axis to plot
%   meandata- m-by-n matrix store the mean of each point
%   ts- time stamps 
%   sigmat- m-by-n matrix made of 0-1, which indicates which position is
%   significant
%   pmat- m-by-n matrix with each point indicates the p value of that
%   posision

meandata = nanmean(data,1);
sigmat= (p_data<pSig);
pmat=p_data;
%plot mean
axes(ax);
clims=[0,1];
meandata=squeeze(meandata);
meandata=flipud(meandata);%so y-axis increase from down to up
imagesc(ts,ts,meandata,clims);
ts_show=get(gca,'XTick');
ind_show=zeros(size(ts_show));
for i=1:length(ts_show)
    ind_show(i)=find(abs(ts_show(i)-ts)==min(abs(ts_show(i)-ts)));
end

ind_show_y=length(ts)+1-ind_show;
set(gca,'YTick',ts(fliplr(ind_show_y)),'YTickLabel',fliplr(ts_show));
hold on;
if strcmp(strbar,'colorbar')
    colormap(jet);
    chandle=colorbar;
    set(chandle,'Position',[0.91,0.1,0.02,0.8]);
    chandle.Label.String= 'Decoding accuracy';
end
%plot significance contour
contour(ts,fliplr(ts),sigmat,1,'--','LineColor','w');
ind0=find(min(abs(ts))==abs(ts));
if abs(ts(ind0))<0.1%which mean is not epoch SVM, zero point exist
    ind0y=length(ts)+1-ind0(1);
    plot([0,0],[ts(1),ts(end)],'k--');
    plot([ts(1),ts(end)],[ts(ind0y),ts(ind0y)],'k--');
end
set(gca,'FontSize',12);

end


