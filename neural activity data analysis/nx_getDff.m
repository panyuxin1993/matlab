

cd 'F:\2P\pyx272_20200113\im_data_reg\result_save';
CurrFolder = pwd;
rootname='pyx272_20200113';
savefolder=rootname;%'1-200trials';%'segNP';
dirmat=strcat(CurrFolder,'\*.mat');
dirs=dir(dirmat);
dircell=struct2cell(dirs);
filenames=dircell(1,:);
file_imaging=cellfun(@(x) contains(x,'Ca'), filenames);
i_file_imaging=find(file_imaging);
load(filenames{1,i_file_imaging});%load imaging data

if ~exist(savefolder)
    mkdir(savefolder);
end

f_cat  = [];
ind_1stFr(1) = 1;
ind_tr_1=1;
ntr = length(SavedCaTrials.f_raw); %cancate all trials together, or  just set the number of trials to use
for i = ind_tr_1:ntr
    f_cat = [f_cat SavedCaTrials.f_raw{i}];
    ind_1stFr(i+1) = size(f_cat,2) + 1;
end
ind_1stFr(i+1) = [];
dff=zeros(size(f_cat));
for roiNo = 1:size(SavedCaTrials.f_raw{1},1) %SavedCaTrials.nROIs may be not true
%%
%substract changing baseline
    totalFr = size(f_cat,2);
    frT = SavedCaTrials.FrameTime;
    span = round(20*1000/frT/2); %every 40s, get a mean and substrate to compensate for baseline change
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
    f = f_cat(roiNo,:) - x + mean(x);
    % figure; histogram(f_cat_sub(roiNo,:),100)   
    %%
    %get f mode as f0
    [N,edges,bin] = histcounts(f,100);
    f0 = edges(N == max(N));
    if length(f0)==2
        f0=mean(f0);
    end
    %get dff
    dff(roiNo,:) = (f - f0)/f0;
    %
    B=figure;
    set(gcf, 'position', [0 0 1500 500]);
    plot(f_cat(roiNo,:),'k-');
    hold on;
    plot(x,'c-');
    plot(f,'g-');
    plot(dff(roiNo,:),'b-');
    legend('f cat','moving f mean','f baseline correction','dff');
    set(gca,'FontName','Arial','FontSize',14);
    saveas(B,[CurrFolder,'\',rootname,'\',rootname,'ROI-',num2str(roiNo),'-cat_f.jpg'],'jpg');
    close all;
    %}
end
save('dff','dff');