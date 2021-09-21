%%
clear
Path_Input=[];
Path_Output=[];

Path_Input{1}='E:\20210917_MAC\20210917QAQ_1';
Path_Output{1}='E:\20210917_MAC\20210917QAQ_1';


%%
clearvars -except ROImaskAll Path_Input Path_Output
for n_tiff=1:length(Path_Input)

close all
clc

cd(Path_Input{n_tiff});
files=dir('*.tif');
n_frame=9;
    warning('off');
%%
for n_file=1:length(files)
    n_FileName=files(n_file).name;
tsStack=TIFFStack(n_FileName);
%%
hh=figure;
im=tsStack(:,:,n_frame);
Im=im(:);
clim=[];
clim(1)=prctile(Im,10);
clim(2)=prctile(Im,90);

imagesc(im,clim);
colormap gray
%%
% plot ROI
if ~exist('ROImaskAll')
    inputchar = 'y';
    ROImaskAll = {};
    k = 1;
    while ~strcmpi(inputchar,'n')
        hROI = imfreehand;
        pos = hROI.getPosition;
        BW = createMask(hROI);
        ROImaskAll{k} = BW;
        k = k + 1;
        inputchar = input('Would you like to add a new ROI or delete(d) the current ROI','s');
        if strcmpi(inputchar,'d')
            k = k - 1;
        end
    end
end

    %%
    % extract ROI data
    nFrame=size(tsStack,3);
    nROIs = length(ROImaskAll);
    ROIdata = zeros(nFrame,nROIs);
    for cROI = 1 : nROIs
        cMask = ROImaskAll{cROI};
        for n = 1 : nFrame
            cData = tsStack(:,:,n);
            WithinData = cData(cMask);
            ROIdata(n,cROI) = mean(WithinData);
        end
    end
    cd(Path_Output{n_tiff});
    save(n_FileName(1:(end-8)),'ROIdata','ROImaskAll');
    cd(Path_Input{n_tiff});

        %% sum all tif data into ROIdata_all
    if n_file==length(files)
        ROIdata_all=zeros(0);
        cd(Path_Output{n_tiff});
        mat_Files=dir(['*_MMStack_' '*.mat']);
        for n_mat=1:length(mat_Files)
            n_matFileName=mat_Files(n_mat).name;
            load(n_matFileName);
            ROIdata_all(length(ROIdata_all)+1:(length(ROIdata_all)+length(ROIdata)),:)=ROIdata;
        end
        matFileName_All=[mat_Files(1).name(1:(end-4)) '_All'];
        save(matFileName_All,'ROIdata_all','ROImaskAll');
    end
end
disp('ROI_read_example Done!')
msgbox({Path_Input{n_tiff}(end-20:end) 'Done'})

warning('on');
close
end
%%
% plot_FP_RF(Path_Output);