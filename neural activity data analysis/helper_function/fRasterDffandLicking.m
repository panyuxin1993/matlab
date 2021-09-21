function [] = fRasterDffandLicking(x_lim,axActivity,neuralActivity,behEvent_aligned,selectedTrialInd,behEventSort,ColLimit,varargin)
%FRASTERDFFANDLICKING plot a block licking raster and corresponding licking
%raster, need to call function fRasterDff
%Input-
%   neuralActivity,behEvent_aligned,selectedTrialInd,behEventSort,ColLimit-
%       see fRasterDff
%   x_lim- indicating the x limit, which is useful especially when licking
%   raster were needed in accompany with activities raster
%   axActivity- ax to plot the activity raster
%   axLicking- ax to plot licking
%   licking_aligned- licking aligned to some behavior events
if isempty(varargin)
    subplot(axActivity);
    fRasterDff(neuralActivity,behEvent_aligned,selectedTrialInd,behEventSort,ColLimit);
    xlim(x_lim);
else
    axLicking=varargin{1};
    licking_aligned=varargin{2};
    subplot(axActivity);
    I=fRasterDff(neuralActivity,behEvent_aligned,selectedTrialInd,behEventSort,ColLimit);
    xlim(x_lim);
    temp=find(selectedTrialInd);
    leftLick=cell(1,length(temp));
    rightLick=cell(1,length(temp));
    for i=1:length(temp)
        leftLick{i}= licking_aligned.leftLick{temp(i)};
        rightLick{i}=licking_aligned.rightLick{temp(i)};
    end
    %raster plot of lickings
    subplot(axLicking);
    for i=1:length(leftLick)%each trial as a row
        for jl=1:length(leftLick{I(i)})
            line([leftLick{I(i)}(jl),leftLick{I(i)}(jl)],[i-0.5,i+0.5],'color','b','linewidth',1);%left lick
            hold on;
        end
        for jr=1:length(rightLick{I(i)})
            line([rightLick{I(i)}(jr),rightLick{I(i)}(jr)],[i-0.5,i+0.5],'color','r','linewidth',1);%right lick
            hold on;
        end
    end
    set(gca,'ytick',[]);%licking raster have corresponding color plot label, so not need
    set(gca,'xtick',[]);
    xlim(x_lim);
    if sum(selectedTrialInd)>0
        ylim([0.5 sum(selectedTrialInd)+0.5]);
    end
    
end
end

