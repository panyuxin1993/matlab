function [fig_A,fig_B] = fActivityInspect(activities_data,activity_form,ind_tr_1,Data_extract,SavedCaTrials,frameNumTime,behEventAlign,masklick,i_selectivity,behEventSort,trial2include,trial2exclude,title_fig,trialTypeStr,chosen_result,chosen_stim,n_trial_show)
%UNTITLED see individual trial activities of single ROI, for different
%trial type, selected number of trials, etc.
%Input-
% activity_form={'dff','spkr','zscored-spkr'},indicating the activity form

ntr = length(SavedCaTrials.f_raw); %use all data, or  just set the number of trials to use
frT = SavedCaTrials.FrameTime;
% align to behavior event
nFrameEachTrial=cellfun(@(x) size(x,2),SavedCaTrials.f_raw);
ind_1stFrame=zeros(1,length(nFrameEachTrial));
ind_1stFrame(1)=1;
ind_1stFrame(2:end)=cumsum(nFrameEachTrial(1:end-1))+1;
ind_1stFrame=ind_1stFrame(ind_tr_1:ind_tr_1+ntr-1);%if number of trials unsed for analysis is not whole but part of trials
if isnan(trial2include)
    indTrial2include=fExcludeTrials(trial2exclude,ind_1stFrame,'logical');
    trial2include=[1,length(ind_1stFrame)];
elseif strcmp(trial2include,'all')
    trial2include=[1,length(ind_1stFrame)];
    indTrial2include=fIncludeTrials(trial2include,ind_1stFrame,'logical');
else
    indTrial2include=fIncludeTrials(trial2include,ind_1stFrame,'logical');
end
frameNum=double(round(frameNumTime*1000/frT));
[behEventFrameIndex,lickingFrameIndex] = fGetBehEventTime( Data_extract, ind_1stFrame, SavedCaTrials.FrameTime ,ind_tr_1);%get behavior event time

%   align to which event(string can be in {'stimOnset', 'go cue','first lick','first left lick','first right lick', 'answer','reward'},
if strcmp(behEventAlign,'delay onset') && strcmp(masklick,'yes')%by default, if align to stim onset, then see delay activity and mask all after go cue as nan
    [ activities_aligned, behEvent_aligned,licking_aligned ] = fAlignDelaySigal( activities_data, behEventFrameIndex,  frameNum );
else
    [ activities_aligned, behEvent_aligned, licking_aligned ] = fAlignSigBehEvent( activities_data, behEventFrameIndex,lickingFrameIndex,behEventAlign,frameNum );%decide which behavior event to align to
end

[trialType,rule] = fGetTrialType( Data_extract,[],i_selectivity,'matrix','left','divideCorErr');%decide trial type, 1d cor/err/miss/vio, each 2d one stimulus, 3d trials


%%
%plot color plot,4 column(correct/error/miss/violation), 7 row(for each
%stimlus)+ mean trace

ColLimit = prctile(activities_data,98);
results_str={'Correct','Error','Miss','Violation'};
titlestr=strcat(title_fig,results_str);
titlestr=strrep(titlestr,'_','\_');
choice_str={'ipsi','contra'};
ind_result = fInd_chosen(chosen_result, results_str);
ind_stim=fInd_chosen(chosen_stim, choice_str);
if size(trialType,2)==6
    color_mean_trace={[0 0 1],[0.2 0 0.8],[0.4 0 0.6],[0.6 0 0.4],[0.8 0 0.2],[1 0 0]};
elseif size(trialType,2)==2
    color_mean_trace={[0 0 1],[1 0 0]};
end
%     if strcmp(rule,'low click rate-right')
%         color_mean_trace=fliplr(color_mean_trace);
%     end
for nStim=1:length(chosen_stim) %for each stimulus
    for  nResult=1:length(chosen_result) %4 column(correct/error/miss/violation),companied with 4 lick raster
        selectedTrialInd=trialType(ind_result(nResult),ind_stim(nStim),:);
        selectedTrialInd=reshape(logical(squeeze(selectedTrialInd)),[],1);
        indTrial2include=reshape(indTrial2include,[],1);
        selectedTrialInd_all=logical(selectedTrialInd.*indTrial2include);
        if isempty(n_trial_show)
            n_fig_rep=1;
            n_trial_show_fig=sum(selectedTrialInd);
        else
            n_fig_rep=ceil(sum(selectedTrialInd)/n_trial_show);
            n_trial_show_fig=n_trial_show;
        end
        selectedTrialInd_all_cell=fPackSelectInd(selectedTrialInd_all,n_trial_show_fig);
        for i_fig=1:n_fig_rep
            if nStim==1 &&nResult==1
                fig_A(i_fig)=figure;
            else
                figure(fig_A(i_fig));
            end
            set(gcf, 'position', [0 100 300*length(chosen_result) 200*(length(chosen_stim)+1)]);%last row, mean trace
            sot=behEvent_aligned.stimOnset(selectedTrialInd_all_cell{i_fig});% stim onset, white
            flt=behEvent_aligned.lickFirst(selectedTrialInd_all_cell{i_fig});% first lick, black
            flt_l=behEvent_aligned.lickFirst_left(selectedTrialInd_all_cell{i_fig});%left lick, magenta dots
            flt_r=behEvent_aligned.lickFirst_right(selectedTrialInd_all_cell{i_fig});%right lick, red dots
            rwt=behEvent_aligned.rewTime(selectedTrialInd_all_cell{i_fig});%reward, cyan dots
            go=behEvent_aligned.go(selectedTrialInd_all_cell{i_fig});%go cue,white
            neuralActivity=activities_aligned(selectedTrialInd_all_cell{i_fig},:);
            temp=find(selectedTrialInd_all_cell{i_fig});
            clear lickTime;
            if strcmp(masklick,'no')
                for i=1:length(temp)
                    lickTime(i).lickleft= licking_aligned.leftLick{temp(i)};%unit, frame number
                    lickTime(i).lickright=licking_aligned.rightLick{temp(i)};
                end               
                main_licking_dir=Data_extract.Action_choice(temp(i));%this trial type, the mainly licking direction
            end
            %decide which behavior event to be sort
            if strcmp(behEventSort,'first lick')
                [B,I]=sort(flt);
            elseif strcmp(behEventSort,'reward')
                [B,I]=sort(rwt);
            elseif strcmp(behEventSort,'go cue')
                [B,I]=sort(go);
            end
            go=go(I);
            sot=sot(I);
            flt=flt(I);
            flt_l=flt_l(I);
            flt_r=flt_r(I);
            rwt=rwt(I);
            figure(fig_A(i_fig));
            
            subplot(length(chosen_stim)+1,2*length(chosen_result),2*nResult-1+2*length(chosen_result)*(nStim-1));
            hold on;
            imagesc(neuralActivity(I,:));
            x_lim=[0,size(neuralActivity,2)];%get(gca,'Xlim');
            plot(sot,[1:sum(selectedTrialInd_all_cell{i_fig})],'w.','markersize',8);
            hold on;
            plot(go,[1:sum(selectedTrialInd_all_cell{i_fig})],'w.','markersize',8);
            hold on;
            plot(flt,[1:sum(selectedTrialInd_all_cell{i_fig})],'k.','markersize',8);
            hold on;
            plot(flt_l,[1:sum(selectedTrialInd_all_cell{i_fig})],'m.','markersize',8);
            hold on;
            plot(flt_r,[1:sum(selectedTrialInd_all_cell{i_fig})],'r.','markersize',8);
            hold on;
            plot(rwt,[1:sum(selectedTrialInd_all_cell{i_fig})],'c.','markersize',8);
            hold on;
            set(gca,'clim',[0 ColLimit]);
            set(gca,'ytick',sum(selectedTrialInd_all_cell{i_fig}),'yticklabel',sum(selectedTrialInd_all_cell{i_fig}),'ydir','normal');
            if nStim==length(chosen_stim)
                set(gca,'xtick',[round(1000/frT*(frameNumTime(1)-floor(frameNumTime(1))))+1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
                %                 xlabel(['time(s) from ',behEventAlign],'FontName','Arial','FontSize',14);
            else
                set(gca,'xtick',[]);
            end
            ColBar = colorbar;
            set(ColBar,'position',[0.91 0.1100 0.02 0.1445]);
            if sum(selectedTrialInd_all_cell{i_fig})>0
                ylim([0.5 sum(selectedTrialInd_all_cell{i_fig})+0.5]);
            end
            title(titlestr{nResult});
            if nResult==1
                if strcmp(Data_extract.rule,'low click rate-right')
                    ylabel([num2str(Data_extract.Stimuli(size(trialType,2)+1-nStim))]);
                else
                    ylabel([num2str(Data_extract.Stimuli(nStim))]); %' clicks/s'
                end
            end
            xlim(x_lim);
            
            subplot(length(chosen_stim)+1,2*length(chosen_result),2*nResult+2*length(chosen_result)*(nStim-1));
            if exist('lickTime','var')
                for i=1:length(lickTime)%each trial as a row
                    for jl=1:length(lickTime(I(i)).lickleft)
                        line([lickTime(I(i)).lickleft(jl),lickTime(I(i)).lickleft(jl)],[i-0.5,i+0.5],'color','b','linewidth',1);%left lick
                        hold on;
                    end
                    for jr=1:length(lickTime(I(i)).lickright)
                        line([lickTime(I(i)).lickright(jr),lickTime(I(i)).lickright(jr)],[i-0.5,i+0.5],'color','r','linewidth',1);%right lick
                        hold on;
                    end
                end
                set(gca,'ytick',[]);%licking raster have corresponding color plot label, so not need
                if nStim==length(chosen_stim)
                    set(gca,'xtick',[round(1000/frT*(frameNumTime(1)-floor(frameNumTime(1))))+1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
                else
                    set(gca,'xtick',[]);
                end
                xlim(x_lim);
                if sum(selectedTrialInd_all_cell{i_fig})>0
                    ylim([0.5 sum(selectedTrialInd_all_cell{i_fig})+0.5]);
                end
            end
            
            %plot mean trace
            subplot(length(chosen_stim)+1,2*length(chosen_result),2*nResult-1+2*length(chosen_result)*length(chosen_stim));
            
            ts=double((-frameNum(1):1:frameNum(2))*frT/1000);
            curve_meanTrace(nStim)=fPlotMean_SE(ts,neuralActivity,color_mean_trace{nStim});
            
            %plot mean trace of licking raster
            if exist('lickTime')
                subplot(length(chosen_stim)+1,2*length(chosen_result),2*nResult+2*length(chosen_result)*length(chosen_stim));
                binSize=5;%frame number
                binStep=1;
                [ lickMat,lickRate ] = fLickRate( lickTime, binSize,sum(frameNum)+1,0, ones(length(lickTime),1), binStep);%r=unit ms, so need transform to frames
                switch main_licking_dir
                    case 0
                        i_lickrate=2;
                    case 1
                        i_lickrate=3;
                    otherwise
                        i_lickrate=1;
                end
                ts=double((-frameNum(1):binStep:frameNum(2))*frT/1000);
                indx=min(length(ts),size(lickRate,1));
                curve_lick_meanTrace(nStim)=plot(ts(1:indx),lickRate(1:indx,i_lickrate)/frT,'color',color_mean_trace{nStim});hold on;
            end
            
            %piled up activities
            if nStim==1 &&nResult==1
                fig_B(i_fig)=figure;
            else
                figure(fig_B(i_fig));
            end
            set(gcf, 'position', [300 100 300*length(chosen_result) 400*(length(chosen_stim))]);%last row, mean trace
            subplot(length(chosen_stim),length(chosen_result),nResult+length(chosen_result)*(nStim-1));

            range=[];
            f_data=neuralActivity(I,:);
            dur_thresh=round(300/frT);%300ms above std
            ind_baseline=1:size(f_data,2);
            if strcmp(behEventAlign,'stim onset')
                ind_baseline=(ind_baseline<round(1000/frT*(frameNumTime(1))));
            elseif strcmp(behEventAlign,'delay onset')
                ind_baseline=(ind_baseline<round(1000/frT*(frameNumTime(1)-0.5)));
            end
            for i_row=1:size(f_data,1)
                current_f=f_data(i_row,:);
                current_f=current_f-min(current_f)+sum(range);
                
                fPlotSigTrace(current_f,ind_baseline,dur_thresh);
                hold on;
                xlim(x_lim);
                range_sum_m1=sum(range);
                range=[range,max(current_f)-min(current_f)];
                if mod(i_row,10)==0
                    text(x_lim(end),sum(range),strcat('Trial-',num2str(i_row)));
                end
                if exist('lickTime','var')         
                    for jl=1:length(lickTime(I(i_row)).lickleft)
                        line([lickTime(I(i_row)).lickleft(jl),lickTime(I(i_row)).lickleft(jl)],[range_sum_m1,range_sum_m1+range(end)],'color','b','linewidth',1);%left lick
                        hold on;
                    end
                    for jr=1:length(lickTime(I(i_row)).lickright)
                        line([lickTime(I(i_row)).lickright(jr),lickTime(I(i_row)).lickright(jr)],[range_sum_m1,range_sum_m1+range(end)],'color','r','linewidth',1);%right lick
                        hold on;
                    end
                    %plot go cue using black lines
                    plot([go(I(i_row)),go(I(i_row))],[range_sum_m1,range_sum_m1+range(end)],'color','k','linewidth',1);%go cue
                end
            end
            plot([round(1000/frT*(frameNumTime(1))),round(1000/frT*(frameNumTime(1)))],[0,sum(range)],'k-');%aligned time point
            if strcmp(behEventAlign,'stim onset')
                plot([round(1000/frT*(frameNumTime(1)+0.5)),round(1000/frT*(frameNumTime(1)+0.5))],[0,sum(range)],'k-');%aligned time point
            elseif strcmp(behEventAlign,'delay onset')
                plot([round(1000/frT*(frameNumTime(1)-0.5)),round(1000/frT*(frameNumTime(1)-0.5))],[0,sum(range)],'k-');%aligned time point
            end
            %plot y scale bar and text
            if length(range)>=10
                scalebar_lengthy=fRound(prctile(range,95));
            else
                scalebar_lengthy=fRound(max(range)/2);
            end
            plot([x_lim(end),x_lim(end)],[0,scalebar_lengthy],'k-','LineWidth',2);
            if strcmp(activity_form,'dff')
                text(x_lim(end),scalebar_lengthy/2,['\it\DeltaF/F ','\rm',num2str(scalebar_lengthy)]);
            elseif strcmp(activity_form,'spkr')
                text(x_lim(end),scalebar_lengthy/2,[num2str(scalebar_lengthy),' spikes/s']);
            end
            set(gca,'Ylim',[0,sum(range)],'ytick',[],'yticklabel',[]);
            set(gca,'xtick',[round(1000/frT*(frameNumTime(1)-floor(frameNumTime(1))))+1:round(1000/frT):size(f_data,2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
            ylabel(chosen_stim{nStim});
            xlabel(['Time (s) from ', behEventAlign]);
            title(titlestr{nResult});
            box off;
        end
    end
end

%% label and general plots
% find common ylim
figure(fig_A(i_fig));
for nResult=length(chosen_result)
    subplot(length(chosen_stim)+1,2*length(chosen_result),2*nResult-1+2*length(chosen_result)*length(chosen_stim));
    flagclearyrange=0;
    if ~exist('yrange','var')
        flagclearyrange=1;
        if ~exist('y_lim','var')
            y_lim=get(gca,'Ylim');
        else
            y_lim2=get(gca,'Ylim');
            y_lim=[min(y_lim(1),y_lim2(1)),max(y_lim(2),y_lim2(2))];
        end
    else
        y_lim=yrange;
    end
end
for nResult=length(chosen_result)
    subplot(length(chosen_stim)+1,2*length(chosen_result),2*nResult-1+2*length(chosen_result)*length(chosen_stim));
    xlim([ts(1),ts(end)]);
    plot([0,0],[y_lim(1),y_lim(2)],'k-');%align to a behavior event
    if strcmp(behEventAlign,'stim onset')
        plot(0.5*ones(1,2),[y_lim(1),y_lim(2)],'k-');%align to a behavior event
    elseif strcmp(behEventAlign,'delay onset')
        plot(-0.5*ones(1,2),[y_lim(1),y_lim(2)],'k-');%align to a behavior event
    end
    xlabel(['Time (s) from ',behEventAlign],'FontSize',12);
    if nResult==1
        switch activity_form
            case 'dff'
                ylabel('\it\DeltaF/F');
            case 'spkr'
                ylabel('spikes/s');
            case 'zscored_spkr'
                ylabel('z-scored firing rate');
        end
    end
    %set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
    set(gca,'Ylim',y_lim,'xtick',[-floor(frameNumTime(1)):1:frameNumTime(2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
    %         set(gca,'FontName','Arial','FontSize',14);
end

% label and general plots
% find common ylim
clear y_lim;
for nResult=length(chosen_result)
    subplot(length(chosen_stim)+1,2*length(chosen_result),2*nResult+2*length(chosen_result)*length(chosen_stim));
    flagclearyrange=0;
    if ~exist('yrange','var')
        flagclearyrange=1;
        if ~exist('y_lim','var')
            y_lim=get(gca,'Ylim');
        else
            y_lim2=get(gca,'Ylim');
            y_lim=[min(y_lim(1),y_lim2(1)),max(y_lim(2),y_lim2(2))];
        end
    else
        y_lim=yrange;
    end
end
for nResult=length(chosen_result)
    subplot(length(chosen_stim)+1,2*length(chosen_result),2*nResult+2*length(chosen_result)*length(chosen_stim));
    xlim([ts(1),ts(end)]);
    plot([0,0],[y_lim(1),y_lim(2)],'k-');%align to a behavior event
    if strcmp(behEventAlign,'stim onset')
        plot(0.5*ones(1,2),[y_lim(1),y_lim(2)],'k-');%align to a behavior event
    elseif strcmp(behEventAlign,'delay onset')
        plot(-0.5*ones(1,2),[y_lim(1),y_lim(2)],'k-');%align to a behavior event
    end
    %         xlabel(['time(s) from ',behEventAlign],'FontSize',12);
    if nResult==1
        ylabel('licks/s');
    end
    box off;
    %set(gca,'xtick',[1:round(1000/frT):size(neuralActivity,2)],'xticklabel',[-frameNumTime(1):1:frameNumTime(2)]);
    set(gca,'Ylim',y_lim,'xtick',[-floor(frameNumTime(1)):1:frameNumTime(2)],'xticklabel',[-floor(frameNumTime(1)):1:frameNumTime(2)]);
    %         set(gca,'FontName','Arial','FontSize',14);
end

subplot(length(chosen_stim)+1,2*length(chosen_result),2*length(chosen_result)*(length(chosen_stim)+1));
if contains(trialTypeStr,'stimuli')
    %         h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','','contra hard','','contra easy'},'Location','best');
    h=legend(curve_meanTrace(:),num2str(Data_extract.Stimuli(:)),'Location','best');
elseif contains(trialTypeStr,'difficulty')
    h=legend(curve_meanTrace(:),{'ipsi easy','ipsi hard','contra hard','contra easy'},'Location','best');
elseif contains(trialTypeStr,'choice')
    h=legend(curve_meanTrace(:),{'ipsi choice','contra choice'},'Location','best');
elseif contains(trialTypeStr,'sensory')
    h=legend(curve_meanTrace(:),{'low click rate','high click rate'},'Location','best');
end
set(h,'box','off');


if flagclearyrange==1
    clear y_lim;
end
end

%% helper function

%fuction for getting a round number, example: 166-->100, 25-->10
function [out]=fRound(in)
    n_digit=floor(log10(in));
    out=round(in/10^n_digit)*10^n_digit;
end
