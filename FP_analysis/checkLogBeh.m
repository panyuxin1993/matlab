rootpath='H:\FP\pyx309_20200904';
load([rootpath,filesep,'2020_09_04_pyx309_FP_Virables.mat'])

FrameInfo = dir([rootpath,'\*.log']);
fileID = fopen([rootpath,'\',FrameInfo.name]);
C=textscan(fileID,'%d %d','HeaderLines',16);
%     ImagingSetup
ImagingSetup=FrameInfo.name(1:end-4);
disp(ImagingSetup);
ImagingSetup(ImagingSetup=='_')='-';
fclose(fileID);
TrialCount=C{1,1};
TrialStart_FrameCount=C{1,2};
fiberstr={'Left','Right'};

timeTrialStartBeh=cellfun(@fTrialStartStr2double,Data_extract.Time_trialStart);% time(s) relative to start of session day
framerate=40;
trialStart=double(C{2});
[ trial_start_frame_final,trial_start_time_final ] = fRefineTrialStartRefBeh( trialStart, framerate ,timeTrialStartBeh);
%% auxiliary fuction
function [timeTrialStart]=fTrialStartStr2double(strTrialStart)
%transform the trial start time from string form to double form
t=strsplit(strTrialStart,'_');
t=str2double(t);
timeTrialStart=t(4)*3600+t(5)*60+t(6)+t(7)/1000; %time(s) relative to trial1
end
function [out]=fBaselineCorrection(in,span)
%span = round(20*1000/frT/2); %every 40s, get a mean and substrate to compensate for baseline change
    x = [];
    for i = 1:length(in)
        %     ind_x = ind_x + 1;
        ind1 = max(i- span,1);
        ind2 = min(i+ span,length(in));
        x(i) = prctile(in(ind1:ind2),5);
    end
    out = in - x + nanmean(x);
end