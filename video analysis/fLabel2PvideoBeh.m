function [outputVideo] = fLabel2PvideoBeh(videofile,behfile,trialInd,framerate)
%FLABEL2PVIDEOBEH label the video file using info from behfile, including
%stim, delay, response, licking epoch, etc
%Input-
%   videofile- filepath of the video
%   behfile- filepath of the behfile
%   trialInd- index of the trial, since behfile include whole session info
%   framerate- frame rate of the video, used to transform beh event time in
%   frame in the video
%  modified from ref:https://blog.csdn.net/u013921430/java/article/details/79283305
load(behfile);%treat with beh info
stimOn=double(SessionResults{1, trialInd}.Time_stimOnset)*framerate/1000;
stimOff=double(SessionResults{1, trialInd}.Time_stimOffset)*framerate/1000;
delayOn=double(SessionResults{1, trialInd}.Time_delayOnset)*framerate/1000;
delayOff=double(SessionResults{1, trialInd}.Time_delayOffset)*framerate/1000;
answer=double(SessionResults{1, trialInd}.Time_answer)*framerate/1000;
%create output file
indtemp=strfind(videofile,filesep);
save_folder=videofile(1:indtemp(end)-1);
rawfilename=videofile(indtemp(end)+1:end);
file_videoOut=[save_folder,filesep,rawfilename,'_labeled_beh.avi'];
% WriterObj=VideoWriter(file_videoOut,'Uncompressed AVI');%���ϳɵ���Ƶ(��������avi��ʽ)���ļ�·��
WriterObj=VideoWriter(file_videoOut);
WriterObj.FrameRate=framerate;
open(WriterObj);
fileextendtemp=strsplit(videofile,'.');
fileextend=fileextendtemp{end};
if strcmp(fileextend,'avi')
    ReaderObj = VideoReader(videofile);
    %     WriterObj.FrameRate=ReaderObj.FrameRate;%use defined fr
    ReaderObj.CurrentTime=0;%The 1st frame in a video is 0
    while hasFrame(ReaderObj)
        frame = readFrame(ReaderObj);  % ��ȡÿһ֡
        imshow(frame);                 % ��ʾÿһ֡
        i=ReaderObj.CurrentTime*ReaderObj.FrameRate+1;
        if (i>=stimOn) && (i<=stimOff)
            text(0.6,0.9,['stimuli delivery'],'FontName','Arial','FontSize',14,'Unit','normalized','Color','y');
        elseif (i>=delayOn) && (i<=delayOff)
            text(0.6,0.9,['delay'],'FontName','Arial','FontSize',14,'Unit','normalized','Color','y');
        elseif (i<answer) && (i>delayOff)
            text(0.6,0.9,['response'],'FontName','Arial','FontSize',14,'Unit','normalized','Color','y');
        elseif i>answer
            text(0.6,0.9,['licking'],'FontName','Arial','FontSize',14,'Unit','normalized','Color','y');
        end
        frameNew=getframe(gca);
        writeVideo(WriterObj,frameNew);
    end
    
else
    Info=imfinfo(videofile);                                      %%��ȡͼƬ��Ϣ���ж��Ƿ�Ϊtif
    tif='tif';
    format=Info.Format;
    if  (strcmp(format ,tif)==0)
        disp('����Ĳ���tifͼ�������ļ�');                %%ȷ�������ͼ����tiffͼ��
    else
        Slice=size(Info,1);                                          %%��ȡͼƬz��֡��
        Width=Info.Width;
        Height=Info.Height;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Image=zeros(Height,Width,Slice);
        for i=1:Slice
            [Image(:,:,i),map]=imread(videofile,i);           %%һ��һ��Ķ���ͼ��
            imshow(Image(:,:,i),map);
            if (i>=stimOn) && (i<=stimOff)
                text(0.6,0.9,['stimuli delivery'],'FontName','Arial','FontSize',14,'Unit','normalized','Color','y');
            elseif (i>=delayOn) && (i<=delayOff)
                text(0.6,0.9,['delay'],'FontName','Arial','FontSize',14,'Unit','normalized','Color','y');
            elseif (i<answer) && (i>delayOff)
                text(0.6,0.9,['response'],'FontName','Arial','FontSize',14,'Unit','normalized','Color','y');
            elseif i>answer
                text(0.6,0.9,['licking'],'FontName','Arial','FontSize',14,'Unit','normalized','Color','y');
            end
            frameNew=getframe(gca);
            writeVideo(WriterObj,frameNew);
        end
    end
end
close(WriterObj);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:Slice
%     J=uint8(Image(:,:,i));                                   %%һ��һ��д��ͼ�� 
%     %%imwrite(J,[num2str(i,'%4d'),'.tif'],'WriteMode','Append');
%     imwrite(J,[num2str(i,'%04d'),'.tif']);
% end


end

