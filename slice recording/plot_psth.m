rootpath='D:\xulab\project\slice recording\data3\2019-4-2\sc\';
cellpath=dir(rootpath);
for n_cell=3:length(cellpath) %first is . ,second is ..
    cellFolder=strcat(rootpath,cellpath(n_cell).name);
    if contains(cellFolder,'abf')%not folder, but abf files
        continue;
    end
% cellFolder='D:\xulab\project\slice recording\XUNL\summary\vglut2-';
    cd(cellFolder);
    AbfFiles = dir(strcat(cellFolder,'\*.abf'));
    for n_abf=1:length(AbfFiles)
        %     strcat(cellFolder,AbfFiles(n_abf).name);
        [d,si,h]=abfload(strcat(cellFolder,'\',AbfFiles(n_abf).name),'channels','a');
        ts=(1:h.sweepLengthInPts)*h.si/10^6;
        xlim=[1.5,2.5];
        nchannel=size(d,2);
        figraw=figure;
        for nsubplot=1:nchannel
            figure(figraw);
            subplot(nchannel,1,nsubplot);
            for i=1:size(d,3)
                plot(ts,d(:,nsubplot,i),'-','Color',[0.7,0.7,0.7],'LineWidth',0.5);
                hold on;
%                 pause(0.1);
            end
            figure(figraw);
            plot(ts,mean(d(:,nsubplot,:),3),'-k','LineWidth',2);
            hold on;
            %         set(gca,'xlim',xlim);
            xlabel('time(s)');
            set(gca,'FontName','Arial','FontSize',14);
            saveas(figraw,[cellpath(n_cell).name,'-',AbfFiles(n_abf).name,'.fig'],'fig');
        end
        close all;
        %     figure(gcf);
        %     subplot(2,1,2);
        %     plot(ts,mean(d(:,2,:),3),'-k','LineWidth',2);
        %     set(gca,'xlim',xlim);
        %     xlabel('time(s)');
        %     set(gca,'FontName','Arial','FontSize',14);
    end
end
