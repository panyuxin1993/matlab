cellFolder='D:\xulab\project\slice recording\XUNL\summary\EPSP example';
cd(cellFolder);
AbfFiles = dir(strcat(cellFolder,'\*.abf'));
for n_abf=1:length(AbfFiles)
    %     strcat(cellFolder,AbfFiles(n_abf).name);
    [d,si,h]=abfload(strcat(cellFolder,'\',AbfFiles(n_abf).name),'channels','a');
    t_stim=find(mean(d(:,2,:),3)>1)*h.si/10^6;%time when stim is on
    xlim=[t_stim(1)-0.05,t_stim(1)+0.1];
    ts=xlim(1):h.si/10^6:xlim(2);
    ts_offset=ts-t_stim(1);
    voltage=d(xlim(1)/(h.si/10^6):xlim(2)/(h.si/10^6),1,:);%choose setted range
    figraw=figure;
    set(gcf,'position',[100,100,300,300]);
    for i=1:size(d,3)
        plot(ts_offset,voltage(:,1,i),'-','Color',[0.8,0.8,0.8],'LineWidth',0.1);%time is offset to stim_onset
        hold on;
    end
    figure(figraw);
    plot(ts_offset,mean(voltage(:,1,:),3),'-k','LineWidth',2);%mean trace
    hold on;
    ylim=get(gca,'ylim');
    plot([0,t_stim(end)-t_stim(1)],[ylim(end),ylim(end)],'-b','LineWidth',5);%light stimuli
%     plot([0.08,0.1],[-77,-77],'-k','LineWidth',3);%scale bar for time
%     plot([0.08,0.08],[-77,-76],'-k','LineWidth',3);%scale bar for voltage
    plot([0.06,0.1],[ylim(end)-0.5,ylim(end)-0.5],'-k','LineWidth',3);%scale bar for time
    plot([0.06,0.06],[ylim(end)-0.5,ylim(end)+0.5],'-k','LineWidth',3);%scale bar for voltage
    set(gca,'xlim',[-0.06,0.11]);
    set(gca,'ylim',[-80,-76]);
    %set(gca,'ylim',[ylim(1),ylim(end)+1]);
    xlabel('time(s)');
    box off;
    set(gca,'FontName','Arial','FontSize',14);
    saveas(figraw,[AbfFiles(n_abf).name,'.png'],'png');
    saveas(figraw,[AbfFiles(n_abf).name,'.fig'],'fig');
end
