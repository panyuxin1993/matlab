cellFolder='D:\xulab\project\slice recording\XUNL\summary of short term plasticity\example';
cd(cellFolder);
AbfFiles = dir(strcat(cellFolder,'\*.abf'));
nexample=length(AbfFiles);   
colorMean={[41 171 226]/255,[41 171 226]/255,[57 181 74]/255,[57 181 74]/255};
figraw=figure;
% set(gcf,'PaperPosition',[1,1,3,3]);
set(gcf,'Position',[1,1,250,250]);
for n_abf=1:nexample
    %     strcat(cellFolder,AbfFiles(n_abf).name);
    [d,si,h]=abfload(strcat(cellFolder,'\',AbfFiles(n_abf).name),'channels','a');
    t_stim=find(mean(d(:,2,:),3)>1)*h.si/10^6;%time when stim is on
    xlim=[t_stim(1)-0.25,t_stim(end)+0.35];
    ts=xlim(1):h.si/10^6:xlim(2);
    ts_offset=ts-t_stim(1);
    voltage=d(xlim(1)/(h.si/10^6):xlim(2)/(h.si/10^6),1,:);%choose data from setted range
%     for i=1:size(d,3)
%         plot(ts_offset,voltage(:,1,i),'-','Color',[0.8,0.8,0.8],'LineWidth',0.1);%time is offset to stim_onset
%         hold on;
%     end
    figure(figraw);
    subplot(nexample,1,n_abf);
    plot(ts_offset,mean(voltage(:,1,:),3),'Color',colorMean{n_abf},'LineWidth',1);%mean trace
    hold on;
    ylim=get(gca,'ylim');
    %find each pulse start and end
    temp2=[t_stim; nan];
    temp1=[0;t_stim];
    temp3=temp2-temp1;
    ind_stim_start=find(temp3>0.01);
    ind_stim_end=ind_stim_start-1;
    ind_stim_end(1)=[];
    ind_stim_end=[ind_stim_end; length(t_stim)];    
    set(gca,'ylim',[-73,-67]);
    for istim=1:length(ind_stim_start)
        plot([t_stim(ind_stim_start(istim))-t_stim(1),t_stim(ind_stim_end(istim))-t_stim(1)],[-67,-67],'-b','LineWidth',3);%light stimuli
        hold on;
    end
    set(gca,'xlim',[xlim(1)-t_stim(1)-0.2,xlim(2)-t_stim(1)+0.2]);
    xl=get(gca,'xlim');
    set(gca,'ylim',[-73,-67]);
    %set(gca,'ylim',[ylim(1),ylim(end)+1]);
    if n_abf==1
        plot([xl(2)-0.2,xl(2)],[-69,-69],'-k','LineWidth',1);%scale bar for time
        plot([xl(2)-0.2,xl(2)-0.2],[-69,-67],'-k','LineWidth',1);%scale bar for voltage
%         text(1.3,ylim(end)-1.5,'0.2s','FontName','Arial','FontSize',10);
%         text(1.1,ylim(end)-0.5,'2mV','FontName','Arial','FontSize',10);
    end
    if n_abf==nexample
        xlabel('time(s)');
    end
%     set(gca,'xtick',[],'xticklabel',[])
set(gca,'xticklabel',[]);
    box off;
    set(gca,'FontName','Arial','FontSize',10);

end
saveas(figraw,['example.png'],'png');
saveas(figraw,['example.fig'],'fig');
saveas(figraw,['example.eps'],'eps');
saveas(figraw,['example.pdf'],'pdf');