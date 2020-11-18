function [ fig ] = fPlotEphysExample( fig,position  )
%FPLOTEPHYSEXAMPLE Summary of this function goes here
%   Detailed explanation goes here
cellFolder='D:\xulab\project\slice recording\data3\examplePSTH\';
celltype={'vglut2+','vgat+'};
% cd(cellFolder);
color_celltype={'F16820','646464'};%
colorMean=fHex2RGB(color_celltype);%{[41 171 226]/255,[57 181 74]/255};

% figure('PaperPosition',[1,1,1,4]);
% set(gcf,'PaperPositionMode','manual','PaperPosition',[1,1,4,2]);
set(fig,'PaperPosition',position);
ratio=[1,1];%[1.3155, 0.4050];
width=0.3347*2*[ratio(1)/(ratio(1)+ratio(2)), ratio(2)/(ratio(1)+ratio(2))];

for i=1:length(celltype)
    AbfFiles = dir(strcat(cellFolder,celltype{i},'\*.abf'));
    nexample=length(AbfFiles);
    for n_abf=nexample:(-1):1
        %     strcat(cellFolder,AbfFiles(n_abf).name);
        [d,si,h]=abfload(strcat(cellFolder,celltype{i},'\',AbfFiles(n_abf).name),'channels','a');
        t_stim=find(mean(d(:,2,:),3)>1)*h.si/10^6;%time when stim is on
        xlim=[t_stim(1)-0.01,t_stim(end)+0.05];
        ts=xlim(1):h.si/10^6:xlim(2);
        ts_offset=ts-t_stim(1);
        voltage=d(xlim(1)/(h.si/10^6):xlim(2)/(h.si/10^6),1,:);%choose data from setted range
        %     for i=1:size(d,3)
        %         plot(ts_offset,voltage(:,1,i),'-','Color',[0.8,0.8,0.8],'LineWidth',0.1);%time is offset to stim_onset
        %         hold on;
        %     end        
        figure(fig);
        subh=subplot(length(celltype),2,n_abf+i*2-2);
%         p=get(gca,'position');
% %         width(1)=p(3)*2/(xlim(2)-xlim(1)+0.4)*(xlim(2)-xlim(1));
% %         width(2)=p(3)*2/(xlim(2)-xlim(1)+0.4)*0.4;
%         if n_abf==2
% %             ax = axes('Parent',gca);
% %             gca.YAxis.Visible = 'off';  
% % set(gca.axes1, 'visible', 'off' );
% set(subh,'Position',[p(1)+(width(1)-width(2))/2-0.05,p(2),width(2),p(4)]);
% % axis off
%         elseif n_abf==1
% %             set(subh,'Position',[p(1),p(2),p(3)+0.1,p(4)]);
%             set(subh,'Position',[p(1),p(2),width(1),p(4)]);
% %             axis off
%         end
        plot(ts_offset,mean(voltage(:,1,:),3),'Color',colorMean{i},'LineWidth',2);%mean trace
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
        %set range for different neurons differently
        if i==1
            set(gca,'ylim',[-70,-60]);
            set(gca,'ytick',[-70,-65,-60],'yticklabel',[-70,-65,-60]);
        else
            set(gca,'ylim',[-73,-63]);
            set(gca,'ytick',[-73,-68,-63],'yticklabel',[-73,-68,-63]);
        end
        if n_abf==2
            set(gca,'ytick',[],'yticklabel',[]);
            set(gca,'xtick',[],'xticklabel',[]);
        end
        ylim=get(gca,'ylim');
        for istim=1:length(ind_stim_start)
            plot([t_stim(ind_stim_start(istim))-t_stim(1),t_stim(ind_stim_end(istim))-t_stim(1)],[ylim(end)-1,ylim(end)-1],'-b','LineWidth',3);%light stimuli
            hold on;
        end
        set(gca,'xlim',[xlim(1)-t_stim(1),xlim(2)-t_stim(1)]);
        xl=get(gca,'xlim');
        yl=get(gca,'ylim');
        %     set(gca,'ylim',[-73,-67]);
        %set(gca,'ylim',[ylim(1),ylim(end)+1]);
        if n_abf<3
            plot([xl(2)-0.01,xl(2)],[yl(2)-3,yl(2)-3],'-k','LineWidth',2);%scale bar for time
            plot([xl(2)-0.01,xl(2)-0.01],[yl(2)-3,yl(2)],'-k','LineWidth',2);%scale bar for voltage
            %         text(1.3,ylim(end)-1.5,'0.2s','FontName','Arial','FontSize',10);
            %         text(1.1,ylim(end)-0.5,'2mV','FontName','Arial','FontSize',10);
            hold on;
        end

        %     set(gca,'xtick',[],'xticklabel',[])
        set(gca,'xticklabel',[]);
        box off;
        set(gca,'FontName','Arial','FontSize',10);
        
    end
end
