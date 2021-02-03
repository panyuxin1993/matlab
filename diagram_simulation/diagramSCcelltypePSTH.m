%plot the guessed SC neural activities with cell type specificity
%brain stem projecting SC neurons
step=0.001;
x=-1:step:2;%aligned to delay onset
[y_bs,y_e]=deal(x);
t_go=1;
t_delay=0;
t_stim=-0.5;
t_lick=1.1;
y_baseline=0;
y_bs(x<=t_go)=y_baseline;
y_bs(x>=t_lick)=exp(-((x(x>=t_lick)-t_lick-0.1)*2).^2);
n_x=(t_lick-t_go)/step-2;
step_temp=(y_bs(x==t_lick)-y_baseline)/n_x;
y_bs(t_go<x&x<t_lick)=0:step_temp:y_bs(x==t_lick);

figure;
subplot(1,3,1);%brain stem projecting SC neurons
plot(x,y_bs,'r-');
title('Brain stem projecting excitatory neuron');

subplot(1,3,2);


for i=1:3
    subplot(1,3,i);
    set(gca,'Xlim',[-1,2]);
    hold on;
    box off;
    y_lim=get(gca,'Ylim');
    p=patch([t_stim,t_stim,t_delay,t_delay],[y_lim,fliplr(y_lim)],[0.5,0.5,0.5]);
    set(p,'FaceAlpha',0.5);
    plot([t_go,t_go],y_lim,'k-');
    xlabel('Time (s) from delay onset');
end
