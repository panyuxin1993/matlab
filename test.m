for i=1:4
plot(dff{1,i}(1,:));
hold on;
end
        xlabel('time(frames)');
        ylabel('dff(%)');
legend('dff baseline correction first','dff motion correction first','dff470','dff410');
saveas(gcf,'dff_comparison.fig','fig');