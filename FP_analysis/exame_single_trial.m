indTrial=136;
clf;
figure(gcf);
plot(frT*(1:size(dff_aligned,2)),dff_aligned(indTrial,:));
line(behEvent_aligned.go(indTrial)*frT*[1,1],ylim(gca),'color',[0,0,0]);
line(behEvent_aligned.stimOnset(indTrial)*frT*[1,1],ylim(gca),'color',[0,1,0]);
line(behEvent_aligned.lickFirst_left(indTrial)*frT*[1,1],ylim(gca),'color',[0,0,1]);
line(behEvent_aligned.lickFirst_right(indTrial)*frT*[1,1],ylim(gca),'color',[1,0,0]);
if isnan(behEvent_aligned.rewTime(indTrial))
    text(0.2,0.8,'Not correct','Unit','Normalized','color',[1,0,0],'FontSize',20);
else
    text(0.2,0.8,'Correct','Unit','Normalized','color',[0,1,0],'FontSize',23);
end
set(gca,'Ylim',[-0.1,0.3]);
    
