for i = 1:length(SessionResults)
    choice(i) = double(SessionResults{i}.Action_choice);
    stim(i) = double(SessionResults{i}.Stim_clickRate);
end
inds_miss = choice==2;
inds_violation = choice == 3;
inds_exclude = (inds_miss | inds_violation);
choice(inds_exclude) = [];
stim(inds_exclude) = [];
stim_log = log2(stim./min(stim));
%%
ffun=fittype('a+b./(1+exp(-k.*(x-c)))','independent','x');
    f =fit(stim_log', choice',ffun,'Startpoint',[0,1,1,1]);
    slope=f.b*f.k/4;
    x = linspace(min(stim_log), max(stim_log), 100);
    yfit =f(x);
%     x=2.^x*20;
%%
 figure; plot(x,yfit)
 hold on;
  stim_types = unique(stim_log);
for i = 1:length(stim_types)
    choice_stim_type(i) = mean(choice(stim_log == stim_types(i)));
end
plot(stim_types, choice_stim_type, 'ro','markersize',15)
stim_types = unique(stim_log);
for i = 1:length(stim_types)
    choice_stim_type(i) = mean(choice(stim_log == stim_types(i)));
end
